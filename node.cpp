#include "node.h"

node::node(){
    associates.set_name("node_associates");
    boundaryPoints.set_name("node_boundaryPoints");
    basisVectors.set_name("node_basisVectors");
    
    time_ricochet=0.0;
    time_coulomb=0.0;
    time_total=0.0;
    last_set_bases=0;
    center_dex=-1;
    
    farthest_associate=0.0;
    
    gg=NULL;
    dice=NULL;
}

node::~node(){}

void node::copy(const node &in){

    gg=in.gg;
    dice=in.dice;
    
    associates.reset();
    boundaryPoints.reset();
    basisVectors.reset();
    
    center_dex=in.center_dex;
    last_set_bases=in.last_set_bases;
    time_ricochet=in.time_ricochet;
    time_coulomb=in.time_coulomb;
    time_total=in.time_total;
    
    int i,j,dim;
    
    farthest_associate=in.farthest_associate;
    
    for(i=0;i<in.associates.get_dim();i++){
        associates.add(in.associates.get_data(i));
    }
    
    for(i=0;i<in.boundaryPoints.get_dim();i++){
        boundaryPoints.add(in.boundaryPoints.get_data(i));
    }
    
    if(gg!=NULL){
        if(gg->is_gp_null()==0){
            dim=gg->get_dim();
            
            basisVectors.set_cols(dim);
            for(i=0;i<dim;i++){
                for(j=0;j<dim;j++){
                    basisVectors.set(i,j,in.basisVectors.get_data(i,j));
                }
            } 
        }
    }
   
}

node::node(const node &in){
    copy(in);
}

node& node::operator=(const node &in){
    if(this==&in) return *this;
    copy(in);
    
    return *this;
}

void node::set_dice(Ran *ddin){
    dice=ddin;
}

void node::evaluate(array_1d<double> &pt, double *chiout, int *dexout){
    evaluate(pt,chiout,dexout,1);
}

void node::evaluateNoAssociate(array_1d<double> &pt, double *chiout, int *dexout){
    evaluate(pt,chiout,dexout,0);
}

void node::evaluate(array_1d<double> &pt, double *chiout, int *dexout, int doAssociate){
    if(gg==NULL){
        printf("WARNING cannot call node::evaluate; gg is NULL\n");
    }
    
    gg->evaluate(pt,chiout,dexout);
    
    double dd;
    if(chiout[0]<gg->get_target() && doAssociate==1){
        associates.add(dexout[0]);
        dd=gg->distance(dexout[0],center_dex);
        if(dd>farthest_associate){
            farthest_associate=dd;
        }
    }
    
}

void node::project_to_unit_sphere(array_1d<double> &in, array_1d<double> &out){
    
    array_1d<double> dir;
    out.reset();
    double norm=0.0;
    int i;
    for(i=0;i<gg->get_dim();i++){
        dir.set(i,in.get_data(i)-gg->get_pt(center_dex,i));
        norm+=power(dir.get_data(i)/(gg->get_max(i)-gg->get_min(i)),2);
    }
    norm=sqrt(norm);
    for(i=0;i<gg->get_dim();i++){
        out.set(i,gg->get_pt(center_dex,i)+dir.get_data(i)/norm);
    }
    
}

void node::add_as_boundary(int dex){

    array_1d<double> vv,vv_norm;
    int i;
    for(i=0;i<gg->get_dim();i++)vv.set(i,gg->get_pt(dex,i));
    project_to_unit_sphere(vv,vv_norm);
    
    gg->add_to_unitSpheres(vv_norm);
    boundaryPoints.add(dex);    
}

int node::bisection(int lowDex, int highDex){
    
    /*
    lowDex and highDex are the indices of the initial highball and lowball poitns
    
    will return the best point it found
    */
    
    int iout;
    double flow,fhigh;
    array_1d<double> lowball,highball;
    lowball.set_name("node_bisection_lowball");
    highball.set_name("node_bisection_highball");
    
    int i;
    
    flow=gg->get_fn(lowDex);
    for(i=0;i<gg->get_dim();i++)lowball.set(i,gg->get_pt(lowDex,i));
    fhigh=gg->get_fn(highDex);
    for(i=0;i<gg->get_dim();i++)highball.set(i,gg->get_pt(highDex,i));
    
    array_1d<double> trial;
    double ftrial;
    int itrial;
    trial.set_name("node_bisection_trial");
    
    double dd=gg->distance(lowball,highball);
    int ct;
    
    iout=lowDex;
    
    while(ct<100 && dd>1.0e-6){
        
        for(i=0;i<gg->get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        
        evaluateNoAssociate(trial,&ftrial,&itrial);
        
        if(ftrial<gg->get_target()){
            if(itrial>=0)iout=itrial;
            flow=ftrial;
            for(i=0;i<gg->get_dim();i++)lowball.set(i,trial.get_data(i));
        }
        else{
            fhigh=ftrial;
            for(i=0;i<gg->get_dim();i++)highball.set(i,trial.get_data(i));
        }
        
        ct++;
        dd*=0.5;
    }
    
    if(iout!=lowDex)add_as_boundary(iout);
    
    return iout;
}

int node::coulomb_search(){
    /*
    will return the index of the best point it found
    */
    int iout=-1;
    double before=double(time(NULL));
    double eps=1.0e-6,tol=1.0;
    
    gg->reset_cache(); //so that results of previous GP interpolation are not carried over
    array_1d<double> av_dir,dir;
    av_dir.set_name("node_coulomb_search_av_dir");
    dir.set_name("node_coulomb_searc_dir");
    
    int i,j;
    for(i=0;i<gg->get_dim();i++)av_dir.set(i,0.0);
    
    /*calculate the negative of the average direction from the node's center to its associates*/
    for(i=0;i<associates.get_dim();i++){
        for(j=0;j<gg->get_dim();j++){
            dir.set(j,gg->get_pt(associates.get_data(i),j)-gg->get_pt(center_dex,j));
        }
        dir.normalize();
        
        for(j=0;j<gg->get_dim();j++){
            av_dir.subtract_val(j,dir.get_data(j));
        }
    }
    
    double nn=av_dir.get_square_norm();
    
    if(associates.get_dim()==0 || nn<eps){
        for(j=0;j<gg->get_dim();j++){
            av_dir.set(j,dice->doub());
        }
    }
    
    av_dir.normalize();
    
    double mu;
    array_1d<double> walker;
    walker.set_name("node_coulomb_search_walker");
    
    if(farthest_associate>0.0){
        nn=20.0*farthest_associate;
    }
    else{
        nn=20.0;
    }
    
    int ct_abort=0,abort_max=30;
    
    /*
    First try to start the search at a point that is much farther away than
    the farthest associate.
    
    Note that we want to start at a point for which the GP predict chisq<chisq_lim
    */
    mu=2.0*chisq_exception;
    while(mu>gg->get_target() && ct_abort<abort_max){
        gg->reset_cache();
        nn*=0.5;
        for(j=0;j<gg->get_dim();j++){
            walker.set(j,gg->get_pt(center_dex,j)+nn*av_dir.get_data(j));
        }
        mu=gg->user_predict(walker);
        ct_abort++;
    }
    
    if(ct_abort==abort_max){
        /*
        If we have failed to find a point for which the GP predicts chisq<chisq_lim,
        try to initialize the search from a random point that is very near to the
        center
        */
        ct_abort=0;
        mu=2.0*chisq_exception;
        while(mu>gg->get_target() && ct_abort<1000){
            ct_abort++;
            for(i=0;i<gg->get_dim();i++)av_dir.set(i,dice->doub());
            av_dir.normalize();
            for(i=0;i<gg->get_dim();i++){
                walker.set(i,gg->get_pt(center_dex,i)+eps*av_dir.get_data(i)*(gg->get_max(i)-gg->get_min(i)));
            }
            gg->reset_cache();
            mu=gg->user_predict(walker);
        }
        
        if(ct_abort==1000){
            /*
            If we still have not been able to find a point for which the GP predict
            chisq<chisq_lim, pick a random point, evaluate chisq at that point and
            abort the searc
            */
        
            for(i=0;i<gg->get_dim();i++){
                walker.set(i,gg->get_pt(center_dex,i)+0.001*av_dir.get_data(i)*(gg->get_max(i)-gg->get_min(i)));
            }
            
            evaluateNoAssociate(walker,&nn,&iout);
            return iout;
            
        }
        
    }
    
    time_coulomb+=double(time(NULL))-before;
    return iout;
}

void node::search(int *out){
    if(gg==NULL){
        printf("WARNING cannot call node search; gpWrapper is null\n");
        exit(1);
    }
    
    if(gg->is_gp_null()==1){
        printf("WARNING cannot call node search; gaussian_process in gpWrapper is null\n");
        exit(1);
    }
    
    if(dice==NULL){
        printf("WARNING cannot call node search; dice is null\n");
        
        exit(1);
    }
}
