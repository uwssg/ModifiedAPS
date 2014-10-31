#include "node.h"

void node::set_names(){
    associates.set_name("node_associates");
    boundaryPoints.set_name("node_boundaryPoints");
    basisVectors.set_name("node_basisVectors");
    basisModel.set_name("node_basisModel");
    range_max.set_name("node_range_max");
    range_min.set_name("node_range_min");
    candidates.set_name("node_candidates");
    geographicCenter.set_name("node_geographicCenter");
    centerCandidates.set_name("node_centerCandidates");
    oldCenters.set_name("node_oldCenters");
}

node::node(){
    
    set_names();
    time_ricochet=0.0;
    time_coulomb=0.0;
    time_search=0.0;
    time_bases=0.0;
    last_nAssociates=0;
    last_nBasisAssociates=0;
    center_dex=-1;
    min_dex=-1;
    last_expanded=0;
    activity=1;
    
    ct_search=0;
    ct_ricochet=0;
    ct_coulomb=0;
    ct_bases=0;
    calls_to_bases=0;
    
    farthest_associate=0.0;
    time_penalty=0.5;
    
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
    basisModel.reset();
    range_min.reset();
    range_max.reset();
    candidates.reset();
    
    center_dex=in.center_dex;
    min_dex=in.min_dex;
    last_nAssociates=in.last_nAssociates;
    last_nBasisAssociates=in.last_nBasisAssociates;
    time_ricochet=in.time_ricochet;
    time_coulomb=in.time_coulomb;
    time_search=in.time_search;
    time_bases=in.time_bases;
    time_penalty=in.time_penalty;
    
    ct_search=in.ct_search;
    ct_ricochet=in.ct_ricochet;
    ct_coulomb=in.ct_coulomb;
    ct_bases=in.ct_bases;
    calls_to_bases=in.calls_to_bases;
    
    last_expanded=in.last_expanded;
    activity=in.activity;
    
    int i,j;
    
    farthest_associate=in.farthest_associate;
    
    centerCandidates.reset();
    for(i=0;i<in.centerCandidates.get_dim();i++){
        centerCandidates.set(i,in.centerCandidates.get_data(i));
    }
    
    oldCenters.reset();
    for(i=0;i<in.oldCenters.get_dim();i++){
        oldCenters.set(i,in.oldCenters.get_data(i));
    }
    
    for(i=0;i<in.geographicCenter.get_dim();i++){
        geographicCenter.set(i,in.geographicCenter.get_data(i));
    }
    
    for(i=0;i<in.range_max.get_dim();i++){
        range_max.set(i,in.range_max.get_data(i));
    }
    
    for(i=0;i<in.range_min.get_dim();i++){
        range_min.set(i,in.range_min.get_data(i));
    }
    
    for(i=0;i<in.associates.get_dim();i++){
        associates.add(in.associates.get_data(i));
    }
    
    for(i=0;i<in.boundaryPoints.get_dim();i++){
        boundaryPoints.add(in.boundaryPoints.get_data(i));
    }
    
    for(i=0;i<in.candidates.get_dim();i++){
        candidates.add(in.candidates.get_data(i));
    }
    
    
    for(i=0;i<in.basisModel.get_dim();i++){
        basisModel.set(i,in.basisModel.get_data(i));
    }
    
    if(gg!=NULL){
        if(gg->is_gp_null()==0){
            
            basisVectors.set_cols(gg->get_dim());
            for(i=0;i<gg->get_dim();i++){
                for(j=0;j<gg->get_dim();j++){
                    basisVectors.set(i,j,in.basisVectors.get_data(i,j));
                }
            } 
        }
    }
   
}

node::node(const node &in){
    set_names();
    copy(in);
}

node& node::operator=(const node &in){
    if(this==&in) return *this;
    copy(in);
    
    return *this;
}

void node::flush_candidates(array_1d<int> &out){
    int i;
    for(i=0;i<candidates.get_dim();i++){
        out.add(candidates.get_data(i));
    }
    candidates.reset();
}

double node::get_farthest_associate(){
    return farthest_associate;
}

int node::get_n_associates(){
    return associates.get_dim();
}

void node::set_gpWrapper(gpWrapper *ggin){
    int i,j;
    basisVectors.reset();
    gg=ggin;
    basisVectors.set_dim(gg->get_dim(),gg->get_dim());
    
    for(i=0;i<gg->get_dim();i++){
        basisModel.set(i,1.0);
        for(j=0;j<gg->get_dim();j++){
            if(i==j)basisVectors.set(i,j,1.0);
            else basisVectors.set(i,j,0.0);
        }
    }

}

void node::set_dice(Ran *ddin){
    dice=ddin;
}

void node::set_center_dex(int ii){
    center_dex=ii;
    
    int i;
    
    if(gg!=NULL){
        if(range_max.get_dim()!=gg->get_dim() ||
           range_min.get_dim()!=gg->get_dim()){
    
            for(i=0;i<gg->get_dim();i++){
                range_max.set(i,gg->get_pt(ii,i));
                range_min.set(i,gg->get_pt(ii,i));
                geographicCenter.set(i,gg->get_pt(ii,i));
            }
    
        }
    }
    
    if(min_dex<0 || gg->get_fn(ii)<gg->get_fn(min_dex)){
        min_dex=ii;
    }
    
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
    if(chiout[0]<gg->get_target() && doAssociate==1 && dexout[0]>=0){
        associates.add(dexout[0]);
        dd=gg->distance(dexout[0],center_dex);
        if(dd>farthest_associate){
            farthest_associate=dd;
        }
    }
    
    if(dexout[0]>=0){
        if(min_dex<0 || gg->get_fn(dexout[0])<gg->get_fn(min_dex)){
            min_dex=dexout[0];
        }
    }
    
    int i,rangeChanged=0;
    double ddCenters,ddTest;
    if(chiout[0]<gg->get_target()){
        for(i=0;i<gg->get_dim();i++){
            if(i>=range_max.get_dim() || pt.get_data(i)>range_max.get_data(i)){
                range_max.set(i,pt.get_data(i));
                rangeChanged=1;
            }
            
            if(i>=range_min.get_dim() || pt.get_data(i)<range_min.get_data(i)){
                range_min.set(i,pt.get_data(i));
                rangeChanged=1;
            }
            
        }
        
        if(rangeChanged==1){
            for(i=0;i<gg->get_dim();i++){
                geographicCenter.set(i,0.5*(range_max.get_data(i)+range_min.get_data(i)));
            }
        }
        
        if(dexout[0]>=0 && chiout[0]<0.9*gg->get_chimin()+0.1*gg->get_target()){
            ddCenters=gg->distance(center_dex,geographicCenter);
            ddTest=gg->distance(dexout[0],geographicCenter);
            
            if(ddTest<ddCenters){
                centerCandidates.add(dexout[0]);
                last_expanded=ct_search;
            }
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

int node::bisection(int low, int high){
    return bisection(low,high,0);
}

int node::bisectionAssociate(int low, int high){
    return bisection(low,high,1);
}

int node::bisection(int lowDex, int highDex, int asAssociates){
    if(lowDex<0 || lowDex>=gg->get_pts()){
        return -1;
    }
    
    if(highDex<0 || highDex>=gg->get_pts()){
        return -1;
    }
    
    array_1d<double> lowball,highball;
    double flow,fhigh;
    
    lowball.set_name("node_bisection(int)_lowball");
    highball.set_name("node_bisection(int)_highball");
    
    int i;
    for(i=0;i<gg->get_dim();i++){
        lowball.set(i,gg->get_pt(lowDex,i));
        highball.set(i,gg->get_pt(highDex,i));
    }
    
    flow=gg->get_fn(lowDex);
    fhigh=gg->get_fn(highDex);
    
    return bisection(lowball,flow,highball,fhigh,asAssociates);
    
}

int node::bisection(array_1d<double> &ll, double l, array_1d<double> &hh, double h){
    return bisection(ll,l,hh,h,0);
}

int node::bisectionAssociate(array_1d<double> &ll, double l,
array_1d<double> &hh, double h){
    return bisection(ll,l,hh,h,1);
}

int node::bisection(array_1d<double> &lowball_in, double flow_in, 
array_1d<double> &highball_in, double fhigh_in, int asAssociates){
    
    /*
    lowDex and highDex are the indices of the initial highball and lowball poitns
    
    will return the best point it found
    */
    
    int i;
    array_1d<double> lowball,highball;
    double flow,fhigh;
    
    lowball.set_name("node_bisection_lowball");
    highball.set_name("node_bisection_highball");
    
    fhigh=fhigh_in;
    flow=flow_in;
    for(i=0;i<gg->get_dim();i++){
        lowball.set(i,lowball_in.get_data(i));
        highball.set(i,highball_in.get_data(i));
    }
    
    
    if(flow>gg->get_target()){
        printf("WARNING in node bisection target %e but flow %e\n",
        gg->get_target(),flow);
        
        exit(1);
    }
    
    if(fhigh<gg->get_target()){
        printf("WARNING in node bisection target %e but fhigh %e\n",
        gg->get_target(),fhigh);
        
        exit(1);
    }
    
    if(gg->get_iWhere()==iCoulomb){
        gg->set_iWhere(iNodeBisect);
    }
    
    int iout;
    double bisection_tolerance=0.01*(gg->get_target()-gg->get_chimin());
    
    array_1d<double> trial;
    double ftrial;
    int itrial;
    trial.set_name("node_bisection_trial");
    
    double dd=gg->distance(lowball,highball);
    int ct=0;
    
    iout=-1;
    double target=gg->get_target();
    
    while(ct<100 && dd>1.0e-6 && target-flow>bisection_tolerance){
        
        for(i=0;i<gg->get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        
        if(asAssociates==0){
            evaluateNoAssociate(trial,&ftrial,&itrial);
        }
        else{
            evaluate(trial,&ftrial,&itrial);
        }
        
        if(ftrial<gg->get_target()){
            //printf("changing lowball\n");
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

    if(iout>=0)add_as_boundary(iout);
    
    return iout;
}

int node::coulomb_search(){
    /*
    will return the index of the best point it found
    */
    
    gg->set_iWhere(iCoulomb);
    
    int iout=-1;
    double before=double(time(NULL));
    int ibefore=gg->get_called();
    
    double eps=1.0e-6,tol=0.01*(gg->get_target()-gg->get_chimin());
    
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
            time_coulomb+=double(time(NULL))-before;
            ct_coulomb+=gg->get_called()-ibefore;
            return iout;
            
        }
        
    }
    
    array_1d<double> velocity,acceleration,newpt;
    velocity.set_name("node_coulomb_velocity");
    acceleration.set_name("node_coulomb_acceleration");
    newpt.set_name("node_coulomb_newpt");
    
    double speed,aa,delta=0.1;
    
    /*initialize with a random velocity pointing (mostly) away from the node's center*/
    for(i=0;i<gg->get_dim();i++){
        velocity.set(i,walker.get_data(i)-gg->get_pt(center_dex,i));
        velocity.add_val(i,1.0e-5*(dice->doub()-0.5));
    }
    speed=velocity.normalize();
    
    while(speed<1.0e-10){
        /*in case the vector between the walker and the center is too small*/
        for(i=0;i<gg->get_dim();i++)velocity.set(i,dice->doub());
        speed=velocity.normalize();
    }
    
    speed=0.001;
    for(i=0;i<gg->get_dim();i++){
        velocity.multiply_val(i,speed);
    }
    
    int istep=0,step_max=100,walker_in_bounds=1;
    double dtv,dta,dt,newmu,dx,distance_to_center;
 
    delta=0.1;
    dx=1.0;
    while(dx>1.0e-6 && delta>1.0e-6 && walker_in_bounds==1 && (gg->get_target()-mu>tol || istep<step_max)){
        istep++;
        
        for(i=0;i<gg->get_dim();i++){
            acceleration.set(i,0.0);
        }
        
        for(i=0;i<associates.get_dim()+1;i++){
            /*loop over points and center, adding repulsive force to acceleration vector*/
            if(i<associates.get_dim()){
                for(j=0;j<gg->get_dim();j++){
                    dir.set(j,gg->get_pt(associates.get_data(i),j)-walker.get_data(j));
                }
            }
            else{
                for(j=0;j<gg->get_dim();j++){
                    dir.set(j,gg->get_pt(center_dex,j)-walker.get_data(j));
                }
            }
            
            nn=dir.normalize();
            if(i>=associates.get_dim()){
                distance_to_center=nn;
            }
            
            for(j=0;j<gg->get_dim();j++){
                acceleration.subtract_val(j,dir.get_data(j)/(nn*nn+eps));
            }
        }//loop over points repelling the walker
        
        aa=sqrt(acceleration.get_square_norm());//this does not actually renormalize the acceleration vector
        speed=sqrt(velocity.get_square_norm());//this does not actually renormalize the velocity vector
        
        dtv=delta*distance_to_center/speed;
        dta=delta*speed/aa;
        
        /*dtv is the time step if we want the walker's motion to be small;
        dta is the time step if we want the walker's acceleration to be small*/
        
        if(isnan(dta) && isnan(dtv)){
            printf("WARNING both coulomb steps are nans\n");
            exit(1);
        }
        else if(isnan(dtv) && !(isnan(dta)))dt=dta;
        else if(isnan(dta) && !(isnan(dtv)))dt=dtv;
        else if(dtv<dta) dt=dtv;
        else if(dta<dtv) dt=dta;
        else{
            dt=eps;
        }

        for(i=0;i<gg->get_dim();i++){
            newpt.set(i,walker.get_data(i)+dt*velocity.get_data(i));
        }

        dx=gg->distance(walker,newpt);
        newmu=gg->user_predict(newpt);
        
        /*
        We do not want the coulomb search to carry us outside of the region
        where the GP predicts chisq<chisq_lim
        */
        
        if(newmu<gg->get_target()){
            mu=newmu;
            for(i=0;i<gg->get_dim();i++){
                walker.add_val(i,dt*velocity.get_data(i));
                velocity.add_val(i,dt*acceleration.get_data(i));
            }
            
            for(i=0;i<gg->get_dim() && walker_in_bounds==1;i++){
                if(walker.get_data(i)>gg->get_max(i) || walker.get_data(i)<gg->get_min(i)){
                    walker_in_bounds=0;
                }
            }
            
            
        }
        else{
            delta*=0.5;
        }
        
        
    }//loop over number of coulomb steps
    
    evaluate(walker,&nn,&iout);
    
    time_coulomb+=double(time(NULL))-before;
    ct_coulomb+=gg->get_called()-ibefore;
    
    return iout;
}

int node::ricochet_driver(int istart, array_1d<double> &vstart, array_1d<double> &vout){
    /*
    returns the index of the point found
    */
    
    if(istart<0)return -1;
    
    array_1d<double> gradient,trial;
    gradient.set_name("node_ricochet_gradient");
    trial.set_name("node_ricochet_trial");
    
    int i,j,k;
    double nn;
    
    //seed some points around istart so that the gradient is accurate
    for(i=0;i<gg->get_dim();i++){
        for(j=0;j<gg->get_dim();j++){
            trial.set(j,gg->get_pt(istart,j));
        }
        
        trial.add_val(i,1.0e-7*(gg->get_max(i)-gg->get_min(i)));
        evaluateNoAssociate(trial,&nn,&j);
        
    }
    
    try{
        gg->actual_gradient(istart,gradient);
    }
    catch(int iex){
        printf("    ricochet gradient calculation failed\n");
        return -1;
    }
    
    double dnorm=sqrt(vstart.get_square_norm());
    
    double vdotg,gnorm=gradient.normalize();
    
    vdotg=0.0;
    for(i=0;i<gg->get_dim();i++){
        vdotg+=vstart.get_data(i)*gradient.get_data(i);
    }
    
    /*reflect the incoming velocity about the gradient of chisquared*/
    array_1d<double> velocity;
    velocity.set_name("node_ricochet_velocity");
    for(i=0;i<gg->get_dim();i++){
        velocity.set(i,vstart.get_data(i)-2.0*vdotg*gradient.get_data(i));
    }
    
    double speed=velocity.normalize(),ss,chibest=2.0*chisq_exception;
    double ftrial=2.0*chisq_exception,flow,fhigh;
    array_1d<double> lowball,highball;
    int iLow,iHigh;
    
    lowball.set_name("ricochet_driver_lowball");
    highball.set_name("ricochet_driver_highball");
    
    /*try to find seeds for bisection along the reflected direction*/
    
    flow=gg->get_fn(istart);
    for(i=0;i<gg->get_dim();i++){
        lowball.set(i,gg->get_pt(istart,i));
    }
  
    
    double startingSpeed=vstart.normalize();
    
    ss=0.01*startingSpeed;
    
    int ct=0;
    while(flow>=gg->get_target() && ct<20){
        for(i=0;i<gg->get_dim();i++){
            lowball.set(i,gg->get_pt(istart,i)-ss*vstart.get_data(i));
        }
        
        evaluateNoAssociate(lowball,&flow,&iLow);
        
        ss+=0.01*startingSpeed;
        
        ct++;
    }
    
    if(ct>=20){
        printf("could not find iLow\n");
        return -1;
    }
    
    ct=0;
    ss=2.0*speed;
    fhigh=-2.0*chisq_exception;
    
    while(fhigh<=gg->get_target() && ct<20){
        
        for(i=0;i<gg->get_dim();i++){
            highball.set(i,lowball.get_data(i)+ss*velocity.get_data(i));
        }
        
        evaluateNoAssociate(highball,&fhigh,&iHigh);
        
        ss*=2.0;
        ct++;
    }
    
    if(ct>=20){
        printf("could not find iHigh %e\n",fhigh);
        return -1;
    }
    
    int iout=-1,ii;
    //iout=bisection(lowball,flow,highball,fhigh);spock
    
    /*implement independent bisection because need different tolerance*/
    for(ii=0;ii<15;ii++){
        for(i=0;i<gg->get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        
        evaluateNoAssociate(trial,&ftrial,&j);
        if(ftrial<=gg->get_target()){
            for(i=0;i<gg->get_dim();i++)lowball.set(i,trial.get_data(i));
            flow=ftrial;
            if(j>=0)iout=j;
           
        }
        else{
            for(i=0;i<gg->get_dim();i++)highball.set(i,trial.get_data(i));
            
        }
        
    }
    
    if(iout>=0)add_as_boundary(iout);
    
    for(i=0;i<gg->get_dim();i++){
        vout.set(i,velocity.get_data(i));
    }
    
    return iout;
}


void node::ricochet_search(int iStart, array_1d<double> &dir){
    
    int ibefore=gg->get_called();
    double before=double(time(NULL));
    
    gg->set_iWhere(iRicochet);
    
    double ftrial;
    array_1d<double> trial;
    trial.set_name("node_ricochet_search_trial");
    
    int iEnd,ii,itrial,i,j;
    double dotproduct;
    array_1d<double> vout;
    vout.set_name("node_ricochet_search_vout");
    array_1d<double> distance_traveled;
    array_1d<int> pts_visited;
    
    distance_traveled.set_name("node_ricochet_search_distance_traveled");
    /*distance_traveled is a running total of how far the ricochet has come*/
    
    pts_visited.set_name("node_ricochet_search_pts_visited");    
    /*pts_visited logs the end points of the individual ricochets*/
    
    int iMedian;
    double medianDistance,nn;
    

    dotproduct=1.0;
    distance_traveled.add(0.0);
    pts_visited.add(iStart);
        
    for(ii=0;ii<10*gg->get_dim() && dir.get_square_norm()>1.0e-20 && dotproduct>0.0; ii++){
        try{
            iEnd=ricochet_driver(iStart,dir,vout);
            
        }
        catch (int iex){
            printf("ending Ricochet because of exception\n");
            ii=10*gg->get_dim()+1;
        }
        
        if(iEnd>=0){    
            i=distance_traveled.get_dim();
            distance_traveled.add(distance_traveled.get_data(i-1)+gg->distance(iStart,iEnd));
            pts_visited.add(iEnd);
            
            dotproduct=0.0;
            for(i=0;i<gg->get_dim();i++){
                /*can this ever be negative...?*/
                dotproduct+=vout.get_data(i)*(gg->get_pt(iEnd,i)-gg->get_pt(center_dex,i));
                dir.set(i,gg->get_pt(iEnd,i)-gg->get_pt(iStart,i));
            }
            
            iStart=iEnd;
        }
        else{
            printf("ending Ricochet because iEnd %d\n",iEnd);
            ii=10*gg->get_dim()+1;
        }

    }//loop on ii
    
    printf("ending Ricochet %d %e %e\n",ii,dir.get_square_norm(),dotproduct);
    printf("points visited %d -- %e\n",
    pts_visited.get_dim(),distance_traveled.get_data(distance_traveled.get_dim()-1));
    printf("number of points %d\n\n",gg->get_whereCt(iRicochet));
    
    array_1d<int> neigh;
    array_1d<double> ddneigh;
    
    neigh.set_name("ricochet_neigh");
    ddneigh.set_name("ricochet_ddneigh");
    
    if(pts_visited.get_dim()>1){
        /*
        First do a compass search in the middle of the last ricochet path
        */
        j=pts_visited.get_dim()-1;
        for(i=0;i<gg->get_dim();i++){
            trial.set(i,0.5*(gg->get_pt(pts_visited.get_data(j),i)+gg->get_pt(pts_visited.get_data(j-1),i)));
        }
            
        evaluateNoAssociate(trial,&ftrial,&itrial);
            
        if(itrial>=0)compass_search(itrial);
            
        /*
        Now do a compass search as near to the middle of the full richochet as possible
        */
        for(i=0;i<pts_visited.get_dim();i++){
            nn=fabs(distance_traveled.get_data(i)-0.5*distance_traveled.get_data(distance_traveled.get_dim()-1));
            if(i==0 || nn<medianDistance){
                medianDistance=nn;
                iMedian=i;
            }
        }
            
        if(iMedian!=distance_traveled.get_dim()-1){
            if(iMedian==0){
                for(i=0;i<gg->get_dim();i++){
                    trial.set(i,0.5*(gg->get_pt(pts_visited.get_data(iMedian),i)+gg->get_pt(pts_visited.get_data(iMedian+1),i)));
                }
            }
            else{
                for(i=0;i<gg->get_dim();i++){
                    trial.set(i,0.5*(gg->get_pt(pts_visited.get_data(iMedian),i)+gg->get_pt(pts_visited.get_data(iMedian-1),i)));
                }
            }
                
            evaluateNoAssociate(trial,&ftrial,&itrial);
                
            if(itrial>=0){
                compass_search(itrial);
            }
        }
        
        /*now find the points nearest to the center of each leg, and save them
        as potential nodes themselves*/
        
        for(i=1;i<pts_visited.get_dim();i++){
            for(j=0;j<gg->get_dim();j++){
                trial.set(j,0.5*(gg->get_pt(pts_visited.get_data(i),j)+gg->get_pt(pts_visited.get_data(i-1),j)));
            }
            
            gg->nn_srch(trial,1,neigh,ddneigh);
            
            if(neigh.get_data(0)>=0){
                candidates.add(neigh.get_data(0));
            }
            
        }
            
    }
    
    time_ricochet+=double(time(NULL))-before;
    ct_ricochet+=gg->get_called()-ibefore;
}

void node::compass_search(int istart){
    /*perform a compass search centered on the point designated by istart*/
    
    if(gg==NULL){
        printf("WARNING cannot compass search; gg is null\n");
        exit(1);
    }
    
    if(basisVectors.get_rows()!=gg->get_dim() || basisVectors.get_cols()!=gg->get_dim()){
        printf("WARNING in compass search dim %d but bases %d by %d\n",
        gg->get_dim(),basisVectors.get_rows(),basisVectors.get_cols());
        
        exit(1);
    }
    
    if(gg->get_fn(istart)>=gg->get_target()){
        return;
    }
    
    if(gg->get_iWhere()!=iRicochet){
        gg->set_iWhere(iCompass);
    }
    
    int idim,i,iHigh;
    double ftrial,sgn,scale,flow;
    array_1d<double> trial,lowball;
    
    trial.set_name("node_compass_trial");
    lowball.set_name("node_compass_lowball");
    
    for(i=0;i<gg->get_dim();i++){
        lowball.set(i,gg->get_pt(istart,i));
    }
    flow=gg->get_fn(istart);
    
    for(idim=0;idim<gg->get_dim();idim++){
        for(sgn=-1.0;sgn<1.5;sgn+=2.0){
            for(i=0;i<gg->get_dim();i++)trial.set(i,gg->get_pt(istart,i));
            
            ftrial=-2.0*chisq_exception;
            scale=0.5;
            while(ftrial<=gg->get_target()){
                //scale*=2.0;
                for(i=0;i<gg->get_dim();i++){
                    trial.add_val(i,basisVectors.get_data(idim,i)*scale*sgn);
                }
                
                evaluateNoAssociate(trial,&ftrial,&iHigh);
            }
            
            bisection(lowball,flow,trial,ftrial);
            
        }//loop over sign (direction along basisVector)
    }//loop over dimension (which basisVector we are bisecting along)
}

//double node::basis_error(int ix, array_1d<double> &dx, array_1d<int> &basis_associates){

int node::perturb_bases(array_2d<double> &bases_in, int ix, array_1d<double> &dx, array_2d<double> &bases_out){
    /*
    This will perturb the ixth basisVector by the small vector dx,
    reconstruct the rest of the basisVectors to be perpendicular to
    each other
    */
    
    if(ix>=gg->get_dim() || ix<0){
        printf("WARNING in basis_error ix %d but dim %d\n",
        ix,gg->get_dim());
        
        exit(1);
    }
    
    double nn,tol=1.0e-12;
    int i,j,jx,kx;
    
    bases_out.set_dim(gg->get_dim(),gg->get_dim());
    
    for(i=0;i<gg->get_dim();i++){
        for(j=0;j<gg->get_dim();j++){
            bases_out.set(i,j,bases_in.get_data(i,j));
        }
    }
    
    /*perturb the ixth basis vector by dx*/
    nn=0.0;
    for(i=0;i<gg->get_dim();i++){
        bases_out.add_val(ix,i,dx.get_data(i));
        nn+=bases_out.get_data(ix,i)*bases_out.get_data(ix,i);
    }
    nn=sqrt(nn);
    for(i=0;i<gg->get_dim();i++)bases_out.divide_val(ix,i,nn);
    
    for(jx=ix+1;jx!=ix;){
        if(jx==gg->get_dim())jx=0;
        
        for(kx=ix;kx!=jx;){
            /*make sure that the jxth basis vector is orthogonal to all of the vectors that
            have already been orthogonalized*/
            nn=0.0;
            for(i=0;i<gg->get_dim();i++)nn+=bases_out.get_data(kx,i)*bases_out.get_data(jx,i);
            for(i=0;i<gg->get_dim();i++)bases_out.subtract_val(jx,i,nn*bases_out.get_data(kx,i));
            
            if(kx<gg->get_dim()-1)kx++;
            else kx=0;
        }
        
        nn=0.0;
        for(i=0;i<gg->get_dim();i++){
            nn+=bases_out.get_data(jx,i)*bases_out.get_data(jx,i);
        }
        
        if(nn<1.0e-20){
            printf("WARNING had to abort perturb_bases because of zero-magnitude vector\n");
            return -1;
        }
        
        
        nn=sqrt(nn);
        for(i=0;i<gg->get_dim();i++){
            bases_out.divide_val(jx,i,nn);
        }
        
        if(jx<gg->get_dim()-1)jx++;
        else jx=0;
    }
    
    double normerr,ortherr,err;
    for(i=0;i<gg->get_dim();i++){
        nn=0.0;
        for(j=0;j<gg->get_dim();j++){
            nn+=bases_out.get_data(i,j)*bases_out.get_data(i,j);
        }
        err=fabs(nn-1.0);
        if(i==0 || err>normerr){
            normerr=err;
        }
        
        for(j=i+1;j<gg->get_dim();j++){
            nn=0.0;
            for(jx=0;jx<gg->get_dim();jx++){
                nn+=bases_out.get_data(i,jx)*bases_out.get_data(j,jx);
            }
            nn=fabs(nn);
            if((i==0 && j==1) || nn>ortherr){
                ortherr=nn;
            }
        }
    }
    
    if(normerr>tol || ortherr>1.0e-6){
        printf("WARNING in basis_err normerr %e ortherr %e\n",normerr,ortherr);
        exit(1);
    }
    
    return 1;
}

void node::project_distance(array_1d<double> &p1, int ip, array_2d<double> &bb, array_1d<double> &dd){
    project_distance(ip,p1,bb,dd);
}

void node::project_distance(int ip, array_1d<double> &p2, array_2d<double> &bb, array_1d<double> &dd){
    if(ip<0 || ip>=gg->get_pts()){
        printf("WARNING want to project distance to %d but %d\n",
        ip,gg->get_pts());
        
        exit(1);
    }

    array_1d<double> p1;
    p1.set_name("node_projection_int_p1");
    int i;
    for(i=0;i<gg->get_dim();i++){
        p1.set(i,gg->get_pt(ip,i));
    }
    project_distance(p1,p2,bb,dd);
}

void node::project_distance(int i1, int i2, array_2d<double> &bb, array_1d<double> &dd){
    if(i1<0 || i2<0 || i1>=gg->get_pts() || i2>=gg->get_pts()){
        printf("WARNING want to project distance between %d and %d but %d\n",
        i1,i2,gg->get_pts());
        
        exit(1);
    }

    array_1d<double> p1,p2;
    p1.set_name("node_projection_p1");
    p2.set_name("node_projection_p2");
    int i;
    for(i=0;i<gg->get_dim();i++){
        p1.set(i,gg->get_pt(i1,i));
        p2.set(i,gg->get_pt(i2,i));
    }
    project_distance(p1,p2,bb,dd);
}

void node::project_distance(array_1d<double> &p1, array_1d<double> &p2, array_2d<double> &_bases,
array_1d<double> &dd){
    
    if(_bases.get_rows()!=gg->get_dim() || _bases.get_cols()!=gg->get_dim()){
        printf("WARNING in project_distance bases is %d by %d but dim %d\n",
        _bases.get_rows(),_bases.get_cols(),gg->get_dim());
        
        exit(1);
    }
    
    dd.reset();
    int i,ix;
    for(ix=0;ix<_bases.get_rows();ix++){
        dd.set(ix,0.0);
        for(i=0;i<gg->get_dim();i++){
            dd.add_val(ix,(p1.get_data(i)-p2.get_data(i))*_bases.get_data(ix,i));
        }
    }
}

double node::apply_model(array_1d<double> &pt){
    return apply_model(pt,basisVectors,basisModel);
}

double node::apply_model(array_1d<double> &pt, array_2d<double> &_bases, array_1d<double> &_model){
    double ans=gg->get_fn(min_dex);
    array_1d<double> dd;
    dd.set_name("node_apply_model_dd");
    project_distance(pt,min_dex,_bases,dd);
    int i;
    for(i=0;i<gg->get_dim();i++){
        ans+=_model.get_data(i)*power(dd.get_data(i),2);
    }
    return ans;
}

double node::basis_error_fn(array_2d<double> &trial_bases, array_1d<int> &basis_associates,
array_1d<double> &trial_model){
    /*
    This is the function that basis_error is trying to minimize (with the constraint that all
    of the model coefficients are positive) using a simplex
    */
    
    int i;
    for(i=0;i<trial_model.get_dim();i++){
        if(trial_model.get_data(i)<0.0)return 2.0*chisq_exception;
    }
    
    array_1d<double> vv;
    vv.set_name("node_basis_error_fn_vv");
    double ans=0.0,chisqModel,ratio;
    int ix,ipt;
    for(ix=0;ix<basis_associates.get_dim();ix++){
        ipt=basis_associates.get_data(ix);
        for(i=0;i<gg->get_dim();i++){
            vv.set(i,gg->get_pt(ipt,i));
        }
        chisqModel=apply_model(vv,trial_bases,trial_model);
        
        ans+=fabs(1.0-chisqModel/gg->get_fn(ipt));
    }
    
    if(isnan(ans)){
        printf("WARNING in basis_error_fn ans is nan\n");
        for(ix=0;ix<basis_associates.get_dim();ix++){
            printf("%e\n",gg->get_fn(basis_associates.get_data(ix)));
        }
        exit(1);
    }
    
    return ans/double(basis_associates.get_dim());;

}

double node::basis_error(array_2d<double> &trial_bases, 
array_1d<int> &basis_associates, array_1d<double> &trial_model){    
    /*
    trial_bases is now made up of a bunch of orthonormal vectors which resulted from
    the small perturbation of the original best_bases.  Now we will see how well a
    multi-dimensional parabola on those bases fits the chisquared data we have observed
    */
    
    if(min_dex<0){
        return 2.0*chisq_exception;
    }
    
    /*this will have to do a simplex to find the best model that is positive definite*/
    
    double ans;
    double alpha=1.0,best=0.5,gamma=2.1;
    double fstar,fstarstar,nn;
    int ih,il,i,j,k;
    array_2d<double> pts;
    array_1d<double> pbar,ff,pstar,pstarstar;
    
    pts.set_name("node_basis_error_pts");
    pbar.set_name("node_basis_error_pbar");
    ff.set_name("node_basis_error_ff");
    pstar.set_name("node_basis_error_pstar");
    pstarstar.set_name("node_basis_error_pstarstar");
    
    for(i=0;i<gg->get_dim()+1;i++){
        for(j=0;j<gg->get_dim();j++){
            pts.set(i,j,1000.0*dice->doub());
        }
        
        nn=basis_error_fn(trial_bases,basis_associates,pts(i)[0])
        ff.set(i,nn);
        if(i==0 || nn<ff.get_data(il)){
            il=i;
        }
        
        if(i==0 || nn>ff.get_data(ih)){
            ih=i;
        }
    }
    
    
    int maxAbort,iterations,last_found;
    double oldmin;
    
    iterations=0;
    last_found=0;
    maxAbort=10*gg->get_dim();
    
    int dim=gg->get_dim();
    
    while(iterations-last_found<maxAbort){
        iterations++;
        oldmin=ff.get_data(il);
        
        for(i=0;i<dim;i++){
            pbar.set(i,0.0);
            for(j=0;j<dim+1;j++){
                if(j!=ih){
                    pbar.add_val(i,pts.get_data(j,i));
                }
            }
            pbar.divide_val(i,double(dim));
        }
        fstar=basis_error_fn(trial_bases,basis_associates,pbar);
        
        if(fstar<ff.get_data(ih) && fstar>ff.get_data(il)){
            ff.set(ih,fstar);
            for(i=0;i<dim;i++){
                pts.set(ih,i,pstar.get_data(i));
            }
        }
        else if(fstar<ff.get_data(il)){
            for(i=0;i<dim;i++){
                pstarstar.set(i,gamma*pstar.get_data(i)+(1.0-gamma)*pbar.get_data(i));
            }
            fstarstar=basis_error_fn(trial_bases,basis_associates,pstarstar);
            
            if(fstarstar<ff.get_data(il)){
                for(i=0;i<dim;i++)pst.set(ih,i,pstarstar.get_data(i));
                ff.set(ih,fstarstar);
            }
            else{
                for(i=0;i<dim;i++)pts.set(ih,i,pstar.get_data(i));
                ff.set(ih,fstar);
            }
        }
        
        for(i=0;i<dim+1;i++){
            if(i==0 || ff.get_data(i)<ff.get_data(il))il=i;
            if(i==0 || ff.get_data(i)>ff.get_data(ih))ih=i;
        }
        
        j=1;
        for(i=0;i<dim+1;i++){
            if(fstar<ff.get_data(i) && i!=ih){
                j=0;
            }
        }
        
        if(j==1){
            for(i=0;i<dim;i++){
                pstarstar.set(i,beta*pts.get_data(ih,i)+(1.0-beta)*pbar.get_data(i));
            }
            fstarstar=basis_error_fn(trial_bases,basis_associates,pstarstar);
            
            if(fstarstar<ff.get_data(ih)){
                for(i=0;i<dim;i++)pts.set(ih,i,pstarstar.get_data(i));
                ff.set(ih,fstarstar);
            }
            else{
                for(i=0;i<dim+1;i++){
                    if(i==0 || ff.get_data(i)<ff.get_data(il)){
                        il=i;
                    }
                }
                for(i=0;i<dim+1;i++){
                    if(i!=il){
                        for(j=0;j<dim;j++){
                            nn=0.5*(pts.get_data(i,j)+pts.get_data(il,j));
                            pts.set(i,j,nn);
                        }
                        nn=basis_error_fn(trial_bases,basis_associates,pts(i)[0]);spock
                    }
                }
            }
        }
    }
    
    return ans/double(basis_associates.get_dim());
}

void node::find_bases(){
    /*find best basis vectors for this node*/
    
    if(min_dex<0){
        return;
    }
    
    gg->set_iWhere(iCompass);
    
    if(dice==NULL){
        printf("WARNING cannot call node::find_bases; dice is null\n");
        exit(1);
    }
    
    if(gg==NULL){
        printf("WARNING cannot call node::find_bases; gg is null\n");
        exit(1);
    }
    
    last_nAssociates=associates.get_dim();
    
    int i,j;
    
    array_1d<int> basis_associates;
    basis_associates.set_name("node_basis_associates");
    
    for(i=0;i<associates.get_dim();i++){
        if(gg->get_fn(associates.get_data(i))-gg->get_fn(min_dex)>1.0e-10){
            basis_associates.add(associates.get_data(i));
        }
    }
    
    if(basis_associates.get_dim()<10){
        return;
    }
    
    /*if(basis_associates.get_dim()<100 || 
      (last_nBasisAssociates>0 && basis_associates.get_dim()<last_nBasisAssociates+1000)){
      
        return;
    }*/
    
    printf("finding bases %d\n",basis_associates.get_dim());
    calls_to_bases++;
    last_nBasisAssociates=basis_associates.get_dim();
    
    double before=double(time(NULL));
    int ibefore=gg->get_called();
    
    array_2d<double> bases_best,bases_trial;
    bases_best.set_name("node_find_bases_bases_best");
    bases_trial.set_name("node_find_bases_bases_trial");
    
    bases_best.set_dim(gg->get_dim(),gg->get_dim());
    
    double Ebest,Etrial,Ebest0,lastEbest;
    for(i=0;i<gg->get_dim();i++){
        for(j=0;j<gg->get_dim();j++){
            bases_best.set(i,j,basisVectors.get_data(i,j));
        }
    }
    
    Ebest0=basis_error(bases_best,basis_associates,basisModel);
    Ebest=Ebest0;
    lastEbest=Ebest0;
    
    array_1d<double> trial_model,best_model;
    trial_model.set_name("find_bases_trial_model");
    best_model.set_name("find_bases_best_model");
    
    int ix,changed_bases=0,aborted=0,max_abort=10*gg->get_dim(),total_aborted=0,total_ct=0;
    double stdevlim=1.0e-5/sqrt(double(gg->get_dim()));
    double stdev=0.1/sqrt(double(gg->get_dim()));
    
    array_1d<double> dx;
    dx.set_name("node_find_bases_dx");
    
    int go_on=1;
    
    while(aborted<max_abort && stdev>stdevlim && Ebest>0.01*Ebest0 && go_on==1){

        ix=-1;
        while(ix>=gg->get_dim() || ix<0){
            ix=dice->int32()%gg->get_dim();
        }
        
        for(i=0;i<gg->get_dim();i++)dx.set(i,normal_deviate(dice,0.0,stdev));
        
        total_ct++;
        
        perturb_bases(bases_best,ix,dx,bases_trial);
        Etrial=basis_error(bases_trial,basis_associates,trial_model);
        
        if(Etrial<Ebest){
            aborted=0;
            changed_bases=1;
            for(i=0;i<gg->get_dim();i++){
                best_model.set(i,trial_model.get_data(i));
                for(j=0;j<gg->get_dim();j++){
                    bases_best.set(i,j,bases_trial.get_data(i,j));
                }
            }
            Ebest=Etrial;
        }
        else{
            aborted++;
            total_aborted++;
        }
        
        if(total_ct%(max_abort/2)==0){
            if(total_aborted<(3*total_ct)/4)stdev=1.5*stdev;
            else if(total_aborted>(3*total_ct)/4)stdev=stdev*0.5;
        }
        
        if(total_ct%1000==0){
            /*
            If Ebest has not changed by more than 10% in the last 1000 calls,
            just stop trying
            */
            if((lastEbest-Ebest)/lastEbest<0.1)go_on=0;
            
            lastEbest=Ebest;
        }
        
    }
    
    printf("done finding bases %d %d \n",aborted,total_ct);
    printf("criteria %e <> %e // %e <> %e\n",stdev,stdevlim,Ebest,0.01*Ebest0);
    printf("changed bases %d time %e\n",changed_bases,double(time(NULL))-before);
    printf("\n");
    
    if(changed_bases==1){
        for(i=0;i<gg->get_dim();i++){
            basisModel.set(i,best_model.get_data(i));
            for(j=0;j<gg->get_dim();j++){
                basisVectors.set(i,j,bases_best.get_data(i,j));
            }
        }
        
        compass_search(min_dex);
    }
    
    time_bases+=double(time(NULL))-before;
    ct_bases+=gg->get_called()-ibefore;
}

double node::volume(){
    if(range_max.get_dim()==0 || range_min.get_dim()==0)return 0.0;
    
    double ans=1.0;
    int i;
    
    for(i=0;i<range_max.get_dim();i++){
        ans*=range_max.get_data(i)-range_min.get_data(i);
    }
    
    return ans;
}

double node::get_max(int idim){
    if(range_max.get_dim()<idim) return 1.0;
    return range_max.get_data(idim);
}

double node::get_min(int idim){
    if(range_max.get_dim()<idim) return 0.0;
    return range_min.get_data(idim);
}

int node::search(){
    
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
    
    if(center_dex<0){
        printf("WARNING cannot call node search; center_dex is %d\n",center_dex);
        exit(1);
    }
    
    double before=double(time(NULL));
    int ibefore=gg->get_called();
    double vstart=volume();
    
    int iCoulomb;
    
    iCoulomb=coulomb_search();
    
    int iLow,iHigh,i,j;
    array_1d<double> dir,trial;
    dir.set_name("node_search_dir");
    trial.set_name("node_search_trial");
    double length,ftrial;
    
    while(iCoulomb<0){
        for(i=0;i<gg->get_dim();i++){
            trial.set(i,gg->get_pt(center_dex,i)+0.01*dice->doub()*(gg->get_max(i)-gg->get_min(i)));
        }
        evaluate(trial,&ftrial,&iCoulomb);
    }
    
    if(gg->get_fn(iCoulomb)<gg->get_target()){
        iLow=iCoulomb;
        for(i=0;i<gg->get_dim();i++){
            dir.set(i,gg->get_pt(iCoulomb,i)-gg->get_pt(center_dex,i));
        
            trial.set(i,gg->get_pt(iCoulomb,i));
        }
        length=0.5*dir.normalize();
        
        ftrial=-2.0*chisq_exception;
        while(ftrial<=gg->get_target()){
            
            //length*=2.0;
            
            for(i=0;i<gg->get_dim();i++){
                trial.add_val(i,length*dir.get_data(i));
            }
            
            evaluateNoAssociate(trial,&ftrial,&iHigh);
        }
    }
    else{
        iLow=center_dex;
        iHigh=iCoulomb;
    }
    
    int iBisection;
    int triedBasesOrRicochet=0;
    
    iBisection=bisectionAssociate(iLow,iHigh);
    
    if((associates.get_dim()>last_nAssociates+10 || last_nAssociates==0) && associates.get_dim()>20){
        try{
            triedBasesOrRicochet=1;
            find_bases();
        }
        catch(int iex){}
    }
    
    /*
    Now, of the two points (the end of the Coulomb search and the end of the bisection), choose the one
    that is closest to chisq_limit as the point where we will begin our ricochet search
    */
    int iStart;
    double nn;
    if(time_ricochet+ct_ricochet*time_penalty<0.5*(time_search+ct_search*time_penalty) && associates.get_dim()>100){
        triedBasesOrRicochet=1;
        if(iBisection<0 || fabs(gg->get_fn(iCoulomb)-gg->get_target())<fabs(gg->get_fn(iBisection)-gg->get_target())){
            iStart=iCoulomb;
        }
        else{
            iStart=iBisection; 
        }
    
    
        if(iStart==iCoulomb){
            for(i=0;i<gg->get_dim();i++){
                dir.set(i,gg->get_pt(iCoulomb,i)-gg->get_pt(center_dex,i));
            }
        }
        else{
            if(iBisection<0 || gg->get_fn(iCoulomb)<gg->get_target()){
                /*
                if the Coulomb point was inside the limit, use the direction from
                the Coulomb point to the bisection point as the initial ricochet
                direction
                */
                for(i=0;i<gg->get_dim();i++){
                    dir.set(i,gg->get_pt(iBisection,i)-gg->get_pt(iCoulomb,i));
                }
            }
            else{
                /*
                otherwise, use the direction from the center to the bisection point
                as the ricochet direction
                */
                for(i=0;i<gg->get_dim();i++){
                    dir.set(i,gg->get_pt(iBisection,i)-gg->get_pt(center_dex,i));
                }
            }
        }
    

        if(iStart>=0){
            ricochet_search(iStart,dir);
            
            //now try ricocheting in the opposite direction relative to the center of the node
            for(i=0;i<gg->get_dim();i++){
                dir.set(i,gg->get_pt(center_dex,i)-gg->get_pt(iStart,i));
            }
            iLow=center_dex;
            iHigh=-1;
            ftrial=-2.0*chisq_exception;
            nn=1.0;
            dir.normalize();
            while(ftrial<gg->get_target()){
                for(i=0;i<gg->get_dim();i++){
                    trial.set(i,gg->get_pt(center_dex,i)+nn*dir.get_data(i));
                }
                evaluateNoAssociate(trial,&ftrial,&iHigh);
                nn*=1.5;
            }
            
            if(iHigh>=0){
                iStart=bisection(iLow,iHigh);
            }
            
            if(iStart>=0){
                printf("    trying backwards ricochet\n");
                ricochet_search(iStart,dir);
            }
            
        }
    }
    
    evaluate(geographicCenter,&nn,&i);
    ct_search+=gg->get_called()-ibefore;
    ibefore=gg->get_called();
    
    double vend=volume();
    
    if(triedBasesOrRicochet==1){
        if(vend>vstart*(1.00001)){
            last_expanded=ct_search;
        }
        else{
            find_bases();
            vend=volume();

            if(!(vend>vstart*1.00001) && ct_ricochet>0){
                activity=0;
            }
            else{
                last_expanded=ct_search;
            }
        }
    }
    
    if(centerCandidates.get_dim()>0){
        recenter();
    }
    
    ct_search+=gg->get_called()-ibefore;
    time_search+=double(time(NULL))-before;
}

void node::recenter(){
    
    //printf("\nrecentering\n");
    
    int i,ibest=-1;
    double dd,ddmin,ddCenter,chiTol;
    
    chiTol=0.9*gg->get_chimin()+0.1*gg->get_target();
    ddCenter=gg->distance(center_dex,geographicCenter);
    
    for(i=0;i<centerCandidates.get_dim();i++){
        if(gg->get_fn(centerCandidates.get_data(i))<chiTol){
            dd=gg->distance(centerCandidates.get_data(i),geographicCenter);
            
            if(dd<ddCenter){
                if(ibest<0 || dd<ddmin){
                    ddmin=dd;
                    ibest=centerCandidates.get_data(i);
                }
            }
        }
    }
    
    if(ibest>=0){
        oldCenters.add(center_dex);
        center_dex=ibest;
        activity=1;
        last_expanded=ct_search;
    }
    
    centerCandidates.reset();
}

int node::get_n_oldCenters(){
    return oldCenters.get_dim();
}

int node::get_oldCenter(int dex){
    return oldCenters.get_data(dex);
}

int node::is_it_active(){
    return activity;
}

int node::get_ct_search(){
    return ct_search;
}

int node::get_ct_ricochet(){
    return ct_ricochet;
}

int node::get_ct_coulomb(){
    return ct_coulomb;
}

int node::get_ct_bases(){
    return ct_bases;
}

double node::get_basis_model(int ix){
    if(basisModel.get_dim()!=gg->get_dim()){
        printf("WARNING basis model has %d but dim %d\n",
        basisModel.get_dim(),gg->get_dim());
        
        exit(1);
    }
    
    if(ix<0 || ix>=basisModel.get_dim()){
        printf("WARNING asked for basis model %d but %d\n",
        ix,basisModel.get_dim());
        
        exit(1);
    }
    
    return basisModel.get_data(ix);
}

double node::get_basis(int ix, int iy){
    if(ix>=gg->get_dim() || iy >=gg->get_dim() || ix<0 || iy<0){
        printf("IN NODE_GET_BASIS %d %d but %d\n",
        ix,iy,gg->get_dim());
        exit(1);
    }
    
    return basisVectors.get_data(ix,iy);
}

void node::set_basis(int ix, int iy, double nn){
    if(ix>=gg->get_dim() || iy >=gg->get_dim() || ix<0 || iy<0){
        printf("IN NODE_SET_BASIS %d %d but %d\n",
        ix,iy,gg->get_dim());
        exit(1);
    }
    
    basisVectors.set(ix,iy,nn);
}

int node::get_calls_to_bases(){
    return calls_to_bases;
}

double node::get_time(){
    return time_search;
}

void node::set_time(double xx){
    time_search=xx;
}

double node::get_time_coulomb(){
    return time_coulomb;
}

double node::get_time_ricochet(){
    return time_ricochet;
}

double node::get_time_bases(){
    return time_bases;
}

Ran* node::get_Ran(){
    return dice;
}

gpWrapper* node::get_gpWrapper(){
    return gg;
}

int node::get_center(){
    return center_dex;
}

void node::set_farthest_associate(double xx){
    farthest_associate=xx;
}

int node::get_n_boundary(){
    return boundaryPoints.get_dim();
}

int node::get_boundary(int dex){
    if(dex>=boundaryPoints.get_dim() || dex<0){
        printf("node::get_boundary asked for %d but only have %d\n",
        dex,boundaryPoints.get_dim());
        
        exit(1);
    }
    return boundaryPoints.get_data(dex);
}

///////////////arrayOfNodes code below//////////

arrayOfNodes::arrayOfNodes(){
    data=NULL;
    ct=0;
    room=0;
}

arrayOfNodes::~arrayOfNodes(){
    if(data!=NULL){
        delete [] data;
    }
}

void arrayOfNodes::add(int cc, Ran *dice, gpWrapper *gg){

    node *buffer;
    int i,j;
    
    if(ct==room){
        if(ct>0){
            buffer=new node[ct];
            for(i=0;i<ct;i++){
                buffer[i].copy(data[i]);
            }
            
            delete [] data;
        }
        
        room+=5;
        data=new node[room];
        
        if(ct>0){
            for(i=0;i<ct;i++){
                data[i].copy(buffer[i]);
            }
            delete [] buffer;
        }
    }
    
    data[ct].set_gpWrapper(gg);
    data[ct].set_center_dex(cc);
    data[ct].set_dice(dice);
    
    ct++;

}

void arrayOfNodes::add(int i, gpWrapper *g, Ran *d){
    add(i,d,g);
}

void arrayOfNodes::add(gpWrapper *g, Ran *d, int i){
    add(i,d,g);
}

void arrayOfNodes::add(gpWrapper *g, int i, Ran *d){
    add(i,d,g);
}

void arrayOfNodes::add(Ran *d, gpWrapper *g, int i){
    add(i,d,g);
}

void arrayOfNodes::add(Ran *d, int i, gpWrapper *g){
    add(i,d,g);
}

int arrayOfNodes::get_dim(){
    return ct;
}

void arrayOfNodes::remove(int ii){
    int i;
    if(ii>=ct) return;
    
    for(i=ii+1;i<ct;i++){
        data[i-1].copy(data[i]);
    }
    ct--;
}

node* arrayOfNodes::operator()(int dex){
    if(dex<0 || dex>=ct){
        printf("WARNING asked arrayOfNodes for dex %d but only have %d (operator)\n",
        dex,ct);
        
        exit(1);
    }
    
    return &data[dex];
}
