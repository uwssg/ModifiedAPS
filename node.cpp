#include "node.h"

void node::set_names(){
    associates.set_name("node_associates");
    boundaryPoints.set_name("node_boundaryPoints");
    basisVectors.set_name("node_basisVectors");
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
    
    center_dex=in.center_dex;
    last_nAssociates=in.last_nAssociates;
    last_nBasisAssociates=in.last_nBasisAssociates;
    time_ricochet=in.time_ricochet;
    time_coulomb=in.time_coulomb;
    time_search=in.time_search;
    time_bases=in.time_bases;
    time_penalty=in.time_penalty;
    
    int i,j;
    
    farthest_associate=in.farthest_associate;
    
    for(i=0;i<in.associates.get_dim();i++){
        associates.add(in.associates.get_data(i));
    }
    
    for(i=0;i<in.boundaryPoints.get_dim();i++){
        boundaryPoints.add(in.boundaryPoints.get_data(i));
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
    
    lowball.set_name("bisection_lowball");
    highball.set_name("bisection_highball");
    
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

int node::bisection(array_1d<double> &lowball, double flow, 
array_1d<double> &highball, double fhigh, int asAssociates){
    
    /*
    lowDex and highDex are the indices of the initial highball and lowball poitns
    
    will return the best point it found
    */
    
    
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
    int i;
    double bisection_tolerance=0.1*gg->get_delta_chisquared();
    
    array_1d<double> trial;
    double ftrial;
    int itrial;
    trial.set_name("node_bisection_trial");
    
    double dd=gg->distance(lowball,highball);
    int ct;
    
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
            time_coulomb+=double(time(NULL))-before;
            time_coulomb+=time_penalty*(gg->get_called()-ibefore);
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
    
    int istep=0,step_max=100;
    double dtv,dta,dt,newmu,dx,distance_to_center;
    
    delta=0.1;
    dx=1.0;
    while(dx>1.0e-20 && delta>1.0e-10 && (gg->get_target()-mu>tol || istep<step_max)){
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
        }
        else{
            delta*=0.5;
        }
        
        
    }//loop over number of coulomb steps
    
    evaluate(walker,&nn,&iout);
    
    time_coulomb+=double(time(NULL))-before;
    time_coulomb+=time_penalty*(gg->get_called()-ibefore);
    return iout;
}

int node::ricochet_driver(int istart, array_1d<double> &vstart, array_1d<double> &vout){
    /*
    returns the index of the point found
    */
    
    double before=double(time(NULL));
    int ibefore=gg->get_called();
    
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
        time_ricochet+=double(time(NULL))-before;
        time_ricochet+=time_penalty*(gg->get_called()-ibefore);
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
    int iLow,iHigh;
    
    /*try to find seeds for bisection along the reflected direction*/
    flow=2.0*chisq_exception;
    ss=2.0*speed;
    int ct=0;
    while(flow>=gg->get_target() && ct<20){
        for(i=0;i<gg->get_dim();i++){
            trial.set(i,gg->get_pt(istart,i)+ss*velocity.get_data(i));
        }
        
        evaluateNoAssociate(trial,&flow,&iLow);
        
        if(flow>=gg->get_target()){
            ss*=0.5;
        }
        
        ct++;
    }
    
    if(ct>=20){
        time_ricochet+=double(time(NULL))-before;
        time_ricochet+=time_penalty*(gg->get_called()-ibefore);
        return -1;
    }
    
    ct=0;
    ss+=0.5*speed;
    fhigh=-2.0*chisq_exception;
    while(fhigh<=gg->get_target() && ct<20){
        for(i=0;i<gg->get_dim();i++){
            trial.set(i,gg->get_pt(istart,i)+ss*velocity.get_data(i));
        }
        
        evaluateNoAssociate(trial,&fhigh,&iHigh);
        
        /*if(fhigh<=gg->get_target()){
            ss*=2.0;
        }*/
        ct++;
    }
    
    if(ct>=20){
        time_ricochet+=double(time(NULL))-before;
        time_ricochet+=time_penalty*(gg->get_called()-ibefore);
        return -1;
    }
    
    int iout;
    iout=bisection(iLow,iHigh);
    
    for(i=0;i<gg->get_dim();i++){
        vout.set(i,velocity.get_data(i));
    }
    
    time_ricochet+=double(time(NULL))-before;
    time_ricochet+=time_penalty*(gg->get_called()-ibefore);
    return iout;
}


void node::ricochet_search(int iStart, array_1d<double> &dir){
    
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
            ii=10*gg->get_dim()+1;
        }

    }//loop on ii
        
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
            
    }
        
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
    double ftrial,sgn,scale;
    array_1d<double> trial;
    trial.set_name("node_compass_trial");
    
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
            
            bisection(istart,iHigh);
            
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


double node::basis_error(array_2d<double> &trial_bases, array_1d<int> &basis_associates){    
    /*
    trial_bases is now made up of a bunch of orthonormal vectors which resulted from
    the small perturbation of the original best_bases.  Now we will see how well a
    multi-dimensional parabola on those bases fits the chisquared data we have observed
    */
    
    int i,j,k,jx;
    double nn;
    
    array_1d<double> matrix,bb;
    matrix.set_name("node_basis_error_matrix");
    bb.set_name("node_basis_error_bb");
    
    array_2d<double> dd;
    dd.set_name("node_basis_error_dd");
    
    matrix.set_dim(gg->get_dim()*gg->get_dim());
    bb.set_dim(gg->get_dim());
    dd.set_dim(basis_associates.get_dim(),gg->get_dim());
    
    for(i=0;i<basis_associates.get_dim();i++){
        k=basis_associates.get_data(i);
        for(j=0;j<gg->get_dim();j++){
            nn=0.0;
            for(jx=0;jx<gg->get_dim();jx++){
                nn+=(gg->get_pt(k,jx)-gg->get_pt(center_dex,jx))*trial_bases.get_data(j,jx);
            }
            dd.set(i,j,nn*nn);
            
            if(isnan(dd.get_data(i,j))){
                printf("WARNING in basis_error dd is nan\n");
                exit(1);
            }
        }
    }
    
    for(i=0;i<gg->get_dim();i++){
        for(j=0;j<gg->get_dim();j++){
            matrix.set(i*gg->get_dim()+j,0.0);
            for(jx=0;jx<basis_associates.get_dim();jx++){
                k=basis_associates.get_data(jx);
                if(gg->get_fn(k)>gg->get_fn(center_dex)){
                    matrix.add_val(i*gg->get_dim()+j,dd.get_data(jx,i)*dd.get_data(jx,j)/power(gg->get_fn(k)-gg->get_fn(center_dex),2));
                }
            }
        }
    }
    
    for(i=0;i<gg->get_dim();i++){
        bb.set(i,0.0);
        for(jx=0;jx<basis_associates.get_dim();jx++){
            k=basis_associates.get_data(jx);
            if(gg->get_fn(k)>gg->get_fn(center_dex)){
                bb.add_val(i,dd.get_data(jx,i)/(gg->get_fn(k)-gg->get_fn(center_dex)));
            }
        }
    }
    
    array_1d<double> trial_model;
    trial_model.set_name("node_basis_error_trial_model");
    
    try{
        naive_gaussian_solver(matrix,bb,trial_model,gg->get_dim());
    }
    catch(int iex){
    }
    
    double ans=0.0;
    for(jx=0;jx<basis_associates.get_dim();jx++){
        k=basis_associates.get_data(jx);
        nn=gg->get_fn(k)-gg->get_fn(center_dex);
        for(i=0;i<gg->get_dim();i++){
            nn-=trial_model.get_data(i)*dd.get_data(jx,i);
        }
        nn=nn/(gg->get_fn(k)-gg->get_fn(center_dex));
        ans+=nn*nn;
    }
    
    if(isnan(ans)){
        printf("WARNING in basis_error ans is nan\n");
        for(jx=0;jx<basis_associates.get_dim();jx++){
	    k=basis_associates.get_data(jx);
	    nn=gg->get_fn(k)-gg->get_fn(center_dex);
	    if(isnan(nn) || nn<1.0e-10){
	        printf("ggfn %e chisq %e %e\n",gg->get_fn(k),gg->get_fn(center_dex),nn);
	    }
	}
     
        return 2.0*chisq_exception;
    
    }
    
    return ans/double(basis_associates.get_dim());
}

void node::find_bases(){
    /*find best basis vectors for this node*/
    
    printf("finding bases\n");
    
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
        if(gg->get_fn(associates.get_data(i))-gg->get_fn(center_dex)>1.0e-10){
            basis_associates.add(associates.get_data(i));
        }
    }
    
    if(basis_associates.get_dim()<100 || 
      (last_nBasisAssociates>0 && basis_associates.get_dim()<last_nBasisAssociates+1000)){
      
        return;
    }
    
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
    
    Ebest0=basis_error(bases_best,basis_associates);
    Ebest=Ebest0;
    lastEbest=Ebest0;
    
    int ix,changed_bases=0,aborted=0,max_abort=10*gg->get_dim(),total_aborted=0,total_ct=0;
    double stdevlim=1.0e-5/sqrt(double(gg->get_dim()));
    double stdev=0.1/sqrt(double(gg->get_dim()));
    
    array_1d<double> dx;
    dx.set_name("node_find_bases_dx");

    while(aborted<max_abort && stdev>stdevlim && Ebest>0.01*Ebest0){

        ix=-1;
        while(ix>=gg->get_dim() || ix<0){
            ix=dice->int32()%gg->get_dim();
        }
        
        for(i=0;i<gg->get_dim();i++)dx.set(i,normal_deviate(dice,0.0,stdev));
        
        total_ct++;
        
        perturb_bases(bases_best,ix,dx,bases_trial);
        Etrial=basis_error(bases_trial,basis_associates);
        
        if(Etrial<Ebest){
            aborted=0;
            changed_bases=1;
            for(i=0;i<gg->get_dim();i++){
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
            if((lastEbest-Ebest)/lastEbest<0.1)stdev=-1.0;
            
            lastEbest=Ebest;
        }
        
    }
    
    if(changed_bases==1){
        for(i=0;i<gg->get_dim();i++){
            for(j=0;j<gg->get_dim();j++){
                basisVectors.set(i,j,bases_best.get_data(i,j));
            }
        }
        
        compass_search(center_dex);
    }
    
    time_bases+=double(time(NULL))-before;
    time_bases+=time_penalty*(gg->get_called()-ibefore);
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
    
    iBisection=bisectionAssociate(iLow,iHigh);
    
    if((associates.get_dim()>last_nAssociates+300 || last_nAssociates==0) && associates.get_dim()>200){
        try{
            find_bases();
        }
        catch(int iex){}
    }
    
    /*
    Now, of the two points (the end of the Coulomb search and the end of the bisection), choose the one
    that is closest to chisq_limit as the point where we will begin our ricochet search
    */
    int iStart;
    if(time_ricochet<0.5*time_search && associates.get_dim()>100){
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
    

        ricochet_search(iStart,dir);
    }
    
    time_search+=double(time(NULL))-before;
    time_search+=time_penalty*(gg->get_called()-ibefore);
}

double node::get_time(){
    return time_search;
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
                buffer->copy(data[i]);
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
    
    data[ct].set_center_dex(cc);
    data[ct].set_gpWrapper(gg);
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
