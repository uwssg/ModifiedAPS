#include "aps.h"

enum{iWIDE,iFOCUS};

straddle_parameter::straddle_parameter(){
    target=-1.0;
}

straddle_parameter::~straddle_parameter(){}

void straddle_parameter::set_target(double tt){
    target=tt;
} 

double straddle_parameter::get_target(){
    return target;
}

double straddle_parameter::operator()(double mu, double sig) const{
    if(target<0.0){
        printf("WARNING target %e in straddle parameter \n",target);
        exit(1);
    }
    
    if(isnan(mu) || isnan(target))return -1.0*chisq_exception;
    
    return sig-fabs(mu-target);
}

gpWrapper::gpWrapper(){
    chimin=-1.0;
    target_asserted=0;
    global_mindex=-1;
    
    chisq=NULL;
    gg=NULL;
    
    good_max.set_name("good_max");
    good_min.set_name("good_min");
    minpt.set_name("minpt");
    good_pts.set_name("good_pts");
}

gpWrapper::~gpWrapper(){}

void gpWrapper::set_gp(gp *gg_in){
    gg=gg_in;
    good_max.set_dim(gg->get_dim());
    good_min.set_dim(gg->get_dim());
    minpt.set_dim(gg->get_dim());
}

void gpWrapper::set_chisq(chisquared *cc){
    chisq=cc;
}

void gpWrapper::set_strad(straddle_parameter *strad_in){
    strad=strad_in;
}

int gpWrapper::in_bounds(array_1d<double> &pt){
    int i;
    for(i=0;i<gg->get_dim();i++){
        if(pt.get_data(i)>chisq->get_max(i) && chisq->get_max(i)>-1.0*chisq_exception){
            return 0;
        }
        
        if(pt.get_data(i)<chisq->get_min(i) && chisq->get_min(i)<chisq_exception){
            return 0;
        }
    }
    
    return 1;
    
    
}

int gpWrapper::is_valid(array_1d<double> &pt){
    double xx;
    return is_valid(pt,&xx);
}

int gpWrapper::is_valid(array_1d<double> &pt, double *chiout){
    
    chiout[0]=-1.0;
    
    if(in_bounds(pt)==0){
        chiout[0]=2.0*chisq_exception;
        return 0;
    }
    
    array_1d<int> neigh;
    array_1d<double> ddneigh;
    
    gg->nn_srch(pt,1,neigh,ddneigh);
    if(ddneigh.get_data(0)<=1.0e-8){
        chiout[0]=gg->get_fn(neigh.get_data(0));
        return 0;
    }
    
    return 1;
}

void gpWrapper::evaluate(array_1d<double> &pt, double *chiout){
    int i;
    evaluate(pt,chiout,&i,-1);
}

void gpWrapper::evaluate(array_1d<double> &pt, double *chiout, int *dex){
    evaluate(pt,chiout,dex,-1);
}

void gpWrapper::evaluate(array_1d<double> &pt, double *chiout, int *dex, int validity){
    
    if(chisq==NULL){
        printf("WARNING chisq is NULL but you just called evaluate\n");
        exit(1);
    }
    
    dex[0]=-1;
    
    if(validity==-1){
        validity=is_valid(pt,chiout);
    }
    
    if(validity==0){
        return;
    }
    else{
        chiout[0]=(*chisq)(pt);
        if(chiout[0]<chisq_exception){
            /*
            Add the point to the Gaussian Process
            */
            add_pt(pt,chiout[0]);
            dex[0]=gg->get_pts()-1;
        }
    }
}

int gpWrapper::add_pt(array_1d<double> &vv, double chitrue){
    
    int i;

    
    gg->add_pt(vv,chitrue);
        
        
    if(chitrue<chimin || chimin<0.0){
            set_chimin(chitrue,vv,gg->get_pts()-1);
    }

    if(chitrue<strad->get_target()){
        good_pts.add(gg->get_pts()-1);
        if(good_pts.get_dim()==0){
            for(i=0;i<gg->get_dim();i++){
                good_max.set(i,vv.get_data(i));
                good_min.set(i,vv.get_data(i));
            }
        }
        else{
            for(i=0;i<gg->get_dim();i++){
                if(vv.get_data(i)<good_min.get_data(i))good_min.set(i,vv.get_data(i));
                if(vv.get_data(i)>good_max.get_data(i))good_max.set(i,vv.get_data(i));
            }
        }
        
    }

}


void gpWrapper::set_chimin(double cc,array_1d<double> &pt, int dex){
    chimin=cc;
    
    int i;
    for(i=0;i<pt.get_dim();i++){
        minpt.set(i,pt.get_data(i));
    }
    
    if(target_asserted==0){
        strad->set_target(cc+delta_chisquared);
    }
    
    global_mindex=dex;

}

int gpWrapper::get_ngood(){
    return good_pts.get_dim();
}

double gpWrapper::get_good_max(int dex){
    return good_max.get_data(dex);
}

double gpWrapper::get_good_min(int dex){
    return good_min.get_data(dex);
}

void gpWrapper::set_good_max(int dex, double nn){
    good_max.set(dex,nn);
}

void gpWrapper::set_good_min(int dex, double nn){
    good_min.set(dex,nn);
}

double gpWrapper::get_chimin(){
    return chimin;
}

double gpWrapper::get_minpt(int dex){
    return minpt.get_data(dex);
}

array_1d<double>* gpWrapper::get_minpt(){
    return &minpt;
}

double gpWrapper::get_delta_chisquared(){
    return delta_chisquared;
}

void gpWrapper::set_delta_chisquared(double nn){
    delta_chisquared=nn;
}

void gpWrapper::assert_target(){
    target_asserted=1;
}

double gpWrapper::call_chisq(array_1d<double> &vv){
    if(chisq==NULL){
        printf("WARNING chisq is NULL but you just called call_chisq\n");
        exit(1);
    }
    
    return (*chisq)(vv);
}

int gpWrapper::get_chisq_dim(){
    if(chisq==NULL){
        printf("WARNING chisq is NULL but you just called get_chisq_dim()\n");
        exit(1);
    }
    
    return chisq->get_dim();
}

int gpWrapper::get_chisq_called(){
    if(chisq==NULL){
        printf("WARNING chisq is NULL but you just called get_chisq_called()\n");
        exit(1);
    }
    
    return chisq->get_called();
}

double gpWrapper::get_chisq_time(){
    if(chisq==NULL){
        printf("WARNING chisq is NULL but you just called get_chisq_time()\n");
    }
    
    return chisq->get_time_spent();
}

int gpWrapper::get_global_mindex(){
    return global_mindex;
}

void gpWrapper::evaluate_ngood(){
    
    good_pts.reset();    
    int i,j;
    
    for(i=0;i<gg->get_pts();i++){
        if(gg->get_fn(i)<strad->get_target()){
            good_pts.add(i);
            if(good_pts.get_dim()==0){
                for(j=0;j<gg->get_dim();j++){
                    good_min.set(j,gg->get_pt(i,j));
                    good_max.set(j,gg->get_pt(i,j));
                }
            }
            else{
                for(j=0;j<gg->get_dim();j++){
                    if(gg->get_pt(i,j)<good_min.get_data(j)){
                        good_min.set(j,gg->get_pt(i,j));
                    }
                    if(gg->get_pt(i,j)>good_max.get_data(j)){
                        good_max.set(j,gg->get_pt(i,j));
                    }
                }
            }
            
        }
    }
}

int gpWrapper::get_good_pt(int i){
    return good_pts.get_data(i);
}

int gpWrapper::get_dim(){
    if(gg==NULL){
        printf("WARNING cannot call gpWrapper::get_dim(); gg is NULL\n");
    }
    return gg->get_dim();
}

int gpWrapper::is_gp_null(){
    if(gg==NULL){
        return 1;
    }
    
    return 0;
}

double gpWrapper::get_target(){
    return strad->get_target();
}
