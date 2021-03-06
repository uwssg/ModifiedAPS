#include "aps.h"

function_wrapper::function_wrapper(){}

function_wrapper::~function_wrapper(){}

int function_wrapper::get_called(){
    return -1;
}

void function_wrapper::evaluate(array_1d<double> &pp, double *cc){
    printf("WARNING calling default function_wrapper::evaluate; should not do that\n");
    exit(1);
}

double function_wrapper::diagnostic_evaluate(array_1d<double> &pp){
    printf("WARNING default diagnostic_evaluate not implemented\n");
    exit(1);
}

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
    unitSpheres=NULL;
    
    sphereSeedData.set_name("sphereSeedData");
    good_max.set_name("good_max");
    good_min.set_name("good_min");
    minpt.set_name("minpt");
    good_pts.set_name("good_pts");
    
    whereCt.set_name("gp_wrapper_whereCt");
    whereFrom.set_name("gp_wrapper_whereFrom");
    
    whereCt.set_name("gp_wrapper_whereCt");
    whereCt.set(iAPS,0);
    whereCt.set(iSimplex,0);
    whereCt.set(iCoulomb,0);
    whereCt.set(iCompass,0);
    whereCt.set(iBisect,0);
    whereCt.set(iNodeBisect,0);
    whereCt.set(iRicochet,0);
    iWhere=iAPS;
    
}

gpWrapper::~gpWrapper(){
    
    if(unitSpheres!=NULL){
        delete unitSpheres;
    }

}

chisquared* gpWrapper::get_chifn(){
    return chisq;
}

void gpWrapper::set_gp(gp *gg_in){
    gg=gg_in;
    good_max.set_dim(gg->get_dim());
    good_min.set_dim(gg->get_dim());
    minpt.set_dim(gg->get_dim());
}

int gpWrapper::get_search_ct_box(){
    if(gg==NULL) return 0;
    return gg->get_search_ct_box();
}

double gpWrapper::get_search_time_box(){
    if(gg==NULL) return 0;
    return gg->get_search_time_box();
}

int gpWrapper::get_biggest_bad_box(){
    if(gg==NULL) return 0;
    return gg->get_biggest_bad_box(get_target());
}

int gpWrapper::get_biggest_box(){
    if(gg==NULL) return 0;
    return gg->get_biggest_box();
}

int gpWrapper::get_smallest_box(){
    if(gg==NULL) return 0;
    return gg->get_smallest_box();
}

int gpWrapper::get_n_small_boxes(){
    if(gg==NULL) return 0;
    return gg->get_n_small_boxes();
}

int gpWrapper::get_n_optimal_boxes(){
    if(gg==NULL) return 0;
    return gg->get_n_optimal_boxes();
}

int gpWrapper::get_nboxes(){
    if(gg==NULL) return 0;
    return gg->get_nboxes();
}

int gpWrapper::get_box_contents(int dex){
    if(gg==NULL) return 0;
    return gg->get_box_contents(dex);
}

int gpWrapper::get_box_contents(int dex, int ii){
    if(gg==NULL){
        printf("WARNING asked for specific box contents, but gg is NULL\n");
        exit(1);
    }
    
    return gg->get_box_contents(dex,ii);
}

double gpWrapper::get_box_max(int dex, int idim){
    if(gg==NULL){
        printf("WARNING asked for box_max but gg is NULL\n");
        exit(1);
    }
    return gg->get_box_max(dex,idim);
}

double gpWrapper::get_box_min(int dex, int idim){
    if(gg==NULL){
        printf("WARNING asked for box_min but gg is NULL\n");
        exit(1);
    }
    return gg->get_box_min(dex,idim);
}

int gpWrapper::get_search_ct(){
    if(gg==NULL){
        return 0;
    }
    
    return gg->get_search_ct();
}

double gpWrapper::get_search_time(){
    if(gg==NULL){
        return 0.0;
    }
    
    return gg->get_search_time();
}


int gpWrapper::get_search_ct_solo(){
    if(gg==NULL){
        return 0;
    }
    
    return gg->get_search_ct_solo();
}

double gpWrapper::get_search_time_solo(){
    if(gg==NULL){
        return 0.0;
    }
    
    return gg->get_search_time_solo();
}

int gpWrapper::set_iWhere(int i){
    iWhere=i;
}

int gpWrapper::get_iWhere(){
    return iWhere;
}

int gpWrapper::get_whereCt(int i){
    return whereCt.get_data(i);
}

void gpWrapper::set_whereFrom(int dex, int i){
    whereFrom.set(dex,i);
}

int gpWrapper::get_whereFrom(int i){
    if(i>=whereFrom.get_dim()){
        return 0;
    }
    return whereFrom.get_data(i);
}

void gpWrapper::set_chisq(chisquared *cc){
    chisq=cc;
}

void gpWrapper::set_strad(straddle_parameter *strad_in){
    strad=strad_in;
}

double gpWrapper::straddle_value(double mu, double sig){
    return (*strad)(mu,sig);
}

double gpWrapper::get_pt(int dex, int i){
    return gg->get_pt(dex,i);
}

double gpWrapper::get_fn(int dex){
    return gg->get_fn(dex);
}

double gpWrapper::get_min(int dex){
    return gg->get_min(dex);
}

double gpWrapper::get_max(int dex){
    return gg->get_max(dex);
}

double gpWrapper::get_length(int dex){
    return gg->get_length(dex);
}

double gpWrapper::distance(int i, int j){
    if(i<0 || j<0 || i>=gg->get_pts() || j>=gg->get_pts()){
        printf("    iWhere %d\n",iWhere);
        printf("    int int\n");
    }
    return gg->distance(i,j);
}

double gpWrapper::distance(array_1d<double> &x, array_1d<double> &y){
    return gg->distance(x,y);
}

double gpWrapper::distance(array_1d<double> &x, int i){
    if(i<0 || i>=gg->get_pts()){
        printf("    iWhere %d\n",iWhere);
        printf("    pt, int\n");
    }

    return gg->distance(x,i);
}

double gpWrapper::distance(int i, array_1d<double> &x){
    if(i<0 || i>=gg->get_pts()){
        printf("    iWhere %d\n",iWhere);
        printf("    int, pt\n");
    }
    return gg->distance(i,x);
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
            
            whereFrom.set(dex[0],iWhere);
            whereCt.add_val(iWhere,1);
            
        }
    }
}

double gpWrapper::diagnostic_evaluate(array_1d<double> &pt){
    double ans=(*chisq)(pt);
    chisq->decrement_called();
    return ans;
}

void gpWrapper::freeze_boxes(){
    if(gg==NULL){
        printf("WARNING cannot freeze boxes; gg is null\n");
        exit(1);
    }
    
    gg->freeze_boxes();
}

void gpWrapper::unfreeze_boxes(){
    if(gg==NULL){
        printf("WARNING cannot unfreeze boxes; gg is null\n");
        exit(1);
    }
    
    gg->unfreeze_boxes();
}

int gpWrapper::get_ct_search(){
    return gg->get_ct_search();
}

void gpWrapper::get_box_quartiles(array_1d<int> &qq){
    if(gg==NULL){
        printf("WARNING cannot get box quartiles; gg is null\n");
        exit(1);
    }
    
    gg->get_box_quartiles(qq);
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

void gpWrapper::add_to_unitSpheres(array_1d<double> &pt){
    
    int i,j;
    array_1d<double> _min,_max;
    
    if(unitSpheres==NULL){
        if(sphereSeedData.get_cols()==0){
            sphereSeedData.set_cols(gg->get_dim());
        }
        
        sphereSeedData.add_row(pt);
        
        if(sphereSeedData.get_rows()==gg->get_dim()){
            for(i=0;i<gg->get_dim();i++){
                _min.set(i,0.0);
                _max.set(i,gg->get_max(i)-gg->get_min(i));
            }
            
            unitSpheres=new kd_tree(sphereSeedData,_min,_max);
            sphereSeedData.reset();
        }
    }
    else{
        unitSpheres->add(pt);
    }
    
}

int gpWrapper::is_unitSpheres_null(){
    if(unitSpheres==NULL){
        return 1;
    }
    
    return 0;
}

void gpWrapper::unitSpheres_nn_srch(array_1d<double> &pt, int kk, array_1d<int> &neigh,
    array_1d<double> &dd){

    
    if(unitSpheres==NULL){
        printf("WARNING cannot call unitSpheres_nn_srch; unitSpheres is null\n");
        exit(1);
    }
    
    unitSpheres->nn_srch(pt,kk,neigh,dd);
    
}

int gpWrapper::get_unitSpheres_pts(){
    if(unitSpheres==NULL){
        return 0;
    }
    
    return unitSpheres->get_pts();
}

void gpWrapper::reset_cache(){
    gg->reset_cache();
}

double gpWrapper::user_predict(array_1d<double> &x, double *s) const{
    int flag;
    double mu;
    mu=gg->user_predict(x,s,0,&flag);
    if(flag==-1){
        printf("Box failure %d\n",iWhere);
    }
    return mu;
}

double gpWrapper::user_predict(array_1d<double> &x) const{
    int flag;
    double mu;
    mu=gg->user_predict(x,0,&flag);
    if(flag==-1){
        printf("Box failure %d\n",iWhere);
    }
    return mu;
}

void gpWrapper::actual_gradient(array_1d<double> &x, array_1d<double> &y){
    gg->actual_gradient(x,y);
}

void gpWrapper::actual_gradient(int x, array_1d<double> &y){
    gg->actual_gradient(x,y);
}

int gpWrapper::get_called(){
    return chisq->get_called();
}

int gpWrapper::get_pts(){
    return gg->get_pts();
}

int gpWrapper::find_box(array_1d<double> &p){
    return gg->find_box(p);
}

int gpWrapper::find_box(int dex){
    return gg->find_box(dex);
}

void gpWrapper::nn_srch(array_1d<double> &pt, int kk, array_1d<int> &neigh,
array_1d<double> &dd){

    if(gg==NULL){
        printf("WARNING gpWrapper cannot do nn_srch; gg is NULL\n");
    }
    
    gg->nn_srch(pt,kk,neigh,dd);

}
