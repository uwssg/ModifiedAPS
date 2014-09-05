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
    
    if(chiout[0]<gg->get_target() && doAssociate==1){
        associates.add(dexout[0]);
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
