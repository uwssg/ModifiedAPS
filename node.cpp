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
    if(gg==NULL){
        printf("WARNING cannot call node::evaluate; gg is NULL\n");
    }
    
    gg->evaluate(pt,chiout,dexout);
    
    if(chiout[0]<gg->get_target()){
        associates.add(dexout[0]);
    }
    
    //spock what should I do about boundaryPoints and the unitSphere?
    //actually, wait, the unitSphere kd_tree is a global aps phenomenon
    //(it needs to keep all of the unitSpheres in one place)
    //
    //we should be fine
}
