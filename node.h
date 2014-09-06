#ifndef NODE_H
#define NODE_H

#include "gp_wrapper.h"

class node{
    
    public:
        
        node();
        node(const node&);
        node& operator=(const node&);
        ~node();
        
        void set_center_dex(int);
        
        void set_gpWrapper(gpWrapper*);
        
        void set_dice(Ran*);
        
        void evaluate(array_1d<double>&,double*,int*);
        void evaluateNoAssociate(array_1d<double>&,double*,int*);
        void evaluate(array_1d<double>&,double*,int*,int);
        
        void project_to_unit_sphere(array_1d<double>&, array_1d<double>&);
        void add_as_boundary(int);
        
        int search();
        
        Ran* get_Ran();
        gpWrapper* get_gpWrapper();
        int get_center();
        
        void copy(const node&);
        
        double get_farthest_associate();
        int get_n_associates();
        
    private:
        
        array_1d<int> associates,boundaryPoints;
        array_2d<double> basisVectors;
        
        int center_dex,last_nBasisAssociates,last_nAssociates;
        double time_ricochet,time_coulomb,time_search,time_bases;
        double farthest_associate;
        
        gpWrapper *gg;
        Ran *dice;
        
        void set_names();
        int bisection(int,int);
        int bisectionAssociate(int,int);
        int bisection(int,int,int);
        int coulomb_search();
        void ricochet_search(int,array_1d<double>&);
        int ricochet_driver(int,array_1d<double>&,array_1d<double>&);
        void compass_search(int);
        void find_bases();
        int perturb_bases(array_2d<double>&,int,array_1d<double>&,array_2d<double>&);
        double basis_error(array_2d<double>&,array_1d<int>&);

};

class arrayOfNodes{

    public:
        arrayOfNodes();
        ~arrayOfNodes();
        
        void add(int,gpWrapper*,Ran*);
        void add(int,Ran*,gpWrapper*);
        void add(Ran*,int,gpWrapper*);
        void add(Ran*,gpWrapper*,int);
        void add(gpWrapper*,Ran*,int);
        void add(gpWrapper*,int,Ran*);
        
        int get_dim();
        
        void remove(int);
        
        node* operator()(int);
        
    private:
        node *data;
        int ct,room;


};


#endif
