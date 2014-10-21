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
        
        void set_farthest_associate(double);
        
        Ran* get_Ran();
        gpWrapper* get_gpWrapper();
        int get_center();
        
        void copy(const node&);
        
        double get_farthest_associate();
        int get_n_associates();
        
        int get_n_boundary();
        int get_boundary(int);
        
        double get_time();
        double get_time_coulomb();
        double get_time_bases();
        double get_time_ricochet();
        int get_ct_search();
        int get_ct_ricochet();
        int get_ct_coulomb();
        int get_ct_bases();
        int get_calls_to_bases();
        
        double get_basis(int,int);
        
        double volume();
        int is_it_active();
        
        double get_max(int);
        double get_min(int);
        
        void flush_candidates(array_1d<int>&);
        
    private:
        
        array_1d<int> associates,boundaryPoints;
        array_2d<double> basisVectors;
        array_1d<double> range_max,range_min,geographicCenter;
        
        array_1d<int> candidates,centerCandidates;
        
        int ct_search,ct_coulomb,ct_bases,ct_ricochet;
        int calls_to_bases;
        int last_expanded,activity;
        int center_dex,last_nBasisAssociates,last_nAssociates;
        double time_ricochet,time_coulomb,time_search,time_bases;
        double farthest_associate,time_penalty;
        
        gpWrapper *gg;
        Ran *dice;
        
        void set_names();
        int bisection(int,int);
        int bisection(array_1d<double>&,double,array_1d<double>&,double,int);
        int bisection(array_1d<double>&,double,array_1d<double>&,double);
        int bisectionAssociate(int,int);
        int bisectionAssociate(array_1d<double>&,double,array_1d<double>&,double);
        int bisection(int,int,int);
        int coulomb_search();
        void ricochet_search(int,array_1d<double>&);
        int ricochet_driver(int,array_1d<double>&,array_1d<double>&);
        void compass_search(int);
        void find_bases();
        int perturb_bases(array_2d<double>&,int,array_1d<double>&,array_2d<double>&);
        double basis_error(array_2d<double>&,array_1d<int>&);
        void recenter();

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
