#ifndef NODE_H
#define NODE_H

#include "gp_wrapper.h"

class node{

    public:
        
        node();
        node(const node&);
        node& operator=(const node&);
        ~node();
        
        void evaluate(array_1d<double>&,double*,int*);
        void evaluateNoAssociate(array_1d<double>&,double*,int*);
        void evaluate(array_1d<double>&,double*,int*,int);
        
        void project_to_unit_sphere(array_1d<double>&, array_1d<double>&);
        void add_as_boundary(int);
        
    private:
        
        void copy(const node&);
        
        array_1d<int> associates,boundaryPoints;
        array_2d<double> basisVectors;
        
        int center_dex,last_set_bases;
        double time_ricochet,time_coulomb,time_total;
        
        gpWrapper *gg;
        Ran *dice;

};

#endif
