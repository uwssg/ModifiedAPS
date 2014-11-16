#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "containers.h"
#include "goto_tools.h"
#include "gp_wrapper.h"

class simplex_minimizer{

public:

    simplex_minimizer();
    simplex_minimizer(double,double,double);
    ~simplex_minimizer();
    
    void set_chisquared(function_wrapper*);
    void set_cost(function_wrapper*);
    void set_dice(Ran*);
    void set_minmax(array_1d<double>&, array_1d<double>&);
    void use_gradient();
    void freeze_temp();
    void unfreeze_temp();
    
    /*
    the array_2d will be the input array of points;
    the array_1d will be the output minimum point
    */
    void find_minimum(array_2d<double>&, array_1d<double>&);

private:
    
    double evaluate(array_1d<double>&);
    double evaluate_cost(array_1d<double>&);
    void cool_off();
    
    void gradient_minimizer();
    
    void find_il();
    
    void initialize();
    
    double _temp,_min_ff,_true_min_ff,_fstar,_fstarstar,_min_temp;
    double _alpha,_beta,_gamma;
    int _il,_ih,_called_cost,_freeze_temp,_use_gradient;
    int _freeze_called,_last_called_gradient;
    int _last_found,_called_evaluate,_abort_max_factor;
    array_1d<double> _transform, _origin,_ff,_pstar,_pstarstar,_min_pt;
    array_1d<double> _last_improved_ff;
    array_1d<double> _min,_max;
    array_2d<double> _pts,_last_improved_pts;
    
    Ran *dice;
    function_wrapper *chisquared, *cost;
    /*
    cost will need to be a function_wrapper sub-class
    that has pointers to all of the aps nodes.
    */

};

#endif
