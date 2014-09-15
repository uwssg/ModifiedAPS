
#ifndef GP_WRAPPER_H
#define GP_WRAPPER_H

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "goto_tools.h"
#include "containers.h"
#include "eigen_wrapper.h"
#include "kd.h"
#include "gaussian_process.h"
#include "chisq.h"

////////////code for diagnostic testing
enum{iAPS,iSimplex,iCoulomb,iCompass,iBisect,iNodeBisect,iRicochet};
///////////////


class straddle_parameter{
    /*
    This class is meant to store both the target value of chisquared_lim and
    the calculation of the S statistic from equation (1) of the paper.
    
    If the user wanted to try a different combination of sigma, mu, and
    chisquared_lim than equation(1), she should alter the code in the 
    operator () of this class.
    */

public:
    ~straddle_parameter();
    straddle_parameter();
    
    /*set the value of chisquared_lim*/
    void set_target(double);
    
    /*return the value of chisquared_lim*/
    double get_target();
    
    /*accept mu and sigma and return S (equation 1 of the paper)*/
    double operator()(double,double) const;


private:
    double target;
};

class gpWrapper{
    /*
    This class will contain both the chisquared function and the
    gaussian_process for ModifiedAPS.  That way, the node objects can
    see the evaluate function
    */
    
    public:
        
        gpWrapper();
        ~gpWrapper();
        
        /*
        The functions evaluate() below are how APS actually calls the chisquared function.
    
        There are two risks involved in just calling the operator() to the provided chisquared
        function:
     
        1) There is no obvious mechanism in APS preventing APS from calling points outside of the
        allowed bounds of parameter space
    
        2) There is no mechanism preventing APS from calling points that are infinitesimally
        close to points that are already sampled and thus choking the Gaussian Process with neighbor
        points that are essentially identical.
    
        evaluate() fixes these problems.  Whenever a point is proposed for chisquared evaluation,
        evaluate first tests that it is in the bounds allowed by the chisquared function.  If not,
        evaluate returns chisquared = 2 x 10^30 and does not add the point to the Gaussian Process
        (if points with absurdly high values of chisquared are kept in the Gaussian Process, our
        attempts to predict chisquared using the Gaussian Process will become invalid near the boundaries
        of parameter space).
    
        If the point passes that test, evaluate() searches for the nearest neighbor to the proposed point.
        If that neighbor is closer than a (normalized) distance of 1.0 x 10^-8 in parameter space,
        evaluate returns the chisquared value of the nearest neighbor and, again, does not add the
        new point to the Gaussian Process.
    
        The arguments to evaluate are
    
        array_1d<double> -- the point to be evaluated
    
        double* -- a pointer to store the chisquared value of the point
    
        int* -- a pointer to store the index of the point, assuming it is added to the Gaussina Process
        (this is set to -1 if the point is not added to the Gaussian Process)
    
        There are some cases in which APS may evaluate the validity of a point before calling evaluate.
        In that case, one can pass the value 1 as a final int and evaluate() will dispense with
        the validity tests described above.
    
        */
        void evaluate(array_1d<double>&,double*,int*,int);
        void evaluate(array_1d<double>&,double*,int*);
        void evaluate(array_1d<double>&,double*);
        
        void set_gp(gp*);
        void set_chisq(chisquared*);
        void set_strad(straddle_parameter*);
        
        /*
        determine whether or not the specified point is within the bounds allowed
        by the chisquard funciton (return 1 if so; return 0 if not)
        */
        int in_bounds(array_1d<double>&);
    
        /*
        Do all of the validity checks required by evaluate(), i.e.
      
        1) is the point inside of the bounds allowed by the chisquared function
    
        2) is the point farther than a normalized parameter space distance of 10^-8
        from its nearest neighbor
    
        If so return 1.  If not, return 0.
    
        The optional double* will store the chisquared value associated with the nearest
        neighbor, if the nearest neighbor is too close.  If the first test failed, this
        pointer will store an absurdly high value (2 x 10^30) 
        */
        int is_valid(array_1d<double>&);
        int is_valid(array_1d<double>&, double*);
    
        /*
        The function add_pt() actually adds a point to the Gaussian Process.
    
        The arguments are
    
        array_1d<double> -- the point to be added
        double -- the chisquared value associated with that point
    
        the function will return the index of the point, once it is added to the
        Gaussian Process
    
        add_pt() will also assess whether or not the point improves upon chisquared_min
        or is a "good" point (i.e. whether chisquared<=chisquared_im)
        */
        int add_pt(array_1d<double>&,double);
        
            
        /*
        set_chimin() will set the minimum value of chisquared, the point at which that minimum occurred
        and the index by which that point is stored in the Gaussian Process
    
        If chisquared_lim is set relative to chisquared_min, this method will also update
        the value of chisquared_lim
        */
        void set_chimin(double,array_1d<double>&,int);
        
        double get_good_max(int);
        double get_good_min(int);
        int get_ngood();
        void set_good_max(int,double);
        void set_good_min(int,double);
        
        double get_chimin();
        double get_minpt(int);
        array_1d<double>* get_minpt();
        
        void set_delta_chisquared(double);
        double get_delta_chisquared();
        
        void assert_target();
        
        double call_chisq(array_1d<double>&);
        int get_chisq_dim();
        double get_chisq_time();
        
        int get_global_mindex();
        
        void evaluate_ngood();
        
        int get_good_pt(int);
        int get_dim();
        int is_gp_null();
        double get_target();
        
        void add_to_unitSpheres(array_1d<double>&);
        void unitSpheres_nn_srch(array_1d<double>&,int,array_1d<int>&,array_1d<double>&);
        int is_unitSpheres_null();
        int get_unitSpheres_pts();
        double get_pt(int,int);
        double get_fn(int);
        double get_min(int);
        double get_max(int);
        
        double distance(int,int);
        double distance(array_1d<double>&,array_1d<double>&);
        double distance(int,array_1d<double>&);
        double distance(array_1d<double>&,int);
        
        void reset_cache();
        double user_predict(array_1d<double>&,double*)const;
        double user_predict(array_1d<double>&)const;
        
        void actual_gradient(int,array_1d<double>&);
        void actual_gradient(array_1d<double>&,array_1d<double>&);
        
        int get_called();
        int get_pts();
        
        int get_search_ct();
        double get_search_time();
        
        int get_search_ct_solo();
        double get_search_time_solo();
        
        int get_search_ct_box();
        double get_search_time_box();
        int get_smallest_box();
        int get_biggest_box();
        
        /////////code for diagnostic testing
        int set_iWhere(int);
        int get_iWhere();
        int get_whereCt(int);
        int get_whereFrom(int);
        void set_whereFrom(int,int);
        //////////
        
    private:
        gp *gg;
        chisquared *chisq;
        straddle_parameter *strad;
        kd_tree *unitSpheres;
        
        array_2d<double> sphereSeedData;
        array_1d<double> good_max,good_min,minpt;
        array_1d<int> good_pts;
        int global_mindex,target_asserted;
        
        double chimin,delta_chisquared;
        
        ///////code for diagnostic testing
        int iWhere;
        array_1d<int> whereCt,whereFrom;
        
        ///////

};


#endif
