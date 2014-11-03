#include "simplex.h"

simplex_minimizer::simplex_minimizer(){
    _alpha=1.0;
    _beta=0.5;
    _gamma=2.1;
    
    initialize();
}

simplex_minimizer::simplex_minimizer(double aa, double bb, double gg){
    _alpha=aa;
    _beta=bb;
    _gamma=gg;
    
    initialize();
}

simplex_minimizer::~simplex_minimizer(){}

void simplex_minimizer::initialize(){
    cost=NULL;
    chisquared=NULL;
    
    _called_cost=0;
    _freeze_temp=0;
    _temp=-0.1;
    _use_gradient=0;
    
    _last_found=0;
    _called_evaluate=0;
    _abort_max_factor=10;
    _min_ff=2.0*chisq_exception;
    
    _transform.set_name("simplex_transform");
    _origin.set_name("simplex_origin");
    _min_pt.set_name("simplex_min_pt");
    _ff.set_name("simplex_ff");
    _pstar.set_name("simplex_pstar");
    _pstarstar.set_name("simplex_pstarstar");
    _pts.set_name("simplex_pts");
}

void simplex_minimizer::set_chisquared(function_wrapper *ff){
    chisquared=ff;
}

void simplex_minimizer::set_cost(function_wrapper *cc){
    cost=cc;
}

void simplex_minimizer::use_gradient(){
    _use_gradient=1;
}

void simplex_minimizer::set_minmax(array_1d<double> &min, array_1d<double> &max){
    if(min.get_dim()!=max.get_dim()){
        printf("WARNING simplex_minimizer::set_minmax min %d max %d\n",
        min.get_dim(),max.get_dim());
        
        exit(1);
    }
    
    int i;
    for(i=0;i<min.get_dim();i++){
        _origin.set(i,min.get_data(i));
        _transform.set(i,max.get_data(i)-min.get_data(i));
    }
    
}

void simplex_minimizer::find_il(){
    if(_ff.get_dim()==0 || 
       _pts.get_rows()==0 || 
       _pts.get_cols()==0 ||
       _pts.get_rows()!=_ff.get_dim()){
   
       printf("WARNING cannot call find_il %d %d %d\n",
       _ff.get_dim(),_pts.get_rows(),_pts.get_cols());
       
       exit(1);
   }
   
   int i;
   for(i=0;i<_ff.get_dim();i++){
       if(i==0 || _ff.get_data(i)<_ff.get_data(_il))_il=i;
       if(i==0 || _ff.get_data(i)>_ff.get_data(_ih))_ih=i;
   }
}

double simplex_minimizer::evaluate(array_1d<double> &pt){
    /*
    pt will be in the transformed dimensions of the simplex
    */
    
    if(chisquared==NULL){
        printf("WARNING cannot call simplex_minimizer::evaluate\n");
        printf("chisquared is NULL\n");
        exit(1);
    }
    
    _called_evaluate++;
    
    int i;
    array_1d<double> vv;
    vv.set_name("simplex_minimizer_vv");
    for(i=0;i<pt.get_dim();i++){
        vv.set(i,_origin.get_data(i)+pt.get_data(i)*_transform.get_data(i));
    }
    
    double fval;
    chisquared->evaluate(vv,&fval);
    
    if(fval<_min_ff){
        _last_found=_called_evaluate;
        _min_ff=fval;
        for(i=0;i<pt.get_dim();i++){
            _min_pt.set(i,vv.get_data(i));
        }
    }
    
    double cval;
    if(cost!=NULL){
        cval=evaluate_cost(vv);
        fval+=cval;
    }
    
    return fval;
}

double simplex_minimizer::evaluate_cost(array_1d<double> &vv){
    /*
    vv will already be transformed into natural coordinates
    */
    
    if(cost==NULL){
        printf("WARNING cannot call simplex_minimizer::evaluate_cost\n");
        printf("cost is NULL\n");
        exit(1);
    }

    if(_temp<-20.0)return 0.0;

    double cval;
    int temp_changed=0;
    cost->evaluate(vv,&cval);

    if(_freeze_temp==0)_called_cost++;

    if(_called_cost%(vv.get_dim()*10)==0 && _freeze_temp==0){
        _temp*=2.0;
        temp_changed=1;
    }
    
    double mu;
    int i;
    if(temp_changed==1){
        _freeze_temp=1;
        if(_pstar.get_dim()>0){
            _fstar=evaluate(_pstar);
        }
        
        if(_pstarstar.get_dim()>0){
            _fstarstar=evaluate(_pstarstar);
        }
        
        for(i=0;i<_pts.get_rows();i++){
            mu=evaluate(_pts(i)[0]);
            _ff.set(i,mu);
        }
        
        find_il();
        
        _freeze_temp=0;
    }
    
    return exp(_temp)*cval;
}

void simplex_minimizer::find_minimum(array_2d<double> &seed, array_1d<double> &min_pt){

    if(seed.get_rows()!=seed.get_cols()+1){
        printf("WARNING you gave simplex minimizer %d points of %d dim\n",
        seed.get_rows(),seed.get_cols());
        printf("Need dim+1 points\n");
        exit(1);
    }
    
    int i;
    if(_origin.get_dim()==0){
        for(i=0;i<seed.get_cols();i++){
            _origin.set(i,0.0);
            _transform.set(i,1.0);
        }
    }
    
    int j;
    _pts.set_cols(seed.get_cols());
    for(i=0;i<seed.get_rows();i++){
        for(j=0;j<seed.get_cols();j++){
            _pts.set(i,j,(seed.get_data(i,j)-_origin.get_data(j))/_transform.get_data(j));
        }
    }
    
    double mu;
    for(i=0;i<seed.get_rows();i++){
        mu=evaluate(_pts(i)[0]);
        _ff.set(i,mu);
    }
    
    find_il();
    
    int abort_max=_abort_max_factor*seed.get_cols();
    int dim=seed.get_cols();
    double spread;
    
    array_1d<double> pbar;
    pbar.set_name("simplex_pbar");
    
    while(_called_evaluate-_last_found<abort_max){
       for(i=0;i<dim;i++){
           pbar.set(i,0.0);
           for(j=0;j<dim+1;j++){
               if(j!=_ih){
                   pbar.add_val(i,_pts.get_data(j,i));
               }
           }
           pbar.divide_val(i,double(dim));
       }
       
       for(i=0;i<dim;i++){
           _pstar.set(i,(1.0+_alpha)*pbar.get_data(i)-_alpha*_pts.get_data(_ih,i));
       }
       _fstar=evaluate(_pstar);
       
       if(_fstar<_ff.get_data(_ih) && _fstar>_ff.get_data(_il)){
           _ff.set(_ih,_fstar);
           for(i=0;i<dim;i++){
               _pts.set(_ih,i,_pstar.get_data(i));
           }
       }
       else if(_fstar<_ff.get_data(_il)){
           for(i=0;i<dim;i++){
               _pstarstar.set(i,_gamma*_pstar.get_data(i)+(1.0-_gamma)*pbar.get_data(i));
           }
           _fstarstar=evaluate(_pstarstar);
           
           if(_fstarstar<_ff.get_data(_il)){
               for(i=0;i<dim;i++){
                   _pts.set(_ih,i,_pstarstar.get_data(i));
               }
               _ff.set(_ih,_fstarstar);
           }
           else{
               for(i=0;i<dim;i++){
                   _pts.set(_ih,i,_pstar.get_data(i));
               }
               _ff.set(_ih,_fstar);
           }
       }
       
       find_il();
       
       j=1;
       for(i=0;i<dim+1;i++){
           if(_fstar<_ff.get_data(i) && i!=_ih){
               j=0;
           }
       }
       
       if(j==1){
           for(i=0;i<dim;i++){
               _pstarstar.set(i,_beta*_pts.get_data(_ih,i)+(1.0-_beta)*pbar.get_data(i));
           }
           _fstarstar=evaluate(_pstarstar);
           
           if(_fstarstar<_ff.get_data(_ih)){
               for(i=0;i<dim;i++)_pts.set(_ih,i,_pstarstar.get_data(i));
               _ff.set(_ih,_fstarstar);
           }
           else{
               find_il();
               for(i=0;i<dim+1;i++){
                   if(i!=_il){
                       for(j=0;j<dim;j++){
                           mu=0.5*(_pts.get_data(i,j)+_pts.get_data(_il,j));
                       }
                       mu=evaluate(_pts(i)[0]);
                       _ff.set(i,mu);
                   }
               }
           }
       }
       
       find_il();
       spread=_ff.get_data(_ih)-_ff.get_data(_il);
       if(spread<1.0 && _use_gradient==1){
           gradient_minimizer();
       }
    }
    
    for(i=0;i<dim;i++){
        min_pt.set(i,_min_pt.get_data(i));
    }
}
