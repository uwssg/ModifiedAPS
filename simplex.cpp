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
    dice=NULL;
    
    _min_temp=-8.0;
    
    _use_gradient=0;
    
    _abort_max_factor=10;
    
    _freeze_temp=-1;
    
    _transform.set_name("simplex_transform");
    _origin.set_name("simplex_origin");
    _min_pt.set_name("simplex_min_pt");
    _ff.set_name("simplex_ff");
    _pstar.set_name("simplex_pstar");
    _pstarstar.set_name("simplex_pstarstar");
    _pts.set_name("simplex_pts");
    _last_improved_ff.set_name("simplex_last_improved_ff");
    _last_improved_pts.set_name("simplex_last_improved_pts");
}

void simplex_minimizer::set_chisquared(function_wrapper *ff){
    chisquared=ff;
}

void simplex_minimizer::set_cost(function_wrapper *cc){
    cost=cc;
}

void simplex_minimizer::set_dice(Ran *dd){
    dice=dd;
}

void simplex_minimizer::use_gradient(){
    _use_gradient=1;
}

void simplex_minimizer::freeze_temp(){
    _freeze_temp=1;
}

void simplex_minimizer::unfreeze_temp(){
    _freeze_temp=0;
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
    
    if(_freeze_called==0)_called_evaluate++;
    
    int i,j;
    array_1d<double> vv;
    vv.set_name("simplex_minimizer_vv");
    for(i=0;i<pt.get_dim();i++){
        vv.set(i,_origin.get_data(i)+pt.get_data(i)*_transform.get_data(i));
    }
    
    double fval,raw;
    chisquared->evaluate(vv,&fval);
    
    if(fval<_true_min_ff){
        _true_min_ff=fval;
        for(i=0;i<pt.get_dim();i++){
            _min_pt.set(i,vv.get_data(i));
        }
    }
    
    raw=fval;
    
    double cval;
    if(cost!=NULL){
        cval=evaluate_cost(vv);
        fval+=cval;
    }
    
   /* printf("    %e %e %d %d -- %d\n",
    fval,raw,_called_evaluate,chisquared->get_called(),
    chisquared->get_called()-_called_evaluate);*/
    
    if(fval<_min_ff){
        _last_found=_called_evaluate;
        _min_ff=fval;
        
        printf("min %e true %e cost %e raw %e\n",
        _min_ff,_true_min_ff,cval,raw);

        if(_ff.get_dim()==pt.get_dim()+1){
            for(i=0;i<_pts.get_rows();i++){
                _last_improved_ff.set(i,_ff.get_data(i));
                for(j=0;j<_pts.get_cols();j++){
                    _last_improved_pts.set(i,j,_pts.get_data(i,j));
                }    
            }
        }
    }
    
    return fval;
}

void simplex_minimizer::cool_off(){
    if(_freeze_temp==1) return;
    _temp-=1.0;
    
    double mu;
    int i;

        
    _freeze_called=1;
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
    _min_ff=_ff.get_data(_il);
        
    _freeze_temp=0;
    _freeze_called=0;
   
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

    if(_temp<_min_temp)return 0.0;

    double cval;
    cost->evaluate(vv,&cval);

    if(_freeze_temp==0)_called_cost++;

    if(_called_cost%(10*vv.get_dim())==0){
        cool_off();
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
    
    
    if(cost!=NULL){
        _temp=0.0;
    }
    else{
        _temp=-1000.0;
    }
    
    if(_freeze_temp<0)_freeze_temp=0;
    
    _freeze_called=0;
    _called_cost=0;
    _last_found=0;
    _called_evaluate=0;
    _last_called_gradient=0;
    
    
    _min_ff=2.0*chisq_exception;
    _true_min_ff=2.0*chisq_exception;
    _min_pt.reset();
    _pstar.reset();
    _pstarstar.reset();
    _last_improved_ff.reset();
    _pts.reset();
    _last_improved_pts.reset();
    
    int i;
    if(_origin.get_dim()==0){
        for(i=0;i<seed.get_cols();i++){
            _origin.set(i,0.0);
            _transform.set(i,1.0);
        }
    }
    
    int j;
    _pts.set_cols(seed.get_cols());
    _last_improved_pts.set_cols(seed.get_cols());
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
    int need_to_thaw,dim=seed.get_cols();
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
                           _pts.set(i,j,mu);
                       }
                       mu=evaluate(_pts(i)[0]);
                       _ff.set(i,mu);
                   }
               }
           }
       }
       
       find_il();
       spread=_ff.get_data(_ih)-_ff.get_data(_il);
       if(spread<0.1*_min_ff && _use_gradient==1 && _called_evaluate>abort_max/2+_last_called_gradient){
           //printf("implementing gradient %e\n",_temp);
           
           if(_freeze_temp==0)need_to_thaw=1;
           else need_to_thaw=0;
           
           _freeze_temp=1;
           
           gradient_minimizer();
           
           if(need_to_thaw==1)_freeze_temp=0;
       }
       
       if(_called_evaluate-_last_found>=abort_max && _temp>_min_temp){
           if(_freeze_temp==0)need_to_thaw=1;
           else need_to_thaw=0;
           
           _freeze_called=1;
           _freeze_temp=1;
           for(i=0;i<dim+1;i++){
               if(i!=_il){
                   for(j=0;j<dim;j++){
                       _pts.set(i,j,_pts.get_data(_il,j)+dice->doub()*0.5);
                   }
                   mu=evaluate(_pts(i)[0]);
                   _ff.set(i,mu);
               }
           }
           _freeze_called=0;
           
           if(need_to_thaw==1){
               _freeze_temp=0;
               _last_found=_called_evaluate;
           }
           find_il();
       }
       //printf("spread %e %e %e\n\n",spread,_temp,_min_ff);
    }
    
    for(i=0;i<dim;i++){
        min_pt.set(i,_min_pt.get_data(i));
    }
    printf("    leaving simplex %d %d %d\n",_called_evaluate,_last_found,abort_max);
    printf("    temp %e\n",_temp);
    
    _freeze_temp=-1;
}

void simplex_minimizer::gradient_minimizer(){
    if(dice==NULL){
        printf("WARNING cannot use simplex gradient; dice is NULL\n");
        exit(1);
    }
    find_il();

    array_1d<double> gradient,trial;
    gradient.set_name("simplex_gradient");
    trial.set_name("simplex_gradient_trial");
    
    int ix,dim,i,j,k;
    double x1,x2,mu1,mu2,mu,dx;
    
    dim=_pts.get_cols();
    
    for(ix=0;ix<dim;ix++){
        for(i=0;i<dim;i++){
            trial.set(i,_pts.get_data(_il,i));
        }
        
        dx=0.1;
        k=0;
        mu1=2.0*chisq_exception;
        while(!(mu1<chisq_exception) && k<5){
            k++;
            x1=_pts.get_data(_il,ix)-dx;
            trial.set(ix,x1);
            mu1=evaluate(trial);
            
            if(!(mu1<chisq_exception)){
                dx*=0.5;
            }
        }
        
        dx=0.1;
        k=0;
        mu2=2.0*chisq_exception;
        while(!(mu2<chisq_exception) && k<5){
            k++;
            x2=_pts.get_data(_il,ix)-dx;
            trial.set(ix,x2);
            mu2=evaluate(trial);
            
            if(!(mu2<chisq_exception)){
                dx*=0.5;
            }
        }
        
        gradient.set(ix,(mu1-mu2)/(x1-x2));
    }
    
    mu2=gradient.normalize();
    array_1d<double> step;
    step.set_name("simplex_gradient_step");
    
    if(_last_improved_ff.get_dim()>0){
        for(i=0;i<dim+1;i++){
            if(i==0 || _last_improved_ff.get_data(i)<_last_improved_ff.get_data(j))j=i;
        }
        
        for(i=0;i<dim;i++){
            step.set(i,_pts.get_data(_il,i)-_last_improved_pts.get_data(j,i));
        }
        mu=step.normalize();
        
        if(!isnan(mu2)){
            for(i=0;i<dim;i++){
                mu1=0.5*(step.get_data(i)+gradient.get_data(i));
                step.set(i,mu1);
            }
        }
    }
    else{
        mu=1.0;//should this be something smaller? in the other case, it is the size of
               //the step from old _il to new _il
        for(i=0;i<dim;i++)step.set(i,gradient.get_data(i));
    }
    
    step.normalize();
    double theta;
    array_1d<double> deviation;
    deviation.set_name("simplex_gradient_deviation");
    
    for(i=0;i<dim+1;i++){
        for(j=0;j<dim;j++){
            _pts.add_val(i,j,mu*step.get_data(j));
        }
        
        if(i!=_il){
            theta=0.0;
            mu1=-1.0;
            while(mu1<0.0 || isnan(mu1)){
                theta=0.0;
                for(j=0;j<dim;j++){
                    deviation.set(j,normal_deviate(dice,0.0,1.0));
                    theta+=deviation.get_data(j)*step.get_data(j);
                }
                
                if(!(_temp>_min_temp)){
                    for(j=0;j<dim;j++){
                        deviation.subtract_val(j,theta*step.get_data(j));
                    }
                }
                mu1=deviation.normalize();
            }
            
            for(j=0;j<dim;j++){
                _pts.add_val(i,j,0.1*deviation.get_data(j));
            }
        }
    }
    
    for(i=0;i<dim+1;i++){
        mu=evaluate(_pts(i)[0]);
        _ff.set(i,mu);
    }
    find_il();
    _last_called_gradient=_called_evaluate;
    
}
