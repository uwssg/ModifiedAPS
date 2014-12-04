#include "aps.h"

gp_cost_function::gp_cost_function(){
    printf("WARNING default constructor not implemented for gp_cost_function\n");
    exit(1);
}

gp_cost_function::gp_cost_function(kd_tree *g, array_1d<double> *ff,
array_1d<double> &min, array_1d<double> &max){
    kptr=g;
    fn=ff;
    
    _covarin.set_name("gp_cost_cached_covarin");
    _neigh.set_name("gp_cost_cached_neigh");
    _min.set_name("gp_cost_min");
    _max.set_name("gp_cost_max");
    
    int i;
    for(i=0;i<min.get_dim();i++){
        _min.set(i,min.get_data(i));
        _max.set(i,max.get_data(i));
    }
}

gp_cost_function::~gp_cost_function(){}

void gp_cost_function::evaluate(array_1d<double> &pt, double *out){
    array_1d<double> dd;
    array_1d<int> neigh;
    double ell,nugget;
    int i,j,npts;
    
    for(i=0;i<pt.get_dim();i++){
        if(pt.get_data(i)>_max.get_data(i)){
            out[0]=0.0;
            return;
        }
        
        if(pt.get_data(i)<_min.get_data(i)){
            out[0]=0.0;
            return;
        }
    }
    
    npts=5;
    nugget=1.0e-4;
    
    dd.set_name("gp_cost_function_dd");
    neigh.set_name("gp_cost_function_neigh");
    
    kptr->nn_srch(pt,npts+1,neigh,dd);
    if(dd.get_data(0)<1.0e-8){
        dd.remove(0);
        neigh.remove(0);
    }
    
    
    ell=dd.get_data(2);
   
    array_2d<double> covar,covarin;
    covar.set_name("gp_cost_covar");
    covarin.set_name("gp_cost_covar");
    
    int calculate=0,found_it=0;
    if(_neigh.get_dim()!=npts)calculate=1;
    else{
        for(i=0;i<npts && calculate==0;i++){
            found_it=0;
            for(j=0;j<npts && found_it==0;j++){
                if(_neigh.get_data(j)==neigh.get_data(i))found_it=1;
            }
            if(found_it==0){
                calculate=1;
            }
        }
    }
    
    calculate=1;
    if(calculate==1){
        covar.set_cols(npts);
        _covarin.set_cols(npts);
        for(i=0;i<npts;i++){
            for(j=i;j<npts;j++){
                covar.set(i,j,exp(-0.5*power(kptr->distance(neigh.get_data(i),neigh.get_data(j))/ell,2)));
                if(i==j){
                    covar.add_val(i,j,nugget);
                }
            }
        }
    
        
        invert_lapack(covar,covarin,0);
        
        for(i=0;i<npts;i++){
            _neigh.set(i,neigh.get_data(i));
            
            for(j=0;j<npts;j++){
                _covarin.set(i,j,covarin.get_data(i,j));
            }
        }
        _ell=ell;
        
    }
    
    /*array_1d<double> provisional_dd,sorted;
    for(i=0;i<npts;i++){
        provisional_dd.set(i,kptr->distance(pt,_neigh.get_data(i)));
    }

    sort_and_check(provisional_dd,sorted,_neigh);
    
    double ddmin=kptr->distance(pt,_neigh.get_data(0));
    double ddtest;
    for(i=0;i<kptr->get_pts();i++){
        if(i!=_neigh.get_data(0)){
            ddtest=kptr->distance(pt,i);
            if(ddtest<ddmin && ddtest>1.0e-8){
                printf("WARNING ddmin %e ddtest %e\n",ddmin,ddtest);
                exit(1);
            }
        }
    }
    
    printf("dd %e\n",ddmin);
    out[0]=-1.0*fn->get_data(_neigh.get_data(0));
    return;*/
        
    array_1d<double> covar_q;
    covar_q.set_name("gp_cost_covar_q");
    for(i=0;i<npts;i++){
        covar_q.set(i,exp(-0.5*power(kptr->distance(_neigh.get_data(i),pt)/_ell,2)));
    }
    
    double fbar=0.0;
    for(i=0;i<npts;i++){
        fbar+=fn->get_data(_neigh.get_data(i));
    }
    fbar=fbar/double(npts);
    
    out[0]=-1.0*fbar;
    for(i=0;i<npts;i++){
        for(j=0;j<npts;j++){
            out[0]-=covar_q.get_data(i)*_covarin.get_data(i,j)*(fn->get_data(_neigh.get_data(j))-fbar);
        }
    }
    
}

node_cost_function::~node_cost_function(){
    if(nodes!=NULL)delete [] nodes;

}

node_cost_function::node_cost_function(){
    printf("WARNING default constructor not implemented for node_cost_function\n");
    exit(1);
}

node_cost_function::node_cost_function(arrayOfNodes &nn,
array_1d<double> &min, array_1d<double> &max){
    _n_nodes=nn.get_dim();
    int i;
    
    _min.set_name("node_cost_function_min");
    _max.set_name("node_cost_function_max");
    
    for(i=0;i<min.get_dim();i++){
        _min.set(i,min.get_data(i));
        _max.set(i,max.get_data(i));
    }
    
    if(_n_nodes>0){
    
        nodes = new node*[_n_nodes];
    
        for(i=0;i<_n_nodes;i++){
            nodes[i]=nn(i);
        }
    }
    else{
        nodes=NULL;
    }
}

void node_cost_function::evaluate(array_1d<double> &pt, double *out){
    
    if(_n_nodes==0){
        out[0]=0.0;
        return;
    }
    
    double cost,costMin;
    int i;
    
    for(i=0;i<pt.get_dim();i++){
        if(pt.get_data(i)<_min.get_data(i)){
            out[0]=0.0;
            return;
        }
        
        if(pt.get_data(i)>_max.get_data(i)){
            out[0]=0.0;
            return;
        }
    }
    
    for(i=0;i<_n_nodes;i++){
        cost=nodes[i]->apply_model(pt);
        if(i==0 || cost<costMin){
            costMin=cost;
        }
    }
    
    out[0]=-1.0*costMin;
}

aps::aps(){
    printf("you called the APS constructor without paramters\n");
    printf("do not do that\n");
    exit(1);
}

aps::~aps(){
     
    int i;
    
    for(i=0;i<gg.get_dim();i++){
        delete [] paramnames[i];
    }
    delete [] paramnames;
    
    delete dice;
}

void aps::set_where(char *word){
    mu_storage.set_where(word);
    sig_storage.set_where(word);
    wide_pts.set_where(word);
}

aps::aps(int dim_in, int kk, double dd, int seed){
    
    ggWrap.set_iWhere(iAPS);
     
    sprintf(outname,"master_output.sav");
    sprintf(timingname,"timing_file.sav");
    
    mu_storage.set_name("aps_mu_storage");
    sig_storage.set_name("aps_sig_storage");
    wide_pts.set_name("aps_wide_pts");
    focus_pts.set_name("aps_focus_pts");
    characteristic_length.set_name("characteristic_length");
    range_max.set_name("range_max");
    range_min.set_name("range_min");
    ddUnitSpheres.set_name("ddUnitSpheres");
    refined_simplex.set_name("refined_simplex");
    simplex_start_pts.set_name("simplex_start_pts");
    
    _aps_wide_contents_buffer.set_name("_aps_wide_contents_buffer");
    _aps_wide_ss_buffer.set_name("_aps_wide_ss_buffer");
    
    write_every=1000;
    n_printed=0;
    do_bisection=1;

    called_focus=0;
    called_wide=0;
    ddNodeRatio=-1.0;
    
    time_penalty=0.5;
    
    n_wide=0;
    n_simplex=0;
    n_simplex_found=0;
        
    last_optimized=0;
    time_optimizing=0.0;
    time_refactoring=0.0;

    mindex_is_candidate=0;
    
    ct_aps=0;
    ct_simplex=0;

    time_aps=0.0;
    time_simplex=0.0;
    time_total=0.0;
    time_cleaning=0.0;
    time_writing=0.0;
    start_time=double(time(NULL));

    gg.set_kk(kk);
    
    ggWrap.set_delta_chisquared(dd);
  
    global_threshold=-2.0*chisq_exception;
    sphere_threshold=-2.0*chisq_exception;
    grat=1.0;
    
    dim=dim_in;
    paramnames=new char*[dim];
    int i;
    for(i=0;i<dim;i++){
        paramnames[i]=new char[letters];
        sprintf(paramnames[i],"p%d",i);
    }
    
    for(i=0;i<dim;i++){
        characteristic_length.set(i,-1.0);
    }
    
    if(seed<-1)seed=int(time(NULL));
    dice=new Ran(seed);
    
    start_timingfile();
    printf("dim is %d dd %d\n",dim,dim_in);
    
}

void aps::disable_bisection(){
    do_bisection=0;
}

void aps::enable_bisection(){
    do_bisection=1;
}

void aps::start_timingfile(){
    FILE *output;
    output=fopen(timingname,"a");
    fprintf(output,"\n# pts_stored calls_to_chisq time_in_chisq ");
    fprintf(output," time_per_chisq total_time time_per_pt -- ");
    
    fprintf(output,"ct_aps time_aps -- ");
    fprintf(output,"ct_simplex time_simplex -- ");
    fprintf(output,"time_optimizing_gp -- time_refactoring -- ");
    fprintf(output,"chisq_min target -- volumes -- called_wide ");
    fprintf(output,"called_focus");
    fprintf(output,"\n");
    
    fclose(output);
}

void aps::set_write_every(int ii){
    write_every=ii;
}

void aps::set_outname(char *word){
    int i;
    for(i=0;i<letters && word[i]!=0;i++)outname[i]=word[i];
    outname[i]=0;
}

void aps::set_timingname(char *word){
    int i;
    for(i=0;i<letters && word[i]!=0;i++)timingname[i]=word[i];
    timingname[i]=0;
    
    start_timingfile();  
}

void aps::set_grat(double nn){
    grat=nn;
}

void aps::assign_covariogram(covariance_function *cc){
    gg.assign_covariogram(cc);
}

void aps::assign_chisquared(chisquared *cc){
    ggWrap.set_chisq(cc);
}

void aps::set_characteristic_length(int dex, double nn){
    characteristic_length.set(dex,nn);
}

void aps::initialize(int npts,array_1d<double> &min, array_1d<double> &max){
    array_2d<double> q;
    initialize(npts,min,max,q);
}

void aps::initialize(int npts, array_1d<double> &min, array_1d<double> &max, 
    array_2d<double> &guesses){

    int nguesses=guesses.get_rows();
    printf("dim %d\n",dim);
    
    set_where("aps_initializer");
    
    int i,j,k;
    
    array_1d<double> ff,vector;
    ff.set_name("aps_initializer_ff");
    vector.set_name("aps_initializer_vector");
    
    array_2d<double> data;
    data.set_name("aps_initializer_data");
    
    ggWrap.set_strad(&strad);
    
    for(i=0;i<nguesses;i++){
        for(j=0;j<dim;j++)vector.set(j,guesses.get_data(i,j));

        ff.set(i,ggWrap.call_chisq(vector));
        
        while(!(ff.get_data(i)<chisq_exception)){
            for(j=0;j<dim;j++){
                vector.set(j,min.get_data(j)+dice->doub()*(max.get_data(j)-min.get_data(j)));
            }

            ff.set(i,ggWrap.call_chisq(vector));
        }
        
        data.add_row(vector);
    }
    
    printf("done guessing\n");
    for(;i<npts;i++){
        //printf("%d\n",i);
        ff.set(i,2.0*chisq_exception);
        while(!(ff.get_data(i)<chisq_exception)){
            
            for(j=0;j<dim;j++)vector.set(j,min.get_data(j)+dice->doub()*(max.get_data(j)-min.get_data(j)));
            ff.set(i,ggWrap.call_chisq(vector));
            
            /*printf("ff %e \n",ff.get_data(i));
            for(j=0;j<dim;j++){
                printf("%e\n",vector.get_data(j));
            }*/
            
        }
        data.add_row(vector);
    }
    printf("done assembling points\n");
    
    
    for(i=0;i<dim;i++){
        range_max.set(i,max.get_data(i));
        range_min.set(i,min.get_data(i));
    }
    
    array_1d<double> ggmin,ggmax;
    
    for(i=0;i<dim;i++){
        ggmin.set(i,0.0);
        if(characteristic_length.get_data(i)<0.0){
            ggmax.set(i,(range_max.get_data(i)-range_min.get_data(i)));
            characteristic_length.set(i,range_max.get_data(i)-range_min.get_data(i));
        }
        else{
            ggmax.set(i,characteristic_length.get_data(i));
        }
    }
    
    printf("calling on gg.initialize\n");
    gg.initialize(data,ff,ggmin,ggmax);
    ggWrap.set_gp(&gg);
    
    if(gg.get_dim()!=dim){
        printf("WARNING gg.get_dim %d dim %d\n",
        gg.get_dim(),dim);
        
        exit(1);
    }
    

    if(ggWrap.get_chisq_dim()!=gg.get_dim() || ggWrap.get_chisq_dim()!=dim){
        printf("WARNING chisq dim %d gg %d dim %d\n",
        ggWrap.get_chisq_dim(),gg.get_dim(),dim);
            
        exit(1);
    }
    
    printf("time to optimize\n");
    
    optimize();
    
    printf("done optimizing\n");
    
    ct_aps=0;
    
    ff.reset();
    data.reset();
    
    
    for(i=0;i<gg.get_pts();i++){
        wide_pts.add(i);
        mu_storage.add(-2.0);
        sig_storage.add(-2.0);
    }
    
    int before_grad=ggWrap.get_called();
        
    double nn;
    for(i=0;i<gg.get_pts();i++){
        if(i==0 || gg.get_fn(i)<nn){
            j=i;
            nn=gg.get_fn(i);
        }
    }
    
    printf("time to set chimin\n");
    if(nn<ggWrap.get_chimin() || ggWrap.get_chimin()<0.0){
        ggWrap.set_chimin(nn,(*gg.get_pt(j)),j);
        mindex_is_candidate=1;
    }
    
    printf("about to write\n");
        
    write_pts();
    
    /*for(i=0;i<gg.get_pts();i++){
        j=is_it_a_candidate(i);
        if(j==1){
            set_as_candidate(i);
            
        }
    }*/
    
    set_where("nowhere");
}

void aps::set_min(array_1d<double> &mn){
    int i;
    if(mn.get_dim()!=dim){
        printf("WARNING setting min with dim %d but should be %d\n",
        mn.get_dim(),dim);
        
        exit(1);
    }
    
    for(i=0;i<dim;i++){
        range_min.set(i,mn.get_data(i));
    }
}

void aps::set_max(array_1d<double> &mx){
    int i;
    if(mx.get_dim()!=dim){
        printf("WARNING setting max with dim %d but should be %d\n",
        mx.get_dim(),dim);
        
        exit(1);
    }
    
    for(i=0;i<dim;i++){
        range_max.set(i,mx.get_data(i));
    }
}

void aps::set_hyper_parameters(array_1d<double> &hh){
    gg.set_hyper_parameters(hh);
}

void aps::resume(){
    resume(outname);
}

void aps::resume(char *filename){
    
    int i,ct=0;
    array_2d<double> data;
    array_1d<double> ff;
    double nn,mu,sig,local_min;
    int ling,i_min;
    char word[500];
    FILE *input=fopen(filename,"r");
    for(i=0;i<dim+5;i++)fscanf(input,"%s",word);
    printf("final word is %s\n",word);
    
    data.set_cols(dim);
    while(fscanf(input,"%le",&nn)>0){
        data.set(ct,0,nn);
        for(i=1;i<dim;i++){
            fscanf(input,"%le",&nn);
            data.set(ct,i,nn);
        }
        
        fscanf(input,"%le",&nn);
        ff.set(ct,nn);
        
        if(ct==0 || nn<local_min){
            local_min=nn;
            i_min=ct;
        }
        
        
        ggWrap.set_whereFrom(ct,ling);
         
        fscanf(input,"%le %le %d",&mu,&sig,&ling);
        if(ling==0){
            
            wide_pts.add(ct);
            mu_storage.add(mu);
            sig_storage.add(sig);
        }
        
        ct++;
    }
    
    
    fclose(input);
    
    array_1d<double> ggmax,ggmin;
    for(i=0;i<dim;i++){
        ggmin.set(i,0.0);
        if(characteristic_length.get_data(i)<0.0){
            ggmax.set(i,range_max.get_data(i)-range_min.get_data(i));
        }
        else{
            ggmax.set(i,characteristic_length.get_data(i));
        }
    }
    
    gg.initialize(data,ff,ggmin,ggmax);
    ggWrap.set_gp(&gg);
    
    printf("initialized gg\n");
    ggWrap.set_chimin(local_min,(*gg.get_pt(i_min)),i_min);
    mindex_is_candidate=1;
    
    assess_node(ggWrap.get_global_mindex());
  
    write_pts();
}

void aps::set_target(double tt){
    ggWrap.assert_target();
    strad.set_target(tt);
}

double aps::get_chimin(){
    return ggWrap.get_chimin();
}

void aps::get_minpt(array_1d<double> &output){
    int i;
    for(i=0;i<gg.get_dim();i++){
        output.set(i,ggWrap.get_minpt(i));
    }
}

int aps::is_it_a_candidate(int dex){

    if(dex>=gg.get_pts() || dex<0){
        printf("WARNING assessing candidacy of %d but total %d\n",dex,gg.get_pts());
        exit(1);
    }
    
    int i,use_it,ic;
    array_1d<double> mid_pt;
    double chitrial;

    for(i=0;i<forbidden_candidates.get_dim();i++){
        if(dex==forbidden_candidates.get_data(i))return 0;
    }
    
    for(i=0;i<known_minima.get_dim();i++){
        if(dex==known_minima.get_data(i))return 0;
    }
    
    if(gg.get_fn(dex)<ggWrap.get_chimin()+grat*(fabs(global_threshold)-ggWrap.get_chimin()) && 
        gg.get_fn(dex)>strad.get_target()){
          
        /*
        * Test the point halfway between the new point (dex) and the midpoint.
        * If chisquared at the test point is more than chimin + 0.75 *(fn(dex) - chimin),
        * then this is a new candidate for function minimization.  If not, this is probably just 
        * a part of the same low chisquared locus associated with chimin.  Do not consider
        * it a candidate for function minimization.
        */
        
        for(ic=0;ic<known_minima.get_dim();ic++){
            for(i=0;i<dim;i++){
                mid_pt.set(i,0.5*(gg.get_pt(known_minima.get_data(ic),i)+gg.get_pt(dex,i)));
            }
        
            ggWrap.evaluate(mid_pt,&chitrial);

            if(chitrial<=gg.get_fn(dex)-0.25*(gg.get_fn(dex)-gg.get_fn(known_minima.get_data(ic)))){
                forbidden_candidates.add(dex);
                return 0;
            }
        }
        return 1;
    }
    else return 0;
}

int aps::find_global_minimum(array_1d<int> &neigh){
    return find_global_minimum(neigh,-1);
}

int aps::find_global_minimum(array_1d<int> &neigh, int limit){
    
    /*
    * Use simplex minimization, seeded by the sampled points whose
    * indices are stored in neigh[] to search for a new local or global
    * minimum in chisquared.
    */
    
    if(neigh.get_dim()!=dim+1){
        printf("WARNING you just called find_global_minimum with an improper number of neighbors\n");
        exit(1);
    }
    
    if(dim!=gg.get_dim()){
        printf("WARNING in find_global_minimum dim %d but gg says %d\n",
        dim,gg.get_dim());
        
        exit(1);
    }
    
    set_where("find_global_minimum");
    
    if(gg.is_kptr_null()==1){
        printf("WARNING gg.kptr is null in find_global_minimum\n");
        exit(1);
    }
    
    array_1d<double> found_minimum;
    array_2d<double> simplex_seeds;
    found_minimum.set_name("aps_found_minimum");
    simplex_seeds.set_name("aps_simplex_seeds");
    int i,j;
    
    int ibefore=ggWrap.get_called();
    double start_min;
    
    simplex_seeds.set_cols(ggWrap.get_dim());
    
    for(i=0;i<ggWrap.get_dim()+1;i++){
        if(i==0 || ggWrap.get_fn(neigh.get_data(i))<start_min){
            start_min=ggWrap.get_fn(neigh.get_data(i));
        }
    }
    
    for(i=0;i<ggWrap.get_dim()+1;i++){
        for(j=0;j<ggWrap.get_dim();j++){
            simplex_seeds.set(i,j,ggWrap.get_pt(neigh.get_data(i),j));
        }
    }
    
    //node_cost_function cost(nodes,range_min,range_max);
    
    array_2d<double> copy_data;
    array_1d<double> copy_ff,k_min,k_max;
    
    copy_data.set_name("aps_copy_data");
    copy_ff.set_name("aps_copy_ff");
    k_min.set_name("aps_k_min");
    k_max.set_name("aps_k_max");
    
    for(i=0;i<ggWrap.get_dim();i++){
        k_min.set(i,0.0);
        k_max.set(i,characteristic_length.get_data(i));
    }
    
    copy_data.set_cols(ggWrap.get_dim());
    for(i=0;i<ggWrap.get_pts();i++){
        copy_ff.set(i,ggWrap.get_fn(i));
        for(j=0;j<ggWrap.get_dim();j++){
            copy_data.set(i,j,ggWrap.get_pt(i,j));
        }
    }
    
    kd_tree copy_kd(copy_data,k_min,k_max);
    
    gp_cost_function cost(&copy_kd, &copy_ff, range_min, range_max);
    
    simplex_minimizer ffmin;
    array_1d<double> min,max;
    min.set_name("find_global_minimum_min");
    max.set_name("find_global_maximum_max");
    
    for(i=0;i<ggWrap.get_dim();i++){
        max.set(i,range_min.get_data(i)+0.1*(range_max.get_data(i)-range_min.get_data(i)));
        min.set(i,range_min.get_data(i));
    }
    
    ffmin.set_minmax(min,max);
    ffmin.set_chisquared(&ggWrap);
    ffmin.set_dice(dice);
    ffmin.use_gradient();
    if(nodes.get_dim()>0){
        ffmin.set_cost(&cost);
    }
    //ffmin.freeze_temp();
    ffmin.find_minimum(simplex_seeds,found_minimum);
    
    array_1d<int> mindex;
    array_1d<double> dd;
    ggWrap.nn_srch(found_minimum,1,mindex,dd);
    if(dd.get_data(0)>1.0e-10){
        printf("WARNING after simplex search dd %e\n",dd.get_data(0));
    }
    

    
    printf("    simplex found %e from %e\n",
    ggWrap.get_fn(mindex.get_data(0)),start_min);
    printf("    in %d\n",ggWrap.get_called()-ibefore);
    
    //if(nodes.get_dim()>0)exit(1);
    
    assess_node(mindex.get_data(0));
    
    return mindex.get_data(0);
}


int aps::get_ct_aps(){
    return ct_aps;
}

int aps::get_ct_simplex(){
    return ct_simplex;
}

int aps::get_called(){
    return ggWrap.get_called();
}

int aps::get_n_pts(){
    return gg.get_pts();
}

int aps::get_n_simplex(){
    return n_simplex;
}

int aps::get_n_simplex_found(){
    return n_simplex_found;
}

int aps::get_n_active_nodes(){
    int i,ans=0;
    for(i=0;i<nodes.get_dim();i++){
        if(nodes(i)->is_it_active()==1)ans++;
    }
    return ans;
}

array_1d<double>* aps::get_pt(int i){
    if(i>=gg.get_pts() || i<0){
        printf("WOAH; asked APS for pt %d but only have %d\n",i,gg.get_pts());
    }
    
    return gg.get_pt(i);
}

int aps::get_nn(array_1d<double> &pt){
    array_1d<int> neigh;
    array_1d<double> ddneigh;
    
    gg.nn_srch(pt,1,neigh,ddneigh);
    
    return neigh.get_data(0);
}

double aps::get_chival(int dex){
    if(dex>=gg.get_pts() || dex<0){
        printf("WARNING asked aps for chi at %d but pts %d\n",
        dex,gg.get_pts());
        
        exit(1);
    }
    
    return gg.get_fn(dex);
}

int aps::get_n_centers(){
    return nodes.get_dim();
}

double aps::get_pt(int dex, array_1d<double> &output){
    if(dex>=gg.get_pts() || dex<0){
        printf("WARNING wanted point %d but total %d\n",dex,gg.get_pts());
    }
    
    gg.get_pt(dex,output);
    
    return gg.get_fn(dex);
}


void aps::get_good_points(array_2d<double> &outpoints){
    outpoints.reset();
    outpoints.set_cols(ggWrap.get_dim());
    int i,j,dex;
    for(i=0;i<ggWrap.get_pts();i++){
        if(ggWrap.get_fn(i)<=ggWrap.get_target()){
            dex=outpoints.get_rows();
            for(j=0;j<ggWrap.get_dim();j++){
                outpoints.set(dex,j,ggWrap.get_pt(i,j));
            }
        }
    }
}

void aps::guess(array_1d<double> &pt){

    double chitrue;
    int ibefore=ggWrap.get_called();
    int actually_added;

    ggWrap.evaluate(pt,&chitrue);

    ct_aps+=ggWrap.get_called()-ibefore;
}

void aps::search(){
    
    double before=double(time(NULL));
 
    double aps_score,simplex_score;
    int i;
    
    //aps_score=ct_aps;
    //simplex_score=ct_simplex;
    
    aps_score=time_aps+ct_aps*time_penalty;
    simplex_score=time_simplex+ct_simplex*time_penalty;
    
    int i_simplex;
    if(simplex_score<aps_score){
        i=ggWrap.get_called();
        simplex_search();
        ct_simplex+=ggWrap.get_called()-i;
        time_simplex+=double(time(NULL))-before;
    }
    
    aps_search();
    
    if(gg.get_pts()>n_printed+write_every){
        write_pts();
    }
      
    time_total+=double(time(NULL))-before;
}

void aps::get_interesting_boxes(array_1d<int> &acceptableBoxes){

    array_1d<double> trial;
    trial.set_name("get_interesting_boxes_trial");
    int ic;
    
    double span,allowedDistance,dd;
    int i,j,use_it;
    for(i=0;i<ggWrap.get_nboxes();i++){
        use_it=1;
        for(j=0;j<ggWrap.get_box_contents(i) && use_it==1;j++){
            if(ggWrap.get_fn(ggWrap.get_box_contents(i,j))<ggWrap.get_target()){
                use_it=0;
            }
        }
        
        if(nodes.get_dim()>0 && use_it==1){
            for(j=0;j<ggWrap.get_dim();j++){
                trial.set(j,0.5*(ggWrap.get_box_max(i,j)+ggWrap.get_box_min(i,j)));
            }
            
            ic=find_nearest_center(trial);
            
            allowedDistance=0.0;
            for(j=0;j<ggWrap.get_dim();j++){
                span=nodes(ic)->get_max(j)-nodes(ic)->get_min(j);
                allowedDistance+=power(span/characteristic_length.get_data(j),2);
            }
            allowedDistance=sqrt(allowedDistance);
            
            dd=ggWrap.distance(ic,trial);
            if(dd<allowedDistance)use_it=0;
        }
        
        
        if(use_it==1){
            acceptableBoxes.add(i);
        }
    
    }


}

double aps::calculate_dchi(int ibox){
    if(nodes.get_dim()==0)return 0.0;
    
    ggWrap.freeze_boxes();

    int i,j,imin;
    double chimin,chi2,dchi,dchibest;
    array_1d<double> midpt;
    midpt.set_name("aps_calculate_dchi_midpt");
    
    /*for(i=0;i<ggWrap.get_box_contents(ibox);i++){
        dchi=ggWrap.get_fn(ggWrap.get_box_contents(ibox,i));
	if(i==0 || dchi<chimin){
	    chimin=dchi;
	    imin=ggWrap.get_box_contents(ibox,i);
	}
    }*/
    
    for(i=0;i<ggWrap.get_dim();i++){
        midpt.set(i,0.5*(ggWrap.get_box_max(ibox,i)+ggWrap.get_box_min(ibox,i)));
    }
    ggWrap.evaluate(midpt,&chimin);

    int ic;
    for(ic=0;ic<nodes.get_dim();ic++){
        i=nodes(ic)->get_center();
	for(j=0;j<ggWrap.get_dim();j++){
	    midpt.set(j,0.5*(ggWrap.get_pt(imin,j)+ggWrap.get_pt(i,j)));
	}
	
	ggWrap.evaluate(midpt,&chi2);
	dchi=chi2-chimin;
	if(ic==0 || dchi<dchibest){
	    dchibest=dchi;
	}
    }
    
    ggWrap.unfreeze_boxes();
    
    return dchibest;
}

int aps::simplex_search(){
    /*
    Select the boxes with no good points in them which are sufficiently
    far from their nearest nodes.
    
    Use a GP to approximate the probability that each of these boxes contains
    zero good points.
    
    Do a simplex minimization of chisquared (directly) seeded with points
    in the box whose probability of containing a good point is the lowest.

    */
    
    ggWrap.set_iWhere(iSimplex);
    
    int i,j,iout;
    array_1d<int> acceptableBoxes,seed;
    acceptableBoxes.set_name("aps_simplex_acceptableBoxes");
    seed.set_name("aps_simplex_seed");
    
    get_interesting_boxes(acceptableBoxes);
    if(acceptableBoxes.get_dim()==0){
        iout=aps_wide();
        n_wide++;
        return iout;
    }
    
    /*do not use any boxes containing smallest seeds from previous simplexes*/
    int ibox,isOkay;
    for(i=0;i<simplex_start_pts.get_dim();i++){
        ibox=ggWrap.find_box(simplex_start_pts.get_data(i));
        isOkay=1;
        for(j=0;j<acceptableBoxes.get_dim() && isOkay==1;j++){
            if(acceptableBoxes.get_data(j)==ibox){
                isOkay=0;
                acceptableBoxes.remove(j);
            }
        }
    }
    
    array_1d<double> trial;
    trial.set_name("simplex_search_trial");
    
    if(acceptableBoxes.get_dim()==0){
        printf("\n\nI'm sorry; there were no acceptable boxes\n\n");
        iout=aps_wide();
        n_wide++;
        return iout;
    }
    
    n_simplex++;
    int dex,iSmallestSeed;
    double smallestSeed,chitrial;
    smallestSeed=2.0*chisq_exception;
    while(seed.get_dim()<ggWrap.get_dim()+1){
        for(i=0;i<ggWrap.get_dim();i++){
            trial.set(i,range_min.get_data(i)
                +dice->doub()*(range_max.get_data(i)-range_min.get_data(i)));
        }
            
        ggWrap.evaluate(trial,&chitrial,&dex);
            
        if(dex>=0){
            if(chitrial<smallestSeed){
                smallestSeed=chitrial;
                iSmallestSeed=dex;
            }
            seed.add(dex);
        }
    } 
    
    j=nodes.get_dim();
    i=find_global_minimum(seed);
    if(nodes.get_dim()>j){
        n_simplex_found++;
    }
    simplex_start_pts.add(iSmallestSeed);
        
    if(i>=0){
        printf("    box found %e from %e\n",ggWrap.get_fn(i),smallestSeed);
    }
    else{
        printf("    box found junk\n");
    }
    
    //if(n_simplex>1)exit(1);
    
}

int aps::aps_wide(){
    
    /*
    Run the simplex search that seeks to maximize S
    
    First, however, select the boxes that have no good points in them
    and are sufficiently far away from their nearest centers.
    
    Evaluate S at the center of each of these boxes.
    
    Maximize S in the box whose center has the maximum S.
    
    After calling simplex_strad, the point at which S is maximized
    will be stored in the global array_1d<double> simplex_best
    
    That is the point at which to evaluate chisquared.
    */
    
    ggWrap.set_iWhere(iAPS);
    
    array_1d<int> acceptableBoxes;
    acceptableBoxes.set_name("aps_wide_acceptableBoxes");
    get_interesting_boxes(acceptableBoxes);

    int i,j,chosenBox,ii,ibox;
    array_1d<double> midpt;
    midpt.set_name("aps_wide_midpt");
    
    double mu,sig,ss,ssbest;
    
    if(acceptableBoxes.get_dim()==0){
        chosenBox=dice->int32()%ggWrap.get_nboxes();
    }
    else{
        for(ii=0;ii<acceptableBoxes.get_dim();ii++){
            ibox=acceptableBoxes.get_data(ii);
            
            if(ibox<_aps_wide_contents_buffer.get_dim() && 
               ggWrap.get_box_contents(ibox)==_aps_wide_contents_buffer.get_data(ibox)){
               
               
               ss=_aps_wide_ss_buffer.get_data(ibox);    
           }
           else{
                if(acceptableBoxes.get_dim()<2000){
                    for(i=0;i<ggWrap.get_dim();i++){
                        midpt.set(i,0.5*(ggWrap.get_box_min(ibox,i)+ggWrap.get_box_max(ibox,i)));
                    }
            
                    mu=ggWrap.user_predict(midpt,&sig);
                }
                else{
                
                    mu=0.0;
                    for(i=0;i<ggWrap.get_box_contents(ibox);i++){
                        mu+=ggWrap.get_fn(ggWrap.get_box_contents(ibox,i));
                    }
                    mu=mu/double(ggWrap.get_box_contents(ibox));
                    sig=0.0;
                    for(i=0;i<ggWrap.get_box_contents(ibox);i++){
                        sig+=power(mu-ggWrap.get_fn(ggWrap.get_box_contents(ibox,i)),2);
                    }
                    sig=sig/double(ggWrap.get_box_contents(ibox));
                    sig=sqrt(sig);
                
                }
                
                ss=ggWrap.straddle_value(mu,sig);
                
                _aps_wide_contents_buffer.set(ibox,ggWrap.get_box_contents(ibox));
                _aps_wide_ss_buffer.set(ibox,ss);
                
            }
            
            if(ii==0 || ss>ssbest){
                chosenBox=ibox;
                ssbest=ss;
            }
            
        }
        
        
    }
    
    /*printf("ABOUT TO CALL simplex_strad %d\n",chosenBox);
    for(i=0;i<ggWrap.get_dim();i++){
        printf("%d -- %e %e -- %e %e\n",
        i,searchMin.get_data(i),searchMax.get_data(i),
        ggWrap.get_box_min(chosenBox,i),ggWrap.get_box_max(chosenBox,i));
    }*/
    
    sig=simplex_strad(chosenBox);
    
    double chitrue;
    int actually_added,ic,use_it;
     
    int iout=-1;
    
    if(simplex_best.get_dim()==gg.get_dim()){
        ggWrap.evaluate(simplex_best,&chitrue,&actually_added);
        
        iout=actually_added;
        
    }
    aps_wide_post_process(actually_added,simplex_mu_best,simplex_sig_best);
    return iout;
}

void aps::aps_wide_post_process(int actually_added, double mu, double sig){    
    
    if(actually_added<0)return;
    
    array_1d<double> pt;
    pt.set_name("aps_wide_post_process_pt");
    
    int i;
    for(i=0;i<ggWrap.get_dim();i++){
        pt.set(i,ggWrap.get_pt(actually_added,i));
    }
    
    double chitrue=ggWrap.get_fn(actually_added);
    int bisect_it=0;
    
    wide_pts.add(actually_added);
    mu_storage.add(mu);
    sig_storage.add(sig);
    
    int ic;
    array_1d<double> unit_v,dd_sphere;
    unit_v.set_name("aps_wide_post_process_unit_v");
    dd_sphere.set_name("aps_wide_post_process_dd_sphere");
    
    array_1d<int> neigh_sphere;
    neigh_sphere.set_name("aps_wide_post_process_neigh_sphere");     
         
    if(do_bisection==1){
        if(sphere_threshold<0.0 && chitrue<global_threshold){
            /*
            If the KD-tree of boundary points projected onto a unit sphere has not
            yet been established, and chisquared is less than the 1/10 quantile of
            points sampled by aps_wide, do bisection
            */
                    
            bisect_it=1;
        }
        else if(sphere_threshold>0.0){
            /*
            If the KD-tree of boundary points projected onto a unit sphere has
            been established, find the nearest low-chisquared center, project
            this point onto a unit sphere about that center, compare the distance
            between that projected point and its nearest projected neighbor
            to the last few hundred such distances.  If that distance is greater
            than the 2/3 quantile of such distances, do bisection.
            */
            ic=find_nearest_center(pt);
            if(ic>=0){
                if(ggWrap.is_unitSpheres_null()==0){
                    project_to_unit_sphere(ic,pt,unit_v);
                    ggWrap.unitSpheres_nn_srch(unit_v,1,neigh_sphere,dd_sphere);
                    ddUnitSpheres.add(dd_sphere.get_data(0));
                   
                    if(dd_sphere.get_data(0)>sphere_threshold){
                        bisect_it=1;
                    }
                }
            }
                
        }
            
            
        if(bisect_it==1){
            bisection(simplex_best,chitrue);
        }
    }
     
    int use_it,inode;
    array_1d<double> midpt;
    midpt.set_name("aps_wide_post_process_midpt");
    
    double chimid;
           
    /*
    If chisquared<chisquared_lim, assess whether or not we have
    discovered a new low-chisquared region.
    */
    if(chitrue<strad.get_target()){
        use_it=1;
        for(ic=0;ic<nodes.get_dim() && use_it==1;ic++){
            inode=nodes(ic)->get_center();
            for(i=0;i<gg.get_dim();i++){
                midpt.set(i,0.5*(simplex_best.get_data(i)+gg.get_pt(inode,i)));
            }
                    
            ggWrap.evaluate(midpt,&chimid,&i);
                    
            if(chimid<strad.get_target())use_it=0;
        }
                
        if(use_it==1){
            assess_node(actually_added);
        }
    }

   
}

double aps::simplex_metric(array_1d<double> &pt, array_1d<double> &min_bound, array_1d<double> &max_bound){
    
    
    if(max_bound.get_dim()!=gg.get_dim() || min_bound.get_dim()!=gg.get_dim()){
        printf("cannot evaluate simplex_metric; have not set min and max yet\n");
        throw -1;
    
    }
    
    int i;
    
    if(ggWrap.in_bounds(pt)==0){
        return 2.0*chisq_exception;
    }
    
    for(i=0;i<gg.get_dim();i++){
        if(pt.get_data(i)>max_bound.get_data(i) || 
            pt.get_data(i)<min_bound.get_data(i)){
            
            return 2.0*chisq_exception;
        
        }
    }
    
   
    double mu,sig,stradval;
    mu=ggWrap.user_predict(pt,&sig);
    
    stradval=strad(mu,sig);
    
    simplex_ct++;
    simplex_calls++;
    
    /*
    If a new local maximum is found for S, store that and reset simplex_ct to zero
    */
    if(stradval>simplex_strad_best){
        simplex_strad_best=stradval;
        
        simplex_mu_best=mu;
        simplex_sig_best=sig;
        
        simplex_ct=0;
        
        for(i=0;i<gg.get_dim();i++){
            simplex_best.set(i,pt.get_data(i));
        }
    }
    
    /*
    Because the Nelder-Mead (1965) simplex is designed to minimize functions,
    we return -S to find the local maxmimum in S
    */
    return -1.0*stradval;
    
}

double aps::simplex_strad(int ibox){
    
   int i,j;
   double nn;
   
   int searches_before=ggWrap.get_ct_search();

   //printf("\nSTARTING APS\n");

   array_1d<double> min_bound,max_bound;
   min_bound.set_name("aps_simplex_strad_min_bound");
   max_bound.set_name("aps_simplex_strad_max_bound");
   
   for(i=0;i<gg.get_dim();i++){
       min_bound.set(i,ggWrap.get_box_min(ibox,i));
       max_bound.set(i,ggWrap.get_box_max(ibox,i));
   }
   
   simplex_best.reset();
   
   simplex_strad_best=-2.0*chisq_exception;
   simplex_ct=0;
   simplex_calls=0;
   ggWrap.reset_cache();

   
   /*
   now do simplex search using simplex_metric
   
   This is essentially the same algorithm as run by find_global_minimum()
   except with different alpha, beta, and gamma parameters and 
   different stopping conditions
   
   Note: this algorithm does not utilize the modified gradient descent algorithm
   by which find_global_minimum tried to avoid getting trapped in narrow valleys
   of chisquared.
   */
   array_2d<double> pts;
   array_1d<double> pbar,ff,pstar,pstarstar;
   double fstar,fstarstar;
   int ih,il;
   double alpha=1.0,beta=0.9,gamma=1.1;
   double mu,sig;
   
   pstar.set_dim(gg.get_dim());
   pstarstar.set_dim(gg.get_dim());
   pts.set_dim(gg.get_dim()+1,gg.get_dim());
   pbar.set_dim(gg.get_dim());
   ff.set_dim(gg.get_dim()+1);
   
   for(i=0;i<gg.get_dim()+1;i++){
       nn=2.0*chisq_exception;
       while(!(nn<chisq_exception)){
           for(j=0;j<gg.get_dim();j++){
               pts.set(i,j,min_bound.get_data(j)+dice->doub()*(max_bound.get_data(j)-min_bound.get_data(j)));
           }
           nn=simplex_metric(*pts(i),min_bound,max_bound);
       }
       
       ff.set(i,nn);
       if(i==0 || ff.get_data(i)<ff.get_data(il))il=i;
       if(i==0 || ff.get_data(i)>ff.get_data(ih))ih=i;
   }
   
   sig=2.0*chisq_exception;
   
   /*
   Stop the simplex search when either:
   
   1) the standard deviation of S values in the simplex is less than delta_chisquared
   
   2) more than 1000 evaluations of chisquared have been made without improving the best
   found value of S
   
   S wil be evaluated by simplex_metric, which will be charged with storing the local maximum
   of S discovered and of keeping track of simplex_ct for convergence purposes
   */
   while(sig>ggWrap.get_delta_chisquared() && simplex_ct<1000){
       for(i=0;i<gg.get_dim();i++){
           pbar.set(i,0.0);
           for(j=0;j<gg.get_dim()+1;j++){
               if(j!=ih){
                   pbar.add_val(i,pts.get_data(j,i));
               }
           }
           pbar.divide_val(i,double(gg.get_dim()));
       }
       
       for(i=0;i<gg.get_dim();i++){
           pstar.set(i,(1.0+alpha)*pbar.get_data(i)-alpha*pts.get_data(ih,i));
       }
       fstar=simplex_metric(pstar,min_bound,max_bound);
       
       if(fstar<ff.get_data(ih) && fstar>ff.get_data(il)){
           ff.set(ih,fstar);
           for(i=0;i<gg.get_dim();i++){
               pts.set(ih,i,pstar.get_data(i));
           }
       }
       else if(fstar<ff.get_data(il)){
           for(i=0;i<gg.get_dim();i++){
               pstarstar.set(i,gamma*pstar.get_data(i)+(1.0-gamma)*pbar.get_data(i));
           }
           fstarstar=simplex_metric(pstarstar,min_bound,max_bound);
           
           if(fstarstar<ff.get_data(il)){
               for(i=0;i<gg.get_dim();i++)pts.set(ih,i,pstarstar.get_data(i));
               ff.set(ih,fstarstar);
           }
           else{
               for(i=0;i<gg.get_dim();i++)pts.set(ih,i,pstar.get_data(i));
               ff.set(ih,fstar);
           }
       }
       
       j=1;
       for(i=0;i<gg.get_dim()+1;i++){
           if(fstar<ff.get_data(i) && i!=ih){
               j=0;
           }
       }
       
       if(j==1){
           for(i=0;i<gg.get_dim();i++){
               pstarstar.set(i,beta*pts.get_data(ih,i)+(1.0-beta)*pbar.get_data(i));
               
           }
           fstarstar=simplex_metric(pstarstar,min_bound,max_bound);
           
           if(fstarstar<ff.get_data(ih)){
               for(i=0;i<gg.get_dim();i++){
                   pts.set(ih,i,pstarstar.get_data(i));
               }
               ff.set(ih,fstarstar);
           }
           else{
               for(i=0;i<gg.get_dim()+1;i++){
                   if(i==0 || ff.get_data(i)<ff.get_data(il)){
                       il=i;
                   }
               }
               for(i=0;i<gg.get_dim()+1;i++){
                   if(i!=il){
                       for(j=0;j<gg.get_dim();j++){
                           mu=0.5*(pts.get_data(i,j)+pts.get_data(il,j));
                           pts.set(i,j,mu);
                       }
                       ff.set(i,simplex_metric(*pts(i),min_bound,max_bound));
                   }
               }
           }
       }
       
       for(i=0;i<gg.get_dim()+1;i++){
           if(i==0 || ff.get_data(i)<ff.get_data(il)){
               il=i;
           }
           
           if(i==0 || ff.get_data(i)>ff.get_data(ih)){
               ih=i;
           }
       }
       
       mu=0.0;
       for(i=0;i<gg.get_dim()+1;i++){
           mu+=ff.get_data(i);
       }  
       mu=mu/double(gg.get_dim()+1);
       
       sig=0.0;
       for(i=0;i<gg.get_dim()+1;i++){
           sig+=power(ff.get_data(i)-mu,2);
       }
       sig=sig/double(gg.get_dim()+1);
       sig=sqrt(sig);
       
       if(ggWrap.get_ct_search()>searches_before+1){
           printf("\nJUST DID APS WIDE WITH %d SEARCHES %d -- %d\n\n",ggWrap.get_ct_search()-searches_before,
           simplex_calls,ibox);
           printf("%d %d\n",ggWrap.get_box_contents(305),ggWrap.get_box_contents(306));
           for(i=0;i<ggWrap.get_dim();i++){
               printf("%d -- %e %e -- %e -- %e %e -- %e %e\n",
               i,min_bound.get_data(i),max_bound.get_data(i),
               max_bound.get_data(i)-min_bound.get_data(i),
               ggWrap.get_box_min(ibox,i),ggWrap.get_box_max(ibox,i),
               ggWrap.get_box_min(306,i),ggWrap.get_box_max(306,i));
           }
           exit(1);

       }
       
   }
   

   
   ggWrap.reset_cache();
   return sig;
}

void aps::random_focus(int ic){
    array_1d<double> trial,rr;
    double chitrue;
    int actually_added,i;
   
    
    actually_added=-1;
    while(actually_added<0){
        for(i=0;i<gg.get_dim();i++){
            rr.set(i,normal_deviate(dice,0.0,1.0));
        }
        rr.normalize();
        for(i=0;i<gg.get_dim();i++){
            trial.set(i,gg.get_pt(nodes(ic)->get_center(),i)+0.1*rr.get_data(i)*(gg.get_max(i)-gg.get_min(i)));
        }
        
        ggWrap.evaluate(trial,&chitrue,&actually_added);
        
    }
    focus_pts.add(actually_added);
    if(do_bisection==1){
        bisection(trial,chitrue);
    }
    
    //printf("found %e\n",chitrue);
    
}

void aps::corner_focus(int ic){
    /*
    ic identifies the center around whic we are focusing
    
    center_dexes.get_data(ic) is the index of the center as stored in the Gaussian Process
    */

    array_1d<double> min,max,trial,sambest,rr_perp,rr,origin;
    array_1d<int> min_dex,max_dex,*extremity;
    double nn,chitrue,stradval,stradmax,mu,sig,mu_chosen,sig_chosen,norm,norm_chosen;
    int i,j,actually_added;
    int ix,idx,ix_chosen,dx_chosen;
    int iy,idy;
    int ict,jct;
    
    min_dex.set_name("corner_focus_min_dex");
    max_dex.set_name("corner_focus_max_dex");
    
    /*
    Find the boudary points that minimize and maximize each parameter
    */
    for(i=0;i<nodes(ic)->get_n_boundary();i++){
        for(j=0;j<gg.get_dim();j++){
            nn=gg.get_pt(nodes(ic)->get_boundary(i),j);
            if(i==0 || nn<min.get_data(j)){
                min.set(j,nn);
                min_dex.set(j,nodes(ic)->get_boundary(i));
            }
            
            if(i==0 || nn>max.get_data(j)){
                max.set(j,nn);
                max_dex.set(j,nodes(ic)->get_boundary(i));
            } 
        }
    }
    
    for(i=0;i<gg.get_dim();i++){
        origin.set(i,gg.get_pt(nodes(ic)->get_center(),i));
    }
    
    for(i=0;i<gg.get_dim();i++){
        while(fabs(max.get_data(i)-min.get_data(i))<1.0e-10){
            max.add_val(i,1.0e-4*(gg.get_max(i)-gg.get_min(i)));
            min.subtract_val(i,1.0e-4*(gg.get_max(i)-gg.get_min(i)));
        }
    }
    
    array_1d<double> length;
    
    for(i=0;i<gg.get_dim();i++)length.set(i,max.get_data(i)-min.get_data(i));
    
    stradmax=-2.0*chisq_exception;
    
    /*idx=0 means we are trying to expand on the minimum points
    idx=1 means we are trying to expand on the maximum points*/
    for(idx=0;idx<2;idx++){
        if(idx==0)extremity=&min_dex;
        else if(idx==1) extremity=&max_dex;
        
        for(ix=0;ix<gg.get_dim();ix++){
            
            /*get the vector pointing from the center to the boundary point in question*/
            for(i=0;i<gg.get_dim();i++){
                rr.set(i,gg.get_pt(extremity->get_data(ix),i)-origin.get_data(i));
            }
            
            /*normalize that vector by the size of the low-chisquared region*/
            nn=0.0;
            for(i=0;i<gg.get_dim();i++){
                nn+=rr.get_data(i)*rr.get_data(i)/power(length.get_data(i),2);
            }
            
            if(nn<1.0e-30 || isnan(nn)){
                printf("WARNNING norm of rr %e\n",nn);
                printf("in corner_focus; before ict loop\n");
                printf("ix %d idx %d\n",ix,idx);
                exit(1);
            }
            
            nn=sqrt(nn);
            for(i=0;i<gg.get_dim();i++){
                rr.divide_val(i,nn);
            }
            
           
            /*propose 20 points that are displaced from the boundary point
            in a direction perpendicular to rr*/
            for(ict=0;ict<20;ict++){
                
                for(i=0;i<gg.get_dim();i++){
                    rr_perp.set(i,normal_deviate(dice,0.0,length.get_data(i)));
                }
                    
                rr_perp.set(ix,0.0);
                
                nn=0.0;
                for(i=0;i<gg.get_dim();i++){
                    nn+=rr_perp.get_data(i)*rr.get_data(i)/power(length.get_data(i),2);
                }
                for(i=0;i<gg.get_dim();i++){
                    rr_perp.subtract_val(i,nn*rr.get_data(i));
                }
                    
                nn=0.0;
                for(i=0;i<gg.get_dim();i++){
                    nn+=rr_perp.get_data(i)*rr_perp.get_data(i)/power(length.get_data(i),2);
                }
                
                if(nn<1.0e-10 || isnan(nn)){
                    printf("WARNING norm of rr_perp is %e\n",nn);
                    printf("in corner focus; inside ict loop; outside jct loop\n");
                    printf("ix %d idx %d\n",ix,idx);
                    exit(1);
                }
                
                nn=sqrt(nn);
                for(i=0;i<gg.get_dim();i++){
                    rr_perp.divide_val(i,nn);
                }
                
                /*try several different distances long the perpendicular displacement;
                with each successive step, step farther away from the boundary point.
                */
                norm=0.1;
                for(jct=0;jct<5;jct++){
                
                    for(i=0;i<gg.get_dim();i++){
                        trial.set(i,gg.get_pt(extremity->get_data(ix),i)+norm*rr_perp.get_data(i));
                    }
                    
                    
                    /*set the dimension whose extremum we are expanding from
                    to remain unchanged; in theory, this means we should be stepping
                    along the surface of the low-chisquared region (that is almost
                    certainly not what is happening in practice; that is why
                    we will follow up with a bisection search)*/
                    if(idx==0){
                        trial.subtract_val(ix,0.1*length.get_data(ix));
                    }
                    else if(idx==1){
                        trial.add_val(ix,0.1*length.get_data(ix));
                    }
                     
                    if(ggWrap.is_valid(trial)==0){
                        stradval=-2.0*chisq_exception;
                    }
                    else{
                        mu=ggWrap.user_predict(trial,&sig);
                        stradval=strad(mu,sig);
                    }
                
                    /*
                    Of these proposed candidate points, select the one which
                    maximizes the S statistic from equation (1) of the paper
                    to actually evaluate chisquared
                    */
                    if(stradval>stradmax){
                        stradmax=stradval;
                
                        for(i=0;i<gg.get_dim();i++){
                            sambest.set(i,trial.get_data(i));
                        }
                    
                        ix_chosen=ix;
                        dx_chosen=idx;
                        mu_chosen=mu;
                        sig_chosen=sig;
                        norm_chosen=norm;
                    }
                    
                    norm+=0.1;
                    
                }//loop over jct
            }//loop over ict (the number of trials proposed from each boundary point)
            
        }//loop over ix (which is the dimension)
        
    }//loop over idx which controls whether this is max or min
    
    
    
    if(stradmax>-1.0*chisq_exception){
        /*
        Evaluate chisquared at the chosen point.
        Perform bisection to make certain that a new boundary point is found.
        */
        ggWrap.evaluate(sambest,&chitrue,&actually_added,1);
        if(actually_added>=0){
            focus_pts.add(actually_added);
            if(do_bisection==1){
                bisection(sambest,chitrue);
            }
        }
    }
    else{
        
        /*in the event that a suitable candidate was not found,
        randomly select a point to sample that is near the center in question*/
        
        random_focus(ic);
    }
}

void aps::aps_focus(){

   int ic,imin=-1;
   double ttmin;
   
   for(ic=0;ic<nodes.get_dim();ic++){
       if(gg.get_fn(nodes(ic)->get_center())<strad.get_target() &&
           nodes(ic)->is_it_active()==1 && 
           (imin==-1 || nodes(ic)->get_time()<ttmin)){
       
           imin=ic;
           ttmin=nodes(ic)->get_time();
       }
   }
   
   nodes(imin)->search();
   
}

int aps::find_nearest_center(array_1d<double> &pt){
    return find_nearest_center(pt,2.0*chisq_exception);
}

int aps::find_nearest_center(array_1d<double> &pt, double chi_in){
    /*
    Find the nearest center with chisq<chi_in
    */
    
    int i,imin,ans;
    double dd,ddmin;
    ans=-1;
    
    for(i=0;i<nodes.get_dim();i++){
        dd=gg.distance(pt,nodes(i)->get_center());
        if((ans<0 || dd<ddmin) && gg.get_fn(nodes(i)->get_center())<chi_in){
            ddmin=dd;
            ans=i;
        }
    }
    
    return ans;
    
}

void aps::project_to_unit_sphere(int ic, array_1d<double> &pt_in, array_1d<double> &pt_out){

    if(ic<0 || ic>=nodes.get_dim()){
        return;
    }
    
    array_1d<double> dir;
    double norm=0.0;
    int i;
    for(i=0;i<gg.get_dim();i++){
        dir.set(i,pt_in.get_data(i)-gg.get_pt(nodes(ic)->get_center(),i));
        
        /*
        gg.get_max-gg.get_min will return the characteristic length if it was set,
        and the range_max-range_min if not.  That is why we use that here
        */
        norm+=power(dir.get_data(i)/(gg.get_max(i)-gg.get_min(i)),2);
    }
    norm=sqrt(norm);
    for(i=0;i<gg.get_dim();i++){
        dir.divide_val(i,norm);
    } 
    
    for(i=0;i<gg.get_dim();i++){
        pt_out.set(i,gg.get_pt(nodes(ic)->get_center(),i)+dir.get_data(i));
    }
    
}

void aps::bisection(array_1d<double> &inpt, double chi_in){
    
    ggWrap.set_iWhere(iBisect);
    
    /*
    If the point is already fairly near to chisquared_lim, or if APS has not yet found
    a low-chisquared region for which chisquared<chisquared_lim, do not do bisection
    */
    if(chi_in<strad.get_target() && 
    strad.get_target()-chi_in<0.1*ggWrap.get_delta_chisquared() || ggWrap.get_chimin()>strad.get_target()){
        return;
    }
    
    array_1d<double> dir_origin,trial,ddneigh;
    array_1d<int> neigh;
    
    double mu,fdir_origin;
    
    double dd,ddmin;
    int origin_dex,i,j,k,use_it_parabola,i_center=-1;
    
    double bisection_tolerance=0.01*ggWrap.get_delta_chisquared();
    
    /*
    The code will work by finding one point with chisquared<chisquared_lim which will
    be stored as the array_1d<double> lowball (and corresponding chisquared value flow)
    and one point with chisquared>chisquared_lim which will be stored as the 
    array_1d<double> highball (and corresponding chisquared value fhigh).
    The code will iteratively step between lowball and highball until it finds a point that is within
    bisection_tolerance of chisquared_lim.
    
    The code will also store an array_1d<double> dir_origin which is the low-chisquared center
    from which the bisection effectively began.  This is for purposes of code at the end of this 
    routine (which is presently commented-out) which had bisection follow up the iterative search
    with one step intentionally just outside of the chisquared=chisquared_lim contour and one
    step intentionally just inside of that contour.  It was hoped that this would give us a better
    characterization of the chisquared function.  The benefit was difficult to prove.
    */
    
   /*
   Finding dir_origin and f_dir_origin;  this point will also be used for lowball,
   if chi_in>chisquared_lim
   
   i_center will be the index (as stored in array_1d<int> center_dexes or array_2d<double> centers
   of the nearest low-chisquared center point) 
   */
    if(ggWrap.get_ngood()==0){
        for(i=0;i<gg.get_dim();i++)dir_origin.set(i,ggWrap.get_minpt(i));
        fdir_origin=ggWrap.get_chimin();
        origin_dex=ggWrap.get_global_mindex();
    }
    else{
        ddmin=chisq_exception;
        
        i_center=find_nearest_center(inpt,chi_in);
        
        if(i_center>=0){
            j=nodes(i_center)->get_center();
            for(i=0;i<gg.get_dim();i++){
                dir_origin.set(i,gg.get_pt(j,i));
            }
            fdir_origin=gg.get_fn(j);
            origin_dex=j;
            
        }
        else{
            fdir_origin=ggWrap.get_chimin();
            origin_dex=ggWrap.get_global_mindex();
            for(i=0;i<gg.get_dim();i++)dir_origin.set(i,ggWrap.get_minpt(i));
        }
    }
    
    array_1d<double> lowball,highball,dir;
    double flow,fhigh,fnearest,rr,new_rr;
    int i_test,i_high=-1,i_low=-1,i_nearest=-1;
    int old_high=-1,old_low=-1,old_nearest=-1;
    
    array_1d<double> mu_to_sort;
    array_1d<int> mu_dex;
    
    /*
    Initially set lowball and highball
    */
    if(chi_in>strad.get_target()){
        for(i=0;i<gg.get_dim();i++){
            lowball.set(i,dir_origin.get_data(i));
            highball.set(i,inpt.get_data(i));
        }
        flow=fdir_origin;
        fhigh=chi_in;
    }
    else{
        /*
        If chi_in<chisquared_lim, then find the direction connecting
        inpt to dir_origin.  Step along this point in larger and larger
        steps until you find a point with chisquared>chisquared_lim.
        Set that to highball and begin bisection.
        */
        for(i=0;i<gg.get_dim();i++){
            dir.set(i,inpt.get_data(i)-dir_origin.get_data(i));
            lowball.set(i,inpt.get_data(i));
        }
        flow=chi_in;
        
        rr=2.0*dir.normalize();
        while(rr<1.0e-20){
            /*
            In case, for some reason, in_pt is dir_origin
            (that should not happen, but one cannot be too careful)
            */
            for(i=0;i<gg.get_dim();i++){
                dir.set(i,dice->doub());
            }
            
            rr=2.0*dir.normalize();
        }
        
        fhigh=-1.0*chisq_exception;
        while(fhigh<strad.get_target()){
            for(i=0;i<gg.get_dim();i++){
                highball.set(i,lowball.get_data(i)+rr*dir.get_data(i));
            }
            
            
            ggWrap.evaluate(highball,&fhigh);
 
            rr*=2.0;
        }
    
    }
    
    /*
    nearest_pt will keep track of the point that is actually nearest to
    chisquared_lim
    
    i_nearest will be the index of that point as stored in the Gaussian Process
    */
    array_1d<double> nearest_pt;
        
    if(strad.get_target()-flow<fhigh-strad.get_target()){
        fnearest=flow;
        for(i=0;i<gg.get_dim();i++){
            nearest_pt.set(i,lowball.get_data(i));
        }
    
    }
    else{
        fnearest=fhigh;
        for(i=0;i<gg.get_dim();i++){
            nearest_pt.set(i,highball.get_data(i));
        }
    }
    
    if(flow>strad.get_target() || fhigh<strad.get_target()){
        /*
        Something went horribly wrong; stop now
        */
        return;
    }
        
    dd=gg.distance(lowball,highball);
    while(dd>1.0e-10 && strad.get_target()-flow>bisection_tolerance){
        for(i=0;i<gg.get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        
        ggWrap.evaluate(trial,&mu,&i_test);
        
        if(mu>strad.get_target()){
            for(i=0;i<gg.get_dim();i++)highball.set(i,trial.get_data(i));
            fhigh=mu;
            old_high=i_high;
            i_high=i_test;
            if(fhigh-strad.get_target() < fabs(fnearest-strad.get_target())){
                fnearest=fhigh;
                old_nearest=i_nearest;
                i_nearest=i_high; 
                for(i=0;i<gg.get_dim();i++){
                    nearest_pt.set(i,highball.get_data(i));
                }
            }
                
        }
        else{
            for(i=0;i<gg.get_dim();i++)lowball.set(i,trial.get_data(i));
            flow=mu; 
            old_low=i_low;
            i_low=i_test;
            if(strad.get_target()-flow < fabs(fnearest-strad.get_target())){
                fnearest=flow;
                old_nearest=i_nearest;
                i_nearest=i_low;  
                for(i=0;i<gg.get_dim();i++){
                    nearest_pt.set(i,lowball.get_data(i));
                }
            }

        }
                        
        dd*=0.5;
        
    }
   
   
    /*
    This is unused code which had bisection sample points that were intentionally
    just inside and just outside of the bound after finding the chisquared=chisquared_lim
    point.  The thought was that this would be useful in constructing Bayesian credible
    limits from APS outputs.  It is not clear that this was helpful, but we are leaving
    the code here for future investigation.
    */
    /*    
    for(i=0;i<dim;i++){
        dir.set(i,nearest_pt.get_data(i)-dir_origin.get_data(i));
    }
    dd=dir.normalize();
    
    for(i=0;i<gg.get_dim();i++){
        trial.set(i,dir_origin.get_data(i)+1.5*dd*dir.get_data(i));
    }
    
    evaluate(trial,&mu,&i_test);
    
    for(i=0;i<gg.get_dim();i++){
        trial.set(i,dir_origin.get_data(i)+0.5*dd*dir.get_data(i));
    }
      
    evaluate(trial,&mu,&i_test);
    */
    
    
    /*add the discovered point to the list of boundary points for the center
    from which we began our bisection search*/
    if(i_center>=0 && i_nearest>=0){
        nodes(i_center)->add_as_boundary(i_nearest); 
    }

    
}

void aps::aps_search(){

    double before=double(time(NULL));
    int ibefore=ggWrap.get_called(),i_wide;
    
    int do_focus=0,i;
    for(i=0;i<nodes.get_dim() && do_focus==0;i++){
        if(gg.get_fn(nodes(i)->get_center())<strad.get_target() && nodes(i)->is_it_active()==1){
            do_focus=1;
        }
    }
    
    int active_nodes=0;
    for(i=0;i<nodes.get_dim();i++){
        if(nodes(i)->is_it_active()==1)active_nodes++;
    }
    
    if(do_focus==1 && called_focus/active_nodes<called_wide){
        aps_focus();
        called_focus+=ggWrap.get_called()-ibefore;
    }
    else{
        i_wide=aps_wide();
        called_wide+=ggWrap.get_called()-ibefore;
    }

    time_aps+=double(time(NULL))-before;
    ct_aps+=ggWrap.get_called()-ibefore;
    set_where("nowhere");
    
}

double aps::distance(int i1, int i2, array_1d<double> &range){
    double dd=0.0;
    int i;
    for(i=0;i<dim;i++){
        if(range.get_data(i)>0.0){
            dd+=power((gg.get_pt(i1,i)-gg.get_pt(i2,i))/range.get_data(i),2);
        }
    }
    
    return sqrt(dd);
}

void aps::simplex_too_few_candidates(array_1d<int> &candidates){
    
    int i,actually_added,i_min;
    double ccmin,cc;
    array_1d<double> trial,step,deviation;
    
    trial.set_name("simplex_too_few_trial");
    step.set_name("simplex_too_few_step");
    deviation.set_name("simplex_too_few_deviation");
    
    for(i=0;i<candidates.get_dim();i++){
        if(i==0 || gg.get_fn(candidates.get_data(i))<ccmin){
            ccmin=gg.get_fn(candidates.get_data(i));
            i_min=candidates.get_data(i);
        }
    }
    
    int ic;
    double stepNorm,devNorm;
    ic=find_nearest_center(*gg.get_pt(i_min));
    
    //printf("    minpt %e\n",gg.get_pt(i_min,0));
    
    stepNorm=0.0;
    for(i=0;i<gg.get_dim();i++){
        step.set(i,gg.get_pt(i_min,i)-gg.get_pt(nodes(ic)->get_center(),i));
        stepNorm+=power(step.get_data(i)/(gg.get_max(i)-gg.get_min(i)),2);
    }
    stepNorm=sqrt(stepNorm);
    for(i=0;i<gg.get_dim();i++){
        step.divide_val(i,stepNorm);
    }
    
    int ct=0;
    double trialNorm;
    while(candidates.get_dim()<gg.get_dim()+1 && ct<(gg.get_dim()+1)*5){
        ct++;
        
        stepNorm=0.0;
        for(i=0;i<gg.get_dim();i++){
            step.set(i,normal_deviate(dice,0.0,gg.get_max(i)-gg.get_min(i)));
            stepNorm+=power(step.get_data(i)/(gg.get_max(i)-gg.get_min(i)),2);
        }
        stepNorm=sqrt(stepNorm);
        
        for(i=0;i<gg.get_dim();i++){
            step.divide_val(i,stepNorm);
        } 
        for(i=0;i<gg.get_dim();i++){
            trial.set(i,gg.get_pt(i_min,i)+0.01*step.get_data(i));
        }
        
        
        ggWrap.evaluate(trial,&cc,&actually_added);
        
        if(actually_added>=0 && is_it_a_candidate(actually_added)){
            candidates.add(actually_added);
        }
        
        
    }
    
    //printf("    after adding have %d candidates\n",candidates.get_dim());
    if(candidates.get_dim()==gg.get_dim()+1){

        find_global_minimum(candidates);
    }
    
}

void aps::refine_center(){
    /*
    Evaluate the known centers of low chisquared regions to see if any of them
    can be improved upon
    */
    
    int ic,ic_chosen=-1;
    double maxchi;
    
    /*choose the center with the highest value of chisquared*/
    for(ic=0;ic<nodes.get_dim();ic++){
        
        if(nodes(ic)->get_n_boundary()>=gg.get_dim()+1 && 
        (ic_chosen<0 || gg.get_fn(nodes(ic)->get_center())>maxchi)){
            
            maxchi=gg.get_fn(nodes(ic)->get_center());
            ic_chosen=ic;
        
        }
        
    }
    
    if(ic_chosen<0){
        return;
    }
    
    array_1d<double> tosort,sorted;
    array_1d<int> inn;
    
    tosort.set_name("refine_tosort");
    sorted.set_name("refine_sorted");
    inn.set_name("refine_inn");
    
    int ii,i;
    
    /*
    sort the boundary points of the chosen center by their distance
    from the chosen center
    */
    for(i=0;i<nodes(ic_chosen)->get_n_boundary();i++){
        ii=nodes(ic_chosen)->get_boundary(i);
        
        inn.set(i,ii);
        tosort.set(i,gg.distance(ii,nodes(ic_chosen)->get_center()));
    
    }
    sort_and_check(tosort,sorted,inn);
    
    array_1d<int> neigh;
    neigh.set_name("refine_neigh");
    
    /*
    Set the simplex to the dim+1 farthest boundary points from the center
    */
    for(i=0;i<gg.get_dim()+1;i++){
        neigh.set(i,inn.get_data(inn.get_dim()-1-i));
    }
    
    /*
    If this simplex is identical to the simplex chosen the last time
    this center was used in refine_center(), do not proceed with 
    the simplex search
    */
    if(refined_simplex.get_rows()>ic_chosen){
        if(compare_int_arr(*refined_simplex(ic_chosen),neigh)==1){
            return;
        }
    }
    
    /*
    Make note that we are starting a simplex search with this seed simplex
    */
    refined_simplex.set_row(ic_chosen,neigh);
    
    /*
    Perform the simplex searc
    */
    find_global_minimum(neigh);
    

}

void aps::calculate_gradient(int center, array_1d<int> &neigh, 
         array_1d<double> &vout){

    array_1d<double> aa,bb,xx;
    
    int i,j;
    
    if(neigh.get_dim()!=dim){
        printf("CANNOT calculate gradient; only %d neighbors\n",
        neigh.get_dim());
        return;
    }
    
    for(i=0;i<gg.get_dim();i++){
        bb.set(i,gg.get_fn(neigh.get_data(i))-gg.get_fn(center));
        for(j=0;j<gg.get_dim();j++){
            aa.set(i*gg.get_dim()+j,(gg.get_pt(i,j)-gg.get_pt(center,j))/(gg.get_max(j)-gg.get_min(j)));
        }
    }
    
    try{
        naive_gaussian_solver(aa,bb,xx,gg.get_dim());
        
        for(i=0;i<gg.get_dim();i++){
            vout.set(i,xx.get_data(i)/(gg.get_max(i)-gg.get_min(i)));
        }
        
    }
    catch(int iex){
        printf("could not get gradient this time\n");
    }
}

void aps::optimize(){
    
    double before=double(time(NULL));
    
    if(wide_pts.get_dim()<=0)return;
    
    gp gg_opt;
    
    gg_opt.set_kk(gg.get_kk());
    array_1d<double> ggmin,ggmax;
    
    int i,j;
    
    for(i=0;i<dim;i++){
        ggmin.set(i,gg.get_min(i));
        ggmax.set(i,gg.get_max(i));
    }
    
    array_1d<int> use_dex;
    array_1d<double> ff_opt;
    array_2d<double> data_opt;
    double rat,roll;
    
    data_opt.set_cols(dim);
    rat=3000.0/double(wide_pts.get_dim());
    for(i=0;i<wide_pts.get_dim();i++){
        roll=dice->doub();
        if(roll<rat){
           
            ff_opt.add(gg.get_fn(wide_pts.get_data(i)));
            data_opt.add_row(*gg.get_pt(wide_pts.get_data(i)));
       }     
      
    }
    
    gg_opt.assign_covariogram(gg.get_covariogram());
    gg_opt.initialize(data_opt,ff_opt,ggmin,ggmax);
    
    use_dex.reset();
    for(i=0;i<data_opt.get_rows();i++)use_dex.set(i,i);
  
    gg_opt.optimize(use_dex);
    
    
    array_1d<double> hh;
    
    gg_opt.get_hyper_parameters(hh);
    set_hyper_parameters(hh);
    
    _aps_wide_contents_buffer.reset();
    _aps_wide_ss_buffer.reset();
    
    last_optimized=wide_pts.get_dim();
    time_optimizing+=double(time(NULL))-before;

    
}

void aps::write_pts(){
    
    set_where("write_pts");
    
    array_1d<double> hyper_params;
    double before=double(time(NULL));
    
    int i,j,k,ii,jj,lling,aps_dex;
    double mu,sig;
    FILE *output;

    double per_chisq,per_total,overhead;
    
    array_1d<int> candidates;
    candidates.set_name("aps_write_pts_candidates");
    
    for(i=0;i<nodes.get_dim();i++){
        nodes(i)->flush_candidates(candidates);
        k=nodes.get_dim();
        for(j=0;j<candidates.get_dim();j++){
            assess_node(candidates.get_data(j));
        }
        
        if(nodes.get_dim()>k){
            for(j=k;j<nodes.get_dim();j++){
                nodes(j)->set_time(nodes(i)->get_time());
                for(ii=0;ii<ggWrap.get_dim();ii++){
                    for(jj=0;jj<ggWrap.get_dim();jj++){
                        nodes(j)->set_basis(ii,jj,nodes(i)->get_basis(ii,jj));
                    }
                }
            }
        }
        
        candidates.reset();
    }
    
    /*how much clock time is spent per call to chisquared just calling chisquared*/
    per_chisq=ggWrap.get_chisq_time()/ggWrap.get_called();
    
    /*how much clock time is spent per call to chisquared on all of the calculations*/
    per_total=(double(time(NULL))-start_time)/ggWrap.get_called();
    
    /*how much extra clock time has been added to each chisquared call
    by all of the extra calculations involved in APS*/
    overhead=per_total-per_chisq;
    
    
    /*
    If we have not optimized hyper parameters yet, or if 100 aps_wide() points have
    been sampled since the last time hyper parameters were optimized AND the overhead
    in clock time per chisquared is less than 1/10 the normal time spent on a chisquared
    evaluation, optimize the Gaussian Process hyper parameters
    */
    if(last_optimized==0 ||
      (wide_pts.get_dim()>last_optimized+100 && 
      overhead<0.1*per_chisq)){
    
        optimize();
    }
    
    /*
    Figure out which points satisfy chisquared<=chisquared_lim
    (since chisquared_lim may have changed because of the discovery of a new
    chisquared_min)
    */
    ggWrap.evaluate_ngood();
    
    
    array_1d<double> volume,vmax,vmin;
    double vol,lvol;
    
    volume.set_name("write_volume");
    vmax.set_name("write_vmax");
    vmin.set_name("write_vmin");
    
    /*
    Find the parameter space volumes of each independent low-chisquared region
    */
    for(i=0;i<nodes.get_dim();i++){
        vmax.reset();
        vmin.reset();
        for(j=0;j<nodes(i)->get_n_boundary();j++){
            for(k=0;k<gg.get_dim();k++){
                if(j==0 || gg.get_pt(nodes(i)->get_boundary(j),k)<vmin.get_data(k)){
                    vmin.set(k,gg.get_pt(nodes(i)->get_boundary(j),k));
                }
                
                if(j==0 || gg.get_pt(nodes(i)->get_boundary(j),k)>vmax.get_data(k)){
                    vmax.set(k,gg.get_pt(nodes(i)->get_boundary(j),k));
                }
            }
        }
        
        if(vmax.get_dim()!=gg.get_dim() || vmin.get_dim()!=gg.get_dim()){
            vol=0.0;
        }
        else{
            lvol=0.0;
            for(j=0;j<gg.get_dim();j++){
                if(vmax.get_data(j)-vmin.get_data(j)>1.0e-20){
                    lvol+=log(vmax.get_data(j)-vmin.get_data(j));
                }
                else lvol-=1.0e10;
            }
            vol=exp(lvol);
        }
        volume.set(i,vol);
        
    }
    
    /*
    Print all of the points sampled by APS
    */
    output=fopen(outname,"w");
    fprintf(output,"# ");
    for(i=0;i<gg.get_dim();i++){
        fprintf(output,"%s ",paramnames[i]);
    }
    fprintf(output,"chisq mu sig ling\n");
    
    int focus_dex=0;
    aps_dex=0;
    for(i=0;i<gg.get_pts();i++){
        lling=ggWrap.get_whereFrom(i);
        if(aps_dex<wide_pts.get_dim() && i==wide_pts.get_data(aps_dex)){

            mu=mu_storage.get_data(aps_dex);
            sig=sig_storage.get_data(aps_dex);
            aps_dex++;
        }
        else if(focus_dex<focus_pts.get_dim() && i==focus_pts.get_data(focus_dex)){

            mu=-2.0;
            sig=-2.0;
            focus_dex++;
        }
        else{

            mu=-2.0;
            sig=-2.0;
        }
    
    
        for(j=0;j<gg.get_dim();j++){
            fprintf(output,"%.18e ",gg.get_pt(i,j));
        }
        fprintf(output,"%.18e %.18e %.18e %d\n",gg.get_fn(i),mu,sig,lling);
    }
    fclose(output);
   
   
    array_1d<double> tosort,sorted;
    tosort.set_name("aps_write_pts_tosort");
    sorted.set_name("aps_write_pts_sorted");
    
    array_1d<int> inn;
    inn.set_name("aps_write_pts_inn");
    
    
    /*
    Find the 1/10 quantile of chisquared among points discovered by aps_wide()
    */
    for(i=0;i<wide_pts.get_dim();i++){
        tosort.set(i,gg.get_fn(wide_pts.get_data(i)));
        inn.set(i,wide_pts.get_data(i));
    }
    
    sort_and_check(tosort,sorted,inn);
    global_threshold=sorted.get_data(tosort.get_dim()/10);
    sorted.reset();
    inn.reset();
    
    if(ddUnitSpheres.get_dim()>0){
        /*
        Find the 2/3 quantile of nearest neighbor distances among the spherical projection of boundary
        points and aps_wide points.  Consider only points discovered since the last call to write_pts
        */
        for(i=0;i<ddUnitSpheres.get_dim();i++)inn.set(i,i);
        sort_and_check(ddUnitSpheres,sorted,inn);
        sphere_threshold=sorted.get_data((2*tosort.get_dim())/3);
        
        ddUnitSpheres.reset();
    }
    
    int ic;
   
    n_printed=gg.get_pts();

    double nn;
    
    i=gg.get_last_refactored()/2;
    if(i<1000)i=1000;
    
    /*
    Consider refactoring the kd_tree in the Gaussian Process object gg.
    
    The hope is that this will make the nearest neighbor search more efficient.
    */
    if(gg.get_pts()>gg.get_last_refactored()+i && gg.get_pts()<20000){
        printf("refactoring\n");
    
        nn=double(time(NULL));
        gg.refactor();
        
        if(ggWrap.get_nboxes()>305){
            printf("%e %e\n",ggWrap.get_box_min(305,3),ggWrap.get_box_max(305,3));
            if(ggWrap.get_box_max(305,3)-ggWrap.get_box_min(305,3)<1.0e-3){
               
                exit(1);
            }
        }
        
        
        _aps_wide_contents_buffer.reset();
        _aps_wide_ss_buffer.reset();
        time_refactoring+=double(time(NULL))-nn;
    }

   
    /*
    Output the timing statistics file
    */
    double time_now = double(time(NULL));
   
    output=fopen(timingname,"a");
    fprintf(output,"%d %d %e %e %e %e -- ",
    gg.get_pts(),ggWrap.get_called(),ggWrap.get_chisq_time(),
    ggWrap.get_chisq_time()/double(ggWrap.get_called()),time_now-start_time,
    (time_now-start_time)/double(ggWrap.get_called()));
    
    fprintf(output,"%d %e -- ",ct_aps,time_aps);
    fprintf(output,"%d %e -- ",ct_simplex,time_simplex);
    
    fprintf(output,"%e -- %e -- ",time_optimizing,time_refactoring);
    
    fprintf(output,"%e %e -- ",ggWrap.get_chimin(),strad.get_target());
    
    for(i=0;i<nodes.get_dim();i++){
        fprintf(output,"%e ",volume.get_data(i));
    }
    
    fprintf(output," -- %d %d ",called_wide,called_focus);
    if(ggWrap.is_unitSpheres_null()==0){
        fprintf(output,"unitSpheres engaged %d ",ggWrap.get_unitSpheres_pts());
    }
    fprintf(output,"-- %d %d",n_simplex,n_simplex_found);
    fprintf(output,"\n");
    
    fclose(output);
    
    printf("\n");
    printf("APS: %d\n",ggWrap.get_whereCt(iAPS));
    printf("Simplex: %d\n",ggWrap.get_whereCt(iSimplex));
    printf("Bisection: %d\n",ggWrap.get_whereCt(iBisect));
    printf("node Bisection: %d\n",ggWrap.get_whereCt(iNodeBisect));
    printf("Coulomb: %d\n",ggWrap.get_whereCt(iCoulomb));
    printf("Compass: %d\n",ggWrap.get_whereCt(iCompass));
    printf("Ricochet: %d\n",ggWrap.get_whereCt(iRicochet));
    printf("nodes: %d\n",nodes.get_dim());
    for(i=0;i<nodes.get_dim();i++){
        printf("node %d associates %d vv %.4e active %d -- bases %d %d\n",
        nodes(i)->get_center(),nodes(i)->get_n_associates(),
        nodes(i)->volume(),nodes(i)->is_it_active(),
        nodes(i)->get_calls_to_bases(),
        nodes(i)->get_ct_bases());
    }
    printf("called time_tot coulomb ricochet bases\n");
    for(i=0;i<nodes.get_dim();i++){
        printf("%d %e %e %e %e\n",nodes(i)->get_ct_search(),
        nodes(i)->get_time(),nodes(i)->get_time_coulomb(),nodes(i)->get_time_ricochet(),
        nodes(i)->get_time_bases());
    }
    printf("search time %e search_ct %d -- %e\n",
    ggWrap.get_search_time(),ggWrap.get_search_ct(),
    ggWrap.get_search_time()/double(ggWrap.get_search_ct()));
    
    printf("search time solo %e ct %d -- %e\n",
    ggWrap.get_search_time_solo(),ggWrap.get_search_ct_solo(),
    ggWrap.get_search_time_solo()/double(ggWrap.get_search_ct_solo()));
    
    array_1d<int> quartiles;
    quartiles.set_name("aps_quartiles");
    ggWrap.get_box_quartiles(quartiles);
    
    printf("search time box %e ct %d -- %e\n",
    ggWrap.get_search_time_box(),ggWrap.get_search_ct_box(),
    ggWrap.get_search_time_box()/double(ggWrap.get_search_ct_box()));
    printf("smallest box %d\n",ggWrap.get_smallest_box());
    printf("biggest box %d\n",ggWrap.get_biggest_box());
    printf("biggest bad box %d\n",ggWrap.get_biggest_bad_box());
    printf("n boxes %d\n",ggWrap.get_nboxes());
    printf("n small boxes %d\n",ggWrap.get_n_small_boxes());
    printf("n optimal boxes %d\n",ggWrap.get_n_optimal_boxes());
    printf("box quartiles %d %d %d\n",
    quartiles.get_data(0),quartiles.get_data(1),quartiles.get_data(2));
    printf("\ntotal time %e\n",time_now-start_time);
    printf("pts %d\n",ggWrap.get_pts());
    printf("n wide %d n box wide %d  found %d\n",n_wide,n_simplex,n_simplex_found);
    printf("\n");
    
    output=fopen("node_dump.sav","w");
    for(i=0;i<nodes.get_dim();i++){
        for(j=0;j<ggWrap.get_dim();j++){
           for(k=0;k<ggWrap.get_dim();k++){
               fprintf(output,"%le ",nodes(i)->get_basis(k,j));
           }
           fprintf(output,"\n");
        }
        fprintf(output,"\n\n");
    }
    fclose(output);
       
    set_where("nowhere");
    time_writing+=double(time(NULL))-before;
    
}

////////////////////node-based code below/////////

void aps::assess_node(int dex){
    if(dex<0 || dex>=gg.get_pts()){
        return;
    }
    
    if(gg.get_fn(dex)>strad.get_target()){
        return;
    }
    
    int i;
    for(i=0;i<nodes.get_dim();i++){
        if(dex==nodes(i)->get_center()){
            return;
        }
    }
    
    int used_because_distance,use_it;
    
    /*used_because_distance allows us to accept node candidates based on their
    distance from the nearest node*/
    used_because_distance=0;
    
    if(ddNodeRatio>1.0e-10){
        for(i=0;i<nodes.get_dim() && used_because_distance==0;i++){
            if(nodes(i)->get_farthest_associate()>0.0 && 
               nodes(i)->get_n_associates()>100){
           
                   used_because_distance=1;   
            }
        }
    }
    
    array_1d<double> midpt;
    midpt.set_name("aps_assess_node_midpt");
    double ftrial,dd,ddmin,wgt;
    int j,itrial,inode,iclosest=-1;
    int iGhost,use_ghost;
    
    use_it=1;
    for(i=0;i<nodes.get_dim() && (use_it==1 || used_because_distance==1);i++){
        inode=nodes(i)->get_center();
        dd=gg.distance(inode,dex);
        if(i==0 || dd<ddmin){
            ddmin=dd;
            iclosest=i;
        }
        
        if(dd<ddNodeRatio*nodes(i)->get_farthest_associate() && ddNodeRatio>1.0e-10){
            /*the point under consideration is too close to an existent node to be
            set based on distance*/
            used_because_distance=0;
        }
        
        //make the same test on old centers of nodes
        for(iGhost=0;iGhost<nodes(i)->get_n_oldCenters();iGhost++){
            dd=gg.distance(nodes(i)->get_oldCenter(iGhost),dex);
            if(dd<ddNodeRatio*nodes(i)->get_farthest_associate() && ddNodeRatio>1.0e-10){
                used_because_distance=0;
            }
        }
        
        use_it=0;
        for(wgt=0.25;wgt<0.8 && use_it==0;wgt+=0.25){
        
            for(j=0;j<gg.get_dim();j++){
                midpt.set(j,wgt*gg.get_pt(inode,j)+(1.0-wgt)*gg.get_pt(dex,j));
            }
        
            ggWrap.evaluate(midpt,&ftrial);
        
            if(ftrial>strad.get_target()){
                use_it=1;
            }
        }
        
        //make the same test on old centers of nodes
        for(iGhost=0;iGhost<nodes(i)->get_n_oldCenters() && use_it==1;iGhost++){
            use_ghost=0;
            for(wgt=0.25;wgt<0.8 && use_ghost==0;wgt+=0.25){
            
                for(j=0;j<gg.get_dim();j++){
                    midpt.set(j,wgt*gg.get_pt(nodes(i)->get_oldCenter(iGhost),j)+(1.0-wgt)*gg.get_pt(dex,j));
                }
                ggWrap.evaluate(midpt,&ftrial);
                if(ftrial>strad.get_target()){
                    use_ghost=1;
                }
            }
            if(use_ghost==0)use_it=0;
        }
    }
    
    if(iclosest>=0){
        if(nodes(iclosest)->get_n_associates()<100)used_because_distance=0;
    }
    
    if(use_it==1 || used_because_distance==1){
        printf("accepting node %d %d -- %e -- %d\n",
        use_it,used_because_distance,gg.get_fn(dex),ggWrap.get_whereCt(iSimplex));
        
        nodes.add(dex,dice,&ggWrap);
        i=nodes.get_dim()-1;
        if(iclosest>=0){
            nodes(i)->set_farthest_associate(nodes(iclosest)->get_farthest_associate());
        }
    }
    
    
    for(i=0;i<nodes.get_dim();i++){
        //printf("about to try node center %d\n",nodes(i)->get_center());
        if(gg.get_fn(nodes(i)->get_center())>strad.get_target()){
            nodes.remove(i);
            i--;
        }
    }
    
}
