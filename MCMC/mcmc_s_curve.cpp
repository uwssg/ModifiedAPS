#include "mcmc.h"
#include <time.h>

double euclideanDistance(array_1d<double>& p1, array_1d<double> &p2){
    
    int i;
    double ans=0.0;
    for(i=0;i<p1.get_dim();i++){
        ans+=power(p1.get_data(i)-p2.get_data(i),2);
    }
    
    return sqrt(ans);

}

main(int iargc, char *argv[]){

int seed=99;

if(iargc>1){
    seed=atoi(argv[1]);
    if(seed<0)seed=int(time(NULL));
}

Ran chaos(seed);
int i,j,dim=22,nchains=5,ncenters=3;

s_curve chifn(dim,ncenters);

array_1d<double> min,max,sig;
array_2d<double> true_centers;

for(i=0;i<dim;i++){
    printf("%e %e\n",chifn.get_real_center(0,i),
    chifn.get_real_center(1,i));
}   

true_centers.set_cols(dim);

for(i=0;i<dim;i++){
    min.set(i,-200.0);
    max.set(i,200.0);
    sig.set(i,10.0);

    for(j=0;j<ncenters;j++){
        true_centers.set(j,i,chifn.get_real_center(j,i));
    }
    
}

printf("time to start declaring stuff\n");
mcmc mcmc_test(dim,nchains,"chains/s_curve_chains",min,max,sig,2.0,&chaos);
mcmc_test.set_chisq(&chifn,1);

printf("done with constructor\n");

mcmc_test.set_statname("chains/s_curve_status.sav");
mcmc_test.set_diagname("chains/s_curve_diagnostic.sav");
mcmc_test.begin_update(1000);
mcmc_test.step_update(1000);
//mcmc_test.cutoff_update(30000);

mcmc_test.generate_random_basis(sig);

mcmc_test.sample(200000);

printf("done sampling\n");

}
