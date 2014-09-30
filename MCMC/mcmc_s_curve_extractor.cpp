#include "mcmc_extractor.h"
#include "../eigen_wrapper.h"
#include "kde.h"
#include "../chisq.h"

main(){

mcmc_extractor extractor;
s_curve chifn(22,3);

int nparams,nchains,i,j;
nparams = 22;
nchains=4;

array_1d<double> vv,vvproj;
char rawName[letters],outname[letters];
FILE *input,*output;
double nn,chival;
int ct,ic;
for(ic=0;ic<nchains;ic++){
    sprintf(rawName,"chains/s_curve_chains_%d.txt",ic+1);
    sprintf(outname,"chains/s_curve_chains_projected_%d.txt",ic+1);
    input=fopen(rawName,"r");
    output=fopen(outname,"w");
    while(fscanf(input,"%d",&ct)>0){
        fscanf(input,"%le",&chival);
        for(i=0;i<nparams;i++){
            fscanf(input,"%le",&nn);
            vv.set(i,nn);
        }
        for(i=0;i<nparams;i++){
            vvproj.set(i,chifn.project_to_basis(i,vv));
        }
        
        fprintf(output,"%d %le ",ct,chival);
        for(i=0;i<nparams;i++)fprintf(output,"%le ",vvproj.get_data(i));
        fprintf(output,"\n");
    }
    
    fclose(input);
    fclose(output);
}


extractor.set_nchains(nchains);
extractor.set_nparams(nparams);

//extractor.set_cutoff(12500);
//control has a total 453,051 samples

extractor.set_chainname("chains/s_curve_chains_projected");

extractor.set_keep_frac(0.5);

extractor.learn_thinby();

printf("independent samples %d\n",extractor.get_nsamples());
printf("thinby %d\nused %d\nkept %d\n",extractor.get_thinby(),extractor.get_total_used(),extractor.get_total_kept());
printf("rows %d\n",extractor.get_total_rows());
printf("best_covar %e\n",extractor.get_best_covar());

extractor.print_samples("processed/scurve_samples.sav");

kde kdeObj;
kdeObj.set_data(extractor.get_samples());

//kdeObj.plot_density(0,0.0001,1,0.0005,0.95,"test_pixelstotal.sav",3);
//kdeObj.plot_boundary(0,0.0001,1,0.0005,0.95,"test_boundarytotal.sav",3);

double dx=0.01;

array_1d<int> ix,iy;

ix.set(0,0);
ix.set(1,0);
ix.set(2,14);
iy.set(0,1);
iy.set(0,2);
iy.set(2,17);

for(i=0;i<3;i++){
 
    sprintf(outname,"processed/scurve_%d_%d.sav",ix.get_data(i),iy.get_data(i));

    kdeObj.plot_boundary(ix.get_data(i),dx,iy.get_data(i),dx,0.95,outname,3);
            
    sprintf(outname,"processed/scurve_scatter_%d_%d_control.sav",ix.get_data(i),iy.get_data(i));
    kdeObj.plot_density(ix.get_data(i),dx,iy.get_data(i),dx,0.95,outname,3);
}
        
extractor.plot_delta("processed/scurve_good_pts_testcontrol.sav",33.93);
extractor.plot_as_aps("processed/scurve_mcmc_like_aps.sav");


array_1d<double> RR,VV,WW,mean,var;

RR.set_name("RR");
VV.set_name("VV");
WW.set_name("WW");

extractor.calculate_r(RR,VV,WW);
extractor.calculate_mean(mean,var);

for(i=0;i<nparams;i++){
    if(fabs(var.get_data(i)/mean.get_data(i))>1.0e-20){
        printf("%d %e -- %e %e\n",i,RR.get_data(i),mean.get_data(i),sqrt(var.get_data(i)));
    }
}

extractor.plot_chimin("processed/scurve_chi_min.sav");
//extractor.show_minpt();

///////////////////////////covariance matrix////////////

/*

array_2d<double> covariance;
array_1d<double> mean;
mean.set_dim(nparams);
int j;

for(i=0;i<nparams;i++)mean.set(i,0.0);
for(i=0;i<extractor.get_nsamples();i++){
    for(j=0;j<nparams;j++)mean.add_val(j,extractor.get_sample(i,j));
}
for(i=0;i<nparams;i++){
    mean.divide_val(i,double(extractor.get_nsamples()));
}

covariance.set_dim(nparams,nparams);
for(i=0;i<nparams;i++){
    for(j=0;j<nparams;j++)covariance.set(i,j,0.0);
}

int k;
for(i=0;i<extractor.get_nsamples();i++){
    for(j=0;j<nparams;j++){
        for(k=0;k<nparams;k++){
            covariance.add_val(j,k,\
                  (extractor.get_sample(i,j)-mean.get_data(j))\
                  *(extractor.get_sample(i,k)-mean.get_data(k)));
        }
    }
}

for(j=0;j<nparams;j++){
    for(k=0;k<nparams;k++){
        covariance.divide_val(j,k,double(extractor.get_nsamples()));
    }
}

printf("\ncovariance matrix\n");
for(i=0;i<nparams;i++){
    for(j=0;j<nparams;j++){
        printf("%.3e ",covariance.get_data(i,j));
    }
    printf("\n");

}

printf("\n");

array_2d<double> e_vectors,evbuff;
array_1d<double> e_values,vbuff;

e_vectors.set_dim(nparams,nparams);

eval_symm(covariance,e_vectors,e_values,nparams-2,nparams,1);
eval_symm(covariance,evbuff,vbuff,2,nparams,-1);

e_values.set(nparams-2,vbuff.get_data(0));
e_values.set(nparams-1,vbuff.get_data(1));

for(i=0;i<nparams;i++){
    e_vectors.set(i,nparams-2,evbuff.get_data(i,0));
    e_vectors.set(i,nparams-1,evbuff.get_data(i,1));
}

printf("\neigen vectors\n");
for(i=0;i<nparams;i++){
    printf("%.3e ",e_values.get_data(i));
}
printf("\n\n");
for(i=0;i<nparams;i++){
    for(j=0;j<nparams;j++)printf("%.3e ",e_vectors.get_data(i,j));
    printf("\n");
}

*/

}
