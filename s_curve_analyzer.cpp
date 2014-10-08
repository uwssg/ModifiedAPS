#include "chisq.h"
#include "gp_wrapper.h"
#include "aps_extractor.h"

main(int iargc, char *argv[]){

int i,j;
array_1d<double> lnchi_hist,xmin,xmax;
array_1d<int> apsHist,simplexHist,coulombHist,compassHist;
array_1d<int> bisectHist,ricochetHist,totalHist;

for(i=0;i<120;i++){
    lnchi_hist.set(i,-2.0+i*0.1);
    apsHist.set(i,0);
    simplexHist.set(i,0);
    coulombHist.set(i,0);
    bisectHist.set(i,0);
    ricochetHist.set(i,0);
    totalHist.set(i,0);
    compassHist.set(i,0);
}

char inputName[letters],word[letters];

int dim,ncenters,seed;
double delta_chi;

dim=22;
ncenters=3;
seed=99;
delta_chi=33.93;

if(iargc>1)ncenters=atoi(argv[1]);
if(iargc>2)dim=atoi(argv[2]);
if(iargc>3)delta_chi=atof(argv[3]);

array_2d<double> data;
array_1d<double> chisq,vv,vvprojected,mu,sig;
array_1d<int> ling;

for(i=0;i<dim;i++){
    xmin.set(i,2.0*chisq_exception);
    xmax.set(i,-2.0*chisq_exception);
}

int hdex;
double nn,chival;
array_1d<int> *hptr;

data.set_name("main_data");
chisq.set_name("main_chisq");
ling.set_name("main_ling");
vv.set_name("main_vv");
mu.set_name("main_mu");
sig.set_name("main_sig");
vvprojected.set_name("main_vvprojected");

s_curve chifn(dim,ncenters);

sprintf(inputName,"outputFiles/s_curve_d%d_c%d_s%d_output.sav",dim,ncenters,seed);
FILE *input,*output;

output=fopen("outputFiles/s_curve_projected.sav","w");
input=fopen(inputName,"r");
for(i=0;i<dim+5;i++){
    fscanf(input,"%s",word);
    fprintf(output,"%s ",word);
}
fprintf(output,"\n");
while(fscanf(input,"%le",&nn)>0){
    vv.set(0,nn);
    for(i=1;i<dim;i++){
        fscanf(input,"%le",&nn);
        vv.set(i,nn);
    }
    
    
    fscanf(input,"%le",&chival);
    chisq.add(chival);
    fscanf(input,"%le",&nn);
    mu.add(nn);
    fscanf(input,"%le",&nn);
    sig.add(nn);
    fscanf(input,"%d",&j);
    ling.add(j);
    
    hdex=get_dex(lnchi_hist,log(chival));
    //printf("log %e hdex %d\n",log(chival),hdex);
    if(j==iAPS){
        hptr=&apsHist;
    }
    else if(j==iSimplex){
        hptr=&simplexHist;
    }
    else if(j==iCoulomb){
        hptr=&coulombHist;
    }
    else if(j==iCompass){
        hptr=&compassHist;
    }
    else if(j==iBisect || j==iNodeBisect){
        hptr=&bisectHist;
    }
    else if(j==iRicochet){
        hptr=&ricochetHist;
    }
    else{
        printf("WARNING ling %d\n",j);
    }
    
    for(i=hdex;i<hptr->get_dim();i++){
        hptr->add_val(i,1);
    }
    for(i=hdex;i<totalHist.get_dim();i++){
        totalHist.add_val(i,1);
    }
    
    if(chival<=delta_chi){
        for(i=0;i<dim;i++){
            if(vv.get_data(i)<xmin.get_data(i)){
                xmin.set(i,vv.get_data(i));
            }
            if(vv.get_data(i)>xmax.get_data(i)){
                xmax.set(i,vv.get_data(i));
            }
        }
    }
    
    
    for(i=0;i<dim;i++){
        vvprojected.set(i,chifn.project_to_basis(i,vv));
    }
    
    for(i=0;i<dim;i++)fprintf(output,"%le ",vvprojected.get_data(i));
    fprintf(output,"%le %le %le %d\n",
        chisq.get_data(chisq.get_dim()-1),
        mu.get_data(mu.get_dim()-1),
        sig.get_data(sig.get_dim()-1),
        ling.get_data(ling.get_dim()-1));

}

fclose(input);
fclose(output);

printf("\n");
for(i=0;i<dim;i++){
    printf("d%d %e to %e\n",
    i,xmin.get_data(i),xmax.get_data(i));
}
printf("\n");

aps_extractor apsExtractor;
apsExtractor.set_filename("outputFiles/s_curve_projected.sav");
apsExtractor.set_delta_chi(delta_chi);

char outname[letters];
array_1d<int> ix,iy;

ix.set(0,0);
ix.set(1,0);
ix.set(2,14);

iy.set(0,1);
iy.set(1,2);
iy.set(2,17);

for(i=0;i<ix.get_dim();i++){

    sprintf(outname,"processedFiles/s_curve_d%d_c%d_s%d_%d_%d_frequentist.sav",dim,ncenters,seed,ix.get_data(i),iy.get_data(i));
    apsExtractor.write_good_points(outname,ix.get_data(i),iy.get_data(i));
        
    //sprintf(outname,"processedFiles/s_curve_d%d_c%d_%d_%d_bayesian.sav",dim,ncenters,ix.get_data(i),iy.get_data(i));
    //apsExtractor.draw_bayesian_bounds(outname,ix.get_data(i),iy.get_data(i),0.95);

}

apsExtractor.write_good_points("processedFiles/s_curve_projected_good.sav");


sprintf(outname,"processedFiles/s_curve_d%d_c%d_s%d_histograms.sav",dim,ncenters,seed);
output=fopen(outname,"w");
fprintf(output,"#lnchi aps simplex coulomb compass bisect ricochet total\n");
for(i=0;i<lnchi_hist.get_dim();i++){
    fprintf(output,"%le %d %d %d %d %d %d %d\n",
    lnchi_hist.get_data(i),
    apsHist.get_data(i),
    simplexHist.get_data(i),
    coulombHist.get_data(i),
    compassHist.get_data(i),
    bisectHist.get_data(i),
    ricochetHist.get_data(i),
    totalHist.get_data(i));
}
fclose(output);

}
