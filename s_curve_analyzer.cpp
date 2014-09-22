#include "chisq.h"
#include "aps_extractor.h"

main(int iargc, char *argv[]){

char inputName[letters],word[letters];

int dim,ncenters;
double delta_chi;

dim=22;
ncenters=3;
delta_chi=33.93;

if(iargc>1)ncenters=atoi(argv[1]);
if(iargc>2)dim=atoi(argv[2]);
if(iargc>3)delta_chi=atof(argv[3]);

array_2d<double> data;
array_1d<double> chisq,vv,vvprojected,mu,sig;
array_1d<int> ling;

int i,j;
double nn;

data.set_name("main_data");
chisq.set_name("main_chisq");
ling.set_name("main_ling");
vv.set_name("main_vv");
mu.set_name("main_mu");
sig.set_name("main_sig");
vvprojected.set_name("main_vvprojected");

s_curve chifn(dim,ncenters);

sprintf(inputName,"outputFiles/s_curve_d%d_c%d_output.sav",dim,ncenters);
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
    
    
    fscanf(input,"%le",&nn);
    chisq.add(nn);
    fscanf(input,"%le",&nn);
    mu.add(nn);
    fscanf(input,"%le",&nn);
    sig.add(nn);
    fscanf(input,"%d",&j);
    ling.add(j);
    
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

    sprintf(outname,"processedFiles/s_curve_d%d_c%d_%d_%d_frequentist.sav",dim,ncenters,ix.get_data(i),iy.get_data(i));
    apsExtractor.write_good_points(outname,ix.get_data(i),iy.get_data(i));
        
    //sprintf(outname,"processedFiles/s_curve_d%d_c%d_%d_%d_bayesian.sav",dim,ncenters,ix.get_data(i),iy.get_data(i));
    //apsExtractor.draw_bayesian_bounds(outname,ix.get_data(i),iy.get_data(i),0.95);

}

apsExtractor.write_good_points("processedFiles/s_curve_projected_good.sav");

}
