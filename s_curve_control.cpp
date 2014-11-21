#include "chisq.h"

int main(){

//d=5 means delta_chisq=11
int dim=22,ncenters=3;

s_curve chifn(dim,ncenters);

double v1=1.0,v2=1.0;
char outname[letters];
FILE *output;

int i,j,k;

printf("centers\n");
for(i=0;i<dim;i++){
    for(j=0;j<ncenters;j++){
        printf("%e ",chifn.get_real_center(j,i));
    }
    printf("\n");
}
printf("\n");


//chifn.build_boundary(11.0);
chifn.build_boundary(33.9);//for 22dof
    
for(i=0;i<dim;i++){
    for(j=i+1;j<dim;j++){

        sprintf(outname,"controlFiles/sNarrowTest_d%d_c%d_%d_%d.sav",dim,ncenters,i,j);
        output=fopen(outname,"w");
        for(k=0;k<chifn.get_n_boundary(i,j);k++){
            fprintf(output,"%e %e %e\n",
            chifn.get_boundary(i,j,k,0),
            chifn.get_boundary(i,j,k,1),
            chifn.get_boundary(i,j,k,2));
        }
        fclose(output);
            
    }
}
    

}
