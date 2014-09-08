#include "chisq.h"

main(){

//d=5 means delta_chisq=11
int dim=22,ncenters=3;

ellipses chifn(dim,ncenters);;

double v1=1.0,v2=1.0;
char outname[letters];
FILE *output;

int i,j,k;


chifn.build_boundary(33.9);//for 22dof
    
for(i=0;i<dim;i++){
    for(j=i+1;j<dim;j++){

        sprintf(outname,"controlFiles/s_d%d_c%d_%d_%d.sav",dim,ncenters,i,j);
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
