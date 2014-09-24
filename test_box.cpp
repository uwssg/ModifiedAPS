#include "box.h"

main(){

Ran chaos(99);

array_2d<double> data;
array_1d<double> min,max;
int rows,cols;

cols=22;
rows=1000;

int i,j;

data.set_cols(cols);
for(i=0;i<cols;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
}

for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
        data.set(i,j,min.get_data(j)+chaos.doub()*(max.get_data(j)-min.get_data(j)));
    }
}

box testBox(&data,20);

array_1d<double> vv;
int ct=0,i_tree,i_dir,i_box;
int betterFit,isOkay;
double before=double(time(NULL));
for(ct=0;ct<100000;ct++){
    for(i=0;i<cols;i++){
        vv.set(i,-200.0+chaos.doub()*400.0);
    }
    
    i_box=testBox.find_box(vv,&i_tree,&i_dir);
    
    isOkay=1;
    for(i=0;i<cols && isOkay==1;i++){
        if(vv.get_data(i)<testBox.get_box_min(i_box,i))isOkay=0;
        if(vv.get_data(i)>testBox.get_box_max(i_box,i))isOkay=0;
    }
    
    if(isOkay==0){
        betterFit=-1;
        for(i=0;i<testBox.get_nboxes() && betterFit==-1;i++){
            betterFit=i;
            for(j=0;j<cols && betterFit==i;j++){
                if(vv.get_data(i)<testBox.get_box_min(i,j))betterFit=-1;
                if(vv.get_data(i)>testBox.get_box_max(i,j))betterFit=-1;
            }
        }
        
        
        if(betterFit>=0){
            printf("FAILED at %d %e\n",ct,double(time(NULL))-before);
            for(i=0;i<cols;i++){
                printf("%e -- %e %e -- %e %e\n",
                vv.get_data(i),testBox.get_box_min(i_box,i),testBox.get_box_max(i_box,i),
                testBox.get_box_min(betterFit,i),testBox.get_box_max(betterFit,i));
            }
            
            exit(1);
        }
    }
    
    data.add_row(vv);
    testBox.add_pt();
    
}


}
