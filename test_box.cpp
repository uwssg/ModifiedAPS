#include "box.h"

main(){

Ran chaos(99);

array_2d<double> data;
array_1d<double> min,max;

data.set_name("test_data");
min.set_name("test_min");
max.set_name("test_max");

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

array_1d<double> vv,base;
vv.set_name("test_vv");
base.set_name("test_base");

for(i=0;i<cols;i++){
    base.set(i,-110.0+chaos.doub()*220.0);
    vv.set(i,base.get_data(i));
}

int idim=chaos.int32()%cols;

int ct=0,i_tree,i_dir,i_box;
int betterFit,isOkay;
double before=double(time(NULL));
for(ct=0;ct<100000;ct++){
    
    if(ct%5000==0){
        for(i=0;i<cols;i++)base.set(i,-110.0+chaos.doub()*220.0);
        for(i=0;i<cols;i++)vv.set(i,base.get_data(i));
        idim=chaos.int32()%cols;
    }
    
    vv.set(idim,base.get_data(idim)+1.0e-6*(ct%5000));
    
    i_box=testBox.find_box(vv,&i_tree,&i_dir);
    i=chaos.int32()%cols;
    vv.set(i,testBox.get_box_max(i_box,i));
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
                if(vv.get_data(j)<testBox.get_box_min(i,j))betterFit=-1;
                if(vv.get_data(j)>testBox.get_box_max(i,j))betterFit=-1;
            }
        }
        
        
        /*if(betterFit>=0){
            printf("FAILED at %d %e\n",ct,double(time(NULL))-before);
            printf("ibox %d betterFit %d\n",i_box,betterFit);
            for(i=0;i<cols;i++){
                printf("%e -- %e %e -- %e %e\n",
                vv.get_data(i),testBox.get_box_min(i_box,i),testBox.get_box_max(i_box,i),
                testBox.get_box_min(betterFit,i),testBox.get_box_max(betterFit,i));
            }
            
            exit(1);
        }*/
    }
    
    data.add_row(vv);
    testBox.add_pt();
    if(ct%1000==0){
        printf("ndata %d \n",data.get_rows());
        printf("smallest box %d\n",testBox.get_smallest_box());
        printf("biggest box %d\n",testBox.get_biggest_box());
        printf("nboxes %d\n",testBox.get_nboxes());
        printf("\n");
    }
    
}


}
