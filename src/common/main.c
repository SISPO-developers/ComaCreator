#include <stdio.h>
#include "cached_array.h"
#include "minunit/minunit.h"


MU_TEST(test_creation) {
    t_float upperCorner[3] = {10.0, 10.0, 10.0};
    t_float lowerCorner[3] = {-10.0, -10.0, -10.0};
    uint chunkSize[3] = {4,4,4};
    uint resolution[3] = {16,16,16}; 

    mu_check(initCachedArrayModule(resolution,
                                    chunkSize,
                                    lowerCorner,
                                    upperCorner)==0);

}

int boundCheck(t_float* vec, uint i){
    return (vec[i]>=-10 && vec[i]<10);
}


MU_TEST(test_insertion){
    t_float vec[3] = {0,0,0};
    int count = 0;
    uint cachedInd[4] = {0,0,0,0};  
    for (int i = -11; i < 12; ++i) {
        for(int j = -11; j<12; ++j){
            for (int k = -11; k < 12; ++k) {
                vec[0] = i; vec[1] = j; vec[2] = k;
                t_float data[7] = {1.0,2.0,3.0,4.0,5.0,6.0};
                if(boundCheck(vec,0) && boundCheck(vec,1) && boundCheck(vec,2)){
                    
                    mu_check(addDataToPos(data,0,7,vec,cachedInd)==0);
                    count += 1;
                }else{
                    mu_check(addDataToPos(data,0,7,vec,cachedInd)==1);
                    //printf("%f %f %f \n",vec[0],vec[1],vec[2]);
                } 
            }       
        }
    }
    //printf("%i \n",count);
    mu_check(count == 20*20*20);
    
}


MU_TEST(test_io_and_insert){
    char fname[] = "test.out";
    t_float vec[3] = {0,0,0};
    t_float data[7];
    t_float data2[7] = {1,1,1,1,1,1,1};
    uint cachedInd[4];
    mu_check(readDataFromPos(data,0,7,vec,cachedInd)==0);
    writeToFile(fname);
    
    mu_check(addDataToPos(data2,0,7,vec,cachedInd)==0);
    mu_check(readDataFromPos(data2,0,7,vec,cachedInd)==0);
    for(uint i = 0; i < 7; ++i) {    
        mu_check(data[i]==(data2[i]-1));
    }

    deleteModule();
    mu_check(readFromFile(fname)==0);


    mu_check(readDataFromPos(data2,0,7,vec,cachedInd)==0);

    for(uint i = 0; i < 7; ++i) {    
        mu_check(data[i]==(data2[i]));
    }
    data2[1] = -1; 
    mu_check(addDataToPos(data2+1,1,2,vec,cachedInd)==0);
    mu_check(readDataFromPos(data2,0,7,vec,cachedInd)==0);

    for(uint i = 0; i < 7; ++i) {    
        if(i==1){
            mu_check(data[i]==(data2[i]+1));
        }else{
            mu_check(data[i]==data2[i]);
        } 
    }
    deleteModule();
}

MU_TEST(test_ptr_switching){
    char fname[] = "test.out";
    readFromFile(fname);
    t_float data[7];
    t_float data2[7];
    uint cachedInd[4];
    t_float vec[3] = {3,3,3};
    mu_check(readDataFromPos(data,0,7,vec,cachedInd)==0);
    storeCurrentPtrTo(0);
    for(uint i = 1; i < 6; ++i) {
        readFromFile(fname);
        storeCurrentPtrTo(i);
    }
    combine(6,0,7);
    switchToPtr(0);
    mu_check(readDataFromPos(data2,0,7,vec,cachedInd)==0);
    
    for(uint i =0; i<CELLSIZE; ++i){
        mu_check(data[i]*6==data2[i]);
    }
}



MU_TEST_SUITE(test_suite) {
    MU_RUN_TEST(test_creation);
    MU_RUN_TEST(test_insertion);
    MU_RUN_TEST(test_io_and_insert);
    MU_RUN_TEST(test_ptr_switching);
}

//Test chunk array
int main (int argc, const char * argv[]) {
    MU_RUN_SUITE(test_suite);
    MU_REPORT();
    return MU_EXIT_CODE;

}
