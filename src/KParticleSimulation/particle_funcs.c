#include "particle_funcs.h"

void save_to_file2(char* fname){
    writeToFile(fname);
}

void load_from_file2(char* fname){
    readFromFile(fname);
}

int add_data_to_cached2(float value, 
                        float* pos, 
                        uint* indStore){
    
    return addDataToPos(&value,4,5,pos,indStore);
}

void normalize(float factor, uint start, uint end){
    normalizeData(factor, start, end);
}
