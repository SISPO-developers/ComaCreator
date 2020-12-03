#include "cached_array.h"

int initCachedArrayModule(uint* resolution,
                            uint* chunkSize,
                            t_float* lowerCorner,
                            t_float* upperCorner){


    debugPrint("%u %u %u \n",resolution[0],resolution[1],resolution[2]); 
    for(uint i=0; i<3; ++i){
        if(chunkSize[i]==0 || resolution[i]==0){
            printf("Resolution and chunkSize must be non-zero\n");
            return 1;
        }
        if(((int)resolution[i] % (int)chunkSize[i]) != 0){
            printf("ERROR, resolution must be divisible by the chunkSize\n");
            return 1;
        }
    }

    cA_ptr = initCachedArray(resolution, chunkSize,
                            lowerCorner,upperCorner);

    if(cA_ptr == NULL){
        return 1;
    }
    return 0;
}

//Initialization
CachedArray* initCachedArray(uint* resolution,
                            uint* chunkSize,
                            t_float* lowerCorner,
                            t_float* upperCorner){
    uint nChunks[3];


    for(uint i=0; i<3; ++i){
        nChunks[i] = resolution[i]/chunkSize[i];
    }

    uint size = resolution[0]*resolution[1]*resolution[2];
    uint cellSize = getCellSize();

    //Allocate the cached array 
    CachedArray* cA = (CachedArray*)malloc(sizeof(*cA)+cellSize*size);
    if(cA == NULL){
        printf("malloc failed\n");
        return NULL;
    }    

    for(uint i=0; i<3; ++i){
        cA->lowerCorner[i] = lowerCorner[i];
        cA->upperCorner[i] = upperCorner[i];
        cA->resolution[i] = resolution[i];
        cA->nChunks[i] = nChunks[i];
        cA->chunkSize[i] = chunkSize[i]; 
    }
    
    memset(cA->data,0,size*cellSize);

    return cA;
}

//Free allocated memory
void deleteModule(){
    free(cA_ptr);
}


//This might the stupidest thing I have done, but let's see
//can't know before testing
uint convertIndx(uint* original, CachedArray* cachedArray){
    uint indx_q[3];
    uint indx_r[3];

    for(uint i=0;i<3;++i){
        indx_q[i] = (uint)original[i]/cachedArray->chunkSize[i];
        indx_r[i] = (uint)original[i] % cachedArray->chunkSize[i]; 
    }

    uint chunkBlock = cachedArray->chunkSize[0]*
                     cachedArray->chunkSize[1]*
                     cachedArray->chunkSize[2]; 

    //number of chunks
    uint a1 =   cachedArray->nChunks[1]*
                cachedArray->nChunks[2];

    uint a2 = cachedArray->nChunks[2];


    uint b1 = cachedArray->chunkSize[1]*
              cachedArray->chunkSize[2];

    uint b2 = cachedArray->chunkSize[2];

    uint indx0 = indx_q[0]*a1+indx_q[1]*a2+indx_q[2];
    uint indx1 = indx_r[0]*b1+indx_r[1]*b2+indx_r[2];

    return indx0*chunkBlock+indx1;
}

void helpDebug(){
    debugPrint("%u %u %u ",cA_ptr->resolution[0],
                        cA_ptr->resolution[1],
                        cA_ptr->resolution[2]);


    debugPrint("%f %f %f ",cA_ptr->lowerCorner[0],
                        cA_ptr->lowerCorner[1],
                        cA_ptr->lowerCorner[2]);


    debugPrint("%f %f %f \n",cA_ptr->upperCorner[0],
                        cA_ptr->upperCorner[1],
                        cA_ptr->upperCorner[2]);

}

int findIndices(uint* cachedInd, t_float* pos, CachedArray* cA){
    int curPos[3];
    uint useCache = 1;
    
    for(uint i=0; i<3; ++i){
        curPos[i] = (pos[i]-cA->lowerCorner[i])*cA->resolution[i]/
                    (cA->upperCorner[i]-cA->lowerCorner[i]);
        
        //debugPrint("%f \n",(pos[i]-cA_ptr->lowerCorner[i])*cA_ptr->resolution[i]);
        if(curPos[i]>=(int)cA->resolution[i] || 
                        pos[i]-cA->lowerCorner[i] < 0){
            return 1;
        }
        if(cachedInd[i] != (uint)curPos[i] && useCache){
            useCache = 0;
        }
        cachedInd[i] = (uint)curPos[i];
    }

    if(useCache==0){
        cachedInd[3] = convertIndx(cachedInd, cA);
    }
    return 0;
}    


int addDataToPos(t_float* data, uint start, uint end,
                t_float* pos, uint* cachedInd){

    int retStatus = findIndices(cachedInd, pos, cA_ptr);
    if(retStatus==1) return 1;

    addTo(cachedInd[3], data, start, end, cA_ptr);
    return 0;
}



int setDataToPos(t_float* data, uint start, uint end,
                t_float* pos, uint* cachedInd){

    int retStatus = findIndices(cachedInd, pos, cA_ptr);
    if(retStatus==1) return 1;

    storeTo(cachedInd[3], data, start, end, cA_ptr);
    return 0;
}



int readDataFromPos(t_float* data, uint start, uint end,
                t_float* pos, uint* cachedInd){

    int retStatus = findIndices(cachedInd, pos, cA_ptr);
    if(retStatus==1) return 1;

    readData(cachedInd[3], data, start, end, cA_ptr);

    return 0;
}

void readData(uint indx, t_float* data, 
                uint start, uint end, CachedArray* cA){
    uint j = 0;
    for(uint i = indx*CELLSIZE+start; 
                    i<indx*CELLSIZE+end; ++i){
        data[j++] = cA->data[i];
    }    
}

int getResolution(uint* resolution){
    if(!cA_ptr) return 1;
    
    for (uint i = 0; i < 3; ++i) {    
        resolution[i] = cA_ptr->resolution[i];
    }

    return 0;
}

void storeCurrentPtrTo(uint i){
    if(max_limit == 0){
        max_limit = 5;
        cA_storage = (CachedArray**) malloc(sizeof(CachedArray*)*max_limit);    
    }else if(max_limit<i+1){
        max_limit = max_limit*2;
        cA_storage = (CachedArray**) realloc(cA_storage, sizeof(CachedArray*)*max_limit);
    }
    cA_storage[i] = cA_ptr;
    
}

void switchToPtr(uint i){
    cA_ptr = cA_storage[i];
}


//Combines data to the cA_storage[0]
void combine(uint i, uint start, uint end){
    if(i==0) return;
    if(i>max_limit){
        printf("There is not that much data available\n");
        return;
    }
    uint* res = cA_storage[0]->resolution;
    for(uint  ji = i-1; ji > 0; --ji) {
        for(uint k=0; k<res[0]*res[1]*res[2]; ++k){
            for (uint jj = start; jj < end; ++jj) {    
                cA_storage[0]->data[k*CELLSIZE+jj] +=
                        cA_storage[ji]->data[k*CELLSIZE+jj];
            }
        }
    }
}

void printDebugStatistics(int cA_indx){
    t_float data[CELLSIZE];
    for(uint i = 0; i<CELLSIZE; ++i){
        data[i] = 0.0;
    }
    t_float* data_ptr;
    uint* res;
    if(cA_indx<0){
        data_ptr = cA_ptr->data;
        res = cA_ptr->resolution;
    }else{
        data_ptr = cA_storage[cA_indx]->data;
        res = cA_ptr->resolution;
    }


    for(uint k=0; k<res[0]*res[1]*res[2]; ++k){
        for (uint jj = 0; jj < CELLSIZE; ++jj) {    
            data[jj] += data_ptr[k*CELLSIZE+jj];
        }
    }
    printf("Debug statistic %u\n",res[0]*res[1]*res[2]);
    for(uint i = 0; i<CELLSIZE; ++i){
        printf("%u %f\n",i,data[i]);
    }
    
}


t_float getCellVolume(){
    t_float x[3];
    for (uint i = 0; i < 3; ++i) {
        x[i] = (cA_ptr->upperCorner[i]-cA_ptr->lowerCorner[i])/cA_ptr->resolution[i];
    }
    return x[0]*x[1]*x[2];
}

void getDataAt(uint* indx, t_float* data, uint start, uint end){
    uint inds[4];
    for (uint i = 0; i < 3; ++i) {    
        inds[i] = indx[i];
    }
    uint tmp_indx = convertIndx(inds, cA_ptr);
    readData(tmp_indx, data, start, end ,cA_ptr);
}

void addTo(uint indx, t_float* data, uint start, uint end, CachedArray* cA){
    uint j = 0;
    for(uint i = indx*CELLSIZE+start; 
                    i<indx*CELLSIZE+end; ++i){
        cA->data[i] += data[j++];
    }    
}


void storeTo(uint indx, t_float* data, uint start, uint end, CachedArray* cA){
    uint j = 0;
    for(uint i = indx*CELLSIZE+start; 
                    i<indx*CELLSIZE+end; ++i){
        cA->data[i] = data[j++];
    }    
}

void writeToFile(char* fname){
    CachedArray* ptr = cA_ptr;
    FILE * file= fopen(fname, "wb");

    uint inds = ptr->resolution[0]*ptr->resolution[1]*ptr->resolution[2];

    if (file != NULL) {
        fwrite(ptr, sizeof(*ptr), 1, file);
        fwrite(ptr->data, sizeof(*(ptr->data)),inds*getCellSize(),file);
        debugPrint("%f %f %f %f %f \n", ptr->data[0],
                                ptr->data[1],
                                ptr->data[2],
                                ptr->data[3],
                                ptr->data[4]);
        fclose(file);
    }
}

uint getCellSize(){
    return sizeof(t_float)*CELLSIZE;
}

int readFromFile(char* fname){
    CachedArray* ptr = malloc(sizeof(*cA_ptr));
    FILE * file= fopen(fname, "rb");
    if (file != NULL) {
        int retStatus = fread(ptr, sizeof(*ptr), 1, file);
        uint res = ptr->resolution[0]*ptr->resolution[1]*ptr->resolution[2];
        CachedArray* retVal = realloc(ptr, 
                        sizeof(*ptr)+sizeof(*(ptr->data))*(res*getCellSize()));
        if(retVal==NULL){
            printf("Couldn't allocate large enough array \n");
            return 1;
        }
        ptr = retVal;

        retStatus = fread(ptr->data,sizeof(t_float),res*getCellSize(),file);
        debugPrint("%f %f %f %f %f\n",ptr->data[0],
                                      ptr->data[1],
                                      ptr->data[2],
                                      ptr->data[3],
                                      ptr->data[4]);
        debugPrint("%u %u %u \n",ptr->resolution[0],
                                ptr->resolution[1],
                                ptr->resolution[2]);
        fclose(file);
    }else{
        printf("cached_array: Error while reading file \n");
    }
    cA_ptr = ptr;
    return 0;
};

void getPos(uint* indx, t_float* vec){
    for(uint i = 0; i < 3; ++i) {
        t_float dx = (cA_ptr->upperCorner[i]-cA_ptr->lowerCorner[i])/
                cA_ptr->resolution[i];
        vec[i] = cA_ptr->lowerCorner[i]+(indx[i]+0.5)*dx;
    }
    //printf("%f %f %f %f \n",vec[0],vec[1],vec[2],dx,cA_ptr->lowerCorner[1],cA_ptr->upperCorner[1], cA_ptr->resolution[1]);
}

int posInCell(t_float* pos, uint* cachedInd){
    int curPos[3];
    for (uint i = 0; i < 3; ++i) { 
        curPos[i] = (pos[i]-cA_ptr->lowerCorner[i])*cA_ptr->resolution[i]/
                    (cA_ptr->upperCorner[i]-cA_ptr->lowerCorner[i]);

        if(cachedInd[i] != (uint)curPos[i]){
            return 0;
        }

    }
    return 1;
}


void normalizeData(t_float factor, uint start, uint end){
    uint* res = cA_ptr->resolution;
    for(uint k=0; k<res[0]*res[1]*res[2]; ++k){
        for (uint jj = start; jj < end; ++jj) {    
            cA_ptr->data[k*CELLSIZE+jj] = (cA_ptr->data[k*CELLSIZE+jj])*factor;
        }
    }
}
