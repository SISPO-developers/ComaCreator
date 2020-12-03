#ifndef CACHED_ARRAY_GUARD_
#define CACHED_ARRAY_GUARD_

#ifndef DEBUG
#define DEBUG 0
#endif

#define debugPrint(fmt, ...) \
            do { if (DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


typedef float t_float;
typedef int t_int;
typedef unsigned int uint;

//g_density, g_velocity, v_density
#define CELLSIZE 7 

typedef struct DataCell DataCell;

//Cached array is an large 3d array
//in which nearby cells are stored close
//to each other in memory, instead of normal 1d way 
//example: chunk 1: (0,0),(0,1),(1,0),(1,1)
//chunk 2 (0,2),(0,3),(1,2),(1,3)
struct CachedArray{
    t_float lowerCorner[3];
    t_float upperCorner[3];
    uint resolution[3];
    uint nChunks[3];
    uint chunkSize[3];
    t_float data[0];
};

typedef struct CachedArray CachedArray;

static CachedArray* cA_ptr = NULL;
static CachedArray** cA_storage = NULL;
static uint max_limit = 0;

int initCachedArrayModule(uint* resolution,
                            uint* chunkSize,
                            t_float* lowerCorner,
                            t_float* upperCorner);


//Initialization
CachedArray* initCachedArray(uint* resolution,
                            uint* chunkSize,
                            t_float* lowerCorner,
                            t_float* upperCorner);


void printDebugStatistics(int indx);

uint getCellSize(void);

void getDataAt(uint* indx, t_float* data, uint start, uint end);

t_float getCellVolume(void);

int readFromFile(char* fname);

void helpDebug(void);

void writeToFile(char* fname);

int getResolution(uint* resolution);

int findIndices(uint* cachedInd, t_float* pos, CachedArray* cA);


void storeTo(uint indx, t_float* data, uint start, uint end, CachedArray* cA);

void addTo(uint indx, t_float* data, uint start, uint end, CachedArray* cA);

void readDataFrom(uint indx, t_float* data, 
                uint start, uint end, uint* cachedInd);

void readData(uint indx, t_float* data, uint start, uint end, CachedArray* cA);

uint convertIndx(uint* original, CachedArray* cachedArray);

int addDataToPos(t_float* data, uint start, uint end, 
                t_float* pos, uint* cachedInd);


int setDataToPos(t_float* data, uint start, uint end, 
                t_float* pos, uint* cachedInd);


int readDataFromPos(t_float* data, uint start, uint end,
                    t_float* pos, uint* cachedInd);

//TODO: go through cA_storage 
void deleteModule(void);

void storeCurrentPtrTo(uint i);

void combine(uint imax, uint start, uint end);

int posInCell(t_float* pos, uint* indx);

//Get middle position
void getPos(uint* dir, t_float* vec);

void switchToPtr(uint i);

void normalizeData(t_float factor, uint start, uint end);
#endif
