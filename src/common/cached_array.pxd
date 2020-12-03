

cdef extern from "cached_array.c":
    
    int initCachedArrayModule(unsigned int* resolution,
                                unsigned int* chunkSize,
                                float* lowerCorner,
                                float* upperCorner);
    float getCellVolume();
    void getDataAt(unsigned int * indx, float* data, unsigned int start, unsigned int end);
    unsigned int getCellSize();
    int readFromFile(char* fname);
    void helpDebug();
    void writeToFile(char* fname);
    int getResolution(unsigned int* resolution);
    void readDataFrom(unsigned int indx, float* data, 
                    unsigned int start, unsigned int end, unsigned int* cachedInd);
    int addDataToPos(float* data, unsigned int start, unsigned int end, 
                    float* pos, unsigned int* cachedInd);
    int readDataFromPos(float* data, unsigned int start, unsigned int end,
                        float* pos, unsigned int* cachedInd);
    int setDataToPos(float* data, unsigned int start, unsigned int end,
                        float* pos, unsigned int* cachedInd);
    
    void deleteModule();
    int posInCell(float* pos, unsigned int* indx);
    void storeCurrentPtrTo(unsigned int i);
    void switchToPtr(unsigned int i);
    void combine(unsigned int imax, unsigned int start, unsigned int end);

    void getPos(unsigned int* indx, float* vec);

    void printDebugStatistics(int indx);
