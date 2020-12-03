import sys
import numpy as np

cimport src.common.cached_array as cA 
cimport numpy as cnp
cimport cython
from cpython cimport array
from libc.math cimport sin, cos, acos, exp, sqrt, fabs, M_PI,pow

#types for convenience
t_double = np.float64  
ctypedef cnp.float64_t t_double_T 

t_float = np.float32
ctypedef cnp.float32_t t_float_T 

t_int = np.int
ctypedef cnp.int_t t_int_T
 

def run(output_path,file_list, logger=None):
    cdef unsigned int nArgs = len(file_list)
    cdef unsigned int i,j,k

    cdef t_double_T threshold = 0.01
    cdef t_double_T dist_from_center = 5.0
   
    filelist = []
    cdef bytes py_bytes
    cdef char* input_file
    if logger is not None:
        logger.info("combining")
    for i in range(0,nArgs):
        py_bytes = file_list[i].encode()
        input_file = py_bytes
        cA.readFromFile(<char*> input_file)
        cA.storeCurrentPtrTo(i)
    cA.combine(nArgs, 0, 4)
    cA.switchToPtr(0)
    cdef cnp.ndarray[unsigned int,ndim=1] resolution = np.array((0,0,0),np.uint32)  
    cdef cnp.ndarray[unsigned int,ndim=1] inds = np.array((0,0,0),np.uint32)
    cdef cnp.ndarray[float,ndim=1] pos = np.array((0,0,0),np.float32)
    cdef cnp.ndarray[float,ndim=1] data = np.array((0,0,0,0),np.float32)
    cdef cnp.ndarray[unsigned int,ndim=1] ind_cache = np.array((0,0,0,0),np.uint32)
    cA.getResolution(<unsigned int*> resolution.data)
    cdef unsigned int i_0,j_0,k_0
    
    for i_0 in range(0,resolution[0]):
        for j_0 in range(0,resolution[1]):
            for k_0 in range(0,resolution[2]):
                inds[0] = i_0
                inds[1] = j_0
                inds[2] = k_0
                cA.getPos(<unsigned int*> inds.data,<float*> pos.data)
                if(sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])<dist_from_center):
                    cA.readDataFromPos(<float*> data.data, 0,4,<float*> pos.data,<unsigned int*> ind_cache.data)
                    if(threshold>data[0]):
                        data[0] = -1.0
                        data[1] = -1.0
                        data[2] = -1.0
                        data[3] = -1.0
                        cA.setDataToPos(<float*> data.data,0, 4, <float*> pos.data, <unsigned int*> ind_cache.data)
    if logger is not None:
        logger.info("combined. Writing to file")
    py_bytes = output_path.encode()
    cdef char* output = py_bytes
    cA.writeToFile(<char*> output)
    cA.printDebugStatistics(0)
