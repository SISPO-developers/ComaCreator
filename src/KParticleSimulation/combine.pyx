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
 

def run(output_path,file_list):
    cdef unsigned int nArgs = len(file_list)
    cdef unsigned int i,j,k

    cdef t_double_T threshold = 0.01
    cdef t_double_T dist_from_center = 0.0
   
    filelist = []
    cdef bytes py_bytes
    cdef char* input_file
    for i in range(0,nArgs):
        py_bytes = file_list[i].encode()
        input_file = py_bytes
        cA.readFromFile(<char*> input_file)
        cA.storeCurrentPtrTo(i)

    cA.combine(nArgs,4,7)

    py_bytes = output_path.encode()
    cdef char* output = py_bytes
    cA.writeToFile(<char*> output)

