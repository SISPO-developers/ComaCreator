
cimport cython
cimport numpy as cnp
cimport src.common.cached_array as cArray


from cpython cimport array
from libc.math cimport sin, cos, acos, exp, sqrt, fabs, M_PI, pow, log


import numpy as np
import json
import tricubic
import cv2

from src.common import meshWrapper as mw
from src.common import misc as misc
from src.common import orbit as orbit


#types for convenience
t_double = np.float64  
ctypedef cnp.float64_t t_double_T 

t_float = np.float32
ctypedef cnp.float32_t t_float_T 

t_int = np.int
ctypedef cnp.int_t t_int_T


##Loads data from the density map (3D) to data_gas and data_dust 
def load_data():

    cdef t_int_T nResolution
    cdef t_float_T corner

    cdef cnp.ndarray[unsigned int, ndim=1] resolution = np.zeros((3),dtype=np.uint32)

    cArray.getResolution(<unsigned int*> resolution.data)

    cdef cnp.ndarray[t_float_T, ndim=3] data_gas = np.zeros((resolution[0], resolution[1], resolution[2]),dtype=t_float)
    cdef cnp.ndarray[t_float_T, ndim=3] data_dust = np.zeros((resolution[0], resolution[1], resolution[2]),dtype=t_float)

    cdef t_int_T j1,j2,j3,j
    cdef cnp.ndarray[unsigned int,ndim=1] inds = np.zeros((3),dtype=np.uint32)
    cdef cnp.ndarray[t_float_T,ndim=1] data_tmp = np.zeros((7),dtype=np.float32)

    cdef t_float_T dVolume = cArray.getCellVolume()

    for j1 in range(0,resolution[0]):
        for j2 in range(0,resolution[1]):
            for j3 in range(0,resolution[2]):
                inds[0] = j1
                inds[1] = j2
                inds[2] = j3
            
                cArray.getDataAt(<unsigned int*> inds.data, <float*> data_tmp.data, 0, 7);
            
                data_gas[j1,j2,j3] = data_tmp[0]
                data_dust[j1,j2,j3] = data_tmp[4]
    
    max_den_gas = np.amax(data_gas[:,:,:])

    for j1 in range(0,resolution[0]):
        for j2 in range(0,resolution[1]):
            for j3 in range(0,resolution[2]):
                if(data_gas[j1,j2,j3]<0):
                    data_gas[j1,j2,j3] = max_den_gas*0.01


    return data_gas, data_dust





#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.cdivision(True)
def run(
        bbox_lcorner,
        bbox_ucorner,
        resolution_0,
        density_file,
        trajectory,
        mesh_file,
        output_dmap,
        aux_info,
        n_processes):
    
    cdef bytes b_input_file = density_file.encode()
    cdef char* fname = b_input_file

    #Load data
    cArray.readFromFile(fname)

    ###Basic unit is kilometers
    #lower corner and upper corner of the bounding box

    #Lower corner of the density map for Blender
    cdef cnp.ndarray[t_float_T, ndim=1] start_pos = np.array((bbox_lcorner[0],
                                                            bbox_lcorner[0],
                                                            bbox_lcorner[0]),dtype=t_float)

    #Upper corner of the density map
    cdef cnp.ndarray[t_float_T, ndim=1] end_pos = np.array((bbox_ucorner[0],
                                                          bbox_ucorner[0],
                                                          bbox_ucorner[0]),dtype=t_float)
    
    cdef t_float_T dmap_length = bbox_ucorner[0]-bbox_lcorner[0]

    cdef t_int_T i, j, k, a
    cdef t_int_T resolution = resolution_0

    
    cdef t_float_T dX = (end_pos[0]-start_pos[0])/resolution
   
    data_gas, data_dust = load_data()
    print(data_gas.shape)
    smoothed_data_gas = tricubic.tricubic(list(data_gas),list(data_gas.shape)) 
    smoothed_data_dust = tricubic.tricubic(list(data_dust), list(data_dust.shape))
     
    cdef cnp.ndarray[cnp.float32_t,ndim=4] smoothed_data = np.zeros((resolution,resolution,resolution,2), dtype=np.float32)

    cdef cnp.ndarray[t_float_T,ndim=1] pos = np.zeros((3),t_float) 
    cdef cnp.ndarray[t_float_T,ndim=1] pos_tmp = np.zeros((3),t_float) 
    cdef cnp.ndarray[t_int_T,ndim=1] inds = np.zeros((3),t_int) 

    cdef t_float_T max_dust_den = 0
    cdef t_float_T max_gas_den = 0

    #Smooth the data
    for i in range(0,resolution):
        for j in range(0,resolution):
            for k in range(0,resolution):
                inds[0] = i
                inds[1] = j
                inds[2] = k

                for a in range(0, 3):
                    pos[a] = start_pos[a] + dX * inds[a] + dX / 2.0
                    pos_tmp[a] = (pos[a]-bbox_lcorner[a])/dmap_length*data_gas.shape[a]

                smoothed_data[i,j,k,0] = max(smoothed_data_gas.ip(list(pos_tmp)),0.0)
                smoothed_data[i,j,k,1] = max(smoothed_data_dust.ip(list(pos_tmp)),0.0)
    
    #For normalization
    max_gas_den = np.max(smoothed_data[:,:,:,0])
    max_dust_den = np.max(smoothed_data[:,:,:,1])
                

    #Create picture
    cdef t_int_T image_columns = np.ceil(sqrt(resolution*1.0));
    cdef t_int_T imageWidth = resolution*image_columns
    cdef t_int_T imageHeight = resolution*image_columns
    cdef t_int_T imageX,imageY,indX,indY

    cdef cnp.ndarray[cnp.float32_t,ndim=3] img = np.zeros((imageHeight,imageWidth,3), dtype=np.float32)

    print(imageWidth,imageHeight,resolution)

    if(max_dust_den == 0.0 or max_gas_den==0.0):
        print("WARNING:",max_dust_den,max_gas_den)

    imageY = 0
    imageX = 0
    for i in range(0, resolution):                
        if(imageX==image_columns):
            imageX=0
            imageY = imageY+1    
        
        for j in range(0,resolution):
            for k in range(0,resolution):
                indX = imageX*resolution+k
                indY = imageY*resolution+j
                if(max_gas_den != 0.0):
                    img[indY,indX,0] = smoothed_data[i,j,k,0]*1.0/max_gas_den
                if(max_dust_den !=0.0):
                    img[indY,indX,1] = smoothed_data[i,j,k,1]*1.0/max_dust_den
                img[indY,indX,2] = 0.0
        imageX=imageX+1

    print("writing to file")
    cv2.imwrite(output_dmap, img)

    #Print auxiliary file
    aux_file = {}
    aux_file["dimensions"] = [bbox_ucorner[0]-bbox_lcorner[0],bbox_ucorner[1]-bbox_lcorner[1],bbox_ucorner[2]-bbox_lcorner[2]]
    aux_file["resolution"] = resolution
    aux_file["tiles_x"] = image_columns
    aux_file["tiles_y"] = image_columns
    aux_file["tiles_file"] = output_dmap
    print(aux_info)
    with open(aux_info, 'w') as outfile:
        json.dump(aux_file, outfile, indent=2)


