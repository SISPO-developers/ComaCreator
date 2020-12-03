import numpy as np
import pathlib
import json
import sys
import os
import math
from libc.math cimport sin, cos, acos, exp, sqrt, fabs, M_PI,pow,floor
cimport src.common.cached_array as cA
from cpython cimport array
cimport numpy as cnp
cimport cython

from pyembree.mesh_construction import TriangleMesh
from pyembree import rtcore_scene as rtcs

from src.common import meshWrapper as mw
from src.common import misc as misc

##Type definitions
t_double = np.float64  
ctypedef cnp.float64_t t_double_T 

t_float = np.float32
ctypedef cnp.float32_t t_float_T 

t_int = np.int32
ctypedef cnp.int_t t_int_T
t_uint = np.uint32
ctypedef cnp.uint32_t t_uint_T

def init_and_run(from_face_in,
                to_face_in,
                mesh_scale,             ##scale of the mesh         [obj space to km]
                fname_mesh,             ##read mesh from this file
                outflow_rates_in,       ##Outflow file
                T_in,                   ##Temperature of the gas
                R_in,                   ##Specific gas constant
                nBounds_in,             ##How many boundaries
                lower_corner_ind_in,
                upper_corner_ind_in,
                sample_points_in,
                bbox_lower_in,
                bbox_upper_in,
                output_to,
                gas_flow_scale,
                rot_mat,
				comet_location,
                chunk_size0):


    cdef t_float_T T = T_in                                     
    cdef t_float_T R = R_in                   
    cdef t_uint_T nBounds = nBounds_in;
    cdef cnp.ndarray[t_uint_T,ndim=1] lower_corner = np.array(lower_corner_ind_in,np.uint32);
    cdef cnp.ndarray[t_uint_T,ndim=1] upper_corner = np.array(upper_corner_ind_in,np.uint32);
    cdef cnp.ndarray[t_uint_T,ndim=1] resolution = np.array(sample_points_in,np.uint32);
    cdef cnp.ndarray[t_uint_T,ndim=1] chunk_sizes = np.array(chunk_size0, np.uint32);
    cdef cnp.ndarray[t_float_T,ndim=1] bbox_lower = np.array(bbox_lower_in,t_float);
    cdef cnp.ndarray[t_float_T,ndim=1] bbox_upper = np.array(bbox_upper_in,t_float); 
    #For paraller runs
    cdef t_uint_T from_face = from_face_in
    cdef t_uint_T to_face = to_face_in
    print(resolution[0],resolution[1],resolution[2])
    cA.initCachedArrayModule(<unsigned int*> resolution.data,
                                <unsigned int*> chunk_sizes.data,
                                <float*> bbox_lower.data,
                                <float*> bbox_upper.data);

    def specFoo(v):
        if np.random.rand(1,t_float) > 0.1:
            return v*0.1
        return v


    mesh = mw.MeshObject(fname_mesh, rot_mat, scale=mesh_scale)
    outflow_rates = None
    if(outflow_rates_in is not None):
        print("Reading outflow_rates_from")
        outflow_rates = np.loadtxt(outflow_rates_in).astype(t_float)
    else:
        print("generating outflow rates")
        n_faces = mesh.n_faces()
        outflow_rates = np.random.rand((n_faces,2),t_float)
        outflow_rates = map(outflow_rates, specFoo)
    outflow_rates = outflow_rates*gas_flow_scale

    main_sim(mesh,
            outflow_rates,
            T, 
            R,
            from_face,
            to_face,
            resolution)

    print(output_to)
    cdef bytes py_bytes = output_to.encode()
    cdef char* output = py_bytes
    cA.writeToFile(<char*> output)
    cA.printDebugStatistics(-1);




##fname: surface of the object
##nDIm_tmp: To how many points the space is divided
##cornder_1_tmp: lower corder of the space
##corner_2_tmp: upper cornder ot the space
##T_tmp: temperature of the gas
##R_tmp: specific gas constant 
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef main_sim(mesh,
              outflow_rates_in,
              t_float_T T, 
              t_float_T R, 
              t_uint_T from_face,
              t_uint_T to_face,
              cnp.ndarray[t_uint_T, ndim=1] resolution):

    cdef cnp.ndarray[t_float_T,ndim=2] normals
    cdef cnp.ndarray[t_float_T,ndim=2] areas
    cdef cnp.ndarray[t_float_T,ndim=2] centroids
    cdef cnp.ndarray[t_float_T,ndim=2] outflow_rates = outflow_rates_in

    normals, areas, centroids = mesh.get_face_information(t_float)

    ##Prepare pyembree
    scene = rtcs.EmbreeScene()
    #embree_mesh = TriangleMesh(scene, mesh.get_vertices(t_float), mesh.get_faces(t_float))
    #TriangleMesh(scene, mesh.get_vertices(t_float), mesh.get_faces(t_float))
    TriangleMesh(scene, mesh.mesh.vertices, mesh.mesh.faces)
    ##pyembree uses these        
    cdef unsigned int RAYBUFFERSIZE = 10000
    cdef cnp.ndarray[t_float_T,ndim=2] ray_dir = np.zeros((RAYBUFFERSIZE,3),dtype=t_float)
    cdef cnp.ndarray[t_float_T,ndim=2] ray_origin = np.zeros((RAYBUFFERSIZE,3),dtype=t_float)  
    cdef cnp.ndarray[t_float_T,ndim=1] distances = np.zeros((RAYBUFFERSIZE),dtype=t_float)
    cdef cnp.ndarray[t_uint_T,ndim=2] inds = np.zeros((RAYBUFFERSIZE,2),dtype=np.uint32)

    
    ##simulation constants
    cdef t_float_T beta = 1.0/(2.0*R* T)
    cdef t_float_T qi
    cdef t_float_T U0
    cdef t_float_T u0 
    cdef t_float_T E
    
    cdef cnp.ndarray[t_float_T,ndim=1] n1 = np.zeros((3),t_float)
    cdef cnp.ndarray[t_float_T,ndim=1] p1 = np.zeros((3),t_float)

    cdef t_float_T normal_offset = 0.000001

    from_face = max(0,from_face)
    to_face = min(normals.shape[0],to_face)

    ##Go though the all the normals
    cdef t_int_T i,j
    for i in range(from_face,to_face):
        if(floor((i-from_face)*100.0/(to_face-from_face))%10==0):
            print("Face index: {} to {}".format(i,to_face))

        ##normal of the face
        ##center of the faces
        for j in range(0,3):
            n1[j] = normals[i,j]
            ##offset is required to prevent pyembree detecting 
            ##the source face
            p1[j] = centroids[i,j]+n1[j]*normal_offset

        E = areas[i]
        qi = outflow_rates[i,0]           ##source strength
        u0 = outflow_rates[i,1]           ##outflow velocity
        

        traverse_points(p1,n1,E,qi,u0,
                        beta,
                        ray_origin,ray_dir,
                        distances,inds,RAYBUFFERSIZE,scene,
                        resolution)





@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef void traverse_points(
                            cnp.ndarray[t_float_T, ndim=1] p1,
                            cnp.ndarray[t_float_T, ndim=1] n1,
                            t_float_T Ei,
                            t_float_T qi,
                            t_float_T u0,
                            t_float_T beta,
                            cnp.ndarray[t_float_T, ndim=2] ray_origin,
                            cnp.ndarray[t_float_T, ndim=2] ray_dir,
                            cnp.ndarray[t_float_T, ndim=1] distances,
                            cnp.ndarray[unsigned int, ndim=2] inds,
                            t_int_T BUFFERSIZE,
                            object scene,
                            cnp.ndarray[t_uint_T, ndim=1] resolution):


    cdef array.array r2_tmp = array.array('f', [1, 2, 3])
    cdef t_float_T[:] r2 = r2_tmp
  
    cdef t_float_T U0
    U0 = u0 * sqrt(beta)
    
   
      
    cdef t_float_T rlen
    cdef cnp.ndarray[t_float_T,ndim=1] pos_tmp = np.empty((3),np.float32)
    cdef unsigned int sIndx = 0

    cdef t_float_T cos_angle
    cdef t_int_T i1
    cdef unsigned int i_0,j_0,k_0,j
    cdef int retVal
    cdef cnp.ndarray[unsigned int,ndim=1] cell_ind = np.empty((3),np.uint32)

    cdef cnp.ndarray[unsigned int,ndim=1] cache_inds = np.array((0,0),np.uint32)
    cdef unsigned int indx = 0  
    #print("resolution",resolution[0],resolution[1],resolution[2])
    for i_0 in range(0,resolution[0]):
        for j_0 in range(0,resolution[1]):
            for k_0 in range(0,resolution[2]):
                cell_ind[0] = i_0
                cell_ind[1] = j_0
                cell_ind[2] = k_0
                cA.getPos(<unsigned int*>cell_ind.data, <float*> pos_tmp.data)
                #print(pos_tmp[0],pos_tmp[1],pos_tmp[2])
                for j in range(0,3):
                    r2[j] = pos_tmp[j]-p1[j]

                cos_angle = n1[0]*r2[0]+n1[1]*r2[1]+n1[2]*r2[2]
                #point behind the face
                if(cos_angle<=0.0):
                    #check if the voxel center is only behind the centroid
                    retVal = cA.posInCell(<float*> p1.data, <unsigned int*> cell_ind.data);
                    if(retVal==1):
                        #shift voxel center so the data will be registered
                        #it's not beautiful solution, but prevents
                        #blocky boundaries
                        for j in range(0,3):
                            pos_tmp[j] = p1[j]+n1[j]*0.05
                            r2[j] = pos_tmp[j]-p1[j]
                    else:
                        continue

                rlen = sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2])
                if(rlen<=0.0001):
                    continue
                #Store values for pyembree
                for j in range(0,3):                
                    ray_dir[indx,j] = r2[j]/rlen
                    ray_origin[indx,j] = p1[j]
                distances[indx] = rlen

                indx=indx+1
                ##If ray tracing buffer is full
                ##evaluate buffer
                if(indx>=BUFFERSIZE):
                    evaluate_buffer(p1,n1,
                                Ei,U0,u0,qi,
                                ray_origin,ray_dir,
                                distances,inds,
                                BUFFERSIZE,scene)
                    indx=0
    ##Evalue the rest of the buffer
    if(indx>0):
        evaluate_buffer(p1,n1,
                        Ei,U0,u0,qi,
                        ray_origin,ray_dir,
                        distances,inds,
                        BUFFERSIZE,scene)

    


##I had to break the function into two because the buffer needs to be evaluated 
##time in order to avoid large chunks of rays that would occupy entire memory
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef void evaluate_buffer(
                            cnp.ndarray[t_float_T, ndim=1] p1,
                            cnp.ndarray[t_float_T, ndim=1] n1,
                            t_float_T Ei,
                            t_float_T U0,
                            t_float_T u0,
                            t_float_T qi,
                            cnp.ndarray[t_float_T, ndim=2] ray_origin,
                            cnp.ndarray[t_float_T, ndim=2] ray_dir,
                            cnp.ndarray[t_float_T, ndim=1] distances,
                            cnp.ndarray[unsigned int, ndim=2] inds,
                            t_int_T BUFFERSIZE,
                            object scene):

    cdef t_float_T rho
    cdef t_float_T tmp

    cdef t_float_T stheta,ctheta
    cdef t_float_T rlen

    #tracer function
    cdef cnp.ndarray ray_status = scene.run(ray_origin, ray_dir, distances)
    
    cdef array.array n2_tmp = array.array('f', [1, 2, 3])
    cdef t_float_T[:] n2 = n2_tmp
    

    #Pass data to the cached_array,
    #ind_cache is for storing indices for temporary use
    cdef cnp.ndarray[t_float_T, ndim=1] p_0 = np.zeros((3),t_float)
    cdef cnp.ndarray[t_float_T,ndim=1] tmp_data = np.zeros((5),t_float)
    cdef cnp.ndarray[unsigned int,ndim=1] ind_cache = np.array((0,0,0,0),np.uint32)
    
    #Avogadros number
    cdef t_float_T Na = 6.02214086e23

    cdef t_uint_T start = 0
    cdef t_uint_T end = 4

    cdef unsigned int j,ray_indx

    #Go throught the ray buffer
    for ray_indx in range(0,BUFFERSIZE):
        #If status code is less than 0 it means that the ray
        #did not hit anything
        if(ray_status[ray_indx]<0):
            for j in range(0,3):
                n2[j] = ray_dir[ray_indx,j] 
            rlen = distances[ray_indx]
            ctheta = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2]
            if(ctheta>=1.0):
                ctheta = 1.0
                stheta = 0.0
            else:
                stheta = sin(acos(ctheta))
            #print(rlen,distances[ray_indx])
            rho = U0/M_PI*ctheta/(rlen*rlen)*Ei*qi*exp(-U0*U0*stheta*stheta)
            #store data to tmp_data that is passed to cached_array
            tmp_data[0] = rho
            tmp = u0*ctheta
            for j in range(0,3):
                tmp_data[1+j] = n2[j]*rho*tmp
                p_0[j] = ray_origin[ray_indx,j]+rlen*ray_dir[ray_indx,j]
            cA.addDataToPos(<float*> tmp_data.data,start,end,<float*> p_0.data,<unsigned int*> ind_cache.data)
