import numpy as np
import time

from libc.math cimport sin, cos, acos, exp, sqrt, fabs, M_PI
from numpy.ctypeslib import ndpointer
from cpython cimport array


from pyembree.mesh_construction import TriangleMesh
from pyembree import rtcore_scene as rtcs


from src.common import meshWrapper as mw
from src.common import misc as misc


cimport cython
cimport numpy as cnp

t_double = np.float64  
ctypedef cnp.float64_t t_double_T 

t_float = np.float32
ctypedef cnp.float32_t t_float_T 

t_int = np.int
ctypedef cnp.int_t t_int_T


cdef extern from "particle_funcs.c":
    int update_particles(double* buffer,
                            double* particle_data, 
                            double h, 
                            long nParticles,
                            double m_comet,
                            double mu_orbit,
                            float* cached_data,
                            double border,
                            unsigned int* indStore);

    void normalize(float factor, 
                    unsigned int start, 
                    unsigned int end);
    void save_to_file2(char* fname);
    void load_from_file2(char* fname);
    int add_data_to_cached2(float value, 
                            float* pos, 
                            unsigned int* indStore);



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef t_int_T generate_particles(t_int_T start_indx, 
                                t_int_T number_of_particles,
                                t_double_T h_init,
                                t_double_T v_init,
                                t_double_T R_mu,
                                t_double_T R_sigma,
                                t_double_T C_d,
                                t_double_T rho_dust,
                                t_double_T molar_mass,
                                cnp.ndarray[t_double_T,ndim=2] normals,
                                cnp.ndarray[t_double_T,ndim=2] area,
                                cnp.ndarray[t_double_T,ndim=2] centroid,
                                t_int_T face_index0,
                                cnp.ndarray[t_double_T,ndim=2] particle_array,
                                sssb_mesh):

    cdef t_int_T i,j,face_index
    sample = None

    for i in range(start_indx,start_indx+number_of_particles):
        if(face_index0 == -1 ):
            sample, face_index = sssb_mesh.get_sample_from_surface()
            sample = sample[0]
        else:
            face_index = face_index0
            sample = sssb_mesh.get_sample_from_face(face_index, np.random.uniform)



        

        for j in range(0,3): 
            particle_array[i,j] = sample[j]+normals[face_index,j]*h_init
            particle_array[i,3+j] = normals[face_index,j]*v_init


        particle_array[i,7] = fabs(np.random.normal(R_mu,R_sigma))
        particle_array[i,6] = (3.0/4.0)/((fabs(particle_array[i,7]))*rho_dust)
        particle_array[i,6] = particle_array[i,6]*0.5*C_d*molar_mass

    
    return face_index

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef void run_sim( output_fname,
                    mesh_fname,
                    dmap_fname,
                    rot_mat,
                    comet_location,
                    t_double_T mu_orbit,
                    t_double_T m_comet,              #kg
                    t_double_T time_step,
                    t_double_T mesh_scale,
                    t_int_T n_particles,
                    cnp.ndarray[t_double_T,ndim=1] omega,
                    t_double_T C_d,
                    t_double_T R_mu,     #km
                    t_double_T R_sigma,   #km
                    t_double_T rho_dust,           #kg/[km^3] 
                    t_double_T molar_mass,    #kg/mol
                    cnp.ndarray[t_double_T,ndim=1] border,
                    t_double_T normalization_factor = 1.0
                    ###include when attitude control is included
                    #cnp.ndarray[t_double_T,ndim=1] dmap_lcorner_t, 
                    #cnp.ndarray[t_double_T,ndim=1] dmap_ucorner_t,
                    #cnp.ndarray[t_int_T,ndim=1] nDims_t, 
                    ):


    print("n_particles",n_particles)
    #print(molar_mass,rho_dust,R_mu,R_sigma)
    #return

    ##load gas density map
    cdef bytes py_bytes = dmap_fname.encode()
    cdef char* fname_in = py_bytes
    load_from_file2(<char*> fname_in)

    #Buffer used to transfer data between python interface and c code
    #change later is necessary
    cdef cnp.ndarray[t_double_T, ndim=1] buffer = np.empty((60),dtype=t_double)

    #cdef cnp.ndarray[t_double_T, ndim=1] lower_corner  = np.empty((3),dtype=t_double)
    #cdef cnp.ndarray[t_double_T, ndim=1] dimensions = np.empty((3),dtype=t_double)
    #cdef cnp.ndarray[t_double_T, ndim=1] sun_loc = np.empty((3),dtype=t_double)
    #cdef cnp.ndarray[t_double_T, ndim=1] omega = np.empty((3),dtype=t_double)

    cdef t_double_T alpha

    ##Comet position
    cdef cnp.ndarray[t_double_T,ndim=1] xyz = np.empty((3),dtype=t_double)
    cdef cnp.ndarray[t_double_T,ndim=1] vel
    cdef cnp.ndarray[t_double_T,ndim=1] n

    vec = comet_location 
    rot = rot_mat
    

   

    sssb_mesh = mw.MeshObject(mesh_fname, rot, scale = mesh_scale)
    cdef cnp.ndarray[t_double_T,ndim=2] normals 
    cdef cnp.ndarray[t_double_T,ndim=2] areas 
    cdef cnp.ndarray[t_double_T,ndim=2] centroids
    normals,areas,centroids = sssb_mesh.get_face_information(t_double)

    ##Particle data
    cdef cnp.ndarray[t_double_T,ndim=2] particle_data = np.zeros((1, 8),dtype=t_double)
    cdef cnp.ndarray[t_float_T,ndim=1] dist = np.array((0,0,0),np.float32)

    cdef cnp.ndarray[t_float_T,ndim=1] cached_data = np.empty((15),np.float32)
    cdef cnp.ndarray[t_float_T,ndim=1] cached_data_prev = np.empty((15),np.float32)

    cdef cnp.ndarray ind_cache = np.array((0,0,0,0),np.uintc)
    cdef cnp.ndarray ind_cache_prev = np.array((0,0,0,0),np.uintc)

    cdef cnp.ndarray[t_float_T,ndim=1] x_vec = np.array((0,0,0),np.float32)
    cdef cnp.ndarray[t_float_T,ndim=1] k_vec = np.array((0,0,0),np.float32)

    cdef t_int_T i, j, i1, i2, particle_number
    cdef t_double_T time 
    
    cdef t_double_T h_init = fabs(np.random.normal(0.025,0.025)) #km
    cdef t_double_T v_init = 0#fabs(np.random.normal(0.0,0.1)) #km/s

    
    buffer[29+15+0] = omega[0]
    buffer[29+15+1] = omega[1]
    buffer[29+15+2] = omega[2]
    cdef t_int_T counter = 0
    cdef t_int_T create_new = 0




    xyz[0] = -vec[0]*0.000001
    xyz[1] = -vec[1]*0.000001  
    xyz[2] = -vec[2]*0.000001  

    print("vec pos:",xyz)

    r2 = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]) 
    n = np.empty((3),t_double)
    for j in range(0,3):
        n[j] = xyz[j]/r2

    xyz[0] = -vec[0]
    xyz[1] = -vec[1] 
    xyz[2] = -vec[2]


    
    scene = rtcs.EmbreeScene()
    TriangleMesh(scene, sssb_mesh.mesh.vertices, sssb_mesh.mesh.faces)

    cdef t_double_T angle
    cdef t_int_T status, normal_index
    cdef t_double_T time_step_mod

    cdef cnp.ndarray[t_float_T,ndim=2] rt_pos = np.empty((1,3),t_float)
    cdef cnp.ndarray[t_float_T,ndim=1] prev_pos = np.empty((3),t_float)


    cdef cnp.ndarray[t_float_T,ndim=2] rt_dir = np.empty((1,3),t_float)
    
    cdef t_int_T particle_count = 0
    cdef t_int_T segments = 6

    cdef t_int_T jet_starters = 0

    cdef t_int_T max_jet_starters = 10
    jet_starter =  np.empty((max_jet_starters),t_float) 

    '''
    #sample jets
    while(jet_starters<max_jet_starters):
        ##Generate location (just use our old script)
        normal_index = generate_particles(0,1,
                            h_init,v_init,
                            R_mu,
                            R_sigma,
                            C_d,
                            rho_dust,
                            molar_mass,
                            normals,areas,centroids,
                            -1,particle_data, sssb_mesh)
        #check position    
        rt_pos[0,0]=particle_data[0,0]
        rt_pos[0,1]=particle_data[0,1]
        rt_pos[0,2]=particle_data[0,2]

        rt_dir[0,0]=n[0]
        rt_dir[0,1]=n[1]
        rt_dir[0,2]=n[2]      
        res = scene.run(rt_pos, rt_dir)

        #is it dark side?
        if(res[0]>=0):
            continue
        jet_starter[jet_starters] = normal_index
        jet_starters=jet_starters+1
    ''' 

    particle_number = 0
    while(particle_number < n_particles):
        #print("particle:",i)
        ind_cache[2] = 0
        h_init = fabs(np.random.normal(0.025,0.025))
        ##Generate particles
        normal_index = generate_particles(0,1,
                            0.001,v_init,
                            R_mu,
                            R_sigma,
                            C_d,
                            rho_dust,
                            molar_mass,
                            normals,areas,centroids,
                            -1,particle_data, sssb_mesh)
                            #np.random.choice(jet_starter,1)[0],
                            #particle_data, sssb_mesh)


        rt_pos[0,0]=particle_data[0,0]
        rt_pos[0,1]=particle_data[0,1]
        rt_pos[0,2]=particle_data[0,2]

        rt_dir[0,0]=n[0]
        rt_dir[0,1]=n[1]
        rt_dir[0,2]=n[2]        
        res = scene.run(rt_pos, rt_dir)
        if(res[0]>=0):
            continue
        #print(normal_index)
        particle_number = particle_number+1
        for j in range(0,3):
            particle_data[0,j] = particle_data[0,j]+normals[normal_index,j]*h_init
            buffer[29+12+j] = xyz[j]*0.001

        #print("{} {} {}".format(rt_pos[0,0],rt_pos[0,1],rt_pos[0,2]))
        
        counter = 0
        create_new = 0
        #print("new")
        #print("out",particle_data[0,0:7],i)
        while True:
            for j in range(0,3):
                prev_pos[j] = particle_data[0,j]

        
            time_step_mod = time_step#1# time_step*(0.5+0.5*np.random.rand())
            counter = counter+1
            status = update_particles(<double*>  buffer.data,
                                <double*>  particle_data.data, 
                                time_step_mod,
                                1, 
                                m_comet,
                                mu_orbit,
                                <float*> cached_data.data,
                                border[0],
                                <unsigned int*> ind_cache.data)    

            if(status<0):
                break;

            for j in range(0,3):
                dist[j] = particle_data[0,j]
                if(fabs(dist[j])>=fabs(border[j])*0.96):
                    #print("reached border",dist[j])
                    create_new = 1
                    break

            if create_new == 1:
                break
            
            particle_count = particle_count+1
   
            for i1 in range(0,3):
                k_vec[i1] = (dist[i1] - prev_pos[i1])/(segments-1.0)
                x_vec[i1] = prev_pos[i1]

            for i2 in range(1,segments):
                for i1 in range(0,3):
                    dist[i1] = x_vec[i1]+k_vec[i1]*i2

                #C_eff approx, change to correct one as soon as possible
                create_new = add_data_to_cached2(1.0/(segments-1.0)*particle_data[0,7]*particle_data[0,7]*1.0e18*0.3*0.3, 
                    <float*> dist.data,  
                    <unsigned int*> ind_cache_prev.data);

            for i1 in range(0,3):
                cached_data_prev[i1] = dist[i1]

            for i1 in range(0,3):
                ind_cache_prev[i1] = ind_cache[i1]

            if(counter>int(20*3600.0/time_step)):
                break
            
    normalize(normalization_factor/particle_count,4,5)

    cdef bytes py_bytes2 = output_fname.encode()
    cdef char* fname_out = py_bytes2
    save_to_file2(<char*> fname_out)
 
