#ifndef PARTICLE_GUARD
#define PARTICLE_GUARD 


#include <math.h>
#include <stdio.h>
#include "cached_array.h"


typedef double t_double;
typedef unsigned int uint;
typedef float t_float;

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

//As the name says, surprise surprise cross prodcut
static void cross_product(t_double* vect_A,
                         t_double* vect_B,
                         t_double* cross_P){
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]; 
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]; 
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]; 
}

//write some better interpolation
static t_double get_gas_information(t_double* p1, 
                                t_double* vgas,  
                                t_float* cached_data,
                                uint* indStore,
                                t_double border){

    t_float pos[3];
    pos[0] = p1[0];
    pos[1] = p1[1];
    pos[2] = p1[2]; 

    for(int i=0;i<3;++i){
        if(fabs(pos[i])>=border*0.98){
            vgas[0] = 0.0;
            vgas[1] = 0.0;
            vgas[2] = 0.0;
            return -1.0;
        }
    }
    
    t_int status = readDataFromPos(cached_data+6,0,4,pos,indStore);
    if(status==1){
        return -1.0;
    }
    t_float ret_value = (t_double)cached_data[6];
    for(int i=0;i<3;++i){
        if(ret_value>0.000000001){
            vgas[i] = (t_double)(cached_data[7+i]/ret_value);
        }else{
            vgas[i] = 0.0;
        }
    //    printf("%f ",vgas[i]);
    }
    
    //printf("%f \n",ret_value);
    //printf("%f %f %f %f \n",ret_value,vgas[0],vgas[1],vgas[2]);
    return ret_value;
}



static int acceleration(t_double* p1, 
                        t_double* v1,  
                        t_double alpha, 
                        t_double* outcome,
                        t_double* buffer,
                        t_double mass_comet,
                        t_double mu_orbit,
                        t_float* cached_data,
                        uint* indStore,
                        t_double border,
                        t_double rad){
        
    t_double* v_gas = buffer;
    t_double* tmp1 = buffer+3;
    t_double* tmp2 = buffer+6;
    t_double* gravity = buffer+9;
    t_double* sun_loc = buffer+12;
    t_double* omega = buffer+15;

    t_double mu = 6.67408e-11*1.0e-9;
    t_double Ngas = get_gas_information(p1,v_gas, cached_data,indStore,border);
    if(Ngas<0.0){
        return -1;
    }
    //v_gas[0] = 0;
    //v_gas[1] = 0;
    //v_gas[2] = 0;
    //printf("in 1\n");

    //coriolis
    cross_product(omega,p1,tmp1);
    cross_product(omega,tmp1,tmp2);
    //printf("omega: %f %f %f \n",omega[0],omega[1],omega[2]);
    //centrifugal
    cross_product(omega,v1,tmp1);
    

    
    //approx vec from comet to sun == vec from particle to sun 
    t_double dist_to_sun = sqrt(sun_loc[0]*sun_loc[0]+sun_loc[1]*sun_loc[1]+sun_loc[2]*sun_loc[2]);
    t_double dist_to_comet = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);


    //Add better gravitation model later. Comet is assumed to have symmetrical
    //gravitational field
    for(int i = 0; i<3; ++i){
        gravity[i]=sun_loc[i]/(pow(dist_to_sun,3))*mu_orbit-p1[i]/(pow(dist_to_comet,3))*mu*mass_comet;
    }

    t_double vel_diff[3];
    for(int i=0;i<3;++i){
        vel_diff[i] = (v_gas[i]-v1[i]);
    }


    t_double vel_diff_abs = sqrt(vel_diff[0]*vel_diff[0]+vel_diff[1]*vel_diff[1]+vel_diff[2]*vel_diff[2]);
    t_double Q_eff = 1; //change to better
    t_float E0 = 3.83e26*0.001*0.001;
    for(int i=0;i<3;++i){
        outcome[i] = (vel_diff[i]*vel_diff_abs*alpha*Ngas+gravity[i]-2*tmp1[i]-tmp2[i])-sun_loc[i]*pow(rad,2)/4.0*Q_eff/(300000.0)*(E0/(4.0*pow(dist_to_sun,3)));
        //outcome[i] = dx[i]-2*tmp1[i]-tmp2[i];
    }
    
    return 0;
}



//4th order runge kutta
static int update_particles(t_double* buffer,
                            t_double* particle_data, 
                            t_double h, 
                            long nParticles,
                            t_double m_comet,
                            t_double mu_orbit,
                            t_float* cached_data,
                            t_double border,
                            uint* indStore){

    t_double* p1 = buffer;
    t_double* v1 = buffer+3;

    t_double* l1 = buffer+6;
    t_double* l2 = buffer+9;
    t_double* l3 = buffer+12;
    t_double* l4 = buffer+15;

    t_double* k1 = buffer+18;
    t_double* k2 = buffer+21;
    t_double* k3 = buffer+24;
    t_double* k4 = buffer+27;
    t_double rad = particle_data[7]; 
    //printf("in 2\n");
     //3 pos + 3 dir + 2 additional
    int particle_data_size = 8;
    int status = 0;
    //Go through all the particles
    for(int indx = 0; indx<nParticles; ++indx){
        //Go through all the time steps during one frame
        for(int j=0; j<1;++j){
            
            t_double* p0 = particle_data + indx*particle_data_size;
            t_double* v0 = particle_data + indx*particle_data_size+3;
            //printf("p1: %f %f %f %f \n",p0[0],p0[1],p0[2], h);
            //printf("v1: %f %f %f \n",v0[0],v0[1],v0[2]);
            for(int k=0;k<3;++k){
                p1[k] = p0[k];
                v1[k] = v0[k];
            }

            t_double alpha = particle_data[indx*particle_data_size+6];
            status = acceleration(p1, v1, alpha, k1, 
                            buffer+29,m_comet,mu_orbit,cached_data,indStore, border, rad);
            if(status<0) return -1;
            for(int k=0;k<3;++k){
                //printf("k1: %f \n",k1[k]);
                k1[k]=k1[k]*h;
                l1[k] = v0[k]*h;
            }            



            for(int k=0;k<3;++k){
                p1[k] = p0[k]+l1[k]*0.5;
                v1[k] = v0[k]+k1[k]*0.5;           
            }

            status =  acceleration(p1, v1, alpha, k2, 
                            buffer+29,m_comet,mu_orbit, cached_data,indStore, border, rad);
            if(status<0) return -1;
            for(int k=0;k<3;++k){
                k2[k]=k2[k]*h;
                l2[k] = (v0[k]+k1[k]*0.5)*h;
            }            

            


            for(int k=0;k<3;++k){
                
                p1[k] = p0[k]+l2[k]*0.5;
                v1[k] = v0[k]+k2[k]*0.5;
            }

            
            status =  acceleration(p1, v1, alpha, k3,
                            buffer+29,m_comet,mu_orbit, cached_data,indStore, border, rad);
            if(status<0) return -1;

            for(int k=0;k<3;++k){
                k3[k]=k3[k]*h;
                l3[k] = (v0[k]+k2[k]*0.5)*h;
            }      

            for(int k=0;k<3;++k){
                p1[k] = p0[k]+l3[k];
                v1[k] = v0[k]+k3[k];    
            }

            status =  acceleration(p1, v1,  alpha, k4, 
                            buffer+29,m_comet,mu_orbit, cached_data,indStore, border, rad);
            if(status<0) return -1;
            for(int k=0;k<3;++k){
                k4[k]=k4[k]*h;
                l4[k] = (v0[k]+k3[k]*0.5)*h;
            }      




            for(int k=0;k<3;++k){
                particle_data[indx*particle_data_size+k] = 
                        p0[k]+1.0/6.0*(l1[k]+2*l2[k]+2*l3[k]+l4[k]);
                particle_data[indx*particle_data_size+3+k] = 
                        v0[k]+1.0/6.0*(k1[k]+2*k2[k]+2*k3[k]+k4[k]);
            }
            
        }
    }
    return 0;
}

void normalize(float factor, uint start, uint end);

void save_to_file2(char* fname);

static void load_from_file2(char* fname);

int add_data_to_cached2(float value, 
                        float* pos, 
                        uint* indStore);

#endif
