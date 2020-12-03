import json
import numpy as np
RAD_CONVERSION = np.pi/180.0
data = {
    "mu_orbit" : 1.32712442099e+20*1.0e-9,     #notice that km   
    "time_step" : 1.0,                              #s
    "mesh_scale" : 1.0,                             #
    "n_particles" : 100000,                         #
    "omega" : [0.0,0.0,np.pi*2/(12.3*3600.0)],              #rad/s           
    "m_comet" : 10.0e13,                            #comet mass [kg]  
    "C_d" : 2.0,                                    #
    "R_mu" : 300.0*1.0e-6*1.0e-3,                   #km                    
    "R_sigma" : 60.0*1.0e-6*1.0e-3,                 #km
    "rho_dust" : 440.0*1.0e9,                        #kg/[km^3]
    "molar_mass" : 18.01528*0.001,                  #kg/mol
    "border" : [20.0,20.0,20.0],                    #km
    "output_to" : "test.out",
    "mesh_fname" : "comet.obj",
    "dmap_fname" : "dmap",
}





with open('defaults.in', 'w') as outfile:
    json.dump(data, outfile, indent=2)
