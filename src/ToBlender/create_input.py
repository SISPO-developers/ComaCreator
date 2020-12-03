import json
import numpy as np
RAD_CONVERSION = np.pi/180.0
data = {
    "bbox_lcorner" : [-20.0,-20.0,-20.0],  
    "bbox_ucorner" : [20.0,20.0,20.0],         
    "resolution" : 400,       
    "dust_particle_field" : "dmap2.in",                       
    "trajectory" : "trajectory",       
    "mesh_file" : "comet.obj",                          
    "output_dmap" : "output.exr",                               
    "aux_info" : "aux.info"
}


with open('defaults.in', 'w') as outfile:
    json.dump(data, outfile)
