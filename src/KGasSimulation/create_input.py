import json
import numpy as np
RAD_CONVERSION = np.pi/180.0
data = {
    "T_gas" : 210.0,
    "R_gas" : 188.92*1.0e-6,
    "from_face_in" : 0,
    "to_face_in" : 9999999999999,
    "mesh_fname" : "comet.obj",     
    "outflow_rates" : "outflow_rates.in", 
    "nBounds" : 3, 
    "lower_corner_ind" : [10,10,10],
    "upper_corner_ind" : [10,10,10],
    "sample_points" : [40,100,400],
    "bbox_lower" : [-20.0,-20.0,-20.0],
    "bbox_upper" : [20.0,20.0,20.0],
    "output_to" : "test.out",
    "gas_flow_scale" : 1.0,
    "chunk_sizes" : [4,4,4]
}

with open('defaults.in', 'w') as outfile:
    json.dump(data, outfile, indent=2)
