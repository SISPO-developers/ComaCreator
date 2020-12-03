import json
import numpy as np
RAD_CONVERSION = np.pi/180.0
data = {
    "gas_emission_base" : 215.87,         #mol/(s km^2)
    "gas_emission_scale" : 1162.4,              #mol/(s km^2)
    "pnoise_scale" : 0.0025,     
    "pnoise_offset" : 0.0, 
    "octaves" : 10,
    "persistence" : 10.9,
    "T_gas" : 210.0,            #K
    "m_gas" : 18.01528*0.001,    #kg/mol
    "mesh_file" : "/media/rokka/main2020/WRK/JetCreator/data/comet.obj",
    "write_output_to" : "/media/rokka/main2020/WRK/JetCreator/data/test.dat"
}


with open('defaults.in', 'w') as outfile:
    json.dump(data, outfile, indent=2)
