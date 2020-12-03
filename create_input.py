import json
import numpy as np
import src.common.codes as JCCODES


RAD_CONVERSION = np.pi/180.0
data = {}
data['comet'] = {}
data["comet"][JCCODES.MESH] = "ROOTPATH/data/comet3.obj"

data["final_output"] = "ROOTPATH/data/out.json"


data['JetCreator'] = {}
data['JetCreator']["ID"] = "myID"
#JCCODES.EMPTYID

data['JetCreator'][JCCODES.PROGRAMID] = {
    "NSDSGenerator" : "src.NSDSGenerator.wrapper",
    "KGasSimulation" : "src.KGasSimulation.wrapper",
    "KParticleSimulation" : "src.KParticleSimulation.wrapper",
    "ToBlender" : "src.ToBlender.wrapper"
}

data["JetCreator"][JCCODES.COMMON] = {
    JCCODES.EXECUTIONID : "THIS_IS_HOLDER_DO_NOT_USE",
    JCCODES.OUTPUTTO : JCCODES.DEFAULT,
    JCCODES.OUTPUTFOLDER : "ROOTPATH/data/",
    JCCODES.NPROCESSES : 3,
    JCCODES.TMPFOLDER : "ROOTPATH/tmp/",
    JCCODES.SPECIALCODE : JCCODES.EMPTY
}


data["JetCreator"][JCCODES.SIMCOMMON] = {
    JCCODES.MESH : JCCODES.DEFAULT, 
    "mesh_scale" : 1.0, 
    "T_gas" : 210.0, 
    "m_gas" : 18.01528*0.001,
    "R_gas" : 4.6152e-04
}


data['JetCreator']["pipeline"] =  [
    "DefaultGasPreprocessor",
    "DefaultGasDensityCreator",
    "DefaultParticleGenerator",
    "DefaultToBlender",
    ]


data["JetCreator"]["DefaultGasPreprocessor"] = {
    JCCODES.PROGRAMID : "NSDSGenerator",
    "scale" : 762.4,
}

data["JetCreator"]["DefaultGasDensityCreator"] = {
    JCCODES.PROGRAMID : "KGasSimulation",
    "outflow_rates" : JCCODES.PREVOUTPUT,
    "nBounds" : 3, 
    "lower_corner_ind" : [10,10,10,10,10],
    "upper_corner_ind" : [10,10,10,10,10],
    "sample_points" : [40,40,40,100,200],
    "bbox_lower" : [-20.0,-20.0,-20.0],
    "bbox_upper" : [20.0,20.0,20.0],
    JCCODES.INCLUDE : ["trajectory"],
    JCCODES.SPECIALCODE : JCCODES.SKIP,
}


data["JetCreator"]["DefaultParticleGenerator"] = {
    JCCODES.PROGRAMID : "KParticleSimulation",
    "gas_field" : JCCODES.PREVOUTPUT,
    JCCODES.INCLUDE : ["trajectory"],
    "n_particles" : 100000,
    JCCODES.SPECIALCODE : JCCODES.SKIP,
}



data["JetCreator"]["DefaultToBlender"] = {
    JCCODES.PROGRAMID : "ToBlender",
    JCCODES.INCLUDE : ["trajectory"],
    JCCODES.OUTPUTTO : JCCODES.FINAL,
    "dust_particle_field" : JCCODES.PREVOUTPUT,
#    JCCODES.SPECIALCODE : JCCODES.SKIP,
}



data['trajectory'] = {}
data['trajectory']["orbit"] = {}
data['trajectory']["date"] = {}
data['trajectory']["encounter_date"] = {}
data['trajectory']["common"] = {}
data['trajectory']["orekit_data_path"] = "ROOTPATH/data/orekit-data.zip"

data['trajectory']["orbit"] = {
    "a" : 518060000000.0,               # semi-major axis [km]
    "e" : 0.64102,                      # eccentricity                   
    "omega" : 12.780 ,  # argument of perapsis
    "Omega" : 50.147 ,    # longitude of ascending node
    "i" : 7.0405 ,      # inclination
    "M" : 303.71 ,     # mean anomaly
    "t0" : 2456879.5,                   # at t0
}

data['trajectory']["date"] = {
    "year"          : 2014,    #simulation day
    "month"         : 8,    #simulation day
    "day"           : 10,    #simulation day
    "hour"          : 0,    #simulation day
    "minutes"        : 0,
    "seconds"       : 0,
}

data['trajectory']["encounter_date"] = {
    "year"          : 2015,    #simulation day
    "month"         : 4,    #simulation day
    "day"           : 21,    #simulation day
    "hour"          : 5,    #simulation day
    "minutes"        : 0,
    "seconds"       : 0,
}

data['trajectory']["att"] = {
    "rotation_rate"         : 0.00,        #Rotation period seconds
    "RA"       : 0.0,   #
    "Dec"      : 0.0
}







with open('./data/defaults.in', 'w') as outfile:
    json.dump(data, outfile, indent=2)
