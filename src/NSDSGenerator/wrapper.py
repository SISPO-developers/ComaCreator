import os,sys,logging
from pathlib import Path

from src.common import codes as JCCODES
from src.common import misc
from src.NSDSGenerator import gen

 
def run(parameters=None, input_file=None, logger=None):
    path = str(Path(__file__).resolve().parent)
    fname = path+"/defaults.in"
    defaults = misc.read_json(fname, "Reading defaults from: {}", logger)
    if type(parameters) is str:
        parameters = misc.read_json(parameters, "Reading input from: {}", logger)
    elif parameters is None:
        logger.error("Give some parameters (dict or string (path))")
    defaults.update(parameters)
    parameters = defaults
    mesh_file = parameters[JCCODES.MESH]
    output_to = parameters[JCCODES.OUTPUTTO]
    gas_emission_base = parameters["gas_emission_base"]
    gas_emission_scale = parameters["gas_emission_scale"]
    pnoise_scale =  parameters["pnoise_scale"]
    pnoise_offset = parameters["pnoise_offset"]
    octaves = parameters["octaves"]
    persistence = parameters["persistence"]
    m_gas = parameters["m_gas"]
    T_gas = parameters["T_gas"]

    gen.execute(mesh_file,
                output_to,
                gas_emission_base,
                gas_emission_scale,
                pnoise_scale,
                pnoise_offset,
                octaves,
                persistence,
                m_gas,
                T_gas,
                logger)
    



if __name__ == "__main__":
    logging.basicConfig()
    logger = logging.getLogger('tcpserver')
    parameters = {}
    if(len(sys.argv)==3):
        parameters[JCCODES.MESH] = sys.argv[0]
        parameters[JCCODES.OUTPUT_TO] = sys.argv[1] 
        parameters[JCCODES.PROCLOC] = "."
    elif(len(sys.argv)==2):
        parameters = sys.argv[1]
    run(parameters,logger=logger)

