import sys,os,math
from src.common import misc
from src.common import codes as JCCODES



def run(parameters,input_file=None,logger=None):
    import convert_to_exr as converter
    path = os.path.dirname(os.path.abspath(__file__))
    fname = path+"/defaults.in"
    defaults = misc.read_json(fname,txt="Reading defaults from: {}")
    if(input_file is not None):
        parameters = misc.read_json(input_file,txt="Reading input from: {}")
    defaults.update(parameters)
    parameters = defaults

    output_file = parameters[JCCODES.OUTPUTFOLDER]+parameters["output_dmap"]
    #aux_info = parameters[JCCODES.OUTPUTFOLDER]+parameters["aux_info"]

    print(parameters)
    converter.run(
                parameters["bbox_lcorner"],      
                parameters["bbox_ucorner"],    
                parameters["resolution"],   
                parameters["dust_particle_field"],     
                parameters["trajectory"],                  
                parameters[JCCODES.MESH],             
                output_file,
                parameters[JCCODES.OUTPUTTO],
                1)






