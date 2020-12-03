import sys,os
from common import misc
from common import codes as JCCODES


def run_gas_sim(from_face, to_face, input_from=None, output_postfix=None):
    import gas_sim
    path = os.path.dirname(os.path.abspath(__file__))
    fname = path+"/defaults.in"
    defaults = misc.read_json(fname,txt="Reading defaults from: {}")
    if(input_from is not None):
        parameters = misc.read_json(input_from,txt="Reading input from: {}")
    defaults.update(parameters)
    parameters = defaults
    mesh_file = parameters[JCCODES.MESH]
    #trajectory_data = parameters[JCCODES.TRAJECTORY]


    if output_postfix is not None:
        parameters["save_output_to"] = parameters["save_output_to"]+output_postfix



    gas_sim.init_and_run(from_face,
                to_face,
                parameters["mesh_scale"],       ##scale of the mesh         [obj space to km]
                parameters[JCCODES.MESH],       ##read mesh from this file  
                parameters["outflow_rates"],    ##Outflow file
                parameters["T_gas"],            ##Temperature of the gas
                parameters["R_gas"],            ##Specific gas constant              
                parameters["nBounds"],                 ##How many boundaries
                parameters["lower_corner_ind"],
                parameters["upper_corner_ind"],
                parameters["sample_points"],
                parameters["bbox_lower"],
                parameters["bbox_upper"],
                parameters["save_output_to"],
                parameters["gas_flow_scale"])







def run_combine():
    import combine
    file_list = []
    output_path = sys.argv[2]
    for i in range(3,len(sys.argv)):
        file_list.append(sys.argv[i])
    combine.run(output_path,file_list)

if __name__ == "__main__":

    if(int(sys.argv[1])==0):
        input_from = sys.argv[2]
        from_face = int(sys.argv[3])
        to_face = int(sys.argv[4])
        output_postfix = sys.argv[5]
        run_gas_sim(from_face,to_face,input_from,output_postfix)
    else:
        run_combine()        









