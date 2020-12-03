import os
import math
import time
from multiprocessing import Process
from src.common import meshWrapper as mw
from src.common import codes as JCCODES
from src.common import misc


def run_gas_sim(from_face, to_face, input_from=None, output_postfix=None):
    import gas_sim
    path = os.path.dirname(os.path.abspath(__file__))
    fname = path+"/defaults.in"
    defaults = misc.read_json(fname,txt="Reading defaults from: {}")
    if input_from is not None:
        parameters = misc.read_json(input_from,txt="Reading input from: {}")
    defaults.update(parameters)
    parameters = defaults

    if output_postfix is not None:
        parameters[JCCODES.OUTPUTTO] = parameters[JCCODES.OUTPUTTO]+output_postfix

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
                parameters[JCCODES.OUTPUTTO],
                parameters["gas_flow_scale"],
                parameters["sim_common"]["rot_mat"],
                parameters["sim_common"]["comet_location"],
                parameters["chunk_sizes"])


def run_combine(output_path, file_list, logger):
    import combine
    combine.run(output_path, file_list, logger)


def run(parameters, input_file, logger):
    n_processes = int(parameters[JCCODES.NPROCESSES])
    #tmp_folder = parameters[JCCODES.TMPFOLDER]
    n_faces = mw.MeshObject(parameters[JCCODES.MESH]).n_faces()
    job_size = math.ceil(n_faces*1.0/n_processes)
    pids = []
    #path = os.path.dirname(os.path.abspath(__file__))

    logger.info("Running {} processes".format(n_processes))
    output_files = []
    start_face = 0

    '''
    for i in range(0,nProcesses):
        end_face = min(start_face+jobSize,nFaces)
        postfix = "_{}".format(i)
        output_files.append("{}{}".format(parameters[JCCODES.OUTPUTTO],postfix))
        ilist = ["python", path+"/run.py", str(0),input_file, str(start_face), str(end_face), postfix]
        print(ilist)
        pids.append(subprocess.Popen(ilist))        
        start_face = end_face
    '''

    for i in range(0, n_processes):
        end_face = min(start_face+job_size, n_faces)
        postfix = "_{}".format(i)
        i_list = [start_face, end_face,input_file,  postfix]
        output_files.append("{}{}".format(parameters[JCCODES.OUTPUTTO], postfix))
        p = Process(target=run_gas_sim, args=(i_list[0], i_list[1], i_list[2], i_list[3]))
        logger.debug(i_list)
        p.start()
        pids.append(p)
        start_face = end_face

    sleep_time = 5.0
    logger.info("waiting processes to finish: start polling")
    misc.wait_processes(pids, sleep_time, logger)
    logger.info("gas simulation finished, combining results")
    pids = []
    arg1 = parameters[JCCODES.OUTPUTTO]
    arg2 = output_files
    p = Process(target=run_combine, args=(arg1,arg2, logger))
    p.start()
    pids.append(p)
    misc.wait_processes(pids, sleep_time, logger)
    logger.info("results combined")



