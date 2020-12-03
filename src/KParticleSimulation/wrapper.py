import os
import numpy as np


from src.common import codes as JCCODES
from src.common import misc as misc
from multiprocessing import Process


def run_program(parameters, job_size, output_file, normalization):
    import particle_sim
    omega = np.array(parameters["omega"])
    border = np.array(parameters["border"])
    particle_sim.run_sim(output_file,
                    parameters[JCCODES.MESH],
                    parameters["gas_field"],
                    parameters["sim_common"]["rot_mat"],
                    parameters["sim_common"]["comet_location"],
                    parameters["mu_orbit"],
                    parameters["m_comet"],
                    parameters["time_step"],
                    parameters["mesh_scale"],
                    job_size,
                    omega,
                    parameters["C_d"],
                    parameters["R_mu"],
                    parameters["R_sigma"],
                    parameters["rho_dust"],
                    parameters["molar_mass"],
                    border,
                    normalization)


def run_combine(output_path, file_list):
    import combine
    combine.run(output_path, file_list)


def run(parameters=None, empty=None, logger=None):
    path = os.path.dirname(os.path.abspath(__file__))
    fname = path+"/defaults.in"
    defaults = misc.read_json(fname, "Reading defaults from: {}", logger)
    if type(parameters) is str:
        parameters = misc.read_json(parameters, "Reading input from: {}", logger)
    elif parameters is None:
        logger.error("Give some parameters (dict or string (path))")
    defaults.update(parameters)
    parameters = defaults

    n_processes = int(parameters[JCCODES.NPROCESSES])
    n_particles = parameters["n_particles"]
    job_size = np.ceil(n_particles*1.0/n_processes)
    pids = []
    logger.info("Running {} processes".format(n_processes))
    output_files = []

    for i in range(0, n_processes):
        postfix = "_{}".format(i)
        output_files.append("{}{}".format(parameters[JCCODES.OUTPUTTO], postfix))
        normalization = 1.0/n_processes
        p = Process(target=run_program, args=(parameters, job_size, output_files[-1], normalization))
        logger.debug(output_files[-1])
        p.start()
        pids.append(p)

    sleep_time = 5.0
    logger.info("waiting processes to finish: start polling")
    misc.wait_processes(pids, sleep_time, logger)
    logger.info("particle simulation finished, combining results")
    arg1 = parameters[JCCODES.OUTPUTTO]
    arg2 = output_files
    p = Process(target=run_combine, args=(arg1, arg2))
    p.start()
    pids.append(p)
    misc.wait_processes(pids, sleep_time, logger)
    logger.info("results combined")







'''
if __name__ == "__main__":
    parameters = {}
    if(len(sys.argv)==3):
        parameters[JCCODES.MESH] = sys.argv[0]
        parameters[JCCODES.OUTPUT_TO] = sys.argv[1] 
        parameters[JCCODES.PROCLOC] = "."
    elif(len(sys.argv)==2):
        parameters = sys.argv[1]
    run(parameters)
'''
