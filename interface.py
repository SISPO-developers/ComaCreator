from pathlib import Path
import importlib
import datetime
import copy
import logging
import os
import numpy as np
import mathutils

import src.common.misc as misc
import src.common.codes as JCCODES
import src.blendertools as bpy_tools

from multiprocessing import Process




class JetCreatorInterface:
    def parse_paths(self, data):
        root = Path(__file__).resolve().parent
        self.root = root
        ROOTPATH = "ROOTPATH"
        if type(data) is dict:
            for key in data:
                if type(data[key]) is str:
                    data[key] = data[key].replace(ROOTPATH, str(root))
                else:
                    self.parse_paths(data[key])
        elif type(data) is list:
            for indx, key in enumerate(data):
                if type(key) is str:
                    data[indx] = data[indx].replace(ROOTPATH, str(root))
                else:
                    self.parse_paths(key)

    # Init
    def __init__(self, input_data=None, defaults=None, logger=None):
        if defaults is None:
            defaults = str(Path(__file__).resolve().parent)+"/data/defaults.in"
        if logger is None:
            logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        # read default values
        data = misc.read_json(defaults, "Reading defaults from: {}", self.logger)
        # read given input and override values from default

        if input_data is not None:
            if type(input_data) is str:
                data2 = misc.read_json(input_data, "Reading input from: {}", self.logger)
                data.update(data2)
            elif type(input_data) is dict:
                data.update(input_data)
                
        self.parse_paths(data)

        self.data = data["JetCreator"]
        self.data_h = data

    @staticmethod
    def fill_input_for_task(task, ID, task_indx, task_history, full_data):
        task[JCCODES.EXECUTIONID] = ID
        if task[JCCODES.MESH] == JCCODES.DEFAULT:
            task[JCCODES.MESH] = full_data["comet"][JCCODES.MESH]

        '''
        # if(task[JCCODES.OUTPUTFOLDER]==JCCODES.DEFAULT):
        #    task[JCCODES.OUTPUTFOLDER] = full_data["JetCreator"][JCCODES.OUTPUTFOLDER]
        '''
        if task[JCCODES.OUTPUTTO] == JCCODES.DEFAULT:
            output = "{}{}.out".format(task[JCCODES.OUTPUTFOLDER], ID)
            task[JCCODES.OUTPUTTO] = output
        elif task[JCCODES.OUTPUTTO] == JCCODES.FINAL:
            task[JCCODES.OUTPUTTO] = full_data["final_output"]

        insert_me = []
        for key in task:
            if task[key] == JCCODES.PREVOUTPUT:
                task[key] = task_history[task_indx - 1][JCCODES.OUTPUTTO]
            if key == JCCODES.INCLUDE:
                for include_me in task[key]:
                    insert_me.append([include_me, full_data[include_me]])
        for entry in insert_me:
            task[entry[0]] = entry[1]

    # Execute entire script
    def execute_pipeline(self):
        task_history = []
        ID = None
        if self.data["ID"] == JCCODES.EMPTYID:
            self.logger.info("ID was not given, generate timestamp ID")
            ID = generate_ID()
        else:
            ID = self.data["ID"]
        self.logger.info("ID: {}".format(ID))

        # Extract data from the input
        for task_indx, line in enumerate(self.data["pipeline"]):
            self.logger.info("Running task: {}".format(line))
            task = copy.deepcopy(self.data[JCCODES.COMMON])
            task.update(self.data[JCCODES.SIMCOMMON])
            task.update(self.data[line])
            ID_gen = "{}_{}".format(ID, task_indx)
            self.fill_input_for_task(task, ID_gen, task_indx, task_history, self.data_h)
            program = self.data[JCCODES.PROGRAMID][task[JCCODES.PROGRAMID]]
            tmp_input = "{}{}.in".format(task[JCCODES.TMPFOLDER], ID_gen)
            misc.save_json(tmp_input, task, "Saving tmp input to: {}", self.logger)
            self.execute_task(task, program, tmp_input, self.logger)
            task_history.append(task)
        return self.data_h["final_output"]


    # Execute task
    def execute_task(self, task, program, input_file, logger):
        logger.debug(program)
        logger.debug(task)
        if task[JCCODES.SPECIALCODE] == JCCODES.SKIP:
            logger.info("Skip code given, Skipping step")
            return
        my_module = importlib.import_module(program)
        my_module.run(task, input_file, logger)


def generate_ID():
    return datetime.datetime.now().strftime('%Y%m%d%H%M%S')




def run(settings, env):
    #modify parameters from sispo to the plugin 
    input_data = {}
    input_data["comet"] = {}
    input_data["sim_common"] = {}
    
    blend_file_path = str(settings["simulation"]["sssb"]["model"]["file"])
    obj_name = settings["simulation"]["sssb"]["model"]["name"]
    mesh_fname = blend_file_path+".obj"
    input_data["comet"][JCCODES.MESH] = mesh_fname

    #Triangulate and also decimate
    comaCreator = settings.get("ComaCreator",None)
    decimateMesh = False
    decimatePercentage = 0.1
    if(comaCreator is not None):
        decimateMesh = comaCreator.get("decimate_mesh",decimateMesh)       
        decimatePercentage = comaCreator.get("decimate_percentage",decimatePercentage)
    

    rot = env.sssb.rot_history[0]

    sssb_axis = np.array(rot.getAxis(env.sssb.rot_conv).toArray())
    sssb_angle = rot.getAngle()
    M = mathutils.Matrix.Rotation(sssb_angle, 4, sssb_axis).to_4x4()
    comet_location = env.sssb.pos_history[0]

    #Process mesh in process because otherwise this will mess Sispo's rendering logic
    p = Process(target=bpy_tools.convert_blend_to_obj, args=(blend_file_path, obj_name, mesh_fname, decimateMesh, decimatePercentage, M)) 
    #wait meshing to finish
    p.start()
    p.join()
    
    java_arr = comet_location.toArray()
    input_data["sim_common"]["rot_mat"] = [list(row) for row in M]
    
    att = settings["simulation"]["sssb"].get("att",None)
    omega = np.array((0.0,0.0,0.0,0.0)) 
        
    if att is not None:
        rotation_rate = att.get("rotation_rate",0.0)
        if(abs(rotation_rate)>0):
            omega = mathutils.Vector((0.0, 0.0, 2.0*np.pi/att["rotation_rate"], 0.0))
            ##Could not get mathutils @ or * to work for some reason
            arr = np.array(M)
            vec = np.array(omega)
            omega = arr.dot(omega) 

    input_data["sim_common"]["omega"] = list(omega[:3])
    input_data["sim_common"]["comet_location"] = [float(java_arr[0]),float(java_arr[1]),float(java_arr[2])]

    myCreator = JetCreatorInterface(input_data)
    aux_info = myCreator.execute_pipeline()

    scenes = None

    filepath = str(Path(__file__).resolve().parent)+"/blend/tail.blend"
    scene_list = ["SssbConstDist","SssbOnly"] 
   

    env.renderer.load_object(filepath, "Tail", scenes=scene_list)

    bpy_tools.prepare_output_for_sispo("Tail",aux_info)



