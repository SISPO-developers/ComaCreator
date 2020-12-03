import src.common.misc as misc

#Export object from the Blender to be in the ComaCreator(obj file)
#Also is capable of decimating and triangulating the meshes 
def convert_blend_to_obj(blend_file_path, obj_name, mesh_fname, decimateMesh, decimatePercentage, M):
    import bpy
    bpy.ops.wm.open_mainfile(filepath=blend_file_path)
    objects = bpy.data.objects
    bpy.ops.object.mode_set(mode='OBJECT')
    for obj in objects:
        if(obj.name == obj_name):
            bpy.ops.object.select_all(action='DESELECT')
            obj.select_set(state=True)
            if decimateMesh:
                modifier=obj.modifiers.new("Decimate",'DECIMATE')
                modifier.ratio=decimatePercentage
                modifier.use_collapse_triangulate=True
            obj.matrix_world = M
            bpy.ops.export_scene.obj(filepath=mesh_fname, use_selection=True, use_materials=False,
                                     axis_forward="-Z", axis_up="Y") 





##Needs refactoring
def prepare_output_for_sispo(model_name, aux_info):
    import bpy
    obj = bpy.data.objects[model_name]

    data = None
    if type(aux_info) is dict:
        data = aux_info
    elif type(aux_info) is str:
        with open(aux_info) as json_file:
            data = misc.read_json(aux_info)

    obj.dimensions = data["dimensions"]

    #Set shader parameters
    mat = bpy.data.materials.get('volumeScatterer')

    tilingNode = mat.node_tree.nodes.get("resolution")
    tilingNode.outputs[0].default_value = data["resolution"]

    tilingNode = mat.node_tree.nodes.get("tiles_x")
    tilingNode.outputs[0].default_value = data["tiles_x"]

    tilingNode = mat.node_tree.nodes.get("tiles_y")
    tilingNode.outputs[0].default_value = data["tiles_y"]

    new_img = bpy.data.images.load(filepath = data["tiles_file"])
    tilingNode = mat.node_tree.nodes.get("tiles_file")
    tilingNode.image = new_img

