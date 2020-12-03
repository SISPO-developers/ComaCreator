import numpy as np
from src.common import meshWrapper as mw
from noise import pnoise3

# If we need to change precision later
t_float = np.float32

# Generates gas source strengths and velocities
# TO DO: study what to do with random generation
def execute(mesh_file,
            output_to,
            gas_emission_base,
            gas_emission_scale,
            pnoise_scale,
            pnoise_offset,
            octaves,
            persistence,
            m_gas,
            T_gas,
            logger,
            chance_to_strong_jets=0.05,
            jets_scale = 3.0,
            render_mode=False):

    N_a = 6.02214086       # Avogadro number [Si unit, mol dropped]
    k_b = 1.38064852       # Boltzman constant [Si unit, mol dropped]

    #Root-Mean-Square velocities of gas particles
    v_th = np.sqrt(3*(k_b*N_a)*T_gas/m_gas)*0.001   #velocity [km/s]

    logger.info("Start gas source generation")
    obj = mw.MeshObject(mesh_file, logger=logger)
    faces, areas, centroids = obj.get_face_information(t_float)
    data = np.zeros((centroids.shape[0], 2), t_float)
    indx = 0
    max_value = 0
    for centroid in centroids:
        x = centroid[0]
        y = centroid[1]
        z = centroid[2]
        v = pnoise3(x * pnoise_scale,
                    y * pnoise_scale,
                    z * pnoise_scale,
                    octaves=octaves,
                    persistence=persistence)
        max_value = max(max_value, v)
        if(np.random.uniform()<chance_to_strong_jets):
            v = v*(1+np.random.normal()*jets_scale)
        data[indx, 0] = v
        indx = indx+1

    min_value = np.amin(data[:, 0])
    data[:, 0] = (data[:, 0]+np.abs(min_value))
    data[:, 0] = data[:, 0]/max_value*gas_emission_scale + gas_emission_base  # surface gas emission rates [mol/(s*km^2)]
    data[:, 1] = v_th                                     # gas source initial velocity [km/s]
    np.savetxt(output_to, data)          # save output
    logger.info("Finished, data saved to: {}".format(output_to))
    # Visualization 
    if(render_mode):
        max_color = np.amax(data[:,0])
        obj.change_to_face_color() 
        colors = np.zeros((centroids.shape[0],3),np.uint8)
        colors[:,0] = data[:,0]/max_color*255
        colors[:,1] = data[:,0]/max_color*255
        colors[:,2] = data[:,0]/max_color*255
        
        obj.mesh.visual.face_colors = colors 
        render(obj);

def render(obj):
    import pyrender, os
    
    m = pyrender.Mesh.from_trimesh(obj.mesh,smooth=False)
    
    scene = pyrender.Scene(ambient_light=[0.1,0.1,0.1,0.1])
    light = pyrender.PointLight(intensity = 500)
    scene.add(m)
    #scene.add(light)
    pyrender.Viewer(scene)
   
