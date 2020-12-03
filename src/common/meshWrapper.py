import trimesh 
import numpy as np
import mathutils
import math

def area(v1, v2, v0): 
    v = np.cross(v1-v0,v2-v0)
    return math.sqrt(v[0]**2+v[1]**2+v[2]**2)*0.5

class MeshObject:

    def __init__(self, fname, rot=None, scale=1, logger=None):
        if(logger is not None):
            logger.info(f"filename {fname}")
        self.mesh = trimesh.load_mesh(fname)

       
        
        ##Trimesh has different coordinate system, so the mesh
        ##needs to be rotated.
        #mat_rot = 
        mat_scale = mathutils.Matrix.Scale(scale,4)
        self.mesh.apply_transform(mat_scale)
        #self.mesh.apply_transform(mat_rot)

        if rot is not None:
            rot2 = np.zeros((4,4))
            rot2[0,:] = rot[0]
            rot2[1,:] = rot[1]
            rot2[2,:] = rot[2]
            rot2[3,:] = rot[3]

            #rot2[:, [2, 1]] = rot2[:, [1, 2]] 
            #rot2[[2, 1], :] = rot2[[1, 2],:] 

            #rot2[2,:] = -rot2[2,:]
            #rot2[:,2] = -rot2[:,2]
            #rot2[2,2] = -rot2[2,2]
            #r1 = mathutils.Matrix.Rotation(math.radians(180.0), 4, 'Z')
            #r2 = mathutils.Matrix.Rotation(math.radians(-90.0), 4, 'Y')

            
            #self.mesh.apply_transform(rot2)
            #self.mesh.apply_transform(r1)

        self.mesh.export('data/test.obj')




    def get_face_information(self, t_float):
        areas = self.mesh.area_faces.astype(t_float).reshape((-1,1))
        face_normals = self.mesh.face_normals
        centroids = self.mesh.triangles_center
        return face_normals.astype(t_float), areas, centroids.astype(t_float)
       
    def get_vertices(self, t_float):
        return self.mesh.vertices.astype(t_float)

    def get_faces(self, t_float):
        return self.mesh.faces.astype(t_float)


    def get_sample_from_surface(self):
         return trimesh.sample.sample_surface_even(self.mesh, 1, radius=None)
       
    def n_faces(self):
        return self.mesh.faces.shape[0]

    def change_to_face_color(self):
        self.mesh.visual = trimesh.visual.color.ColorVisuals()

    def get_sample_from_face(self, index, rand):
        while(True):
            v1=self.mesh.vertices[self.mesh.faces[index,0]]
            v2=self.mesh.vertices[self.mesh.faces[index,1]]
            v3=self.mesh.vertices[self.mesh.faces[index,2]]
            v2 = v2-v1
            v3 = v3-v1
            x = v2*rand()+v3*rand()
            collective_area = area(x,v2,np.zeros(3))+area(x,v3,np.zeros(3))+area(x,v2,v3)
            triangle_area = area(v2,v3,np.zeros(3))    
            #print(collective_area,triangle_area)    
            if(triangle_area*(1-0.001)<collective_area and triangle_area*(1+0.001)>collective_area):
                return v1+x
                


