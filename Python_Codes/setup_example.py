from dolfin import *
from mshr import *
import numpy as np

class Example(object):
    pass

class Example3D1(Example):
    def __init__(self):
        

        self.T = 12
        self.N = 500
        filename = "cube3.xdmf" 
  
        gdim = 3
        
        ###Mr:
        mymesh = Mesh()
        xdmf_mesh = XDMFFile(filename)
        xdmf_mesh.read(mymesh)    
        self.mesh = mymesh
        ###       
        class Exits(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary 
        self.exits = Exits()
        
        #spot on right
        self.rho_0 = Expression(('0.5*exp(-(pow(x[0]+0.9375,2) + pow(x[1]-0.0,2) + pow(x[2]-1.0,2))/1) \
        '), degree=2)  
        
        #spot on top
        '''
        self.rho_0 = Expression(('0.5*exp(-(pow(x[0]+5.0,2) + pow(x[1]-0.5,2) + pow(x[2]-0.0,2))/1) \
        '), degree=2)  
        '''
 	
 	#another spot on top
        '''
        self.rho_0 = Expression(('0.5*exp(-(pow(x[0]+5.0,2) + pow(x[1]+0.75,2) + pow(x[2]+0.125,2))/0.1) \
        '), degree=2)  
        '''
        
        self.ag_pos_0 = np.array([])
        


