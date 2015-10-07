import MDAnalysis
#from MDAnalysis.analysis import *
import numpy as np
import math

class h_bond_analysis:
    
    def __init__(self,psf_path, dcd_path, ndx_path):
        self.u = MDAnalysis.Universe(psf_path,dcd_path)
        self.ndx = self.ndx_file_read_to_list(ndx_path)
        self.max_oo_dist = 0.0
        self.max_oho_angle = 0.0
        self.max_oh_dist = 0.0
        self.hmat = self.abc_to_hmatrix(*self.u.dimensions)
        self.hinv = np.linalg.inv(self.hmat)
    def set_parameters(self, max_oo_dist, max_oh_dist, max_oho_angle):
        self.max_oh_dist = max_oh_dist
        self.max_oo_dist = max_oo_dist
        self.max_oho_angle = max_oho_angle
    def start_h_bond_analysis(self,oxygen_list, hydrogen_list):
        
        self.calc_si_all(self.hinv)
        
        for i in xrange(oxygen_list.n_atoms):
            i_index = oxygen_list.atoms.indices[i]
            for j in xrange(i+1,oxygen_list.n_atoms):             
                j_index = oxygen_list.atoms.indices[j]
                dist = self.calc_dist(i_index, j_index) 
                if dist <= self.max_oo_dist:
                    for h in xrange(hydrogen_list.n_atoms):
                        h_index = hydrogen_list.atoms.indices[h]
                        angle = self.calc_angle(i_index,h_index,j_index)
                        if angle >= self.max_oho_angle:
                            dist_ih = self.calc_dist(i_index,h_index)
                            dist_jh = self.calc_dist(j_index,h_index)
                            if dist_ih <= self.max_oh_dist or dist_jh <= self.max_oh_dist:
                                print dist, angle, dist_ih, dist_jh
                                print i_index, j_index, h_index                                         
        
    #form H matrix  (given that it's a triclinic cell) 
    # based on code by Guido https://github.com/mdtraj/mdtraj/issues/908
    def abc_to_hmatrix(self, a, b, c, alpha, beta, gamma, degrees=True):
        if degrees:
                alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
        result = np.zeros((3, 3))

        a = np.array((a, 0, 0))
        b = b * np.array((math.cos(gamma), math.sin(gamma), 0))
        bracket = (math.cos(beta) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
        c = c * np.array((math.cos(alpha), bracket, math.sqrt(math.sin(beta) ** 2 - bracket ** 2)))

        result[:, 0] = a
        result[:, 1] = b
        result[:, 2] = c

        return result
    #Calculates the reduced s_i for all atoms 
    def calc_si_all(self, hinv):
        self.s_coords = np.zeros((self.u.atoms.n_atoms,3))
        for i in xrange(self.u.atoms.n_atoms):
           self.s_coords[i] = np.dot(hinv, self.u.atoms.coordinates()[i]) 
    #Calculate inter atom distance 
    def calc_dist(self, i, j):
        s_ij = self.s_coords[i] - self.s_coords[j]
        s_ij = s_ij - np.round(s_ij)
        return np.linalg.norm(np.dot(self.hmat,s_ij))
    #Calculate angle in deg
    def calc_angle(self,i,j,k):
        s_hj = self.s_coords[i] - self.s_coords[j]
        s_hj = s_hj - np.round(s_hj)
        r_hj = np.dot(self.hmat,s_hj)
        
        s_hk = self.s_coords[k] - self.s_coords[j]
        s_hk = s_hk - np.round(s_hk)
        r_hk = np.dot(self.hmat,s_hk)
        
        angle = np.arccos(np.dot(r_hk, r_hj)/(np.linalg.norm(r_hk)*np.linalg.norm(r_hj)))
        angle = np.rad2deg(angle)
        return angle
    ##Reads a NDX file and stores each atom group in list
    def ndx_file_read_to_list(self, file_path):
        ndx={}
        with open(file_path) as f:
            for line in f:
                
                if line.startswith('['):
                    
                    name_index_str = line.split()[1]
                else:
                    id_list = line.split()
                    id_list = [int(i)-1 for i in id_list]
                    if id_list != []:
                        ndx[name_index_str] = id_list 
        return ndx
 
test = h_bond_analysis('/home/rasmus/project/rasmus/input.psf','/home/rasmus/project/rasmus/IOHMD-A-prod.dcd', '/home/rasmus/project/rasmus/input.ndx' )       
test.set_parameters(3.0,1.05, 155)
hydrogen_top = test.u.atoms[test.ndx['hydrogen_top']]
oxygen_A = test.u.atoms[test.ndx['oxygen_A']]
test.start_h_bond_analysis(oxygen_A,hydrogen_top)