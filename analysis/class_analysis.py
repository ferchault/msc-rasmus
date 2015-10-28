import MDAnalysis
import numpy as np
import math


class Analysis:
    def __init__(self, psf_path, dcd_path, ndx_path, trajectory_name):
        self.u = MDAnalysis.Universe(psf_path, dcd_path)
        self.ndx = self.ndx_file_read_to_list(ndx_path)
        self.max_oo_dist = 0.0
        self.max_oho_angle = 0.0
        self.max_oh_dist = 0.0
        self.hmat = np.zeros((3, 3))
        self.hinv = np.zeros((3, 3))
        self.s_coords = np.zeros((len(self.u.atoms), 3))
        self.output_file = ''
        self.trajectory_name = trajectory_name

    def set_parameters(self, max_oo_dist, max_oh_dist, max_oho_angle):
        self.max_oh_dist = max_oh_dist
        self.max_oo_dist = max_oo_dist
        self.max_oho_angle = max_oho_angle

    def start_h_bond_analysis(self, oxygen_list_name, hydrogen_list_name, start_timestep, end_timestep, output_file_path, append=False):
        hydrogen_list = self.u.atoms[hydrogen_list_name]
        oxygen_list = self.u.atoms[oxygen_list_name]

        if append:
            output_file = open(output_file_path, 'a')
        else:
            output_file = open(output_file_path, 'w')
            self.write_line_to_file(output_file, ['trajectory', 'frame', 'donor_oxygen_index',
                                                  'acceptor_oxygen_index', 'hydrogen_index',
                                                  'OH_distance_donor', 'OH_distance_acceptor',
                                                  'OO_distance', 'OHO_angle'])

        for frame in self.u.trajectory[start_timestep:end_timestep]:
            self.h_bond_analysis_frame(oxygen_list, hydrogen_list, output_file)

        output_file.close()

    def h_bond_analysis_frame(self, oxygen_list, hydrogen_list, output_file):

        self.hmat = self.abc_to_hmatrix(*self.u.dimensions)
        self.hinv = np.linalg.inv(self.hmat)
        self.calc_si_all(self.hinv)

        for i in xrange(oxygen_list.n_atoms):
            i_index = oxygen_list.atoms.indices[i]
            for j in xrange(i+1, oxygen_list.n_atoms):
                j_index = oxygen_list.atoms.indices[j]
                dist_oo = self.calc_dist(i_index, j_index)
                if dist_oo <= self.max_oo_dist:
                    for h in xrange(hydrogen_list.n_atoms):
                        h_index = hydrogen_list.atoms.indices[h]
                        angle_oho = self.calc_angle(i_index, h_index, j_index)
                        if angle_oho >= self.max_oho_angle:
                            dist_ih = self.calc_dist(i_index, h_index)
                            dist_jh = self.calc_dist(j_index, h_index)
                            if dist_ih <= self.max_oh_dist or dist_jh <= self.max_oh_dist:
                                if dist_ih < dist_jh:
                                    donor = i_index
                                    acceptor = j_index
                                    donor_dist = dist_ih
                                    acceptor_dist = dist_jh
                                else:
                                    donor = j_index
                                    acceptor = i_index
                                    donor_dist = dist_jh
                                    acceptor_dist = dist_ih
                                self.write_line_to_file(output_file, [self.trajectory_name, self.u.trajectory.frame,
                                                         donor, acceptor, h_index, donor_dist, acceptor_dist,
                                                         dist_oo, angle_oho])

    # form H matrix  (given that it's a triclinic cell)
    # based on code by Guido https://github.com/mdtraj/mdtraj/issues/908
    @staticmethod
    def abc_to_hmatrix(a, b, c, alpha, beta, gamma, degrees=True):
        if degrees:
                alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
        result = np.zeros((3, 3))

        a = np.array((a, 0, 0))
        b = b * np.array((math.cos(gamma), math.sin(gamma), 0))
        bracket = (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
        c = c * np.array((math.cos(beta), bracket, math.sqrt(math.sin(beta) ** 2 - bracket ** 2)))

        result[:, 0] = a
        result[:, 1] = b
        result[:, 2] = c

        return result

    # Calculates the reduced s_i for all atoms
    def calc_si_all(self, hinv):
        for i in xrange(self.u.atoms.n_atoms):
            self.s_coords[i] = np.dot(hinv, self.u.atoms.coordinates()[i])

    # Calculate inter atom distance
    def calc_dist(self, i, j):
        s_ij = self.s_coords[i] - self.s_coords[j]
        s_ij = s_ij - np.round(s_ij)
        return np.linalg.norm(np.dot(self.hmat, s_ij))

    # Calculate angle in deg
    def calc_angle(self, i, j, k):
        s_hj = self.s_coords[i] - self.s_coords[j]
        s_hj = s_hj - np.round(s_hj)
        r_hj = np.dot(self.hmat, s_hj)
        
        s_hk = self.s_coords[k] - self.s_coords[j]
        s_hk = s_hk - np.round(s_hk)
        r_hk = np.dot(self.hmat, s_hk)
        
        angle = np.arccos(np.dot(r_hk, r_hj)/(np.linalg.norm(r_hk)*np.linalg.norm(r_hj)))
        angle = np.rad2deg(angle)
        return angle

    # Reads a NDX file and stores each atom group in list
    @staticmethod
    def ndx_file_read_to_list(file_path):
        ndx = {}
        name_index_str = []
        with open(file_path) as f:
            for line in f:
                if line.startswith('['):
                    name_index_str = line.split()[1]
                else:
                    id_list = line.split()
                    id_list = [int(i)-1 for i in id_list]
                    if id_list:
                        ndx[name_index_str] = id_list 
        return ndx

    def return_id_list_from_name(self, ndx_name):
        return self.ndx[ndx_name]

    @staticmethod
    def write_line_to_file(data_file, line_data_list):
        for i in line_data_list:
            data_file.write(str(i))
            data_file.write(' ')
        data_file.write('\n')

    def get_water_oxygen(self, water_id_list):
        water = self.u.atoms[water_id_list]
        water = water.select_atoms("type OW")
        return [int(i) for i in water.indices]

    def get_water_hydrogen(self, water_id_list):
        water = self.u.atoms[water_id_list]
        water = water.select_atoms("type HW")
        return [int(i) for i in water.indices]


    def start_interface_h_analysis(self, plane_atoms_ids, point_atoms_ids, plane, start_timestep, end_timestep, output_file_path, append=False):
        plane_atoms = self.u.atoms[plane_atoms_ids]
        point_atoms = self.u.atoms[point_atoms_ids]

        if append:
            output_file = open(output_file_path, 'a')
        else:
            output_file = open(output_file_path, 'w')
            self.write_line_to_file(output_file, ['trajectory', 'frame', 'topbottom', 'hydrogen_index',
                                                 'OH_plane_distance', ])

        for frame in self.u.trajectory[start_timestep:end_timestep]:
            self.interface_h_analysis_frame(plane_atoms, point_atoms, plane, output_file)

        output_file.close()

    def interface_h_analysis_frame(self, plane_atoms, point_atoms, plane, output_file):
        #Tranposes the coordinate matrix from N by 3 to 3 by N for use with SVD methods
        matrix = np.transpose(plane_atoms.atoms.coordinates())
        center = matrix.mean(axis=1)
        matrix -= center[:, np.newaxis]

        normal = np.linalg.svd(matrix)[0][:, -1]

        if plane == "b":
            normal = -normal

        for patom in point_atoms:
            plane_to_point = patom.position - center
            plane_to_point_dist = np.dot(normal, plane_to_point)

            i_index = patom.id
            self.write_line_to_file(output_file, [self.trajectory_name, self.u.trajectory.frame,
                                                  plane, i_index, plane_to_point_dist])

        # for i in xrange(len(point_atoms)):
        #     plane_to_point = point_atoms.atoms.coordinates()[i] - center
        #     plane_to_point_dist = np.dot(normal, plane_to_point)

        #     i_index = point_atoms.atoms.indices[i]
        #     self.write_line_to_file(output_file, [self.trajectory_name, self.u.trajectory.frame,
        #                                           plane, i_index, plane_to_point_dist])



