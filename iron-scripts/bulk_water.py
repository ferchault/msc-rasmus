__author__ = 'rasmus'
from analysis.class_analysis import *
##Calculates the atom density from a given plane
base_path = "/home/rasmus/ownCloud/data/iron-term-trajectory/"
analysis_obj_A = Analysis(base_path + 'IO3-index.xyz', base_path + 'IO3-A.dcd', base_path + 'input.ndx', 'A')
analysis_obj_B = Analysis(base_path + 'IO3-index.xyz', base_path + 'IO3-B.dcd', base_path + 'input.ndx', 'B')

oxygen_list_A = analysis_obj_A.return_id_list_from_name("oxygen_B")
oxygen_list_G = analysis_obj_A.return_id_list_from_name("oxygen_F")

water_id_list = analysis_obj_A.return_id_list_from_name("water")

oxygen_water_list_A = []
oxygen_water_list_A += analysis_obj_A.get_water_oxygen(water_id_list)

analysis_obj_A.set_min_plane_dist(9.5)
analysis_obj_B.set_min_plane_dist(9.5)

analysis_obj_A.start_interface_h_analysis(oxygen_list_A, oxygen_water_list_A, "t", 300, analysis_obj_A.u.trajectory.n_frames, base_path + "bulk_oxygen_7.out")
analysis_obj_B.start_interface_h_analysis(oxygen_list_A, oxygen_water_list_A, "t", 300, analysis_obj_A.u.trajectory.n_frames, base_path + "bulk_oxygen_7.out", True)
analysis_obj_A.start_interface_h_analysis(oxygen_list_G, oxygen_water_list_A, "b", 300, analysis_obj_A.u.trajectory.n_frames, base_path + "bulk_oxygen_7.out", True)
analysis_obj_B.start_interface_h_analysis(oxygen_list_G, oxygen_water_list_A, "b", 300, analysis_obj_A.u.trajectory.n_frames, base_path + "bulk_oxygen_7.out", True)
