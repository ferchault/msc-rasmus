#finds water oxygen that is further than 5 to 7 angstroms away from the oxygen hermatite interface
__author__ = 'rasmus'
from analysis.class_analysis import *

data_directory = "/home/rasmus/ownCloud/data/trajectory/"
analysis_obj_A = Analysis(data_directory +'input.psf', data_directory +'IOHMD-A-prod.dcd', data_directory +'input.ndx', 'A')
analysis_obj_B = Analysis(data_directory +'input.psf', data_directory +'IOHMD-B-prod.dcd', data_directory +'input.ndx', 'B')

oxygen_list_A = analysis_obj_A.return_id_list_from_name("oxygen_A")
oxygen_list_G = analysis_obj_A.return_id_list_from_name("oxygen_G")

water_id_list = analysis_obj_A.return_id_list_from_name("water")

oxygen_water_list_A = []
oxygen_water_list_A += analysis_obj_A.get_water_oxygen(water_id_list)

analysis_obj_A.set_min_plane_dist(5.0)
analysis_obj_B.set_min_plane_dist(5.0)


analysis_obj_A.start_interface_h_analysis(oxygen_list_A, oxygen_water_list_A, "t", 0, analysis_obj_A.u.trajectory.n_frames, "bulk_oxygen_5.out")
analysis_obj_B.start_interface_h_analysis(oxygen_list_A, oxygen_water_list_A, "t", 0, analysis_obj_B.u.trajectory.n_frames, "bulk_oxygen_5.out", True)
analysis_obj_A.start_interface_h_analysis(oxygen_list_G, oxygen_water_list_A, "b", 0, analysis_obj_A.u.trajectory.n_frames, "bulk_oxygen_5.out", True)
analysis_obj_B.start_interface_h_analysis(oxygen_list_G, oxygen_water_list_A, "b", 0, analysis_obj_B.u.trajectory.n_frames, "bulk_oxygen_5.out", True)

analysis_obj_A.set_min_plane_dist(7.0)
analysis_obj_B.set_min_plane_dist(7.0)

analysis_obj_A.start_interface_h_analysis(oxygen_list_A, oxygen_water_list_A, "t", 0, analysis_obj_A.u.trajectory.n_frames, "bulk_oxygen_7.out")
analysis_obj_B.start_interface_h_analysis(oxygen_list_A, oxygen_water_list_A, "t", 0, analysis_obj_B.u.trajectory.n_frames, "bulk_oxygen_7.out", True)
analysis_obj_A.start_interface_h_analysis(oxygen_list_G, oxygen_water_list_A, "b", 0, analysis_obj_A.u.trajectory.n_frames, "bulk_oxygen_7.out", True)
analysis_obj_B.start_interface_h_analysis(oxygen_list_G, oxygen_water_list_A, "b", 0, analysis_obj_B.u.trajectory.n_frames, "bulk_oxygen_7.out", True)
