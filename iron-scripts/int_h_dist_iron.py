__author__ = 'rasmus'
from analysis.class_analysis import *
import scipy, numpy as np

data_directory = "/home/rasmus/ownCloud/data/iron-term-trajectory/"

analysis_obj_A = Analysis(data_directory + 'IO3-index.xyz', data_directory + 'IO3-A.dcd', data_directory +  'input.ndx', 'A')
analysis_obj_B = Analysis(data_directory + 'IO3-index.xyz', data_directory + 'IO3-B.dcd', data_directory +  'input.ndx', 'B')

#form 3 by N of coordinates
oxygen_list = analysis_obj_A.return_id_list_from_name("oxygen_B")
hydrogen_list = analysis_obj_A.return_id_list_from_name("hydrogen_top")
analysis_obj_A.set_max_plane_dist(12)
analysis_obj_B.set_max_plane_dist(12)
analysis_obj_A.set_min_plane_dist(0.0)
analysis_obj_B.set_min_plane_dist(0.0)

analysis_obj_A.start_interface_h_analysis(oxygen_list, hydrogen_list, "t", 300, analysis_obj_A.u.trajectory.n_frames, data_directory + "hplanedist.out")
analysis_obj_B.start_interface_h_analysis(oxygen_list, hydrogen_list, "t", 300, analysis_obj_B.u.trajectory.n_frames, data_directory + "hplanedist.out", append=True)


oxygen_list = analysis_obj_A.return_id_list_from_name("oxygen_F")
hydrogen_list = analysis_obj_A.return_id_list_from_name("hydrogen_bottom")

analysis_obj_A.start_interface_h_analysis(oxygen_list, hydrogen_list, "b", 300, analysis_obj_A.u.trajectory.n_frames, data_directory + "hplanedist.out", append=True)
analysis_obj_B.start_interface_h_analysis(oxygen_list, hydrogen_list, "b", 300, analysis_obj_B.u.trajectory.n_frames, data_directory + "hplanedist.out", append=True)
