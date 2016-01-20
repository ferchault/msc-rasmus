__author__ = 'rasmus'
from analysis.class_analysis import *
import scipy, numpy as np

data_directory = "/home/rasmus/Dropbox/Education/UCL/fourth Year/Project/data/"

analysis_obj_A = Analysis(data_directory + 'input.psf', data_directory + 'IOHMD-A-prod.dcd', data_directory +  'input.ndx', 'A')
analysis_obj_B = Analysis(data_directory + 'input.psf', data_directory + 'IOHMD-B-prod.dcd', data_directory +  'input.ndx', 'B')

#form 3 by N of coordinates
oxygen_list = analysis_obj_A.return_id_list_from_name("oxygen_A")
hydrogen_list = analysis_obj_A.return_id_list_from_name("hydrogen_top")



analysis_obj_A.start_interface_h_analysis(oxygen_list, hydrogen_list, "t", 0, analysis_obj_A.u.trajectory.n_frames, data_directory + "hplanedist.out")
analysis_obj_B.start_interface_h_analysis(oxygen_list, hydrogen_list, "t", 0, analysis_obj_B.u.trajectory.n_frames, data_directory + "hplanedist.out", append=True)


oxygen_list = analysis_obj_A.return_id_list_from_name("oxygen_G")
hydrogen_list = analysis_obj_A.return_id_list_from_name("hydrogen_bottom")

analysis_obj_A.start_interface_h_analysis(oxygen_list, hydrogen_list, "b", 0, analysis_obj_A.u.trajectory.n_frames, data_directory + "hplanedist.out", append=True)
analysis_obj_B.start_interface_h_analysis(oxygen_list, hydrogen_list, "b", 0, analysis_obj_B.u.trajectory.n_frames, data_directory + "hplanedist.out", append=True)
