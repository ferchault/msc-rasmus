__author__ = 'rasmus'

from analysis.class_analysis import *
import scipy, numpy as np

data_directory = "/home/rasmus/ownCloud/data/iron-term-trajectory/"

analysis_obj_A = Analysis(data_directory + 'IO3-index.xyz', data_directory + 'IO3-A.dcd', data_directory +  'input.ndx', 'A')
analysis_obj_B = Analysis(data_directory + 'IO3-index.xyz', data_directory + 'IO3-B.dcd', data_directory +  'input.ndx', 'B')

oxygen_list_t = analysis_obj_A.return_id_list_from_name("oxygen_A")
oxygen_list_b = analysis_obj_A.return_id_list_from_name("oxygen_G")

hydrogen_list_t = analysis_obj_A.return_id_list_from_name("hydrogen_top")
hydrogen_list_b = analysis_obj_A.return_id_list_from_name("hydrogen_bottom")

analysis_obj_A.create_nearest_atom_list(hydrogen_list_t, oxygen_list_t , 't' , 300, 0, data_directory + "hyd-nb.out")
analysis_obj_B.create_nearest_atom_list(hydrogen_list_t, oxygen_list_t , 't' , 300, 0, data_directory + "hyd-nb.out", True)

analysis_obj_A.create_nearest_atom_list(hydrogen_list_b, oxygen_list_b , 'b' , 300, 0, data_directory + "hyd-nb.out", True)
analysis_obj_B.create_nearest_atom_list(hydrogen_list_b, oxygen_list_b , 'b' , 300, 0, data_directory + "hyd-nb.out", True)