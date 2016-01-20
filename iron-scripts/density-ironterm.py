__author__ = 'rasmus'
from analysis.class_analysis import *
##Calculates the atom density from a given plane
base_path = "/home/rasmus/ownCloud/data/iron-term-trajectory/"
analysis_obj_A = Analysis(base_path + 'IO3-index.xyz', base_path + 'IO3-A.dcd', base_path + 'input.ndx', 'A')
analysis_obj_B = Analysis(base_path + 'IO3-index.xyz', base_path + 'IO3-B.dcd', base_path + 'input.ndx', 'B')

analysis_obj_A.set_max_plane_dist(12)
analysis_obj_B.set_max_plane_dist(12)

#the second top layer of oxygen
layer_top = analysis_obj_A.return_id_list_from_name("oxygen_B")
layer_bot = analysis_obj_A.return_id_list_from_name("oxygen_F")

hydrogen_top = analysis_obj_A.return_id_list_from_name("hydrogen_top")
hydrogen_bot = analysis_obj_A.return_id_list_from_name("hydrogen_bottom")

analysis_obj_A.start_interface_h_analysis(layer_top, hydrogen_top, "t", 300, 0, base_path+"hydrogen_dist.out")
analysis_obj_A.start_interface_h_analysis(layer_bot, hydrogen_bot, "b", 300, 0, base_path+"hydrogen_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_top, hydrogen_top, "t", 300, 0, base_path+"hydrogen_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_bot, hydrogen_bot, "b", 300, 0, base_path+"hydrogen_dist.out", True)

oxygen_top_int = analysis_obj_A.return_id_list_from_name("oxygen_A")
oxygen_bot_int =analysis_obj_A.return_id_list_from_name("oxygen_G")
iron_top = analysis_obj_A.return_id_list_from_name("iron_A")
iron_bot = analysis_obj_A.return_id_list_from_name("iron_F")

water_id_list = analysis_obj_A.return_id_list_from_name("water")
oxygen_water_list = []
oxygen_water_list += analysis_obj_A.get_water_oxygen(water_id_list)
hydrogen_water_list = []
hydrogen_water_list += analysis_obj_A.get_water_hydrogen(water_id_list)

analysis_obj_A.start_interface_h_analysis(layer_top, oxygen_top_int, "t", 300, 0, base_path+"oxygen_int_dist.out")
analysis_obj_A.start_interface_h_analysis(layer_bot, oxygen_bot_int, "b", 300, 0, base_path+"oxygen_int_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_top, oxygen_top_int, "t", 300, 0, base_path+"oxygen_int_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_bot, oxygen_bot_int, "b", 300, 0, base_path+"oxygen_int_dist.out", True)

analysis_obj_A.start_interface_h_analysis(layer_top, iron_top, "t", 300, 0, base_path+"iron_dist.out")
analysis_obj_A.start_interface_h_analysis(layer_bot, iron_bot, "b", 300, 0, base_path+"iron_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_top, iron_top, "t", 300, 0, base_path+"iron_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_bot, iron_bot, "b", 300, 0, base_path+"iron_dist.out", True)

analysis_obj_A.set_min_plane_dist(0.0)
analysis_obj_B.set_min_plane_dist(0.0)

analysis_obj_A.start_interface_h_analysis(layer_top, hydrogen_water_list, "t", 300, 0, base_path+"hydrogen_water_dist.out")
analysis_obj_A.start_interface_h_analysis(layer_bot, hydrogen_water_list, "b", 300, 0, base_path+"hydrogen_water_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_top, hydrogen_water_list, "t", 300, 0, base_path+"hydrogen_water_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_bot, hydrogen_water_list, "b", 300, 0, base_path+"hydrogen_water_dist.out", True)


analysis_obj_A.start_interface_h_analysis(layer_top, oxygen_water_list, "t", 300, 0, base_path+"oxygen_water_dist.out")
analysis_obj_A.start_interface_h_analysis(layer_bot, oxygen_water_list, "b", 300, 0, base_path+"oxygen_water_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_top, oxygen_water_list, "t", 300, 0, base_path+"oxygen_water_dist.out", True)
analysis_obj_B.start_interface_h_analysis(layer_bot, oxygen_water_list, "b", 300, 0, base_path+"oxygen_water_dist.out", True)

