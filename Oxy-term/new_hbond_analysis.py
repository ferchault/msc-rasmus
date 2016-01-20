__author__ = 'rasmus'
from analysis.class_analysis import *

#Uses critera for the angle at 160 and distance 3.5 to create
#an updated hbond database with filled in hbonds (as in filled in up to delta)

data_directory = "/home/rasmus/Dropbox/Education/UCL/fourth Year/Project/data/"
hbond_file = open(data_directory + "rasmus-hbond-db.txt", 'r')

analysis_obj_A = Analysis(data_directory +'input.psf', data_directory +'IOHMD-A-prod.dcd', data_directory +'input.ndx', 'A')
analysis_obj_B = Analysis(data_directory +'input.psf', data_directory +'IOHMD-B-prod.dcd', data_directory +'input.ndx', 'B')

#create oxygen list for water oxygen and A and G oxygens
oxygen_list = analysis_obj_A.return_id_list_from_name("oxygen_A")
oxygen_list += analysis_obj_A.return_id_list_from_name("oxygen_G")

water_id_list = analysis_obj_A.return_id_list_from_name("water")

oxygen_list += analysis_obj_A.get_water_oxygen(water_id_list)

analysis_obj_A.set_parameters(3.5, 160)

Hbond_db_structure = analysis_obj_A.hbond_analysis_build_structure(data_directory + "rasmus-hbond-db.txt", oxygen_list)
Hbond_db_filled_structure = analysis_obj_A.hbond_get_filled_hbond(Hbond_db_structure, 60)
analysis_obj_A.hbond_analysis_write_filled_hbond(Hbond_db_filled_structure, "new_hb_filled_db.out")

analysis_obj_B.set_parameters(3.5, 160)

Hbond_db_structure = analysis_obj_B.hbond_analysis_build_structure(data_directory + "rasmus-hbond-db.txt", oxygen_list)
Hbond_db_filled_structure = analysis_obj_B.hbond_get_filled_hbond(Hbond_db_structure, 60)
analysis_obj_B.hbond_analysis_write_filled_hbond(Hbond_db_filled_structure, "new_hb_filled_db.out", append=True)