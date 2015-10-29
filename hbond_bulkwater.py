#Will calcualte number of hbonds in bulk water
__author__ = 'rasmus'
from analysis.class_analysis import *
from analysis.class_oxy_struck import *
from datetime import datetime

start_time = datetime.now()
data_directory = "/home/rasmus/Dropbox/Education/UCL/fourth Year/Project/data/"
bulk_oxygen_file = open(data_directory + "bulk_oxygen_5.out", 'r')
hbond_file = open(data_directory + "rasmus-hbond-db.txt", 'r')

analysis_obj_A = Analysis(data_directory +'input.psf', data_directory +'IOHMD-A-prod.dcd', data_directory +'input.ndx', 'A')
analysis_obj_B = Analysis(data_directory +'input.psf', data_directory +'IOHMD-B-prod.dcd', data_directory +'input.ndx', 'B')

analysis_obj_A.set_parameters(3.5, 160)
trajectory, bulk_water = analysis_obj_A.hbond_bulkwater_analysis_build_structure(data_directory + "bulk_oxygen_5.out", data_directory + "rasmus-hbond-db.txt")

#analysis_obj_A.hbond_bulkwater_analysis_simple(trajectory, bulk_water)

analysis_obj_A.hbond_bulkwater_analysis_population(trajectory, bulk_water, 5)

print datetime.now() - start_time
