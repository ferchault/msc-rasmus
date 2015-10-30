#Will calcualte number of hbonds in bulk water
__author__ = 'rasmus'
from analysis.class_analysis import *
from analysis.class_oxy_struck import *
from datetime import datetime
import matplotlib.pyplot as plt

start_time = datetime.now()
data_directory = "/home/rasmus/Dropbox/Education/UCL/fourth Year/Project/data/"
bulk_oxygen_file = open(data_directory + "bulk_oxygen_5.out", 'r')
hbond_file = open(data_directory + "rasmus-hbond-db.txt", 'r')

analysis_obj_A = Analysis(data_directory +'input.psf', data_directory +'IOHMD-A-prod.dcd', data_directory +'input.ndx', 'A')
analysis_obj_B = Analysis(data_directory +'input.psf', data_directory +'IOHMD-B-prod.dcd', data_directory +'input.ndx', 'B')

analysis_obj_A.set_parameters(3.5, 160)
trajectory, bulk_water = analysis_obj_A.hbond_bulkwater_analysis_build_structure(data_directory + "bulk_oxygen_5.out", data_directory + "rasmus-hbond-db.txt")

# analysis_obj_A.hbond_bulkwater_analysis_simple(trajectory, bulk_water)

mean = []
delta = []

plt.figure("fig1")
for i in xrange(0, 1):
    lifetimes, m = analysis_obj_A.hbond_bulkwater_analysis_population(trajectory, bulk_water, 50)
    delta.append(math.pow(2,i))
    mean.append(m)
    bin_number = lifetimes[-1] - lifetimes[0]
    plt.hist(lifetimes, bins=bin_number, histtype='step', label=str(math.pow(2, i)))

print mean
print delta

plt.legend()
plt.show()

plt.figure("fig2")
plt.scatter(delta, mean)
plt.show()

print datetime.now() - start_time
