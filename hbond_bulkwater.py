#Will calcualte number of hbonds in bulk water
__author__ = 'rasmus'
from analysis.class_analysis import *
from analysis.class_oxy_struck import *
from datetime import datetime
import matplotlib.pyplot as plt

max_oo_dist = 3.5
min_oho_angle = 160
start_time = datetime.now()
data_directory = "/home/rasmus/Dropbox/Education/UCL/fourth Year/Project/data/"
bulk_oxygen_file = open(data_directory + "bulk_oxygen_5.out", 'r')
hbond_file = open(data_directory + "rasmus-hbond-db.txt", 'r')

analysis_obj_A = Analysis(data_directory +'input.psf', data_directory +'IOHMD-A-prod.dcd', data_directory +'input.ndx', 'A')
analysis_obj_B = Analysis(data_directory +'input.psf', data_directory +'IOHMD-B-prod.dcd', data_directory +'input.ndx', 'B')

fig = plt.figure("fig2-A-5")
# fighb = plt.figure("fighb-B-7")

ax = plt.subplot(111)
# axhb = plt.subplot(111)
for i in xrange(-1,2):
    for j in xrange(-1,2):

        current_max_oo = max_oo_dist + i*0.05*max_oo_dist
        current_min_oho = min_oho_angle + j*0.05*min_oho_angle
        analysis_obj_A.set_parameters(current_max_oo, current_min_oho)

        trajectory_A_5, bulk_water_A_5 = analysis_obj_A.hbond_bulkwater_analysis_build_structure_bulkwater(data_directory + "bulk_oxygen_5.out", data_directory + "rasmus-hbond-db.txt")

        mean = []
        delta = []
        avg_hb_list = []
        std_list = []

        # plt.figure("fig1-A-5-" + str(current_max_oo) + "-" + str(current_min_oho))
        for k in xrange(0, 6):
            lifetimes, m, std, avg_hb = analysis_obj_A.hbond_bulkwater_analysis_population(trajectory_A_5, bulk_water_A_5, 20 + k*10)
            delta.append(20 + k * 10)
            mean.append(m)
            std_list.append(std)
            avg_hb_list.append(avg_hb)
            # bin_number = (lifetimes[-1] - lifetimes[0])/2
            # plt.hist(lifetimes, bins=bin_number, histtype='step', label=str(30+i*10))

        print "A-5", current_max_oo , current_min_oho
        print mean
        print std_list
        print delta
        print avg_hb_list

        # plt.legend()
        # plt.savefig("fig1-A-5-" + str(current_max_oo) + "-" + str(current_min_oho) + '.png')
        # axhb.plot(delta, avg_hb_list, label=str(str(current_max_oo) + "-" + str(current_min_oho)), linestyle='none', marker='o')
        ax.errorbar(delta, mean, yerr=std_list, label=str(str(current_max_oo) + "-" + str(current_min_oho)), linestyle='none', marker='o')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

ax.set_xlim(25,75)
ax.set_xlabel('Delta/frames')
ax.set_ylabel('Mean life time / fs')
ax.set_title('A trajectory - 5A from interface')
fig.savefig("fig2-A-5.png")


print datetime.now() - start_time
