# Calculates the average hbond lifetime for Hbonds in plane, from plane to water , water to plane and water
from datetime import datetime
from analysis.class_analysis import *
from analysis.class_oxy_struck import *
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
def split_array_by_delta(array, delta):
    steps = array[1:] - array[:-1]
    seperators = np.where(steps > delta)
    return np.split(array, seperators[0]+1)

def calc_avg_lifetime(hbond_db, donor_list, acceptor_list):
    avg_lifetime = 0
    hbond_list = dict()
    hbond_list['A'] = []
    hbond_list['B'] = []

    for traj in hbond_db:
        if traj == 'A':
            nframes = 4421
        elif traj == 'B':
            nframes = 4288
        for donor in hbond_db[traj]:
            if donor in donor_list:
                for frame in hbond_db[traj][donor][60: nframes - 60 ]:
                        if frame.has_hbond:
                            for hb in frame.Hbonds:
                                if hb.acceptor in acceptor_list:
                                    found = False
                                    for hb_list in hbond_list[traj]:

                                        if hb == hb_list.hbond:
                                            hb_list.framelist.append(frame.frame)
                                            found = True
                                    if not found:
                                        hbond_list[traj].append(HbondStrcture(hb, frame.frame))

    unique_hb = []
    for traj in hbond_list:
        for hb in hbond_list[traj]:
            hb.framelist = np.array(hb.framelist)
            parts = split_array_by_delta(hb.framelist, 1)
            unique_hb.append(parts)

    lifetimes = []
    for outer_array in unique_hb:
        for inner_array in outer_array:
            lifetimes.append(inner_array[-1] - inner_array[0])

    lifetimes.sort()
    lifetimes = np.array(lifetimes)

    avg_lifetime = np.mean(lifetimes)*0.5
    median_lifetime = np.median(lifetimes)*0.5

    print "median", median_lifetime
    print "average" , avg_lifetime
    print "max" ,np.max(lifetimes)*0.5
    print "min" ,np.min(lifetimes)*0.5
    print "25 percentile", np.percentile(lifetimes,25)*0.5
    print "75 percentile" ,np.percentile(lifetimes,75)*0.5
    return lifetimes*0.5


def load_hbond_database(hbond_file, full_oxygen_list):
    first = True

    hb_structure_db = dict()
    hb_structure_db['A'] = dict()
    hb_structure_db['B'] = dict()
    n_frames = 4421

    for donor in full_oxygen_list:
        array = []
        for i in xrange(n_frames+1):
            array.append(OxyStruck(donor, i))
        hb_structure_db['A'][donor] = array
        hb_structure_db['B'][donor] = array


    for line in hbond_file:
        if first:
            first = False
            continue
        parts = line.strip().split()
        traj = parts[0]
        frame = int(parts[1])
        donor = int(parts[2])
        acceptor = int(parts[3])
        hydrogen = int(parts[4])
        hb_structure_db[traj][donor][frame].has_hbond = True
        hb_structure_db[traj][donor][frame].Hbonds.append(Hbond(donor, acceptor, hydrogen))
    return hb_structure_db


start_time = datetime.now()

##hbond definitions taken into account in database file - min oho angle 160, max oo dist 3.5 AA and delta of 60
data_directory = "/home/rasmus/ownCloud/UCL/fourth Year/Project/data/"
hbond_file = open(data_directory + "hbond_surface_analysis/filled_hb_db.out", 'r')

# loads the analysis object just to get the oxygen lists easily
analysis_obj_A = Analysis(data_directory + 'trajectory/input.psf', data_directory + 'trajectory/IOHMD-A-prod.dcd', data_directory + 'trajectory/input.ndx',
                        'A')
analysis_obj_B = Analysis(data_directory + 'trajectory/input.psf', data_directory + 'trajectory/IOHMD-B-prod.dcd', data_directory + 'trajectory/input.ndx',
                        'B')

interface_oxygen = analysis_obj_A.return_id_list_from_name("oxygen_A")
interface_oxygen += analysis_obj_A.return_id_list_from_name("oxygen_G")

water_id_list = analysis_obj_A.return_id_list_from_name("water")
water_oxygen = analysis_obj_A.get_water_oxygen(water_id_list)

all_oxygen = interface_oxygen + water_oxygen

hbond_db = load_hbond_database(hbond_file, all_oxygen)

# calc_avg_lifetime, takes list of donor and acceptor oxygens
print 'i i'
lifetime_data_i_i = calc_avg_lifetime(hbond_db, interface_oxygen, interface_oxygen)
print "i w"
lifetime_data_i_w = calc_avg_lifetime(hbond_db, interface_oxygen, water_oxygen)
print "w i"
lifetime_data_w_i = calc_avg_lifetime(hbond_db, water_oxygen, interface_oxygen)
print "w w (all)"
lifetime_data_w_w = calc_avg_lifetime(hbond_db, water_oxygen, water_oxygen)


analysis_obj_A.set_parameters(3.5 ,160)
analysis_obj_B.set_parameters(3.5, 160)
trajectory_A_5, bulk_water_A_5 = analysis_obj_A.hbond_bulkwater_analysis_build_structure_bulkwater(data_directory + "hbond_surface_analysis/bulk_oxygen_5.out", data_directory + "hbond_surface_analysis/filled_hb_db.out")
trajectory_B_5, bulk_water_B_5 = analysis_obj_A.hbond_bulkwater_analysis_build_structure_bulkwater(data_directory + "hbond_surface_analysis/bulk_oxygen_5.out", data_directory + "hbond_surface_analysis/filled_hb_db.out")

lifetimes1, m, std, avg_hb = analysis_obj_A.hbond_bulkwater_analysis_population(trajectory_A_5, bulk_water_A_5, 60)

lifetimes2, m, std, avg_hb = analysis_obj_B.hbond_bulkwater_analysis_population(trajectory_B_5, bulk_water_B_5, 60)

lifetimes = np.concatenate((lifetimes1, lifetimes2))*0.5


boxes=[]

boxes.append(lifetime_data_i_i)
boxes.append(lifetime_data_i_w)
boxes.append(lifetime_data_w_i)
boxes.append(lifetime_data_w_w)
boxes.append(lifetimes)

ax =plt.axes()

plt.boxplot(boxes, whis=[0,100], showmeans=True)
# axe = plt.boxplot(lifetime_data_i_i , whis=[0, 100])
# plt.boxplot(lifetime_data_i_w , whis=[0, 100], ax=axe)
# plt.boxplot(lifetime_data_w_i , whis=[0, 100], ax=axe)
# plt.boxplot(lifetime_data_w_w , whis=[0, 100], ax=axe)
# plt.boxplot(lifetimes , whis=[0, 100], ax=axe)

loc = plticker.MultipleLocator(base=50.0)

ax.yaxis.set_major_locator(loc)
ax.set_xticklabels(["i-i", "i-w", "w-i", "w-w all", "w-w bulk 5AA"])
plt.show()
print datetime.now() - start_time
