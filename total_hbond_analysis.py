# Calculates the average hbond lifetime for Hbonds in plane, from plane to water , water to plane and water
from datetime import datetime
from analysis.class_analysis import *
from analysis.class_oxy_struck import *

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
    return avg_lifetime


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
data_directory = "/home/rasmus/Dropbox/Education/UCL/fourth Year/Project/data/"
hbond_file = open(data_directory + "filled_hb_db.out", 'r')

# loads the analysis object just to get the oxygen lists easily
analysis_obj = Analysis(data_directory + 'input.psf', data_directory + 'IOHMD-A-prod.dcd', data_directory + 'input.ndx',
                        'A')

interface_oxygen = analysis_obj.return_id_list_from_name("oxygen_A")
interface_oxygen += analysis_obj.return_id_list_from_name("oxygen_G")

water_id_list = analysis_obj.return_id_list_from_name("water")
water_oxygen = analysis_obj.get_water_oxygen(water_id_list)

all_oxygen = interface_oxygen + water_oxygen

hbond_db = load_hbond_database(hbond_file, all_oxygen)

# calc_avg_lifetime, takes list of donor and acceptor oxygens
lifetime_avg_i_i = calc_avg_lifetime(hbond_db, interface_oxygen, interface_oxygen)
lifetime_avg_i_w = calc_avg_lifetime(hbond_db, interface_oxygen, water_oxygen)
lifetime_avg_w_i = calc_avg_lifetime(hbond_db, water_oxygen, interface_oxygen)
lifetime_avg_w_w = calc_avg_lifetime(hbond_db, water_oxygen, water_oxygen)

print lifetime_avg_i_i
print lifetime_avg_i_w
print lifetime_avg_w_i
print lifetime_avg_w_w

print datetime.now() - start_time
