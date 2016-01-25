__author__ = 'rasmus'
from analysis.class_analysis import *
base_path = "/home/rasmus/ownCloud/data/iron-term-trajectory/"
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
            nframes = 4570
        elif traj == 'B':
            nframes = 1904
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
    n_frames = 4570

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


filled_hb_file = open(base_path + "hb_filled_db.out")
bulk_water_db = open(base_path + "bulk_oxygen_7.out")

analysis_obj_A = Analysis(base_path + 'IO3-index.xyz', base_path + 'IO3-A.dcd', base_path + 'input.ndx', 'A')

#the second top layer of oxygen
layer_top = analysis_obj_A.return_id_list_from_name("oxygen_A")
layer_bot = analysis_obj_A.return_id_list_from_name("oxygen_G")

interface_oxygen = layer_top + layer_bot

water_id_list = analysis_obj_A.return_id_list_from_name("water")
water_oxygen = analysis_obj_A.get_water_oxygen(water_id_list)

all_oxygen = interface_oxygen + water_oxygen

hbond_db = load_hbond_database(filled_hb_file, all_oxygen)

# calc_avg_lifetime, takes list of donor and acceptor oxygens
print 'i i'
lifetime_data_i_i = calc_avg_lifetime(hbond_db, interface_oxygen, interface_oxygen)
print "i w"
lifetime_data_i_w = calc_avg_lifetime(hbond_db, interface_oxygen, water_oxygen)
print "w i"
lifetime_data_w_i = calc_avg_lifetime(hbond_db, water_oxygen, interface_oxygen)
print "w w (all)"
lifetime_data_w_w = calc_avg_lifetime(hbond_db, water_oxygen, water_oxygen)


def calc_avg_lifetime_bulk(hbond_db_file, bulk_water_db):

    for line in hbond_db_file:
        lineparts = line.strip().split()


    avg_lifetime = 0
    hbond_list = dict()
    hbond_list['A'] = []
    hbond_list['B'] = []

    for traj in hbond_db:
        if traj == 'A':
            nframes = 4570
        elif traj == 'B':
            nframes = 1904
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

print "w w (bulk)"
lifetime_data_w_w_bulk = calc_avg_lifetime_bulk(hbond_db, bulk_water_db, water_oxygen)