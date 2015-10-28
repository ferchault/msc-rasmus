#Will calcualte number of hbonds in bulk water
__author__ = 'rasmus'
from analysis.class_analysis import *
from analysis.class_oxy_struck import *

data_directory = "../data/"
bulk_oxygen_file = open(data_directory + "bulk_oxygen_5.out", 'r')
hbond_file = open(data_directory + "rasmus-hbond-db.txt", 'r')

analysis_obj_A = Analysis('../data/input.psf', '../data/IOHMD-A-prod.dcd', '../data/input.ndx', 'A')
analysis_obj_B = Analysis('../data/input.psf', '../data/IOHMD-B-prod.dcd', '../data/input.ndx', 'B')

n_frames_A = analysis_obj_A.u.trajectory.n_frames
n_frames_B = analysis_obj_B.u.trajectory.n_frames

min_OHO_angle = 160
max_OO_dist = 3.5

bulk_oxygen_file.next()
hbond_file.next()

#creates list of unique oxygen ids that appear in the trajectory
oxygen_bulk_unique_A = set()
oxygen_bulk_unique_B = set()

oxygen_bulk_A = []
oxygen_bulk_B = []
for line in bulk_oxygen_file:
    frame = line.split()[1]
    oxy_id = line.split()[3]
    if line.split()[0] == 'A':
        oxygen_bulk_unique_A.add(int(line.split()[3]))
    elif line.split()[0] == 'B':
        oxygen_bulk_unique_B.add(int(line.split()[3]))

#form data structure
trajectory_A = {}
trajectory_B = {}

for oxy_id in oxygen_bulk_unique_A:
    array = []
    for i in xrange(n_frames_A):
        array.append(OxyStruck(oxy_id, i))
    trajectory_A[oxy_id] = array

for oxy_id in oxygen_bulk_unique_B:
    array = []
    for i in xrange(n_frames_B):
        array.append(OxyStruck(oxy_id, i))
    trajectory_B[oxy_id] = array

for line in hbond_file:
    line_parts = line.split()
    traj_name = line_parts[0]
    frame = int(line_parts[1])
    donor_id = int(line_parts[2])
    hyd_id = int(line_parts[4])
    acceptor_id = int(line_parts[3])
    oo_dist = float(line_parts[7])
    oho_angle = float(line_parts[8])

    if oho_angle > min_OHO_angle and oo_dist < max_OO_dist:
        2


bulk_oxygen_file.close()
hbond_file.close()