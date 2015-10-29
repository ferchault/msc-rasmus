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

n_frames_A = analysis_obj_A.u.trajectory.n_frames
n_frames_B = analysis_obj_B.u.trajectory.n_frames

min_OHO_angle = 160
max_OO_dist = 3.5
allowed_hb_break = 5

bulk_oxygen_file.next()
hbond_file.next()

#creates list of unique oxygen ids that appear in the trajectory
oxygen_bulk_unique_A = set()
oxygen_bulk_unique_B = set()

oxygen_bulk_A = [[] for i in range(n_frames_A)]
oxygen_bulk_B = [[] for i in range(n_frames_B)]

for line in bulk_oxygen_file:
    frame = int(line.split()[1])
    oxy_id = int(line.split()[3])
    if line.split()[0] == 'A':
        oxygen_bulk_A[frame].append(oxy_id)
        oxygen_bulk_unique_A.add(int(line.split()[3]))
    elif line.split()[0] == 'B':
        oxygen_bulk_B[frame].append(oxy_id)
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

    if traj_name == 'A':
        if donor_id in oxygen_bulk_A[frame]:
            trajectory_A[donor_id][frame].in_bulk = True
            if oo_dist <= max_OO_dist and oho_angle >= min_OHO_angle:
                trajectory_A[donor_id][frame].has_hbond = True
                trajectory_A[donor_id][frame].Hbonds.append(HbondStruck(donor_id, acceptor_id, hyd_id))
    elif traj_name == 'B':
        if donor_id in oxygen_bulk_B[frame]:
            trajectory_B[donor_id][frame].in_bulk = True
            if oo_dist <= max_OO_dist and oho_angle >= min_OHO_angle:
                trajectory_B[donor_id][frame].has_hbond = True
                trajectory_B[donor_id][frame].Hbonds.append(HbondStruck(donor_id, acceptor_id, hyd_id))


number_hb_A = 0
number_hb_break_A = 0
number_oxygen_bulk_A = 0

for frame in oxygen_bulk_A:
    number_oxygen_bulk_A += len(frame)

for donor in trajectory_A.itervalues():
    for c_frame in donor:
        if c_frame.in_bulk and c_frame.has_hbond:
            number_hb_A += len(c_frame.Hbonds)
            for hb in c_frame.Hbonds:
                if hb not in donor[c_frame.frame-1].Hbonds and c_frame.frame != 0:
                    number_hb_break_A += 1
print "A"
print number_hb_A
print number_oxygen_bulk_A
print number_hb_break_A

avg_hb_A = float(number_hb_A) / float(number_oxygen_bulk_A)

avg_lt_A = (avg_hb_A * (number_oxygen_bulk_A/n_frames_A) * (n_frames_B*0.5)) / (float(number_hb_break_A))

print avg_hb_A
print avg_lt_A

number_hb_B = 0
number_hb_break_B = 0
number_oxygen_bulk_B = 0

for frame in oxygen_bulk_A:
    number_oxygen_bulk_B += len(frame)

for donor in trajectory_B.itervalues():
    for c_frame in donor:
        if c_frame.in_bulk and c_frame.has_hbond:
            number_hb_B += len(c_frame.Hbonds)
            for hb in c_frame.Hbonds:
                if hb not in donor[c_frame.frame-1].Hbonds and c_frame.frame != 0:
                    number_hb_break_B += 1
print "B"
print number_hb_B
print number_oxygen_bulk_B
print number_hb_break_B

avg_hb_B = float(number_hb_B) / float(number_oxygen_bulk_B)

avg_lt_B = (avg_hb_B * (number_oxygen_bulk_B/n_frames_B) * (n_frames_B*0.5)) / (float(number_hb_break_B))

print avg_hb_B
print avg_lt_B

bulk_oxygen_file.close()
hbond_file.close()

print
print datetime.now() - start_time