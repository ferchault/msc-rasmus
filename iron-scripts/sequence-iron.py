__author__ = 'rasmus'
#Calculates the sequence patterns... will have 3+ states ...
#Oxygen site no hydrogen 0
#oxygen site 1 hydrogen in plane 1 (in plane + below plane)
#--..-- out plane 1
#_---- ..-- more than 1 hydrogen

base_path = "/home/rasmus/ownCloud/data/iron-term-trajectory/"

h_position_file = open(base_path + "h_position.out", 'r')
first = True
hinterface = dict()
#create data structure
for line in h_position_file:
    if first:
        first = False
        continue
    parts = line.strip().split(',')
    traj = parts[0]
    frame = int(parts[1])
    surface = parts[2]
    in_out_below = int(parts[3])
    hyd_id = int(parts[4])
    nearest_o_id = int(parts[5])

    if traj not in hinterface:
        hinterface[traj] = dict()
    if surface not in hinterface[traj]:
        hinterface[traj][surface] = dict()
    if frame not in hinterface[traj][surface]:
        hinterface[traj][surface][frame] = dict()
    if nearest_o_id not in hinterface[traj][surface][frame]:
        hinterface[traj][surface][frame][nearest_o_id] = dict()
    if hyd_id not in hinterface[traj][surface][frame][nearest_o_id]:
        hinterface[traj][surface][frame][nearest_o_id][hyd_id] = in_out_below

# What we want is some sequence for the 12 oxygen sites "001 010 102 210"
# 0 == No hydrogen at site
# 1 == one hydrogen in or below plane
# 2 == one hydrogen out of plane
# 3 == 2 or more

found_patterns = dict()

for traj, maxframe in zip('A B'.split(), (4570, 1904)):
    for surface in 't b'.split():
        for frame in range(300,maxframe):

