__author__ = 'rasmus'
import numpy as np
import operator
#Calculates the sequence patterns... will have 3+ states ...
#Oxygen site no hydrogen 0
#oxygen site 1 hydrogen in plane 1 (in plane + below plane)
#--..-- out plane 1
#_---- ..-- more than 1 hydrogen

base_path = "C:\Users\Rasmus\ownCloud\data\iron-term-trajectory/"
#base_path = "/home/rasmus/ownCloud/data/iron-term-trajectory/"

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


def permute_turple(turple, permutation):
    permutation = [x-1 for x in permutation ]
    new_turple = [turple[x] for x in permutation]
    return new_turple


def build_key(a, b, c, d):
    res = ''
    for _ in (a, b, c, d):
        res += reduce(lambda x, y: x+y, map(str, _))
    return res
# What we want is some sequence for the 12 oxygen sites "001 010 102 210"
# 0 == No hydrogen at site
# 1 == one hydrogen in or below plane
# 2 == one hydrogen out of plane
# 3 == 2 or more

# Z positive interface
oxygen_idx_top = ((76,70,73),(138,132,135),(45,39,42),(107,101,104))
oxygen_idx_bot = ((57,55,56),(119,117,118),(26,24,25),(88,86,87))
found_patterns = dict()

for traj, maxframe in zip('A B'.split(), (4570, 1904)):
    for surface in 't b'.split():
        if surface=='t':
            flat_oxy_list = sum(oxygen_idx_top,())
        elif surface=='b':
            flat_oxy_list = sum(oxygen_idx_bot,())
        for frame in range(300,maxframe):
            sequence = ""
            for oxy_id in flat_oxy_list:
                if oxy_id in hinterface[traj][surface][frame]:
                    if len(hinterface[traj][surface][frame][oxy_id]) == 1:
                        for hyd_id in hinterface[traj][surface][frame][oxy_id]:
                            if hinterface[traj][surface][frame][oxy_id][hyd_id] == 0 or hinterface[traj][surface][frame][oxy_id][hyd_id] == 1:
                                sequence += "1"
                            else:
                                sequence += "2"
                    else:
                        sequence += "3"
                else:
                    sequence += "0"

            a = sequence[0:3]
            b = sequence[3:6]
            c = sequence[6:9]
            d = sequence[9:12]
            key1 = build_key(a,b,c,d)
            key2 = build_key(b,a,d,c)
            key3 = build_key(c,d,a,b)
            key4 = build_key(d,c,b,a)

            found = None
            for x in (key1,key2,key3,key4):
                if x in found_patterns:
                    found = x

            if found is None:
                found_patterns[key1] = 1
            else:
                found_patterns[found] += 1

print found_patterns

print sorted(found_patterns.items(), key=operator.itemgetter(1))

