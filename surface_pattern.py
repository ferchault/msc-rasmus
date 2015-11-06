#!/usr/bin/env python
import numpy as np

data_directory = "/home/rasmus/Dropbox/Education/UCL/fourth Year/Project/data/"

# read hinterface files
hinterface_file = data_directory + "hplanedist.out"
first = True

hinterface_data = dict()
# hinterface_data[traj][t/b][frame][hydrogen] = dist
for line in open(hinterface_file):
	if first:
		first = False
		continue
	parts = line.strip().split()
	frame = int(parts[1])
	surface = parts[2]
	traj = parts[0]
	H_idx = int(parts[3])
	dist = float(parts[4])

	if traj not in hinterface_data:
		hinterface_data[traj] = dict()
	if surface not in hinterface_data[traj]:
		hinterface_data[traj][surface] = dict()
	if frame not in hinterface_data[traj][surface]:
		hinterface_data[traj][surface][frame] = dict()
	if H_idx not in hinterface_data[traj][surface][frame]:
		hinterface_data[traj][surface][frame][H_idx] = dict()
	hinterface_data[traj][surface][frame][H_idx] = dist

# read indices
# per unit cell:
#    1
# 0
#    2
# per surface:
# a b		^
# c d   x-> y
oxygen_idx_top = ((80, 74, 77), (146, 140, 143), (47, 41, 44), (113, 107, 110))
oxygen_idx_bot = ((59, 57, 58), (125, 123, 124), (26, 24, 25), (92, 90, 91))

lookup = dict(((92,8), (26,2), (93,9), (91,7), (27,3), (125,11), (25,1), (59,5), (126,12), (124,10), (58,4), (60,6), (81,18), (111,20), (45,14), (48,15), (42,13), (114,21), (108,19), (78,17), (75,16), (147,24), (144,23), (141,22)))
hydrogen_idx_top = [[lookup[_+1]-1 for _ in x] for x in oxygen_idx_top]
hydrogen_idx_bot = [[lookup[_+1]-1 for _ in x] for x in oxygen_idx_bot]

hs = dict((('t', hydrogen_idx_top), ('b', hydrogen_idx_bot)))
os = dict((('t', oxygen_idx_top), ('b', oxygen_idx_bot)))

# build pattern
threshold = .5

found_patterns = dict()

def build_key(a, b, c, d):
	res = ''
	for _ in (a, b, c, d):
		res += reduce(lambda x, y: x+y, map(str, _))
	return res

for traj, maxframe in zip('A B'.split(), (4421, 4288)):
	for surface in 't b'.split():
		for frame in range(0, maxframe):
			try:
				abcd = map(lambda x: [hinterface_data[traj][surface][frame][_] for _ in x], hs[surface])
			except:
				print traj, surface, frame
				raise
			a, b, c, d = map(lambda x: list((np.array(x) > threshold)*1), abcd)
			key1 = build_key(a, b, c, d)
			key2 = build_key(d, c, b, a)
			key3 = build_key(b, a, d, c)
			key4 = build_key(c, d, a, b)
			found = None
			for x in (key1, key2, key3, key4):
				if x in found_patterns:
					found = x
					break
			if found is None:
				found_patterns[key1] = 1
			else:
				found_patterns[found] += 1

#Read the hb database
hb_data_file = open(data_directory + "filled_hb_db.out")
first = True
for line in hb_data_file:
    if first:
        first = False
        continue
    parts = line.strip().split()
    traj = parts[0]
    frame = int(parts[1])
    donor = int(parts[2])
    acceptor = int(parts[3])
    hydrogen = int(parts[4])
