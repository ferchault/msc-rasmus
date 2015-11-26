#!/usr/bin/env python
import pandas as pd
__author__ = 'Guido von Rudorff'

# input files
basepath = '/Volumes/ALFA/IOHMD-A-prod-22ce9a183b38cedfca608118c1fc99f9/msc-copy/'
hbond_db = basepath + 'new_hb_filled_db.txt'
h_dist = basepath + 'total.txt'

# determine without hbonds
cols = 'trajectory frame donor_oxygen_index acceptor_oxygen_index hydrogen_index OH_distance_donor OH_distance_acceptor OO_distance OHO_angle nn'.split()
hbonds = pd.read_csv(hbond_db, sep=' ', names=cols, comment='#')
cols = 'trajectory frame topbottom hydrogen_index OH_plane_distance nn'.split()
h_dist = pd.read_csv(h_dist, sep=' ', names=cols, comment='#')

# merge
blended = pd.merge(h_dist, hbonds, 'left', 'trajectory frame hydrogen_index'.split())
print 'no H bond in plane', blended.donor_oxygen_index[blended.OH_plane_distance <= 0.5].isnull().sum() / float(len(blended))
print 'no H bond out plane', blended.donor_oxygen_index[blended.OH_plane_distance >= 0.5].isnull().sum() / float(len(blended))

oxygen_idx_top = (80, 74, 77, 146, 140, 143, 47, 41, 44, 113, 107, 110)
oxygen_idx_bot = (59, 57, 58, 125, 123, 124, 26, 24, 25, 92, 90, 91)

rawpatterns = {
	'110100100100': 1425, '010010101100': 117, '001011101100': 34, '001011101101': 10, '110100110100': 2338,
	'110110101100': 2, '111010101101': 22, '111010101100': 367, '110100110000': 4, '110100111110': 3,
	'111100100100': 21, '001010101100': 83, '110100101100': 77, '000101101101': 10, '110100100110': 2200,
	'010101101100': 6, '010101101101': 1828, '111110100100': 4, '010101100101': 2398, '011010101100': 3687,
	'011010101101': 16, '110110100100': 204, '110100110110': 2282, '110010101100': 106, '110100111100': 120,
	'011011101101': 19, '011011101100': 4, '110100101110': 27}

in_plane = 0
out_plane = 0
for pattern in rawpatterns.iterkeys():
	t_in = len([_ for _ in pattern if _ == '0'])
	t_out = len([_ for _ in pattern if _ == '1'])
	in_plane += t_in * rawpatterns[pattern]
	out_plane += t_out * rawpatterns[pattern]

print 'total in/out', in_plane / float(in_plane+out_plane), out_plane / float(in_plane+out_plane)

######
patterns = [None, ]
data = open(basepath + 'd_mat.out').readlines()[1:]
weighted = 0
total = 0
for line in data:
	parts = line.split()
	if parts[0] != patterns[-1] and parts[0] in patterns:
		continue
	if parts[0] not in patterns:
		patterns.append(parts[0])
	if parts[0][int(parts[1])] != '0':
		continue
	vals = sum(map(float, parts[2:]))
	weighted += vals * rawpatterns[parts[0]]
	total += rawpatterns[parts[0]]
print 'in plane w H bond to H', weighted/total
######
patterns = [None, ]
data = open(basepath + 'd_out_mat.out').readlines()[1:]
weighted = 0
total = 0
for line in data:
	parts = line.split()
	if parts[0] != patterns[-1] and parts[0] in patterns:
		continue
	if parts[0] not in patterns:
		patterns.append(parts[0])
	if parts[0][int(parts[1])] != '1':
		continue
	vals = sum(map(float, parts[2:]))
	weighted += vals * rawpatterns[parts[0]]
	total += rawpatterns[parts[0]]
print 'out plane w H bond to W', weighted/total
########
patterns = [None, ]
data = open(basepath + 'a_mat.out').readlines()[1:]
weighted = 0
total = 0
for line in data:
	parts = line.split()
	if parts[0] != patterns[-1] and parts[0] in patterns:
		continue
	if parts[0] not in patterns:
		patterns.append(parts[0])
	if parts[0][int(parts[1])] != '0':
		continue
	vals = float(parts[2])
	weighted += vals * rawpatterns[parts[0]]
	total += rawpatterns[parts[0]]
print 'in plane w H bond from H', weighted/total
########
patterns = [None, ]
data = open(basepath + 'a_mat.out').readlines()[1:]
weighted = 0
total = 0
for line in data:
	parts = line.split()
	if parts[0] != patterns[-1] and parts[0] in patterns:
		continue
	if parts[0] not in patterns:
		patterns.append(parts[0])
	if parts[0][int(parts[1])] != '0':
		continue
	vals = float(parts[3])
	weighted += vals * rawpatterns[parts[0]]
	total += rawpatterns[parts[0]]
print 'in plane w H bond from W', weighted/total

########
patterns = [None, ]
data = open(basepath + 'a_mat.out').readlines()[1:]
weighted = 0
total = 0
for line in data:
	parts = line.split()
	if parts[0] != patterns[-1] and parts[0] in patterns:
		continue
	if parts[0] not in patterns:
		patterns.append(parts[0])
	if parts[0][int(parts[1])] != '1':
		continue
	vals = float(parts[2])
	weighted += vals * rawpatterns[parts[0]]
	total += rawpatterns[parts[0]]
print 'out plane w H bond from H', weighted/total
########
patterns = [None, ]
data = open(basepath + 'a_mat.out').readlines()[1:]
weighted = 0
total = 0
for line in data:
	parts = line.split()
	if parts[0] != patterns[-1] and parts[0] in patterns:
		continue
	if parts[0] not in patterns:
		patterns.append(parts[0])
	if parts[0][int(parts[1])] != '1':
		continue
	vals = float(parts[3])
	weighted += vals * rawpatterns[parts[0]]
	total += rawpatterns[parts[0]]
print 'out plane w H bond from W', weighted/total