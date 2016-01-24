#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
import operator

import time

start = time.time()
data_directory = "C:\Users\Rasmus\ownCloud\data\hbond_surface_analysis/"
# data_directory = "/home/rasmus/Dropbox/Education/UCL/fourth Year/Project/data/"

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


lookup = dict(((92,8), (26,2), (93,9), (91,7), (27,3), (125,11), (25,1), (59,5), (126,12), (124,10), (58,4), (60,6), (81,18), (111,20), (45,14), (48,15), (42,13),
			   (114,21), (108,19), (78,17), (75,16), (147,24), (144,23), (141,22)))
hydrogen_idx_top = [[lookup[_+1]-1 for _ in x] for x in oxygen_idx_top]
hydrogen_idx_bot = [[lookup[_+1]-1 for _ in x] for x in oxygen_idx_bot]

hs = dict((('t', hydrogen_idx_top), ('b', hydrogen_idx_bot)))
os = dict((('t', oxygen_idx_top), ('b', oxygen_idx_bot)))

# build pattern
threshold = .5

found_patterns = dict()

def permute_turple(turple, permutation):
	permutation = [x-1 for x in permutation ]

	new_turple = [turple[x] for x in permutation]
	return new_turple

def build_key(a, b, c, d):
	res = ''
	for _ in (a, b, c, d):
		res += reduce(lambda x, y: x+y, map(str, _))
	return res
#4421, 4288
sequence_array = []
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
			oxylist = ()
			if surface == 't':
				oxylist= oxygen_idx_top
			elif surface == 'b':
				oxylist= oxygen_idx_bot
			for x in (key1, key2, key3, key4):
				if x in found_patterns:
					found = x

					if found == key2:
						oxylist = permute_turple(oxylist, [4, 3 , 2 ,1])
					elif found == key3:
						oxylist = permute_turple(oxylist, [2, 1 , 4 ,3])
					elif found == key4:
						oxylist = permute_turple(oxylist, [3, 4 , 1 ,2])
					#flatten turple
					oxylist = sum(oxylist, ())
					sequence_array.append((traj, frame, surface, str(found), oxylist))
					break
			if found is None:
				found_patterns[key1] = 1
				oxylist = sum(oxylist, ())
				sequence_array.append((traj, frame, surface, str(key1), oxylist))
			else:
				found_patterns[found] += 1

# 110 100 100 100 sequence pattern for 1426 times
#Build neighbour list
nb_list = np.zeros((12,12))

nb_list[:,0] = [0 , 1, 1 , 0 , 1 , 1 , 0 ,0 ,1 ,0 ,1 , 0]
nb_list[:,1] = [1 , 0, 1 , 1 , 0 , 0 , 0 ,0 ,1 ,1 ,0 , 1]
nb_list[:,2] = [1 , 1, 0 , 1 , 0 , 0 , 1 ,1 ,0 ,0 ,1 , 0]

for j in range(0,3):
	a, b ,c ,d = [list(nb_list[i:i+3,j]) for i in xrange(0, len(nb_list), 3)]
	nb_list[:,j+3] = b + a + d + c
	nb_list[:,j+6] = c + d + a + b
	nb_list[:,j+9] = d + c + b + a

#Read the hb database
hb_data_file = open(data_directory + "filled_hb_db.out")
first = True
hb_db = dict()
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

	if traj not in hb_db:
		hb_db[traj] = dict()
	if frame not in hb_db[traj]:
		hb_db[traj][frame] = dict()
	if donor not in hb_db[traj][frame]:
		hb_db[traj][frame][donor] = dict()
	if acceptor not in hb_db[traj][frame][donor]:
		hb_db[traj][frame][donor][acceptor] = list()
	hb_db[traj][frame][donor][acceptor].append(hydrogen)

acceptor_matrix = dict()
donor_matrix = dict()
donor_out_matrix = dict()

matrix_a = np.zeros((12,2))
matrix_d_out = np.zeros((12,1))
matrix_d_in = np.zeros((12,6))

for x in sequence_array:
	tray = x[0]
	frame = x[1]
	key = x[3]
	id_list = x[4]

	if key not in donor_matrix:
		donor_matrix[key] = dict()

	if frame not in donor_matrix[key]:
		donor_matrix[key][frame] = np.zeros((12,6))

	if key not in donor_out_matrix:
		donor_out_matrix[key] = dict()
	if frame not in donor_out_matrix[key]:
		donor_out_matrix[key][frame] = np.zeros((12,1))

	if key not in acceptor_matrix:
		acceptor_matrix[key] = dict()
	if frame not in acceptor_matrix[key]:
		acceptor_matrix[key][frame] = np.zeros((12,2))

	for i in xrange(len(id_list)):
		id = id_list[i]
		donors = hb_db[tray][frame].keys()
		for donor in donors:
			acceptors = hb_db[tray][frame][donor].keys()
			for acceptor in acceptors:
				#acceptor being donated to onto interface
				if acceptor == id and donor not in id_list:
					acceptor_matrix[key][frame][i][1] += 1
				#Acceptor being donated to from interface
				elif acceptor == id and donor in id_list:
					acceptor_matrix[key][frame][i][0] += 1

				if donor == id and acceptor not in id_list:
					donor_out_matrix[key][frame][i][0] += 1

				elif donor == id and acceptor in id_list:
					masked = np.ma.array(id_list, mask = nb_list[:,i]-1).compressed()
					for nb_index in xrange(len(masked)):
						if masked[nb_index] == acceptor:
							donor_matrix[key][frame][i][nb_index] += 1

def write_line_to_file(data_file, line_data_list):
	for i in line_data_list:
		data_file.write(str(i))
		data_file.write(' ')
	data_file.write('\n')

a_file = open(data_directory + "a_mat.out", 'w')
d_file = open(data_directory + "d_mat.out", 'w')
d_out_file = open(data_directory + "d_out_mat.out", 'w')

total_acceptor_matrix = dict()
write_line_to_file(a_file, ["seq" , "row", "in plane", "out plane"])
write_line_to_file(d_out_file, ["seq" , "row", "count"])
write_line_to_file(d_file, ["seq" , "row", "nb1", "nb2" , "nb3", "nb4", "nb5", "nb6"])

for key in acceptor_matrix:
	total_acceptor_matrix[key] = reduce(lambda x, y: x+y, acceptor_matrix[key].itervalues())/found_patterns[key]
	for row_num in xrange(0, 12):
		write_line_to_file(a_file, [key, row_num, total_acceptor_matrix[key][row_num][0], total_acceptor_matrix[key][row_num][1]])

total_donor_out_matrix = dict()
for key in donor_out_matrix:
	total_donor_out_matrix[key] = reduce(lambda x, y: x+y, donor_out_matrix[key].itervalues())/found_patterns[key]
	for row_num in xrange(0, 12):
		row = total_donor_out_matrix[key][row_num]
		write_line_to_file(d_out_file, [key, row_num, row[0]])


total_donor_matrix = dict()
for key in donor_matrix:
	total_donor_matrix[key] = reduce(lambda x, y: x+y, donor_matrix[key].itervalues())/found_patterns[key]
	for row_num in xrange(0, 12):
		row = total_donor_matrix[key][row_num]
		write_line_to_file(d_file, [key, row_num, row[0], row[1], row[2], row[3], row[4], row[5]])


print time.time() - start

