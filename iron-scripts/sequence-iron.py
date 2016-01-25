__author__ = 'rasmus'
import numpy as np
import operator
#Calculates the sequence patterns... will have 3+ states ...
#Oxygen site no hydrogen 0
#oxygen site 1 hydrogen in plane 1 (in plane + below plane)
#--..-- out plane 1
#_---- ..-- more than 1 hydrogen

#base_path = "C:\Users\Rasmus\ownCloud\data\iron-term-trajectory/"
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
sequence_array = []

for traj, maxframe in zip('A B'.split(), (4570, 1904)):
    for surface in 't b'.split():
        for frame in range(300,maxframe):
            sequence = ""
            if surface == 't':
                oxy_list = oxygen_idx_top
                flat_oxy_list = sum(oxy_list,())
            elif surface == 'b':
                oxy_list = oxygen_idx_bot
                flat_oxy_list = sum(oxy_list,())
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
                    if found == key2:
                        oxy_list = permute_turple(oxy_list, [2, 1 , 4 ,3])
                    elif found == key3:
                        oxy_list = permute_turple(oxy_list, [3, 4 , 1 ,2])
                    elif found == key4:
                        oxy_list = permute_turple(oxy_list, [4, 3 , 2 ,1])

                    oxy_list = sum(oxy_list,())
                    sequence_array.append((traj,frame,surface,str(found),oxy_list))
                    break
            if found is None:
                found_patterns[key1] = 1
                oxy_list = sum(oxy_list,())
                sequence_array.append((traj,frame,surface, str(key1), oxy_list))
            else:
                found_patterns[found] += 1


nb_list = np.zeros((12,12))

nb_list[:,0] = [0 , 1, 1 , 0 , 1 , 1 , 0 ,0 ,1 ,0 ,1 , 0]
nb_list[:,1] = [1 , 0, 1 , 1 , 0 , 0 , 0 ,0 ,1 ,1 ,0 , 1]
nb_list[:,2] = [1 , 1, 0 , 1 , 0 , 0 , 1 ,1 ,0 ,0 ,1 , 0]

for j in range(0,3):
    a, b ,c ,d = [list(nb_list[i:i+3,j]) for i in xrange(0, len(nb_list), 3)]
    nb_list[:,j+3] = b + a + d + c
    nb_list[:,j+6] = c + d + a + b
    nb_list[:,j+9] = d + c + b + a

hb_file = open(base_path + "hb_filled_db.out", "r")
first = True
hb_db = dict()

for line in hb_file:
    if first:
        first = False
        continue
    lineparts = line.strip().split()
    traj = lineparts[0]
    frame = int(lineparts[1])
    donor = int(lineparts[2])
    acceptor = int(lineparts[3])
    hyd = int(lineparts[4])
    if traj not in hb_db:
        hb_db[traj] = dict()
    if frame not in hb_db[traj]:
        hb_db[traj][frame] = dict()
    if donor not in hb_db[traj][frame]:
        hb_db[traj][frame][donor] = dict()
    if acceptor not in hb_db[traj][frame][donor]:
        hb_db[traj][frame][donor][acceptor] = list()
    hb_db[traj][frame][donor][acceptor].append(hyd)

acceptor_matrix = dict()
donor_matrix = dict()
donor_out_matrix = dict()

matrix_a = np.zeros((12,2))
matrix_d_out = np.zeros((12,1))
matrix_d_in = np.zeros((12,6))

for x in sequence_array:
    traj = x[0]
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
        donors = hb_db[traj][frame].keys()
        for donor in donors:
            acceptors = hb_db[traj][frame][donor].keys()
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

a_file = open(base_path + "a_mat.out", 'w')
d_file = open(base_path + "d_mat.out", 'w')
d_out_file = open(base_path + "d_out_mat.out", 'w')

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

# prints sequence array to file
sequence_file = open(base_path + "sequence.out" , 'w')
for lineparts in sequence_array:

    write_line_to_file(sequence_file, [lineparts[0],lineparts[2],lineparts[1],lineparts[3]])
