__author__ = 'rasmus'
#creates the final datafile for where the hydrogens are located

data_directory = "/home/rasmus/ownCloud/data/iron-term-trajectory/"

h_dist_db = open(data_directory + 'hplanedist.out', 'r')
h_nb_db = open(data_directory + 'hyd-nb.out', 'r')
h_position = open(data_directory + 'h_position.out', 'w')
#The minimum of the hydrogen density is at 2.7xx.. the peak for the oxygen dist is at 2.xx giving a in plane width of 1AA.

below_p = 1.7
above_p = 2.7

num_lines_1 = sum(1 for line in open(data_directory + 'hplanedist.out', 'r'))
num_lines_2 = sum(1 for line in open(data_directory + 'hyd-nb.out', 'r'))

if num_lines_1 != num_lines_2 :
    print "error"

h_dist_db_list = []
h_nb_db_list = []

for line in h_dist_db:
    h_dist_db_list.append(line)
for line in h_nb_db:
    h_nb_db_list.append(line)

h_position.write("traj, frame, surface, in-out-below, id, nearest_o \n")

for i in xrange(1,num_lines_1):
    dist_list_parts = h_dist_db_list[i].split()
    nb_list_parts = h_nb_db_list[i].split()

    location = 1

    if float(dist_list_parts[4]) > above_p:
        location = 2
    elif float(dist_list_parts[4]) < below_p:
        location = 0

    h_position.write(nb_list_parts[0] + "," + nb_list_parts[1] + "," + nb_list_parts[2] + "," + str(location) + "," + nb_list_parts[3] + "," + nb_list_parts[4] + "\n")