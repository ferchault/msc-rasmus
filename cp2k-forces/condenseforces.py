__author__ = 'rasmus'
import sys
import numpy as np
from math import log 
log10 = np.log(10) 

in_force_path = sys.argv[1]
out_force_path = sys.argv[2]

first = True
def coord_to_xyz(coord):
      return (int(coord)-1)%3
f = open(in_force_path)
out = open(out_force_path, "w")

force_abcd = dict()
force_abcd[0] = dict()
force_abcd[1] = dict()
force_abcd[2] = dict()
force_abcd[3] = dict()

for i in [0,1,2,3]:
    for j in [0,1,2]:
        force_abcd[i][j]=""

prev = 0
for line in f:
    line_parts = line.strip().split()
    cur = line_parts[0:13]
    if cur != prev:
        if first:
            first = False
        elif not first:
            info_string = ",".join(prev)
            for abcd in [0,1,2,3]:
                force_part = "," + prev[abcd]
                xyz = ""
                count = 0
                for xyz in [0,1,2]:
                    try:
                        force_part += "," + force_abcd[abcd][xyz]
                    except:
                        force_part += "," + ""

                    if force_abcd[abcd][xyz] == "":
                        count = count + 1
                    force_abcd[abcd][xyz] = ""
                line_out = info_string + force_part + "\n"
                if count != 3:
                    out.write(line_out)
    coord = line_parts[13]
    force = line_parts[14]
    cart_coord = coord_to_xyz(coord)
    atom_abcd = (int(coord)-1) / 3
    #force_abcd[atom_abcd][cart_coord] = str(int(round(log(float(force),10))))
    force_abcd[atom_abcd][cart_coord] = force
    prev = cur


info_string = ",".join(prev)
for abcd in [0,1,2,3]:
    force_part = "," + prev[abcd]
    xyz = ""
    count = 0 
    for xyz in [0,1,2]:
        if force_abcd[abcd][xyz] == "":
            count = count + 1
        force_part += "," + force_abcd[abcd][xyz]
    line_out = info_string + force_part + "\n"
    if count != 3:
        out.write(line_out)

