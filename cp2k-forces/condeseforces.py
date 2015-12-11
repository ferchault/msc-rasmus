__author__ = 'rasmus'

base_path = "/home/rasmus/ownCloud/UCL/fourth Year/Project/data/cp2k/monomer_ethelyne/"
in_force_path = base_path + "temp.forces"
out_force_path = base_path + "tmp_forces.txt"

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

prev = 0
for line in f:
    line_parts = line.strip().split(",")
    cur = line_parts[0:13]
    if line_parts[0] == "aatom":
        continue
    if cur != prev:
        if first:
            first = False
        elif not first:
            info_string = ",".join(prev)
            for abcd in [0,1,2,3]:
                force_part = "," + prev[abcd]
                xyz = ""
                for xyz in [0,1,2]:
                    force_part += "," + force_abcd[abcd][xyz]
                    force_abcd[abcd][xyz] = "0.0"
                line_out = info_string + force_part + "\n"
                out.write(line_out)
    coord = line_parts[13]
    force = line_parts[14]
    cart_coord = coord_to_xyz(coord)
    atom_abcd = (int(coord)-1) / 3
    force_abcd[atom_abcd][cart_coord] = force
    prev = cur


info_string = ",".join(prev)
for abcd in [0,1,2,3]:
    force_part = "," + prev[abcd]
    xyz = ""
    for xyz in [0,1,2]:
        force_part += "," + force_abcd[abcd][xyz]
    line_out = info_string + force_part + "\n"
    out.write(line_out)