base_path = "/home/rasmus/ownCloud/UCL/fourth Year/Project/data/cp2k/monomer_ethelyne/"
file_path = base_path + "data.txt"
output_file_basis_path = base_path + "basis_txt"
output_file_forces_path =base_path + "forces.txt"

input_file = open(file_path, "r")
output_file_forces = open(output_file_forces_path, "w")
output_file_basis = open(output_file_basis_path, "w")

prefix_forces = "forces"
prefix_orb = "info"

output_file_basis.write("set,n, orb_name\n")
output_file_forces.write("aatom, batom, catom, datom, aset, bset, cset, dset, ma, mb, mc, md, spin_channel, coord, force\n")
for line in input_file:
    if line != "":
        line_parts = line.strip().split()
        if len(line_parts) != 0:
            if line_parts[0] == prefix_forces:
                for part in line_parts[1:-1]:
                    output_file_forces.write(part)
                    output_file_forces.write(",")
                output_file_forces.write(line_parts[-1])
                output_file_forces.write("\n")
            elif line_parts[0] == prefix_orb:
                for part in line_parts[1:-1]:
                    output_file_basis.write(part)
                    output_file_basis.write(",")
                output_file_basis.write(line_parts[-1])
                output_file_basis.write("\n")

