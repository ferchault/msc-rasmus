from config import base_path

file_path = base_path + "data.txt"

output_file_basis_path = base_path + "basis_txt"
output_file_forces_path = base_path + "forces.txt"
output_file_total_force_path = base_path + "total_force.txt"

input_file = open(file_path, "r")
output_file_forces = open(output_file_forces_path, "w")
output_file_basis = open(output_file_basis_path, "w")
output_file_total_force = open(output_file_total_force_path, "w")

prefix_forces = "forces"
prefix_orb = "info"
prefix_total_force = "total_force"

output_file_basis.write("set,n, orb_name\n")
output_file_forces.write("aatom, batom, catom, datom, aset, bset, cset, dset, ma, mb, mc, md, spin_channel, coord, force\n")
output_file_total_force.write("kind,1,2,3,4,5,6,7,8,9,10,11,12\n")

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
            elif line_parts[0] == prefix_total_force:
                for part in line_parts[1:-1]:
                    output_file_total_force.write(part)
                    output_file_total_force.write(",")
                output_file_total_force.write(line_parts[-1])
                output_file_total_force.write("\n")


