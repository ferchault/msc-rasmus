#This tests whether our understanding of the force components is correct
#Compares the total hartree fock exchange to the individual components
#force(kind) % fock_4c(atom_kind_id , xyz)  =? sum of individual components
import pandas as pd

base_path = "C:/Users/Rasmus/ownCloud/UCL/4year/Project/data/cp2k/monomer_ethelyne/"
data = "joined-data.csv"
total_force_data = "total_force.txt"

forces = pd.read_csv(base_path+data)
total_force_4c = pd.read_csv(base_path+total_force_data)

total_summed_forces = forces.groupby(["acton"], as_index=False).sum()

for i in xrange(1, len(total_summed_forces["acton"])+1):
    print "element " + str(i)
    print "fx_summed fx_4c "
    print -total_summed_forces["fx"][i-1] , total_force_4c["fx"][i-1]

    print "fy_summed fy_4c "
    print -total_summed_forces["fy"][i-1] , total_force_4c["fy"][i-1]

    print "fy_summed fz_4c "
    print -total_summed_forces["fz"][i-1] , total_force_4c["fz"][i-1]

    print "abs diff x " + str((-total_summed_forces["fx"][i-1] - total_force_4c["fx"][i-1]))
    print "rel diff x " + str((-total_summed_forces["fx"][i-1] - total_force_4c["fx"][i-1])/abs(total_force_4c["fx"][i-1]))

    print "abs diff y " + str((-total_summed_forces["fy"][i-1] - total_force_4c["fy"][i-1]))
    print "rel diff y " + str((-total_summed_forces["fy"][i-1] - total_force_4c["fy"][i-1])/abs(total_force_4c["fy"][i-1]))

    print "abs diff z " + str((-total_summed_forces["fz"][i-1] - total_force_4c["fz"][i-1]))
    print "rel diff z " + str((-total_summed_forces["fz"][i-1] - total_force_4c["fz"][i-1])/abs(total_force_4c["fz"][i-1]))