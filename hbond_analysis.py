__author__ = 'rasmus'
from analysis.class_analysis import *

test = Analysis('../rasmus/input.psf', '../rasmus/IOHMD-A-prod.dcd', '../rasmus/input.ndx', 'A')
test.set_parameters(4.0, 1.05, 145)

oxygen_list = test.return_id_list_from_name("oxygen_A")
oxygen_list += test.return_id_list_from_name("oxygen_G")

hydrogen_list = test.return_id_list_from_name("hydrogen_top")
hydrogen_list += test.return_id_list_from_name("hydrogen_bottom")

water_id_list = test.return_id_list_from_name("water")

oxygen_list += test.get_water_oxygen(water_id_list)
hydrogen_list += test.get_water_hydrogen(water_id_list)
test.start_h_bond_analysis(oxygen_list, hydrogen_list, 0, test.u.trajectory.n_frames, 'run.out')

test2 = Analysis('../rasmus/input.psf', '../rasmus/IOHMD-B-prod.dcd', '../rasmus/input.ndx', 'B')
test2.set_parameters(4.0, 1.05, 145)

test2.start_h_bond_analysis(oxygen_list, hydrogen_list, 0, test2.u.trajectory.n_frames, 'run.out', True)
