__author__ = 'rasmus'
from analysis.class_analysis import *

base_path = "/home/rasmus/ownCloud/data/iron-term-trajectory/"

traj_A = Analysis(base_path + 'IO3-index.xyz', base_path + 'IO3-A.dcd', base_path + 'input.ndx', 'A')
traj_A.set_parameters(3.5, 160)

oxygen_list = traj_A.return_id_list_from_name("oxygen_A")
oxygen_list += traj_A.return_id_list_from_name("oxygen_G")

hydrogen_list = traj_A.return_id_list_from_name("hydrogen_top")
hydrogen_list += traj_A.return_id_list_from_name("hydrogen_bottom")

water_id_list = traj_A.return_id_list_from_name("water")

oxygen_list += traj_A.get_water_oxygen(water_id_list)
hydrogen_list += traj_A.get_water_hydrogen(water_id_list)
traj_A.start_h_bond_analysis(oxygen_list, hydrogen_list, 0, traj_A.u.trajectory.n_frames, 'run.out')

traj_B = Analysis(base_path + 'IO3-index.xyz', base_path + 'IO3-B.dcd', base_path + 'input.ndx', 'B')
traj_B.set_parameters(3.5, 160)

traj_B.start_h_bond_analysis(oxygen_list, hydrogen_list, 0, traj_B.u.trajectory.n_frames, 'run.out', True)
