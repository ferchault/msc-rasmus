__author__ = 'rasmus'
from analysis.class_analysis import *
import scipy, numpy as np

analysis_obj = Analysis('../rasmus/input.psf', '../rasmus/IOHMD-A-prod.dcd', '../rasmus/input.ndx', 'A')

#form 3 by N of coordinates
oxygen_list = analysis_obj.return_id_list_from_name("oxygen_A")
hydrogen_list = analysis_obj.return_id_list_from_name("hydrogen_top")

analysis_obj.start_interface_h_analysis(oxygen_list, hydrogen_list, "t", 0, analysis_obj.u.trajectory.n_frames, "test.out")