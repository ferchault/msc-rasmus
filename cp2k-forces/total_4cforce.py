__author__ = 'rasmus'
##sums the cp2k force array in x y z components
from config import base_path

force_4c_totals = open(base_path + "total_forces.txt")

force_list = []
for line in force_4c_totals:
    parts = line.strip().split()
    force_list += parts[1:]

