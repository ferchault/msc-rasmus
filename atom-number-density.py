__author__ = 'rasmus'
import numpy as np
base_path = '/home/rasmus/ownCloud/data/oxygenterm-tray/'

oxygen = open(base_path + 'oxygen_int_dist.out', 'r')
oxygen_water = open(base_path + 'oxygen_water_dist.out', 'r')
hydrogen =open(base_path + 'hydrogen_dist.out','r')
hydrogen_water = open(base_path + 'hydrogen_water_dist.out','r')
iron = open(base_path + 'iron_dist.out', 'r')


#OXYGEN INTERFACE
oxygen_dist = []
AB_frame = dict()

AB_frame['A'] = 0
AB_frame['B'] = 0

first = 1
for line in oxygen:
    if first:
        first = False
        continue
    line_parts = line.split()
    if line_parts[0] == 'A':
        AB_frame['A'] = max(float(AB_frame['A']) , float(line_parts[1]))
    elif line_parts[0] =='B':
        AB_frame['B'] = max(float(AB_frame['B']) , float(line_parts[1]))


    oxygen_dist.append(float(line_parts[4]))

oxygen_dist = np.array(oxygen_dist)
total_frames = AB_frame['A'] + AB_frame['B']
oxygen_hist = np.histogram(oxygen_dist, bins=100, range=(0.0,12.0), density=False)

#HYDROGEN INTERFACE
hydrogen_dist = []
AB_frame = dict()

AB_frame['A'] = 0
AB_frame['B'] = 0

first = 1
for line in hydrogen:
    if first:
        first = False
        continue
    line_parts = line.split()
    if line_parts[0] == 'A':
        AB_frame['A'] = max(float(AB_frame['A']) , float(line_parts[1]))
    elif line_parts[0] =='B':
        AB_frame['B'] = max(float(AB_frame['B']) , float(line_parts[1]))


    hydrogen_dist.append(float(line_parts[4]))

hydrogen_dist = np.array(hydrogen_dist)
total_frames = AB_frame['A'] + AB_frame['B']
hydrogen_hist = np.histogram(hydrogen_dist, bins=100, range=(0.0,12.0), density=False)

#IRON
iron_dist = []
AB_frame = dict()

AB_frame['A'] = 0
AB_frame['B'] = 0

first = 1
for line in iron:
    if first:
        first = False
        continue
    line_parts = line.split()
    if line_parts[0] == 'A':
        AB_frame['A'] = max(float(AB_frame['A']) , float(line_parts[1]))
    elif line_parts[0] =='B':
        AB_frame['B'] = max(float(AB_frame['B']) , float(line_parts[1]))


    iron_dist.append(float(line_parts[4]))

iron_dist = np.array(iron_dist)
total_frames = AB_frame['A'] + AB_frame['B']
iron_hist = np.histogram(iron_dist, bins=100, range=(0.0,12.0), density=False)

#OXYGENWATER
oxygen_water_dist = []
AB_frame = dict()

AB_frame['A'] = 0
AB_frame['B'] = 0

first = 1
for line in oxygen_water:
    if first:
        first = False
        continue
    line_parts = line.split()
    if line_parts[0] == 'A':
        AB_frame['A'] = max(float(AB_frame['A']) , float(line_parts[1]))
    elif line_parts[0] =='B':
        AB_frame['B'] = max(float(AB_frame['B']) , float(line_parts[1]))


    oxygen_water_dist.append(float(line_parts[4]))

oxygen_water_dist = np.array(oxygen_water_dist)
total_frames = AB_frame['A'] + AB_frame['B']
oxygen_water_hist = np.histogram(oxygen_water_dist, bins=100, range=(0.0,12.0), density=False)

#HYDROGENWATER
hydrogen_water_dist = []
AB_frame = dict()

AB_frame['A'] = 0
AB_frame['B'] = 0

first = 1
for line in hydrogen_water:
    if first:
        first = False
        continue
    line_parts = line.split()
    if line_parts[0] == 'A':
        AB_frame['A'] = max(float(AB_frame['A']) , float(line_parts[1]))
    elif line_parts[0] =='B':
        AB_frame['B'] = max(float(AB_frame['B']) , float(line_parts[1]))


    hydrogen_water_dist.append(float(line_parts[4]))

hydrogen_water_dist = np.array(hydrogen_water_dist)
total_frames = AB_frame['A'] + AB_frame['B']
hydrogen_water_hist = np.histogram(hydrogen_water_dist, bins=119, range=(0.0,12.0), density=False)

np.savetxt(base_path + "hist_hydrogen_water.csv",  zip(hydrogen_water_hist[1],hydrogen_water_hist[0]/(4*2*total_frames)), delimiter=',')
np.savetxt(base_path + "hist_hydrogen.csv",  zip(hydrogen_hist[1],hydrogen_hist[0]/(4*2*total_frames)), delimiter=',')
np.savetxt(base_path + "hist_oxygen.csv",  zip(oxygen_hist[1],oxygen_hist[0]/(4*2*total_frames)), delimiter=',')
np.savetxt(base_path + "hist_oxygen_water.csv",  zip(oxygen_water_hist[1],oxygen_water_hist[0]/(4*2*total_frames)), delimiter=',')
np.savetxt(base_path + "hist_iron.csv",  zip(iron_hist[1],iron_hist[0]/(4*2*total_frames)), delimiter=',')

