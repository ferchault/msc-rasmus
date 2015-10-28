__author__ = 'guido'
""" Converts the output of hinterface_analysis.py into a XYZ file for visualisation purposes."""

# system modules
import sys
# third-party modules
import MDAnalysis as mda 
# custom modules
import analysis.class_analysis as ca

# parameters
if len(sys.argv) != 6:
	print 'Usage: %s INFILE NDXFILE PSFFILE DCDFILE THRESHOLD' % sys.argv[0]
	exit(1)
infile = sys.argv[1]
ndxfile = sys.argv[2]
psffile = sys.argv[3]
dcdfile = sys.argv[4]
threshold = float(sys.argv[5])

# coordinates
u = mda.Universe(psffile, dcdfile)
ndxdata = ca.Analysis.ndx_file_read_to_list(ndxfile)
coordinates = dict()
for group in ndxdata:
	if group.startswith('termination_'):
		O_idx, H_idx = ndxdata[group]
		coordinates[H_idx] = u.atoms.positions[O_idx]

# convert data
first = True
lastframe = -1
for line in open(infile):
	if first:
		first = False
		continue
	parts = line.strip().split()
	frame = int(parts[1])
	H_idx = int(parts[3])
	dist = float(parts[4])
	if dist >= threshold:
		element = 'O' # out of plane
	else:
		element = 'P' # in plane

	if frame != lastframe:
		print len(coordinates)/2
		print 'Comment line'
		lastframe = frame

	print element, coordinates[H_idx][0], coordinates[H_idx][1], int(dist >= threshold)
