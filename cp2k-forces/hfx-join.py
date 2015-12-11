#!/usr/bin/env python
import pandas as pd

basepath = '/home/rasmus/ownCloud/UCL/fourth Year/Project/data/cp2k/monomer_ethelyne/'

forces = basepath + 'tmp_forces.txt'
elements = basepath + 'elements.txt'
basis = basepath + 'basis.txt'

elements = pd.read_csv(elements)
basis = pd.read_csv(basis)
forces = pd.read_csv(forces, names='aatom batom catom datom aset bset cset dset ma mb mc md spin_channel acton fx fy fz'.split())

basis = basis.drop('basis', 1)
basis = basis.drop('set', 1)

m1 = forces
for col in 'abcd':
	m1 = pd.merge(m1, elements, 'left', left_on=(col + 'atom'), right_on='atom')
	m1 = m1.drop('atom', 1)
	m1.rename(columns={'element': (col + 'element')}, inplace=True)
print m1.columns

m2 = m1
for col in 'abcd':
	m2 = pd.merge(m2, basis, 'left', left_on=(col + 'element', col + 'set'), right_on=('element', 'iset'))
	m2 = m2.drop(col + 'set', 1)
	m2 = m2.drop('iset', 1)
	m2 = m2.drop('element', 1)
	m2.rename(columns={'n': (col + 'n'), 'orb_name': (col + 'l')}, inplace=True)

m2.to_csv(basepath + "joined-data.csv")
#m2 = pd.merge(m1, elements, 'left', left_on='batom', right_on='atom')
#m1 = pd.merge(m2, elements, 'left', left_on='catom', right_on='atom')
#m2 = pd.merge(m1, elements, 'left', left_on='datom', right_on='atom')
#print m2