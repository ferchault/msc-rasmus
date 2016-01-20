#!/usr/bin/env python
import numpy as np
import itertools as it
#import scipy.optimize as sco

__author__ = 'Guido von Rudorff'


class IsingEstimator(object):
	""" Implements a log-likelihood estimator for an arbitrary n-dimensional Ising model."""

	def __init__(self, adjacency_matrix):
		self._adj = np.array(adjacency_matrix)
		# symmetric matrix with zeros on the diagonal
		assert (np.all(self._adj == self._adj.T))
		assert (np.linalg.norm(np.diag(self._adj)) == 0)
		self._site_count = self._adj.shape[0]
		self._patterns = list()
		self._cached = True
		self._frequencies = list()

	def feed_pattern(self, pattern, frequency):
		assert (set(pattern) <= {-1, 1})
		assert (len(pattern) == self._site_count)
		frequency = float(frequency)
		assert (frequency > 0)
		self._patterns.append(pattern)
		try:
			self._frequencies.append(frequency)
		except:
			self._frequencies = list(self._frequencies).append(frequency)
		self._cached = False

	def update_cache(self):
		if self._cached:
			return
		self._build_pattern_cache()
		self._build_global_cache()
		self._frequencies = np.array(self._frequencies)

	def clear(self):
		self._patterns = []
		self._cached = False
		self._frequencies = []

	def _calc_A_from_pattern(self, pattern):
		A = 0
		for site in range(self._site_count):
			A += pattern[site] * self._adj[site, :] * pattern
		return np.sum(A) / 2

	def _build_pattern_cache(self):
		self._pattern_A_cache = np.array([self._calc_A_from_pattern(_) for _ in self._patterns])
		self._pattern_B_cache = np.array([np.sum(_) for _ in self._patterns])

	def _build_global_cache(self):
		self._global_A_cache = []
		self._global_B_cache = []
		for pattern in it.product((-1, 1), repeat=self._site_count):
			self._global_A_cache.append(self._calc_A_from_pattern(pattern))
			self._global_B_cache.append(np.sum(pattern))
		self._global_A_cache = np.array(self._global_A_cache)
		self._global_B_cache = np.array(self._global_B_cache)

	def _calc_z(self, alpha, beta):
		return np.sum(np.exp(-alpha * self._global_A_cache - beta * self._global_B_cache))

	def logprobabilities(self, alpha, beta):
		return (-alpha * self._pattern_A_cache - beta * self._pattern_B_cache) - np.log(self._calc_z(alpha, beta))

	def loglikelihood(self, alpha, beta):
		return np.sum(self._frequencies * self.logprobabilities(alpha, beta))

	def optimize_loglikelihood(self, alphabeta):
		alpha, beta = map(lambda x: np.sqrt(float(x) * float(x)), alphabeta)
		return -self.loglikelihood(alpha, beta)

	def gradloglikelihood(self, alpha, beta):
		base = np.exp(-alpha * self._global_A_cache - beta * self._global_B_cache)
		q = np.sum(self._global_A_cache * base)
		qprime = np.sum(self._global_B_cache * base)
		z = self._calc_z(alpha, beta)
		N = np.sum(self._frequencies)
		grad_a = q * N / z - np.sum(self._frequencies * self._pattern_A_cache)
		grad_b = qprime * N / z - np.sum(self._frequencies * self._pattern_B_cache)
		return np.array((grad_a, grad_b))

	def optimize_gradloglikelihood(self, alphabeta):
		alpha, beta = map(lambda x: np.sqrt(float(x) * float(x)), alphabeta)
		return -self.gradloglikelihood(alpha, beta)


def testIsingEstimator():
	adjmat = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])
	ie = IsingEstimator(adjmat)
	ie.feed_pattern((1, 1, 1), 1)
	ie.feed_pattern((-1, -1, 1), 2)
	ie.update_cache()

	# test calculation of site interaction (A) and caching thereof
	assert (np.all(ie._pattern_A_cache == np.array([2, 0])))

	# test calculation of site sum (B) and caching thereof
	assert (np.all(ie._pattern_B_cache == np.array([3, -1])))

	# test partition sum and gradient thereof
	delta = 0.000001
	assert (abs(ie._calc_z(1., 1.) - 31.7011610664) < 0.000001)
	grad_a_finite_differences = (ie._calc_z(1. + delta / 2, 1.) - ie._calc_z(1. - delta / 2, 1.)) / delta
	grad_b_finite_differences = (ie._calc_z(1., 1. + delta / 2) - ie._calc_z(1., 1. - delta / 2)) / delta
	grad_a_analytic = -np.sum(ie._global_A_cache * np.exp(- ie._global_A_cache - ie._global_B_cache))
	grad_b_analytic = -np.sum(ie._global_B_cache * np.exp(- ie._global_A_cache - ie._global_B_cache))
	assert (abs(grad_a_analytic - grad_a_finite_differences) < 0.0000001)
	assert (abs(grad_b_analytic - grad_b_finite_differences) < 0.0000001)

	# test probabilities
	assert (abs(np.exp(ie.logprobabilities(1., 1.)[0]) - 0.00021255) < 0.0000001)
	assert (abs(np.exp(ie.logprobabilities(1., 1.)[1]) - 0.08574707) < 0.0000001)

	# test loglikelihood and gradient thereof
	assert (abs(ie.loglikelihood(1., 1.)) - 13.3690599207 < 0.000001)
	assert (abs(ie.gradloglikelihood(1., 1.)[0] + 5.80026440034) < 0.0000001)
	assert (abs(ie.gradloglikelihood(1., 1.)[1] + 3.85819419613) < 0.0000001)


def hematite_isomorphisms(iterable):
	iterable = list(iterable)
	assert(len(iterable) == 12)
	a, b, c, d = zip(*[iter(iterable)]*3)
	isomorphisms = (a + b + c + d, b + a + d + c, c + d + a + b, d + c + b + a)
	if type(a[0]) == str:
		return map(''.join, isomorphisms)
	if type(a[0]) == int:
		return isomorphisms
	raise NotImplementedError('Unknown type %s' % type(a[0]))

def test_hematite_isomorphisms():
	# string parameters
	t1 = hematite_isomorphisms('abcdefghijkl')
	assert(t1 == ['abcdefghijkl', 'defabcjklghi', 'ghijklabcdef', 'jklghidefabc'])

	# integer parameters
	t2 = hematite_isomorphisms([1, 2, 3, 4, 5, 6, 7, 8, 9, 0, -1, -2])
	assert(t2 == ((1, 2, 3, 4, 5, 6, 7, 8, 9, 0, -1, -2), (4, 5, 6, 1, 2, 3, 0, -1, -2, 7, 8, 9), (7, 8, 9, 0, -1, -2, 1, 2, 3, 4, 5, 6), (0, -1, -2, 7, 8, 9, 4, 5, 6, 1, 2, 3)))

def hexagonal_move_up(iterable):
		return [iterable[_] for _ in (5, 0, 10, 2, 3, 7, 11, 6, 4, 8, 9, 1)]

def hexagonal_move_down(iterable):
		return [iterable[_] for _ in (4, 8, 0, 1, 11, 3, 10, 2, 6, 7, 5, 9)]

def hexagonal_unique_transformations(iterable):
	iterable = tuple(iterable)
	assert(len(iterable) == 12)
	def _use_mapping(base, mapping):
		return tuple([base[_-1] for _ in mapping])
	mappings = [
		# mirroring
		(1, 5, 6, 4, 2, 3, 10, 8, 9, 7, 11, 12),
		(4, 2, 3, 1, 5, 6, 7, 11, 12, 10, 8, 9),
		(1, 2, 9, 10, 11, 6, 7, 8, 4, 3, 5, 12),
		(7, 8, 3, 4, 5, 12, 1, 2, 9, 10, 11, 6),
		(10, 2, 12, 4, 8, 6, 7, 5, 9, 1, 11, 3),
		(1, 11, 3, 7, 5, 9, 4, 8, 6, 10, 2, 12),
		# rotation
		(12, 11, 7, 3, 2, 4, 9, 8, 10, 6, 5, 1),
		(12, 2, 10, 9, 5, 7, 6, 8, 4, 3, 11, 1),
		(1, 5, 9, 7, 11, 3, 10, 8, 6, 4, 2, 12),
	]
	matches = [iterable]
	for mapping in mappings:
		base = _use_mapping(iterable[:], mapping)
		while base != iterable:
			matches.append(base)
			base = _use_mapping(base, mapping)

	return set(matches)

def hexagonal_ising_isomorphisms(iterable):
	iterable = list(iterable)
	assert(len(iterable) == 12)

	matches = []
	base = iterable[:]
	for shift_up in range(6):
		base2 = base[:]
		for shift_down in range(6):
			matches += hexagonal_unique_transformations(base2)
			base2 = hexagonal_move_down(base2)
		base = hexagonal_move_up(base)

	return set(matches)

def test_hexagonal_isomorphisms():
	# test lattice vector shifts
	assert(hexagonal_move_up(range(12)) == [5, 0, 10, 2, 3, 7, 11, 6, 4, 8, 9, 1])
	assert(hexagonal_move_down(range(12)) == [4, 8, 0, 1, 11, 3, 10, 2, 6, 7, 5, 9])
	assert(len(hexagonal_ising_isomorphisms(range(12))) == 144)

# run tests
test_hematite_isomorphisms()
test_hexagonal_isomorphisms()
testIsingEstimator()

# build hematite adjacency matrix
## BEGIN CODE FROM RASMUS
nb_list = np.zeros((12, 12))

nb_list[:, 0] = [0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0]
nb_list[:, 1] = [1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1]
nb_list[:, 2] = [1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0]

for j in range(0, 3):
	a, b, c, d = [list(nb_list[i:i + 3, j]) for i in xrange(0, len(nb_list), 3)]
	nb_list[:, j + 3] = b + a + d + c
	nb_list[:, j + 6] = c + d + a + b
	nb_list[:, j + 9] = d + c + b + a
## END CODE FROM RASMUS

# previous analysis results
rawpatterns = {
	'110100100100': 1425, '010010101100': 117, '001011101100': 34, '001011101101': 10, '110100110100': 2338,
	'110110101100': 2, '111010101101': 22, '111010101100': 367, '110100110000': 4, '110100111110': 3,
	'111100100100': 21, '001010101100': 83, '110100101100': 77, '000101101101': 10, '110100100110': 2200,
	'010101101100': 6, '010101101101': 1828, '111110100100': 4, '010101100101': 2398, '011010101100': 3687,
	'011010101101': 16, '110110100100': 204, '110100110110': 2282, '110010101100': 106, '110100111100': 120,
	'011011101101': 19, '011011101100': 4, '110100101110': 27}

# build hematite symmetry analogues
hematitepatterns = dict()
hematitelookup = dict()
for pattern in rawpatterns:
	variants = set()
	for elem in hematite_isomorphisms(pattern):
		variants.update(hexagonal_ising_isomorphisms(elem))
	for variant in variants:
		hematitepatterns[variant] = float(rawpatterns[pattern])/len(variants)
		hematitelookup[''.join(variant)] = pattern
print 'The %d configurations observed map to %d configuration in the hexagonal Ising model' % (len(rawpatterns), len(hematitepatterns))

# prepare estimator
ie = IsingEstimator(nb_list)
for pattern, count in hematitepatterns.iteritems():
	# convert visual to Ising representation
	pattern = map(lambda _: 1 if _ == '1' else -1, pattern)

	ie.feed_pattern(pattern, count)
ie.update_cache()

# estimate
alpha, beta = 1., 1.
for i in range(1000):
	curval = ie.loglikelihood(alpha, beta)
	curgrad = ie.gradloglikelihood(alpha, beta)
	steplength = 0.000001
	nextval = ie.loglikelihood(alpha+curgrad[0]*steplength, beta+curgrad[1]*steplength)
	diff = curval-nextval
	alpha += curgrad[0]*steplength
	beta += curgrad[1]*steplength
	alpha = abs(alpha)
	beta = abs(beta)
	if diff > -0.002:
		steplength /= 10

# output
print 'Optimal Parameters'
print 'alpha', alpha, 'beta', beta
print 'gradient here', ie.gradloglikelihood(alpha, beta)

print

# print 'probabilities and frequencies for observed patterns'
# print '  pattern, freq, log(probability)'
# lgps = ie.logprobabilities(alpha, beta)
# pentries = ie._patterns
# pfreq = ie._frequencies
# for p, freq, ps in sorted(zip(pentries, pfreq, lgps), key=lambda x: x[1], reverse=True):
# 	strpattern = ''.join(map(lambda _: str(_), np.maximum(p, 0)))
# 	print strpattern, freq, ps

print '50 most probable configurations of all 2^12 mapped back to their raw hematite pattern (if it exists)'
ie.clear()
for pattern in it.product((-1, 1), repeat=12):
	ie.feed_pattern(pattern, 1)
ie.update_cache()
lgps = ie.logprobabilities(alpha, beta)
fullpatternlist = [''.join(['1' if _ == 1 else '0' for _ in pattern]) for pattern in ie._patterns]

# print hematitematches
ordered = sorted(zip(lgps, fullpatternlist), key=lambda x: x[0], reverse=True)
printed = list()
print 'best/worst probability', max(lgps), min(lgps)
print '# pattern, log(probability), raw hematite pattern, hematite count'
pvals = []
hemcount = []
for ps, pattern in ordered:
	if pattern in hematitelookup:
		ref = hematitelookup[pattern]
		if ref not in printed:
			print pattern, ps, ref, rawpatterns[ref]
			pvals.append(ps)
			hemcount.append(rawpatterns[ref])
			printed.append(ref)
	else:
		variants = set()
		for elem in hematite_isomorphisms(pattern):
			variants.update(hexagonal_ising_isomorphisms(elem))
		already_printed = False
		for variant in variants:
			if ''.join(variant) in printed:
				already_printed = True
				break
		if not already_printed:
			print pattern, ps, 'NOHEMATITE', 'NOCOUNT'
			pvals.append(ps)
			hemcount.append(0)
			printed.append(pattern)

import matplotlib.pyplot as plt
plt.plot(range(len(pvals)), pvals, '-', color='black', label='Ising Model')
hemhas_idx = []
hemhas_p = []
hemhas_idx_common = []
hemhas_p_common = []
for idx, combo in enumerate(zip(pvals, hemcount)):
	p, h = combo
	if h > 0:
		hemhas_idx.append(idx)
		hemhas_p.append(p)
	if h > 1000:
		hemhas_idx_common.append(idx)
		hemhas_p_common.append(p)
plt.scatter(hemhas_idx, hemhas_p, color='orange', label='Hematite (all)', s=100)
plt.scatter(hemhas_idx_common, hemhas_p_common, color='red', label='Hematite (frequent)', s=50)
plt.xlabel('Index of configuration')
plt.ylabel('Ln (Probability)')
plt.legend()
plt.show()