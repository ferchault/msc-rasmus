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

# take preprocessed list of patterns found and build model
patterns = {
	'110100100100': 1425, '010010101100': 117, '001011101100': 34, '001011101101': 10, '110100110100': 2338,
	'110110101100': 2, '111010101101': 22, '111010101100': 367, '110100110000': 4, '110100111110': 3,
	'111100100100': 21, '001010101100': 83, '110100101100': 77, '000101101101': 10, '110100100110': 2200,
	'010101101100': 6, '010101101101': 1828, '111110100100': 4, '010101100101': 2398, '011010101100': 3687,
	'011010101101': 16, '110110100100': 204, '110100110110': 2282, '110010101100': 106, '110100111100': 120,
	'011011101101': 19, '011011101100': 4, '110100101110': 27}

hematitematches = []
# start estimation
ie = IsingEstimator(nb_list)
for pattern, count in patterns.iteritems():
	# build equivalent patterns
	a, b, c, d = map(''.join, zip(*[iter(pattern)]*3))
	p1 = ''.join((a, b, c, d))
	assert(p1 == pattern)
	p2 = ''.join((d, c, b, a))
	p3 = ''.join((b, a, d, c))
	p4 = ''.join((c, d, a, b))
	ps = set((p1, p2, p3, p4))
	for elem in ps:
		hematitematches.append(elem)

	# convert from visual to Ising representation
	for elem in ps:
		elem = map(lambda _: 1 if _ == '1' else -1, elem)
		ie.feed_pattern(elem, max(1, int(count/float(len(ps)))))

ie.update_cache()
#print sco.minimize(ie.optimize_loglikelihood, (1., 2.), jac=ie.optimize_gradloglikelihood)
#print ie.loglikelihood( 2.63644577e-01,  -1.39931567e-08), ie.gradloglikelihood( 2.63644577e-01,  -1.39931567e-08)
alpha = 1.
beta = 1.

for i in range(1000):
	curval = ie.loglikelihood(alpha, beta)
	curgrad = ie.gradloglikelihood(alpha, beta)
	steplength = 0.000001
	nextval = ie.loglikelihood(alpha+curgrad[0]*steplength, beta+curgrad[1]*steplength)
	diff = curval-nextval
	#print diff, curgrad
	alpha+=curgrad[0]*steplength
	beta+=curgrad[1]*steplength
	alpha = abs(alpha)
	beta = abs(beta)
	if diff > -0.002:
		steplength /= 10

print 'Optimal Parameters'
print 'alpha', alpha, 'beta', beta

print

print 'probabilities and frequencies for observed patterns'
print '  pattern, freq, log(probability)'
lgps = ie.logprobabilities(alpha, beta)
pentries = ie._patterns
pfreq = ie._frequencies
for p, freq, ps in sorted(zip(pentries, pfreq, lgps), key=lambda x: x[1], reverse=True):
	strpattern = ''.join(map(lambda _: str(_), np.maximum(p, 0)))
	print strpattern, freq, ps


print '50 most probable configurations (roughly top 1%)'
z = ie._calc_z(alpha, beta)
ie.clear()
fullpatternlist = []
for pattern in it.product((-1, 1), repeat=12):
	strpattern = ''.join(map(lambda _: str(_), np.maximum(pattern, 0)))
	a, b, c, d = map(''.join, zip(*[iter(strpattern)]*3))

	# ignore equivalent (in the hematite case) patterns
	if ''.join((d, c, b, a)) in fullpatternlist:
		continue
	if ''.join((b, a, d, c)) in fullpatternlist:
		continue
	if ''.join((c, d, a, b)) in fullpatternlist:
		continue

	ie.feed_pattern(pattern, 1)
	fullpatternlist.append(strpattern)
ie.update_cache()
print 'unique pattern count', len(fullpatternlist)
lgps = ie.logprobabilities(alpha, beta)
assert(len(lgps) == len(fullpatternlist))

#import matplotlib.pyplot as plt
#plt.plot(range(len(lgps)), sorted(lgps))
#plt.show()

#print hematitematches
ordered = sorted(zip(lgps, fullpatternlist), key=lambda x: x[0], reverse=True)[:50]
print 'best/worst probability', max(lgps), min(lgps)
for ps, pattern in ordered:
	print pattern, ps, pattern in hematitematches
