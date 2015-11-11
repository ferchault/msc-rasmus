#!/usr/bin/env python
import numpy as np
import itertools as it
import scipy.optimize as sco

__author__ = 'Guido von Rudorff'
""" Implements a log-likelihood estimator for an arbitrary n-dimensional Ising model."""


class IsingEstimator(object):
	def __init__(self, adjacency_matrix):
		self._adj = np.array(adjacency_matrix)
		# symmetric matrix with zeros on the diagonal
		assert(np.all(self._adj == self._adj.T))
		assert(np.linalg.norm(np.diag(self._adj)) == 0)
		self._site_count = self._adj.shape[0]
		self._patterns = list()
		self._cached = True
		self._frequencies = list()

	def feed_pattern(self, pattern, frequency):
		assert(set(pattern) <= {-1, 1})
		assert(len(pattern) == self._site_count)
		assert(frequency > 0)
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

	def _calc_A_from_pattern(self, pattern):
		A = 0
		for site in range(self._site_count):
			A += pattern[site] * self._adj[site, :] * pattern
		return np.sum(A)/2

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
		alpha, beta = map(lambda x: np.sqrt(float(x)*float(x)), alphabeta)
		return -self.loglikelihood(alpha, beta)

	def gradloglikelihood(self, alpha, beta):
		base = np.exp(-alpha * self._global_A_cache - beta * self._global_B_cache)
		q = np.sum(self._global_A_cache * base)
		qprime = np.sum(self._global_B_cache * base)
		z = self._calc_z(alpha, beta)
		N = np.sum(self._frequencies)
		grad_a = q*N/z - np.sum(self._frequencies * self._pattern_A_cache)
		grad_b = qprime*N/z - np.sum(self._frequencies * self._pattern_B_cache)
		return np.array((grad_a, grad_b))

	def optimize_gradloglikelihood(self, alphabeta):
		alpha, beta = map(lambda x: np.sqrt(float(x)*float(x)), alphabeta)
		return -self.gradloglikelihood(alpha, beta)

def testIsingEstimator():
	adjmat = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]])
	ie = IsingEstimator(adjmat)
	ie.feed_pattern((1, 1, 1), 1)
	ie.feed_pattern((-1, -1, 1), 2)
	ie.update_cache()

	# test calculation of site interaction (A) and caching thereof
	assert(np.all(ie._pattern_A_cache == np.array([2, 0])))

	# test calculation of site sum (B) and caching thereof
	assert(np.all(ie._pattern_B_cache == np.array([3, -1])))

	# test partition sum and gradient thereof
	delta = 0.000001
	assert(abs(ie._calc_z(1., 1.)-31.7011610664) < 0.000001)
	grad_a_finite_differences = (ie._calc_z(1. + delta/2, 1.)-ie._calc_z(1.-delta/2, 1.))/delta
	grad_b_finite_differences = (ie._calc_z(1., 1.+delta/2)-ie._calc_z(1., 1.-delta/2))/delta
	grad_a_analytic = -np.sum(ie._global_A_cache*np.exp( - ie._global_A_cache-ie._global_B_cache))
	grad_b_analytic = -np.sum(ie._global_B_cache*np.exp( - ie._global_A_cache-ie._global_B_cache))
	assert(abs(grad_a_analytic-grad_a_finite_differences) < 0.0000001)
	assert(abs(grad_b_analytic-grad_b_finite_differences) < 0.0000001)

	# test probabilities
	assert(abs(np.exp(ie.logprobabilities(1., 1.)[0])-0.00021255) < 0.0000001)
	assert(abs(np.exp(ie.logprobabilities(1., 1.)[1])-0.08574707) < 0.0000001)

	# test loglikelihood and gradient thereof
	assert(abs(ie.loglikelihood(1., 1.))-13.3690599207 < 0.000001)
	assert(abs(ie.gradloglikelihood(1., 1.)[0]+5.80026440034) < 0.0000001)
	assert(abs(ie.gradloglikelihood(1., 1.)[1]+3.85819419613) < 0.0000001)

	#print sco.minimize(ie.optimize_loglikelihood, (1., 2.), jac=ie.optimize_gradloglikelihood)

testIsingEstimator()


