import numpy as np
import scipy.linalg as linalg
from abc import ABC, abstractmethod

class SubstitutionModel(ABC):
	"""Base class for substitution models"""
	@abstractmethod
	def get_states(self):
		"""Gets all possible states."""
		pass

	@abstractmethod
	def get_pi(self):
		"""Returns equilibrium state frequencies."""
		pass

	@abstractmethod
	def get_q(self):
		"""Returns the normalized instantaneous rate matrix Q."""
		pass

	def get_p(self, t):
		"""Returns the transition probability matrix P(t)."""
		q = self.get_q()
		return linalg.expm(q * t)

	def __init__(self):
		pass

class SiFit3(SubstitutionModel):
	"""SiFit substitution model for ternary data from Zafar et al (2017).

	Parameters:
	lambda_d -- rate of deleteion (default 1.0)
	lambda_l -- rate of loss of heterozygousity (default 1.0)
	"""
	STATES = ['0', '1', '2']

	lambda_d = 1.0
	lambda_l = 1.0

	def get_states(self):
		return self.STATES

	def get_pi(self):
		d = self.lambda_d
		l = self.lambda_l
		x = 1 + 2 / (d + l) + 1 / d
		pi_0 = 1 / x
		pi_1 = 2 / (x * (d + l))
		pi_2 = 1 / (x * d)
		return np.array([pi_0, pi_1, pi_2])

	def get_beta(self):
		d = self.lambda_d
		l = self.lambda_l
		pi = self.get_pi()
		diag = np.array([1, d + l, d])
		return 1 / np.matmul(pi, diag)

	def get_q(self):
		d = self.lambda_d
		l = self.lambda_l
		q = np.matrix([
			[-1, 1, 0],
			[(d + l)/2, -(d + l), (d + l)/2],
			[0, d, -d]])
		beta = self.get_beta()
		return beta * q

	def __init__(self, lambda_d, lambda_l):
		if lambda_d > 0 and lambda_l > 0:
			self.lambda_d = lambda_d
			self.lambda_l = lambda_l
		else:
			raise Exception("SiFit3.init() lambda_d and lambda_l parameters must be greater than zero.")
