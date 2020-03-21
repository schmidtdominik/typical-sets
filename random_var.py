import math
from typing import List
import numpy as np
import itertools

class RandomVariable():
    def __init__(self, distribution: List):
        self.values = np.arange(0, len(distribution))
        self.distribution = np.array(distribution)
        assert np.abs(np.sum(self.distribution)-1) < 1e-3

    @property
    def entropy(self):
        return -np.sum(self.distribution * np.log(self.distribution)/np.log(2))

    def p(self, seq):
        return np.product(self.distribution[np.array(seq)])

    def get_seq(self, n):
        return np.random.choice(self.values, p=self.distribution, size=n)

    def get_all_sequences(self, n):
        return itertools.product(*[self.values for i in range(n)])

    def partition_weakly_typical_set(self, seq_set, epsilon):
        wts = set()
        comp = set()
        for seq in seq_set:
            empirical_entropy = (-1/len(seq) * math.log(self.p(seq))/math.log(2))
            true_entropy = self.entropy
            if abs(empirical_entropy - true_entropy) < epsilon:
                wts.add(seq)
            else:
                comp.add(seq)

        return wts, comp