from math import log, floor, sqrt, ceil
import random
import json
import numpy


def py2_sample(population, k):
    """Chooses k unique random elements from a population sequence.

    This implementation replicates the behavior of Python2's random.sample
    in Python3.
    """

    if not isinstance(population, (list, tuple)):
        # Convert sets and dicts to list (Python2 allowed sampling from sets)
        population = list(population)

    n = len(population)
    if not 0 <= k <= n:
        raise ValueError("Sample larger than population or is negative")

    result = [None] * k
    setsize = 21
    if k > 5:
        setsize += 4 ** ceil(log(k * 3, 4))

    if n <= setsize:
        pool = list(population)
        for i in range(k):
            j = int(random.random() * (n - i))
            result[i] = pool[j]
            pool[j] = pool[n - i - 1]
    else:
        try:
            selected = set()
            for i in range(k):
                j = int(random.random() * n)
                while j in selected:
                    j = int(random.random() * n)
                selected.add(j)
                result[i] = population[j]
        except (KeyError, IndexError):
            raise ValueError("Sample larger than population or is negative")

    return result


class PRNG(object):
    def __init__(self, K, delta, c, np = None):
        self.K = float(K)
        self.K_int = int(K)
        self.delta = delta
        self.c = c

        self.S = self.c * log(self.K/self.delta) * sqrt(self.K)
        self.cdf, self.Z = self._gen_rsd_cdf(K, self.S, self.delta)
        # seed:
        self.state = 1
        # nprand?
        self.np_rand = numpy.random.RandomState(1)
        self.np = np

    def _gen_rsd_cdf(self,K, S, delta):
        pivot = int(floor(K/S))
        val1 =  [S/K * 1/d for d in range(1, pivot)]
        val2 =  [S/K * log(S/delta)]
        val3 =  [0 for d in range(pivot, K)]
        tau = val1 + val2 + val3
        rho = [1.0/K] + [1.0/(d*(d-1)) for d in range(2, K+1)]
        Z = sum(rho) + sum(tau)
        mu = [(rho[d] + tau[d])/Z for d in range(K)]
        cdf = numpy.cumsum(mu)
        return cdf, Z

    def get_S(self):
        return self.S
    def get_state(self):
        return self.state

    def set_seed(self, seed):
        self.state = seed
    def get_src_blocks_wrap(self, seed=None):
        if seed:
            self.state = seed
        if self.np:
            self.np_rand.seed(self.state)
            p = self.np_rand.rand()
            d = self._sample_d(p)
            nums = self.np_rand.randint(0, self.K_int, d)
        else:
            random.seed(self.state)
            p = random.random()
            d = self._sample_d(p)
            nums = py2_sample(range(self.K_int), d)
        return d, nums
    def _sample_d(self,p):
        for ix, v in enumerate(self.cdf):
            if v > p:
                return ix + 1
        return ix + 1

    def debug(self):
        return json.dumps({'K': self.K, 'delta': self.delta, 'c' : self.c, 'S' : self.S,
                       #'cdf': self.cdf,
                           'Z': self.Z , 'K_prime': self.K * self.Z})
