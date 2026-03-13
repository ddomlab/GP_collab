from graphdot.kernel.marginalized.starting_probability import StartingProbability, Uniform
from graphdot.codegen.cpptool import cpptype

import numpy as np



class Additive_p(StartingProbability):
    def __init__(self, **kwargs):
        self._names = list(kwargs.keys())
        for name, kernel in kwargs.items():
            setattr(self, name, kernel)

    @property
    def probabilities(self):
        return [getattr(self, name) for name in self._names]

    def __call__(self, nodes):
        results = [p(nodes) for p in self.probabilities]
        p_total = np.sum([r[0] for r in results], axis=0)
        d_p_list = [r[1] for r in results if r[1].size > 0]
        d_p_total = np.vstack(d_p_list) if d_p_list else np.empty((0, len(nodes)))
        return p_total, d_p_total

    def gen_expr(self):
        f_exprs = []
        all_jacs = []
        for name in self._names:
            p = getattr(self, name)
            scope = f'self.{name}.'
            f, j = p.gen_expr() 
            f_exprs.append(f"{scope}{f}")
            all_jacs.extend(j)
        return " + ".join([f"({e})" for e in f_exprs]), all_jacs

    @property
    def theta(self):
        return np.concatenate([np.atleast_1d(p.theta) for p in self.probabilities])

    @theta.setter
    def theta(self, value):
        start = 0
        for p in self.probabilities:
            n = len(np.atleast_1d(p.theta))
            p.theta = value[start:start + n]
            start += n

    @property
    def bounds(self):
        return np.concatenate([np.atleast_1d(p.bounds) for p in self.probabilities])

# Apply the decorator MANUALLY here
Additive_p = cpptype(Additive_p)