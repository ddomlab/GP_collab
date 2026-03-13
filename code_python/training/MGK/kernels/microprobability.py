from graphdot.kernel.marginalized.starting_probability import StartingProbability, Uniform
from graphdot.codegen.cpptool import cpptype

import numpy as np



@cpptype
class Additive_p(StartingProbability):
    def __init__(self, **kwargs):
        self._names = list(kwargs.keys())
        # Do NOT store self.probabilities as an attribute here
        for name, kernel in kwargs.items():
            setattr(self, name, kernel)

    @property
    def probabilities(self):
        # This allows your __call__, theta, and bounds to still work
        # without confusing the @cpptype serializer
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
            
            # Since Uniform.gen_expr() returns 'p', 
            # this correctly creates 'self.AtomicNumber.p'
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