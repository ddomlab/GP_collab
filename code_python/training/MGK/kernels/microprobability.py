from graphdot.kernel.marginalized.starting_probability import StartingProbability, Uniform
from graphdot.codegen.cpptool import cpptype

import numpy as np



@cpptype
class Additive_p(StartingProbability):
    def __init__(self, **kwargs):
        self._names = list(kwargs.keys())
        # Based on your JSON, kwargs will contain {'AtomicNumber': <Const_p object>}
        for name, kernel in kwargs.items():
            setattr(self, name, kernel)

    def gen_expr(self):
        f_exprs = []
        all_jacs = []
        for name in self._names:
            # Match the attribute name to the C++ scope
            p = getattr(self, name)
            scope = f'self.{name}.'
            try:
                f, j = p.gen_expr('n', theta_scope=scope)
            except TypeError:
                f, j = p.gen_expr()
            
            f_exprs.append(f)
            if j is not None:
                all_jacs.extend(j)
        
        return " + ".join([f"({e})" for e in f_exprs]), all_jacs

    @property
    def theta(self):
        return np.concatenate([np.atleast_1d(getattr(self, n).theta) for n in self._names])

    @theta.setter
    def theta(self, value):
        start = 0
        for n in self._names:
            p = getattr(self, n)
            size = len(np.atleast_1d(p.theta))
            p.theta = value[start:start + size]
            start += size

    @property
    def bounds(self):
        return np.concatenate([np.atleast_1d(getattr(self, n).bounds) for n in self._names])
    