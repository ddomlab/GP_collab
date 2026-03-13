from graphdot.kernel.marginalized.starting_probability import StartingProbability, Uniform
from graphdot.codegen.cpptool import cpptype

import numpy as np


@cpptype(_placeholder=np.int32)
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
        exprs = []
        jacs = []
        for name in self._names:
            p = getattr(self, name)
            f, j = p.gen_expr()  # get C++ variable (e.g., "p") and jacobian
            exprs.append(f"({f})")  # just wrap in parentheses
            jacs.extend(j)
        return " + ".join(exprs), jacs

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
    





def ConstantStartingProbability(c, c_bounds='fixed'):
    r"""Creates a starting probability that assigns the same constant
    probability to each node.

    Parameters
    ----------
    c: float > 0
        The constant value assigned to all nodes.
    """

    class ConstantStarting(StartingProbability):
        @property
        def name(self):
            return 'ConstantStarting'

        def __init__(self, c, c_bounds):
            self.c = float(c)
            self.c_bounds = c_bounds
            self._assert_bounds('c', c_bounds)

        def __call__(self, nodes):
            n_nodes = len(nodes)
            # Return constant array and gradient
            p = np.full(n_nodes, self.c, dtype=np.float32)
            d_p = np.ones((1, n_nodes), dtype=np.float32)  # gradient w.r.t c
            return p, d_p

        def __repr__(self):
            return f'{self.name}({self.c})'

        def gen_expr(self):
            # C++ expression: constant value and derivative 1.0
            return f'c', ['1.0f']

        @property
        def theta(self):
            return np.array([self.c], dtype=np.float32)

        @theta.setter
        def theta(self, seq):
            self.c = seq[0]

        @property
        def bounds(self):
            return (self.c_bounds,)

        def _assert_bounds(self, hyp, bounds):
            if not ((isinstance(bounds, tuple) and len(bounds) == 2)
                    or bounds == 'fixed'):
                raise ValueError(
                    f'Bounds for hyperparameter {hyp} of probability {self.name} '
                    f'must be a 2-tuple or "fixed": {bounds} provided.'
                )

    return ConstantStarting(c, c_bounds)


