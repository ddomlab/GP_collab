from graphdot.kernel.marginalized.starting_probability import StartingProbability, Uniform
from graphdot.codegen.cpptool import cpptype

import numpy as np



@cpptype(p=np.float32)
class Additive_p(StartingProbability):
    def __init__(self, **kwargs):
        """
        kwargs: Dictionary of feature names and their corresponding probability kernels.
        Example: Additive_p(element=AssignProb, charge=AssignProb)
        """
        self._names = list(kwargs.keys())
        self.probabilities = list(kwargs.values())
        
        # Crucial: Set attributes so that C++ 'self.name.p' resolution works
        for name, kernel in kwargs.items():
            setattr(self, name, kernel)

    def __call__(self, nodes):
        # r[0] is the probability array p, r[1] is the gradient matrix d_p
        results = [p(nodes) for p in self.probabilities]
        
        # Sum the probabilities across all kernels
        p_total = np.sum([r[0] for r in results], axis=0)
        
        # Stack the gradients vertically (each row is a hyperparameter)
        # Handle cases where d_p might be empty for fixed parameters
        d_p_list = [r[1] for r in results if r[1].size > 0]
        d_p_total = np.vstack(d_p_list) if d_p_list else np.empty((0, len(nodes)))
        
        return p_total, d_p_total

    def gen_expr(self):
        f_exprs = []
        all_jacs = []
        
        for i, p in enumerate(self.probabilities):
            # We use 'n' as the node variable to match the CUDA template 'operator()(N const &n)'
            scope = f'self.{self._names[i]}.'
            
            # Use an adaptive call to handle different microkernel signatures
            try:
                # Try the older signature: gen_expr(x, theta_scope)
                f, j = p.gen_expr('n', theta_scope=scope)
            except TypeError:
                # Fallback to the modern signature: gen_expr()
                f, j = p.gen_expr()
            
            f_exprs.append(f)
            if j is not None:
                all_jacs.extend(j)
        
        combined_f = " + ".join([f"({e})" for e in f_exprs])
        return combined_f, all_jacs

    @property
    def theta(self):
        # Concatenate log-scale hyperparameters
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
        # Concatenate 2D bounds arrays
        return np.concatenate([np.atleast_1d(p.bounds) for p in self.probabilities])
    