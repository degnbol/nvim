import nanobind

ctorus: nanobind.nb_func
"""
ctorus(n_samples: int, dim: int, sample: str = 'random') -> numpy.ndarray[dtype=float64]

Generate random i.i.d. points on a d-torus in R^2d or as a grid

:param n_samples: The number of points to be generated.
:type n_samples: integer
:param dim: The dimension of the torus on which points would be generated in R^2*dim.
:type dim: integer
:param sample: The sample type. Available values are: `"random"` and `"grid"`. Default value is `"random"`.
:type sample: string
:returns: the generated points on a torus.

The shape of returned numpy array is:

If sample is 'random': (n_samples, 2*dim).

If sample is 'grid': (⌊n_samples**(1./dim)⌋**dim, 2*dim), where shape[0] is rounded down to the closest perfect 'dim'th power.
        
"""
sphere: nanobind.nb_func
"""
sphere(n_samples: int, ambient_dim: int, radius: float = 1.0, sample: str = 'random') -> numpy.ndarray[dtype=float64]

Generate random i.i.d. points uniformly on a (d-1)-sphere in R^d

:param n_samples: The number of points to be generated.
:type n_samples: integer
:param ambient_dim: The ambient dimension d.
:type ambient_dim: integer
:param radius: The radius. Default value is `1.`.
:type radius: float
:param sample: The sample type. Default and only available value is `"random"`.
:type sample: string
:returns: the generated points on a sphere.
          
"""
