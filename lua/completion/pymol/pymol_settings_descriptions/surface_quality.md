This controls how well PyMOL draws surfaces. Lower values, like 0, are
rough surface approximations. These low values are good for
speed\--especially for larger surfaces. For rendering of publication
quality photos, and truer representations of the biological surface, set
the value higher\--to something like 2, 3 or 4. In practice typical
values are 1, 2 and 3.

This value, I believe, changes how finely the polygons are sampled for
surface representation. Settings of 2 or 3 are far more computationally
intensive to compute and show than are the default and low values. For
example, using a value of 3 took my computer about 2 minutes just to
prepare the surface for showing in the PyMOL GUI (this does not include
any ray-tracing or rendering). Lastly, ray tracing surfaces with high
quality settings takes much longer.
