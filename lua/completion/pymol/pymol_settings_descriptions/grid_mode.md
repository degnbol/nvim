**grid_mode** partitions the screen into a grid and displays each
molecule in one grid location. Each molecule rotates, zooms, etc in that
grid. A possibly very useful option for PyMOL. Each grid area is
assigned a [grid_slot](/index.php/Grid_slot "Grid slot") number. This
allows you to assign objects to certain grids; very helpful. This can be
useful when comparing homologous structures and you want to view and
rotate the aligned structures side by side. This is similar to a
multiple split screen view of the proteins.

**Hint:** When using grid_mode with many molecules, it\'s sometimes good
to align their centers of mass. This puts them all squarely in the
middle of their grid element. The **alignto** command from
[cealign](/index.php/Cealign "Cealign") can do this for you.
