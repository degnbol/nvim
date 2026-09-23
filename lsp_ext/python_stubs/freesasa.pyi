# patch_stubs: freesasa 2.2.1 f196cd89
"""The :py:mod:`freesasa` python module wraps the FreeSASA `C API`_"""

import _cython_3_0_5
from typing import ClassVar

LeeRichards: str
ShrakeRupley: str
__reduce_cython__: _cython_3_0_5.cython_function_or_method
__setstate_cython__: _cython_3_0_5.cython_function_or_method
__test__: dict
apolar: str
calc: _cython_3_0_5.cython_function_or_method
"""
Calculate SASA of Structure

Args:
    structure: :py:class:`.Structure` to be used
    parameters: :py:class:`.Parameters` to use (if not specified defaults are used)

Returns:
    :py:class:`.Result`: The results

Raises:
    Exception: something went wrong in calculation (see C library error messages)
"""
calcBioPDB: _cython_3_0_5.cython_function_or_method
"""
Calc SASA from `BioPython` PDB structure.

Usage::

    result, sasa_classes, residue_areas = calcBioPDB(structure, ...)

Experimental, not thorougly tested yet

Args:
    bioPDBStructure: A `Bio.PDB` structure
    parameters: A :py:class:`.Parameters` object (uses default if none specified)
    classifier: A :py:class:`.Classifier` object (uses default if none specified)
    options (dict): Options supported are 'hetatm', 'skip-unknown' and 'halt-at-unknown'
        (uses :py:attr:`.Structure.defaultOptions` if none specified

Returns:
    A :py:class:`.Result` object, a dictionary with classes
    defined by the classifier and associated areas,
    and a dictionary of the type returned by :py:meth:`.Result.residueAreas`.

Raises:
    Exception: if unknown atom is encountered and the option
        'halt-at-unknown' is active. Passes on exceptions from
        :py:func:`.calc()`, :py:func:`.classifyResults()` and
        :py:func:`.structureFromBioPDB()`.
"""
calcCoord: _cython_3_0_5.cython_function_or_method
"""
Calculate SASA for a set of coordinates and radii

Args:
    coord (list): array of size 3*N with atomic coordinates
       `(x1, y1, z1,  x2, y2, z2, ..., x_N, y_N, z_N)`.
    radii (list): array of size N with atomic radii `(r_1, r_2, ..., r_N)`.
    parameters: :py:class:`.Parameters` to use (if not specified, defaults are used)
Raises:
    AssertionError: mismatched array-sizes
    Exception: Out of memory
    Exception: something went wrong in calculation (see C library error messages)
"""
classifyResults: _cython_3_0_5.cython_function_or_method
"""
Break SASA result down into classes.

Args:
    result: :py:class:`.Result` from SASA calculation.
    structure: :py:class:`Structure` used in calculation.
    classifier: :py:class:`.Classifier` to use (if not specified default is used).

Returns:
    dict: Dictionary with names of classes as keys and their SASA values as values.

Raises:
    Exception: Problems with classification, see C library error messages
        (or Python exceptions if run with derived classifier).
"""
debug: int
getVerbosity: _cython_3_0_5.cython_function_or_method
"""
Get global verbosity

Returns:
    int: Verbosity :py:const:`.silent`, :py:const:`.nowarnings`
    or :py:const:`.normal`
"""
normal: int
nowarnings: int
polar: str
selectArea: _cython_3_0_5.cython_function_or_method
"""
Sum SASA result over a selection of atoms

Args:
    commands (list): A list of commands with selections using Pymol
        syntax, e.g. ``"s1, resn ala+arg"`` or ``"s2, chain A and resi 1-5"``.
        See `select-syntax`_.
    structure: A :py:class:`.Structure`.
    result: :py:class:`.Result` from sasa calculation on structure.

Returns:
    dict: Dictionary with names of selections (``"s1"``, ``"s2"``, ...) as
    keys, and the corresponding SASA values as values.

Raises:
    Exception: Parser failed (typically syntax error), see
       library error messages.
"""
setVerbosity: _cython_3_0_5.cython_function_or_method
"""
Set global verbosity

Args:
    verbosity (int): Can have values :py:const:`.silent`, :py:const:`.nowarnings`
        or :py:const:`.normal`
Raises:
    AssertionError: if verbosity has illegal value
"""
silent: int
structureArray: _cython_3_0_5.cython_function_or_method
"""
Create array of structures from PDB file.

Split PDB file into several structures by either by treating
chains separately and/or by treating each MODEL as a separate
structure and/or grouping chains.

Args:
    fileName (str): The PDB file.
    options (dict): Specification for how to read the PDB-file
        (see :py:attr:`.Structure.defaultStructureArrayOptions` for
        options and default value).
    classifier: :py:class:`.Classifier` to assign atoms radii, default is used
        if none specified.

Returns:
    list: An array of :py:class:`.Structure`

Raises:
    AssertionError: if `fileName` is None
    AssertionError: if an option value is not recognized
    AssertionError: if neither of the options `'separate-chains'`
        and `'separate-models'` are specified.
    IOError: if can't open file
    Exception: if there are problems parsing the input
"""
structureFromBioPDB: _cython_3_0_5.cython_function_or_method
"""
Create a freesasa structure from a Bio.PDB structure

Experimental, not thorougly tested yet.
Structures generated this way will not preserve whitespace in residue numbers, etc,
as in :py:class:`.Structure`.

Args:
    bioPDBStructure: a `Bio.PDB` structure
    classifier: an optional :py:class:`.Classifier` to specify atomic radii
    options (dict): Options supported are `'hetatm'`, `'skip-unknown'` and `'halt-at-unknown'`

Returns:
    :py:class:`.Structure`: The structure

Raises:
    Exception: if option 'halt-at-unknown' is selected and
        unknown atoms are encountered. Passes on exceptions from
        :py:meth:`.Structure.addAtom()` and
        :py:meth:`.Structure.setRadiiWithClassifier()`.
"""

class Classifier:
    """
    Assigns class and radius to atom by residue and atom name.

    Subclasses derived from :py:class:`.Classifier` can be used to define custom
    atomic radii and/or classes. Can also be initialized from
    config-files_ with a custom classifier.

    If initialized without arguments the default classifier is used.

    Derived classifiers must set the member :py:attr:`.purePython` to ``True``

    Residue names should be of the format ``"ALA"``, ``"ARG"``, etc.
    Atom names should be of the format ``"CA"``, ``"N"``, etc.
    """
    purePython: ClassVar[bool] = ...
    def __init__(self, *args, **kwargs) -> None:
        """
        Constructor.

        If no file is provided the default classifier is used.

        Args:
            fileName (str): Name of file with classifier configuration.

        Raises:
            IOError:   Problem opening/reading file
            Exception: Problem parsing provided configuration or
                       initializing defaults
        """
        ...
    def classify(self, *args, **kwargs):
        """
        Class of atom.

        Depending on the configuration these classes can be
        anything, but typically they will be ``"Polar"`` and ``"Apolar"``.
        Unrecognized atoms will get the class ``"Unknown"``.

        Args:
            residueName (str): Residue name (`"ALA"`, `"ARG"`,...).
            atomName (str): Atom name (`"CA"`, `"C"`,...).

        Returns:
            str: Class name
        """
        ...
    @staticmethod
    def getStandardClassifier(*args, **kwargs):
        """
        Get a standard classifier (ProtOr, OONS or NACCESS)

        Args:
            type (str): The type, can have values ``'protor'``, ``'oons'`` or ``'naccess'``

        Returns:
            :py:class:`.Classifier`: The requested classifier

        Raises:
            Exception: If type not recognized
        """
        ...
    def radius(self, *args, **kwargs):
        """
        Radius of atom.

        This allows the classifier to be used to calculate the atomic
        radii used in calculations. Unknown atoms will get a negative
        radius.

        Args:
            residueName (str): Residue name (`"ALA"`, `"ARG"`, ...).
            atomName (str): Atom name (`"CA"`, `"C"`, ...).

        Returns:
            float: The radius in Å.
        """
        ...
    def __reduce__(self): ...

class Parameters:
    """
    Stores parameter values to be used by calculation.

    Default parameters are ::

        Parameters.defaultParameters = {
            'algorithm'    : LeeRichards,
            'probe-radius' : freesasa_default_parameters.probe_radius,
            'n-points'     : freesasa_default_parameters.shrake_rupley_n_points,
            'n-slices'     : freesasa_default_parameters.lee_richards_n_slices,
            'n-threads'    : freesasa_default_parameters.n_threads
        }

    Attributes:
        defaultParamers (dict): The default parameters
    """
    defaultParameters: ClassVar[dict] = ...
    def __init__(self, *args, **kwargs) -> None:
        """
        Initializes Parameters object.

        Args:
              param (dict): optional argument to specify parameter-values,
                  see :py:attr:`.Parameters.defaultParameters`.

        Raises:
              AssertionError: Invalid parameter values supplied
        """
        ...
    def algorithm(self, *args, **kwargs):
        """
        Get algorithm.

        Returns:
            str: Name of algorithm
        """
        ...
    def nPoints(self, *args, **kwargs):
        """
        Get number of test points in Shrake & Rupley algorithm.

        Returns:
            int: Number of points.
        """
        ...
    def nSlices(self, *args, **kwargs):
        """
        Get the number of slices per atom in Lee & Richards algorithm.

        Returns:
            int: Number of slices.
        """
        ...
    def nThreads(self, *args, **kwargs):
        """
        Get the number of threads to use in calculations.

        Returns:
            int: Number of threads.
        """
        ...
    def probeRadius(self, *args, **kwargs):
        """
        Get probe radius.

        Returns:
             float: Probe radius in Å
        """
        ...
    def setAlgorithm(self, *args, **kwargs):
        """
        Set algorithm.

        Args:
             alg (str): algorithm name, only allowed values are
             :py:data:`freesasa.ShrakeRupley` and :py:data:`freesasa.LeeRichards`

        Raises:
             AssertionError: unknown algorithm specified
        """
        ...
    def setNPoints(self, *args, **kwargs):
        """
        Set number of test points in Shrake & Rupley algorithm.

        Args:
            n (int): Number of points (> 0).

        Raises:
            AssertionError: n <= 0.
        """
        ...
    def setNSlices(self, *args, **kwargs):
        """
        Set the number of slices per atom in Lee & Richards algorithm.

        Args:
            n (int): Number of slices (> 0)

        Raises:
            AssertionError: n <= 0
        """
        ...
    def setNThreads(self, *args, **kwargs):
        """
        Set the number of threads to use in calculations.

        Args:
            n (int): Number of points (> 0)

        Raises:
            AssertionError: n <= 0
        """
        ...
    def setProbeRadius(self, *args, **kwargs):
        """
        Set probe radius.

        Args:
            r (float): probe radius in Å (>= 0)

        Raises:
            AssertionError: r < 0
        """
        ...
    def __reduce__(self): ...

class ResidueArea:
    """
    Stores absolute and relative areas for a residue

    Attributes:
        residueType (str): Type of Residue
        residueNumber (str): Residue number
        hasRelativeAreas (bool): False if there was noe reference area to calculate relative areas from

        total (float): Total SASA of residue
        polar (float): Polar SASA
        apolar (float): Apolar SASA
        mainChain (float): Main chain SASA
        sideChain (float): Side chain SASA

        relativeTotal (float): Relative total SASA
        relativePolar (float): Relative polar SASA
        relativeApolar (float): Relative Apolar SASA
        relativeMainChain (float): Relative main chain SASA
        relativeSideChain (float): Relative side chain SASA
    """
    apolar: ClassVar[int] = ...
    hasRelativeAreas: ClassVar[bool] = ...
    mainChain: ClassVar[int] = ...
    polar: ClassVar[int] = ...
    relativeApolar: ClassVar[int] = ...
    relativeMainChain: ClassVar[int] = ...
    relativePolar: ClassVar[int] = ...
    relativeSideChain: ClassVar[int] = ...
    relativeTotal: ClassVar[int] = ...
    residueNumber: ClassVar[str] = ...
    residueType: ClassVar[str] = ...
    sideChain: ClassVar[int] = ...
    total: ClassVar[int] = ...

class Result:
    """
    Stores results from SASA calculation.

    The type of object returned by :py:func:`freesasa.calc()`,
    not intended to be used outside of that context.
    """
    def __init__(self, *args, **kwargs) -> None: ...
    def atomArea(self, *args, **kwargs):
        """
        SASA for a given atom.

        Args:
            i (int): index of atom.

        Returns:
            float: SASA of atom i in Å^2.

        Raise:
            AssertionError: If no results have been associated
                      with the object or if index is out of bounds
        """
        ...
    def nAtoms(self, *args, **kwargs):
        """
        Number of atoms in the results.

        Returns:
            int: Number of atoms.
        """
        ...
    def residueAreas(self, *args, **kwargs):
        """
        Get SASA for all residues including relative areas if available for the
        classifier used.

        Returns dictionary of results where first dimension is chain label and
        the second dimension residue number. I.e. ``result["A"]["5"]`` gives the
        :py:class:`.ResidueArea` of residue number 5 in chain A.

        Relative areas are normalized to 1, but can be > 1 for
        residues in unusual conformations or at the ends of chains.

        Returns:
            dictionary

        Raise:
            AssertionError: If no results or structure has been associated
                 with the object.
        """
        ...
    def totalArea(self, *args, **kwargs):
        """
        Total SASA.

        Returns:
            The total area in Å^2.
        Raises:
            AssertionError: If no results have been associated with the object.
        """
        ...
    def write_pdb(self, *args, **kwargs): ...
    def __reduce__(self): ...

class Structure:
    """
    Represents a protein structure, including its atomic radii.

    Initialized from PDB-file. Calculates atomic radii using default
    classifier, or custom one provided as argument to initalizer.

    Default options are ::

        Structure.defaultOptions = {
            'hetatm' : False,          # False: skip HETATM
                                       # True: include HETATM

            'hydrogen' : False,        # False: ignore hydrogens
                                       # True: include hydrogens

            'join-models' : False,     # False: Only use the first MODEL
                                       # True: Include all MODELs

            'skip-unknown' : False,    # False: Guess radius for unknown atoms
                                       #     based on element
                                       # True: Skip unknown atoms

            'halt-at-unknown' : False  # False: set radius for unknown atoms,
                                       #    that can not be guessed to 0.
                                       # True: Throw exception on unknown atoms.
        }

    Attributes:
        defaultOptions:  Default options for reading structure from PDB.
    """
    defaultOptions: ClassVar[dict] = ...
    defaultStructureArrayOptions: ClassVar[dict] = ...
    def __init__(self, *args, **kwargs) -> None:
        """
        Constructor

        If a PDB file is provided, the structure will be constructed
        based on the file. If not, this simply initializes an empty
        structure with the given classifier and options. Atoms will then
        have to be added manually using `:py:meth:`.Structure.addAtom()`.

        Args:
            fileName (str): PDB file (if `None` empty structure generated).
            classifier: An optional :py:class:`.Classifier` to calculate atomic
                radii, uses default if none provided.
                This classifier will also be used in calls to :py:meth:`.Structure.addAtom()`
                but only if it's the default classifier, one of the standard
                classifiers from :py:meth:`.Classifier.getStandardClassifier()`,
                or defined by a config-file (i.e. if it uses the underlying
                C API).
            options (dict): specify which atoms and models to include, default is
                :py:attr:`.Structure.defaultOptions`

        Raises:
            IOError: Problem opening/reading file.
            Exception: Problem parsing PDB file or calculating
                atomic radii.
            Exception: If option 'halt-at-unknown' selected and
                unknown atom encountered.
        """
        ...
    def addAtom(self, *args, **kwargs):
        """
        Add atom to structure.

        This function is meant to be used if the structure was not
        initialized from a PDB. The options and classifier passed to
        the constructor for the :py:class:`.Structure` will be used
        (see the documentation of the constructor for restrictions).
        The radii set by the classifier can be overriden by calling
        :py:meth:`.Structure.setRadiiWithClassifier()` afterwards.

        There are no restraints on string lengths for the arguments, but
        the atom won't be added if the classifier doesn't
        recognize the atom and also cannot deduce its element from the
        atom name.

        Args:
            atomName (str): atom name (e.g. `"CA"`)
            residueName (str): residue name (e.g. `"ALA"`)
            residueNumber (str or int): residue number (e.g. `'12'`)
                or integer. Some PDBs have residue-numbers that aren't
                regular numbers. Therefore treated as a string primarily.
            chainLabel (str): 1-character string with chain label (e.g. 'A')
                x,y,z (float): coordinates

        Raises:
            Exception: Residue-number invalid
            AssertionError:
        """
        ...
    def addAtoms(self, *args, **kwargs):
        """
        Add multiple atoms to structure.

        Args:
            atomNames (list): list of atom name (e.g. `["CA"]`)
            residueNames (list): list of residue name (e.g. `["ALA"]`)
            residueNumbers (list): list of residue number (e.g. `['12']`)
                or integer. Some PDBs have residue-numbers that aren't
                regular numbers. Therefore treated as a string primarily.
            chainLabels (list): list of 1-character string with chain label (e.g. ['A'])
                xs,ys,zs (list): list of coordinates

        Raises:
            AssertionError: inconsistent size of input args
        """
        ...
    def atomName(self, *args, **kwargs):
        """
        Get atom name

        Args:
            i (int): Atom index.

        Returns:
            str: Atom name as 4-character string.

        Raises:
            AssertionError: if index out of range or Structure not properly initialized.
        """
        ...
    def chainLabel(self, *args, **kwargs):
        """
        Get chain label for given atom.

        Args:
            i (int): Atom index.

        Returns:
            str: Chain label as 1-character string.

        Raises:
            AssertionError: if index out of range or Structure not properly initialized
        """
        ...
    def coord(self, *args, **kwargs):
        """
        Get coordinates of given atom.

        Args:
            i (int): Atom index.

        Returns:
            list: array of x, y, and z coordinates

        Raises:
            AssertionError: if index out of range or Structure not properly initialized
        """
        ...
    def nAtoms(self, *args, **kwargs):
        """
        Number of atoms.

        Returns:
            int: Number of atoms

        Raises:
            AssertionError: if not properly initialized
        """
        ...
    def radius(self, *args, **kwargs):
        """
        Radius of atom.

        Args:
            i (int): Index of atom.

        Returns:
            float: Radius in Å.

        Raises:
            AssertionError: if index out of bounds, object not properly initalized.
        """
        ...
    def residueName(self, *args, **kwargs):
        """
        Get residue name of given atom.

        Args:
            i (int): Atom index.

        Returns:
            str: Residue name as 3-character string.

        Raises:
            AssertionError: if index out of range or Structure not properly initialized
        """
        ...
    def residueNumber(self, *args, **kwargs):
        """
        Get residue number for given atom.

        Residue number will include the insertion code if there is one.

        Args:
            i (int): Atom index.

        Returns:
            str: Residue number as 5-character string (last character is either whitespace or insertion code)

        Raises:
            AssertionError: if index out of range or Structure not properly initialized
        """
        ...
    def setRadii(self, *args, **kwargs):
        """
        Set atomic radii from an array

        Args:
            radiusArray (list): Array of atomic radii in Ångström, should
                have nAtoms() elements.
        Raises:
            AssertionError: if radiusArray has wrong dimension, structure
                not properly initialized, or if the array contains
                negative radii (not properly classified?)
        """
        ...
    def setRadiiWithClassifier(self, *args, **kwargs):
        """
        Assign radii to atoms in structure using a classifier.

        Args:
            classifier: A :py:class:`.Classifier` to use to calculate radii.

        Raises:
            AssertionError: if structure not properly initialized
        """
        ...
    def setRadius(self, *args, **kwargs):
        """
        Set radius for a given atom

        Args:
            atomIndex (int): Index of atom
            radius (float): Value of radius

        Raises:
            AssertionError: if index out of bounds, radius
                negative, or structure not properly initialized
        """
        ...
    def __reduce__(self): ...
