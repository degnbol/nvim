"""
Module containing classes and functions for working with chemical reactions.
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['CartesianProductStrategy', 'ChemicalReaction', 'Compute2DCoordsForReaction', 'CreateDifferenceFingerprintForReaction', 'CreateStructuralFingerprintForReaction', 'EnumerateLibrary', 'EnumerateLibraryBase', 'EnumerateLibraryCanSerialize', 'EnumerationParams', 'EnumerationStrategyBase', 'EvenSamplePairsStrategy', 'FingerprintType', 'GetChemDrawRxnAdjustParams', 'GetDefaultAdjustParams', 'HasAgentTemplateSubstructMatch', 'HasProductTemplateSubstructMatch', 'HasReactantTemplateSubstructMatch', 'HasReactionAtomMapping', 'HasReactionSubstructMatch', 'IsReactionTemplateMoleculeAgent', 'MOL_SPTR_VECT', 'MatchOnlyAtRgroupsAdjustParams', 'MrvBlockIsReaction', 'MrvFileIsReaction', 'PreprocessReaction', 'RandomSampleAllBBsStrategy', 'RandomSampleStrategy', 'ReactionFingerprintParams', 'ReactionFromMolecule', 'ReactionFromMrvBlock', 'ReactionFromMrvFile', 'ReactionFromPNGFile', 'ReactionFromPNGString', 'ReactionFromRxnBlock', 'ReactionFromRxnFile', 'ReactionFromSmarts', 'ReactionFromSmiles', 'ReactionMetadataToPNGFile', 'ReactionMetadataToPNGString', 'ReactionToCXSmarts', 'ReactionToCXSmiles', 'ReactionToMolecule', 'ReactionToMrvBlock', 'ReactionToMrvFile', 'ReactionToRxnBlock', 'ReactionToSmarts', 'ReactionToSmiles', 'ReactionToV3KRxnBlock', 'ReactionsFromCDXMLBlock', 'ReactionsFromCDXMLFile', 'ReduceProductToSideChains', 'RemoveMappingNumbersFromReactions', 'SANITIZE_ADJUST_REACTANTS', 'SANITIZE_ALL', 'SANITIZE_ATOM_MAPS', 'SANITIZE_MERGEHS', 'SANITIZE_NONE', 'SANITIZE_RGROUP_NAMES', 'SanitizeFlags', 'SanitizeRxn', 'SanitizeRxnAsMols', 'UpdateProductsStereochemistry', 'VectMolVect']
class CartesianProductStrategy(EnumerationStrategyBase):
    """
    CartesianProductStrategy produces a standard walk through all possible
    reagent combinations:
    
    (0,0,0), (1,0,0), (2,0,0) ...
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __copy__(self) -> EnumerationStrategyBase:
        """
        """
    def __init__(self) -> None:
        ...
class ChemicalReaction(Boost.Python.instance):
    """
    A class for storing and applying chemical reactions.
    
    Sample Usage:
      >>> from rdkit import Chem
      >>> from rdkit.Chem import rdChemReactions
      >>> rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
      >>> reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('CNC'))
      >>> products = rxn.RunReactants(reacts)
      >>> len(products)
      1
      >>> len(products[0])
      1
      >>> Chem.MolToSmiles(products[0][0])
      'CN(C)C=O'
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 40
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def AddRecursiveQueriesToReaction(reaction: ChemicalReaction, queries: dict = {}, propName: str = 'molFileValue', getLabels: bool = False) -> typing.Any:
        """
            adds recursive queries and returns reactant labels
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddAgentTemplate(self, mol: Mol) -> int:
        """
            adds a agent (a Molecule)
        
        """
    def AddProductTemplate(self, mol: Mol) -> int:
        """
            adds a product (a Molecule)
        
        """
    def AddReactantTemplate(self, mol: Mol) -> int:
        """
            adds a reactant (a Molecule) to the reaction
        
        """
    def ClearComputedProps(self) -> None:
        """
            Removes all computed properties from the reaction.
            
            
        
        """
    def ClearProp(self, key: str) -> None:
        """
            Removes a property from the reaction.
            
              ARGUMENTS:
                - key: the name of the property to clear (a string).
            
        
        """
    def GetAgentTemplate(self, which: int) -> rdkit.Chem.Mol:
        """
            returns one of our agent templates
        
        """
    def GetAgents(self) -> MOL_SPTR_VECT:
        """
            get the agent templates
        
        """
    def GetBoolProp(self, key: str) -> typing.Any:
        """
            Returns the Bool value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a bool
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    def GetDoubleProp(self, key: str) -> typing.Any:
        """
            Returns the double value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    def GetIntProp(self, key: str) -> typing.Any:
        """
            Returns the integer value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    def GetNumAgentTemplates(self) -> int:
        """
            returns the number of agents this reaction expects
        
        """
    def GetNumProductTemplates(self) -> int:
        """
            returns the number of products this reaction generates
        
        """
    def GetNumReactantTemplates(self) -> int:
        """
            returns the number of reactants this reaction expects
        
        """
    def GetProductTemplate(self, which: int) -> rdkit.Chem.Mol:
        """
            returns one of our product templates
        
        """
    def GetProducts(self) -> MOL_SPTR_VECT:
        """
            get the product templates
        
        """
    def GetProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Returns a tuple with all property names for this reaction.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to 0.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to 0.
            
              RETURNS: a tuple of strings
            
        
        """
    def GetPropsAsDict(self, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict:
        """
            Returns a dictionary populated with the reaction's properties.
             n.b. Some properties are not able to be converted to python types.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to False.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to False.
            
              RETURNS: a dictionary
            
        
        """
    def GetReactantTemplate(self, which: int) -> rdkit.Chem.Mol:
        """
            returns one of our reactant templates
        
        """
    def GetReactants(self) -> MOL_SPTR_VECT:
        """
            get the reactant templates
        
        """
    def GetReactingAtoms(self, mappedAtomsOnly: bool = False) -> typing.Any:
        """
            returns a sequence of sequences with the atoms that change in the reaction
        
        """
    def GetSubstructParams(self) -> rdkit.Chem.SubstructMatchParameters:
        """
            get the parameter object controlling the substructure matching
        
        """
    def GetUnsignedProp(self, key: str) -> typing.Any:
        """
            Returns the unsigned int value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an unsigned integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    def HasProp(self, key: str) -> int:
        """
            Queries a molecule to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            
        
        """
    def Initialize(self, silent: bool = False) -> None:
        """
            initializes the reaction so that it can be used
        
        """
    def IsInitialized(self) -> bool:
        """
            checks if the reaction is ready for use
        
        """
    def IsMoleculeAgent(self, mol: Mol) -> bool:
        """
            returns whether or not the molecule has a substructure match to one of the agents.
        
        """
    def IsMoleculeProduct(self, mol: Mol) -> bool:
        """
            returns whether or not the molecule has a substructure match to one of the products.
        
        """
    def IsMoleculeReactant(self, mol: Mol) -> bool:
        """
            returns whether or not the molecule has a substructure match to one of the reactants.
        
        """
    def RemoveAgentTemplates(self, targetList: typing.Any = None) -> None:
        """
            Removes agents from reaction. If targetList is provide the agents will be transferred to that list.
        
        """
    def RemoveUnmappedProductTemplates(self, thresholdUnmappedAtoms: float = 0.2, moveToAgentTemplates: bool = True, targetList: typing.Any = None) -> None:
        """
            Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from product templates to the agent templates or to a given targetList
        
        """
    def RemoveUnmappedReactantTemplates(self, thresholdUnmappedAtoms: float = 0.2, moveToAgentTemplates: bool = True, targetList: typing.Any = None) -> None:
        """
            Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from reactant templates to the agent templates or to a given targetList
        
        """
    def RunReactant(self, reactant: typing.Any, reactionIdx: int) -> typing.Any:
        """
            apply the reaction to a single reactant
        
        """
    def RunReactantInPlace(self, reactant: Mol, removeUnmatchedAtoms: bool = True) -> bool:
        """
            apply the reaction to a single reactant in place. The reactant itself is modified. This can only be used for single reactant - single product reactions.
        
        """
    @typing.overload
    def RunReactants(self, reactants: tuple, maxProducts: int = 1000) -> typing.Any:
        """
            apply the reaction to a sequence of reactant molecules and return the products as a tuple of tuples.  If maxProducts is not zero, stop the reaction when maxProducts have been generated [default=1000]
        
        """
    @typing.overload
    def RunReactants(self, reactants: list, maxProducts: int = 1000) -> typing.Any:
        """
            apply the reaction to a sequence of reactant molecules and return the products as a tuple of tuples.  If maxProducts is not zero, stop the reaction when maxProducts have been generated [default=1000]
        
        """
    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None:
        """
            Sets a boolean valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a bool.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None:
        """
            Sets a double valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a double.
                - computed: (optional) marks the property as being computed.
                            Defaults to 0.
            
            
        
        """
    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None:
        """
            Sets an integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (an unsigned number).
                - value: the property value as an integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    def SetProp(self, key: str, val: str, computed: bool = False) -> None:
        """
            Sets a molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None:
        """
            Sets an unsigned integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as an unsigned integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    @typing.overload
    def ToBinary(self) -> typing.Any:
        """
            Returns a binary string representation of the reaction.
        
        """
    @typing.overload
    def ToBinary(self, propertyFlags: int) -> typing.Any:
        """
            Returns a binary string representation of the reaction.
        
        """
    def Validate(self, silent: bool = False) -> tuple:
        """
            checks the reaction for potential problems, returns (numWarnings,numErrors)
        
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, binStr: str) -> None:
        ...
    @typing.overload
    def __init__(self, other: ChemicalReaction) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def _getImplicitPropertiesFlag(self) -> bool:
        """
            EXPERT USER: returns whether or not the reaction can have implicit properties
        
        """
    def _setImplicitPropertiesFlag(self, val: bool) -> None:
        """
            EXPERT USER: indicates that the reaction can have implicit properties
        
        """
class EnumerateLibrary(EnumerateLibraryBase):
    """
    EnumerateLibrary
    This class allows easy enumeration of reactions.  Simply provide a reaction
    and a set of reagents and you are off the races.
    
    Note that this functionality should be considered beta and that the API may
    change in a future release.
    
    EnumerateLibrary follows the python enumerator protocol, for example:
    
    library = EnumerateLibrary(rxn, bbs)
    for products in library:
       ... do something with the product
    
    It is useful to sanitize reactions before hand:
    
    SanitizeRxn(rxn)
    library = EnumerateLibrary(rxn, bbs)
    
    If ChemDraw style reaction semantics are prefereed, you can apply
    the ChemDraw parameters:
    
    SanitizeRxn(rxn, params=GetChemDrawRxnAdjustParams())
    
    For one, this enforces only matching RGroups and assumes all atoms
    have fully satisfied valences.
    
    Each product has the same output as applying a set of reagents to
    the libraries reaction.
    
    This can be a bit confusing as each product can have multiple molecules
    generated.  The returned data structure is as follows:
    
       [ [products1], [products2],... ]
    Where products1 are the molecule products for the reactions first product
    template and products2 are the molecule products for the second product
    template.  Since each reactant can match more than once, there may be
    multiple product molecules for each template.
    
    for products in library:
        for results_for_product_template in products:
            for mol in results_for_product_template:
                Chem.MolToSmiles(mol) # finally have a molecule!
    
    For sufficiently large libraries, using this iteration strategy is not
    recommended as the library may contain more products than atoms in the
    universe.  To help with this, you can supply an enumeration strategy.
    The default strategy is a CartesianProductStrategy which enumerates
    everything.  RandomSampleStrategy randomly samples the products but
    this strategy never terminates, however, python supplies itertools:
    
    import itertools
    library = EnumerateLibrary(rxn, bbs, rdChemReactions.RandomSampleStrategy())
    for result in itertools.islice(library, 1000):
        # do something with the first 1000 samples
    
    for result in itertools.islice(library, 1000):
        # do something with the next 1000 samples
    
    Libraries are also serializable, including their current state:
    
    s = library.Serialize()
    library2 = EnumerateLibrary()
    library2.InitFromString(s)
    for result in itertools.islice(libary2, 1000):
        # do something with the next 1000 samples
    """
    __instance_size__: typing.ClassVar[int] = 384
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetReagents(self) -> VectMolVect:
        """
            Return the reagents used in this library.  These are the subset of the input reagents that are compatible with the reaction so may be smaller than the input reagent sets.
        
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, rxn: ChemicalReaction, reagents: list, params: EnumerationParams) -> None:
        ...
    @typing.overload
    def __init__(self, rxn: ChemicalReaction, reagents: tuple, params: EnumerationParams) -> None:
        ...
    @typing.overload
    def __init__(self, rxn: ChemicalReaction, reagents: list, enumerator: EnumerationStrategyBase, params: EnumerationParams) -> None:
        ...
    @typing.overload
    def __init__(self, rxn: ChemicalReaction, reagents: tuple, enumerator: EnumerationStrategyBase, params: EnumerationParams) -> None:
        ...
class EnumerateLibraryBase(Boost.Python.instance):
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetEnumerator(self) -> EnumerationStrategyBase:
        """
            Returns the enumation strategy for the current library
        
        """
    def GetPosition(self) -> VectSizeT:
        """
            Returns the current enumeration position into the reagent vectors, as returned by GetReagents().  They do not necessarily refer to the input reagent sets as they only refer to reagents compatible with the reaction.
        
        """
    def GetReaction(self) -> ChemicalReaction:
        """
            Returns the chemical reaction for this library
        
        """
    def GetState(self) -> str:
        """
            Returns the current enumeration state (position) of the library.
            This position can be used to restart the library from a known position
        
        """
    def InitFromString(self, data: str) -> None:
        """
            Inititialize the library from a binary string
        
        """
    def ResetState(self) -> None:
        """
            Returns the current enumeration state (position) of the library to the start.
        
        """
    def Serialize(self) -> typing.Any:
        """
            Serialize the library to a binary string.
            Note that the position in the library is serialized as well.  Care should
            be taken when serializing.  See GetState/SetState for position manipulation.
        
        """
    def SetState(self, state: str) -> None:
        """
            Sets the enumeration state (position) of the library.
        
        """
    def __bool__(self) -> bool:
        """
        """
    def __iter__(self) -> typing.Any:
        """
        """
    def __next__(self) -> typing.Any:
        """
            Return the next molecule from the enumeration.
        
        """
    def __nonzero__(self) -> bool:
        """
        """
    def next(self) -> typing.Any:
        """
            Return the next molecule from the enumeration.
        
        """
    def nextSmiles(self) -> VectorOfStringVectors:
        """
            Return the next smiles string from the enumeration.
        
        """
class EnumerationParams(Boost.Python.instance):
    """
    EnumerationParams
    Controls some aspects of how the enumeration is performed.
    Options:
      reagentMaxMatchCount [ default Infinite ]
        This specifies how many times the reactant template can match a reagent.
    
      sanePartialProducts [default false]
        If true, forces all products of the reagent plus the product templates
         pass chemical sanitization.  Note that if the product template itself
         does not pass sanitization, then none of the products will.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def reagentMaxMatchCount(self) -> int:
        """default: 2147483647"""
    @reagentMaxMatchCount.setter
    def reagentMaxMatchCount(self, value: int) -> None: ...
    @property
    def sanePartialProducts(self) -> bool:
        """default: False"""
    @sanePartialProducts.setter
    def sanePartialProducts(self, value: bool) -> None: ...
class EnumerationStrategyBase(Boost.Python.instance):
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetNumPermutations(self) -> int:
        """
            Returns the total number of results for this enumeration strategy.
            Note that some strategies are effectively infinite.
        
        """
    def GetPosition(self) -> VectSizeT:
        """
            Return the current indices into the arrays of reagents, as returned by GetReagents().  They do not necessarily refer to the input reagent sets as they only refer to reagents compatible with the reaction.
        
        """
    def Initialize(self, rxn: ChemicalReaction, ob: list) -> None:
        """
        """
    def Skip(self, skipCount: int) -> bool:
        """
            Skip the next Nth results. note: this may be an expensive operation
            depending on the enumeration strategy used. It is recommended to use
            the enumerator state to advance to a known position
        
        """
    def Type(self) -> str:
        """
            Returns the enumeration strategy type as a string.
        
        """
    def __bool__(self) -> bool:
        """
        """
    @typing.overload
    def __copy__(self) -> EnumerationStrategyBase:
        """
        """
    @typing.overload
    def __copy__(self) -> None:
        """
        """
    @typing.overload
    def __next__(self) -> VectSizeT:
        """
            Return the next indices into the arrays of reagents
        
        """
    @typing.overload
    def __next__(self) -> None:
        """
        """
    def __nonzero__(self) -> bool:
        """
        """
    @typing.overload
    def next(self) -> VectSizeT:
        """
            Return the next indices into the arrays of reagents
        
        """
    @typing.overload
    def next(self) -> None:
        """
        """
class EvenSamplePairsStrategy(EnumerationStrategyBase):
    """
    Randomly sample Pairs evenly from a collection of building blocks
    This is a good strategy for choosing a relatively small selection
    of building blocks from a larger set.  As the amount of work needed
    to retrieve the next evenly sample building block grows with the
    number of samples, this method performs progressively worse as the
    number of samples gets larger.
    See EnumerationStrategyBase for more details.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def Stats(self) -> str:
        """
            Return the statistics log of the pairs used in the current enumeration.
        
        """
    def __copy__(self) -> EnumerationStrategyBase:
        """
        """
    def __init__(self) -> None:
        ...
class FingerprintType(Boost.Python.enum):
    AtomPairFP: typing.ClassVar[FingerprintType]  # value = rdkit.Chem.rdChemReactions.FingerprintType.AtomPairFP
    MorganFP: typing.ClassVar[FingerprintType]  # value = rdkit.Chem.rdChemReactions.FingerprintType.MorganFP
    PatternFP: typing.ClassVar[FingerprintType]  # value = rdkit.Chem.rdChemReactions.FingerprintType.PatternFP
    RDKitFP: typing.ClassVar[FingerprintType]  # value = rdkit.Chem.rdChemReactions.FingerprintType.RDKitFP
    TopologicalTorsion: typing.ClassVar[FingerprintType]  # value = rdkit.Chem.rdChemReactions.FingerprintType.TopologicalTorsion
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'AtomPairFP': rdkit.Chem.rdChemReactions.FingerprintType.AtomPairFP, 'TopologicalTorsion': rdkit.Chem.rdChemReactions.FingerprintType.TopologicalTorsion, 'MorganFP': rdkit.Chem.rdChemReactions.FingerprintType.MorganFP, 'RDKitFP': rdkit.Chem.rdChemReactions.FingerprintType.RDKitFP, 'PatternFP': rdkit.Chem.rdChemReactions.FingerprintType.PatternFP}
    values: typing.ClassVar[dict]  # value = {1: rdkit.Chem.rdChemReactions.FingerprintType.AtomPairFP, 2: rdkit.Chem.rdChemReactions.FingerprintType.TopologicalTorsion, 3: rdkit.Chem.rdChemReactions.FingerprintType.MorganFP, 4: rdkit.Chem.rdChemReactions.FingerprintType.RDKitFP, 5: rdkit.Chem.rdChemReactions.FingerprintType.PatternFP}
class MOL_SPTR_VECT(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class RandomSampleAllBBsStrategy(EnumerationStrategyBase):
    """
    RandomSampleAllBBsStrategy randomly samples from the reagent sets
    with the constraint that all building blocks are samples as early as possible.
    Note that this strategy never halts and can produce duplicates.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __copy__(self) -> EnumerationStrategyBase:
        """
        """
    def __init__(self) -> None:
        ...
class RandomSampleStrategy(EnumerationStrategyBase):
    """
    RandomSampleStrategy simply randomly samples from the reagent sets.
    Note that this strategy never halts and can produce duplicates.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __copy__(self) -> EnumerationStrategyBase:
        """
        """
    def __init__(self) -> None:
        ...
class ReactionFingerprintParams(Boost.Python.instance):
    """
    A class for storing parameters to manipulate the calculation of fingerprints of chemical reactions.
    """
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, includeAgents: bool, bitRatioAgents: float, nonAgentWeight: int, agentWeight: int, fpSize: int, fpType: FingerprintType) -> None:
        ...
    @property
    def agentWeight(self) -> int:
        """default: 1"""
    @agentWeight.setter
    def agentWeight(self, value: int) -> None: ...
    @property
    def bitRatioAgents(self) -> float:
        """default: 0.2"""
    @bitRatioAgents.setter
    def bitRatioAgents(self, value: float) -> None: ...
    @property
    def fpSize(self) -> int:
        """default: 2048"""
    @fpSize.setter
    def fpSize(self, value: int) -> None: ...
    @property
    def fpType(*args, **kwargs):
        ...
    @fpType.setter
    def fpType(*args, **kwargs):
        ...
    @property
    def includeAgents(self) -> bool:
        """default: False"""
    @includeAgents.setter
    def includeAgents(self, value: bool) -> None: ...
    @property
    def nonAgentWeight(self) -> int:
        """default: 10"""
    @nonAgentWeight.setter
    def nonAgentWeight(self, value: int) -> None: ...
class SanitizeFlags(Boost.Python.enum):
    SANITIZE_ADJUST_REACTANTS: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS
    SANITIZE_ALL: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL
    SANITIZE_ATOM_MAPS: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS
    SANITIZE_MERGEHS: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS
    SANITIZE_NONE: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE
    SANITIZE_RGROUP_NAMES: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'SANITIZE_NONE': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE, 'SANITIZE_ATOM_MAPS': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS, 'SANITIZE_RGROUP_NAMES': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES, 'SANITIZE_ADJUST_REACTANTS': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS, 'SANITIZE_MERGEHS': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS, 'SANITIZE_ALL': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE, 2: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS, 1: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES, 4: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS, 8: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS, 4294967295: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL}
class VectMolVect(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
def Compute2DCoordsForReaction(reaction: ChemicalReaction, spacing: float = 1.0, updateProps: bool = True, canonOrient: bool = True, nFlipsPerSample: int = 0, nSample: int = 0, sampleSeed: int = 0, permuteDeg4Nodes: bool = False, bondLength: float = -1.0) -> None:
    """
        Compute 2D coordinates for a reaction. 
          ARGUMENTS: 
             - reaction - the reaction of interest
             - spacing - the amount of space left between components of the reaction
             - canonOrient - orient the reactants and products in a canonical way
             - updateProps - if set, properties such as conjugation and
                hybridization will be calculated for the reactant and product
                templates before generating coordinates. This should result in
                better depictions, but can lead to errors in some cases.
             - nFlipsPerSample - number of rotatable bonds that are
                        flipped at random at a time.
             - nSample - Number of random samplings of rotatable bonds.
             - sampleSeed - seed for the random sampling process.
             - permuteDeg4Nodes - allow permutation of bonds at a degree 4
                         node during the sampling process 
             - bondLength - change the default bond length for depiction
        
    
    """
def CreateDifferenceFingerprintForReaction(reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ...) -> UIntSparseIntVect:
    """
        construct a difference fingerprint for a ChemicalReaction by subtracting the reactant fingerprint from the product fingerprint
    
    """
def CreateStructuralFingerprintForReaction(reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ...) -> ExplicitBitVect:
    """
        construct a structural fingerprint for a ChemicalReaction by concatenating the reactant fingerprint and the product fingerprint
    
    """
def EnumerateLibraryCanSerialize() -> bool:
    """
        Returns True if the EnumerateLibrary is serializable (requires boost serialization
    
    """
def GetChemDrawRxnAdjustParams() -> rdkit.Chem.AdjustQueryParameters:
    """
        (deprecated, see MatchOnlyAtRgroupsAdjustParams)
        	Returns the chemdraw style adjustment parameters for reactant templates
    
    """
def GetDefaultAdjustParams() -> rdkit.Chem.AdjustQueryParameters:
    """
        Returns the default adjustment parameters for reactant templates
    
    """
def HasAgentTemplateSubstructMatch(reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
        tests if the agents of a queryReaction are the same as those of a reaction
    
    """
def HasProductTemplateSubstructMatch(reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
        tests if the products of a queryReaction are substructures of the products of a reaction
    
    """
def HasReactantTemplateSubstructMatch(reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
        tests if the reactants of a queryReaction are substructures of the reactants of a reaction
    
    """
def HasReactionAtomMapping(rxn: ChemicalReaction) -> bool:
    """
        tests if a reaction obtains any atom mapping
    
    """
def HasReactionSubstructMatch(reaction: ChemicalReaction, queryReaction: ChemicalReaction, includeAgents: bool = False) -> bool:
    """
        tests if the queryReaction is a substructure of a reaction
    
    """
def IsReactionTemplateMoleculeAgent(molecule: Mol, agentThreshold: float) -> bool:
    """
        tests if a molecule can be classified as an agent depending on the ratio of mapped atoms and a give threshold
    
    """
def MatchOnlyAtRgroupsAdjustParams() -> rdkit.Chem.AdjustQueryParameters:
    """
        Only match at the specified rgroup locations in the reactant templates
    
    """
def MrvBlockIsReaction(mrvData: str) -> bool:
    """
        returns whether or not an MRV block contains reaction data
    
    """
def MrvFileIsReaction(filename: str) -> bool:
    """
        returns whether or not an MRV file contains reaction data
    
    """
def PreprocessReaction(reaction: ChemicalReaction, queries: dict = {}, propName: str = 'molFileValue') -> typing.Any:
    """
        A function for preprocessing reactions with more specific queries.
        Queries are indicated by labels on atoms (molFileAlias property by default)
        When these labels are found, more specific queries are placed on the atoms.
        By default, the available quieries come from 
          FilterCatalog.GetFlattenedFunctionalGroupHierarchy(True)n
        Sample Usage:
          >>> from rdkit import Chem, RDConfig
          >>> from rdkit.Chem import MolFromSmiles, AllChem
          >>> from rdkit.Chem.rdChemReactions import PreprocessReaction
          >>> import os
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','boronic1.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          >>> nWarn
          0
          >>> nError
          0
          >>> nReacts
          2
          >>> nProds
          1
          >>> reactantLabels
          (((0, 'halogen.bromine.aromatic'),), ((1, 'boronicacid'),))
        
        If there are functional group labels in the input reaction (via atoms with molFileValue properties),
        the corresponding atoms will have queries added to them so that they only match such things. We can
        see this here:
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> r1 = rxn.GetReactantTemplate(0)
          >>> m1 = Chem.MolFromSmiles('CCBr')
          >>> m2 = Chem.MolFromSmiles('c1ccccc1Br')
          
        These both match because the reaction file itself just has R1-Br:
          >>> m1.HasSubstructMatch(r1)
          True
          >>> m2.HasSubstructMatch(r1)
          True
        
        After preprocessing, we only match the aromatic Br:
          >>> d = PreprocessReaction(rxn)
          >>> m1.HasSubstructMatch(r1)
          False
          >>> m2.HasSubstructMatch(r1)
          True
        
        We also support or queries in the values field (separated by commas):
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','azide_reaction.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> reactantLabels = PreprocessReaction(rxn)[-1]
          >>> reactantLabels
          (((1, 'azide'),), ((1, 'carboxylicacid,acidchloride'),))
          >>> m1 = Chem.MolFromSmiles('CC(=O)O')
          >>> m2 = Chem.MolFromSmiles('CC(=O)Cl')
          >>> m3 = Chem.MolFromSmiles('CC(=O)N')
          >>> r2 = rxn.GetReactantTemplate(1)
          >>> m1.HasSubstructMatch(r2)
          True
          >>> m2.HasSubstructMatch(r2)
          True
          >>> m3.HasSubstructMatch(r2)
          False
        
        unrecognized final group types are returned as None:
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value1.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          Traceback (most recent call last):
            ...
          KeyError: 'boromicacid'
        
        One unrecognized group type in a comma-separated list makes the whole thing fail:
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value2.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          Traceback (most recent call last):
            ...
          KeyError: 'carboxylicacid,acidchlroide'
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value3.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          Traceback (most recent call last):
            ...
          KeyError: 'carboxyliccaid,acidchloride'
          >>> rxn = rdChemReactions.ChemicalReaction()
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          >>> reactantLabels
          ()
          >>> reactantLabels == ()
          True
        
    
    """
def ReactionFromMolecule(mol: Mol) -> ChemicalReaction:
    """
        construct a ChemicalReaction from an molecule if the RXN role property of the molecule is set
    
    """
@typing.overload
def ReactionFromMrvBlock(rxnblock: typing.Any, sanitize: bool = False, removeHs: bool = False) -> ChemicalReaction:
    """
        construct a ChemicalReaction from a string in Marvin (mrv) format
    
    """
@typing.overload
def ReactionFromMrvBlock(rxnblock: typing.Any, sanitize: bool = False, removeHs: bool = False) -> ChemicalReaction:
    """
        construct a ChemicalReaction from a string in Marvin (mrv) format
    
    """
def ReactionFromMrvFile(filename: str, sanitize: bool = False, removeHs: bool = False) -> ChemicalReaction:
    """
        construct a ChemicalReaction from an Marvin (mrv) rxn file
    
    """
def ReactionFromPNGFile(fname: str) -> ChemicalReaction:
    """
        construct a ChemicalReaction from metadata in a PNG file
    
    """
def ReactionFromPNGString(data: str) -> ChemicalReaction:
    """
        construct a ChemicalReaction from an string with PNG data
    
    """
def ReactionFromRxnBlock(rxnblock: str, sanitize: bool = False, removeHs: bool = False, strictParsing: bool = True) -> ChemicalReaction:
    """
        construct a ChemicalReaction from a string in MDL rxn format
    
    """
def ReactionFromRxnFile(filename: str, sanitize: bool = False, removeHs: bool = False, strictParsing: bool = True) -> ChemicalReaction:
    """
        construct a ChemicalReaction from an MDL rxn file
    
    """
def ReactionFromSmarts(SMARTS: str, replacements: dict = {}, useSmiles: bool = False) -> ChemicalReaction:
    """
        construct a ChemicalReaction from a reaction SMARTS string. 
        see the documentation for rdkit.Chem.MolFromSmiles for an explanation
        of the replacements argument.
    
    """
def ReactionFromSmiles(SMILES: str, replacements: dict = {}) -> ChemicalReaction:
    """
        construct a ChemicalReaction from a reaction SMILES string. 
        see the documentation for rdkit.Chem.MolFromSmiles for an explanation
        of the replacements argument.
    
    """
def ReactionMetadataToPNGFile(mol: ChemicalReaction, filename: typing.Any, includePkl: bool = True, includeSmiles: bool = True, includeSmarts: bool = False, includeMol: bool = False) -> typing.Any:
    """
        Reads the contents of a PNG file and adds metadata about a reaction to it. The modified file contents are returned.
    
    """
def ReactionMetadataToPNGString(mol: ChemicalReaction, pngdata: typing.Any, includePkl: bool = True, includeSmiles: bool = True, includeSmarts: bool = False, includeRxn: bool = False) -> typing.Any:
    """
        Adds metadata about a reaction to the PNG string passed in.The modified string is returned.
    
    """
@typing.overload
def ReactionToCXSmarts(reaction: ChemicalReaction) -> str:
    """
        construct a reaction SMARTS string for a ChemicalReaction
    
    """
@typing.overload
def ReactionToCXSmarts(reaction: ChemicalReaction, params: SmilesWriteParams, flags: int = ...) -> str:
    """
        construct a reaction CXSMARTS string for a ChemicalReaction
    
    """
@typing.overload
def ReactionToCXSmiles(reaction: ChemicalReaction, canonical: bool = True) -> str:
    """
        construct a reaction SMILES string for a ChemicalReaction
    
    """
@typing.overload
def ReactionToCXSmiles(reaction: ChemicalReaction, params: SmilesWriteParams, flags: int = ...) -> str:
    """
        construct a reaction CXSMILES string for a ChemicalReaction
    
    """
def ReactionToMolecule(reaction: ChemicalReaction) -> rdkit.Chem.Mol:
    """
        construct a molecule for a ChemicalReaction with RXN role property set
    
    """
def ReactionToMrvBlock(reaction: ChemicalReaction, prettyPrint: bool = False) -> str:
    """
        construct a string in Marvin (MRV) rxn format for a ChemicalReaction
    
    """
def ReactionToMrvFile(reaction: ChemicalReaction, filename: str, prettyPrint: bool = False) -> None:
    """
        write a Marvin (MRV) rxn file for a ChemicalReaction
    
    """
def ReactionToRxnBlock(reaction: ChemicalReaction, separateAgents: bool = False, forceV3000: bool = False) -> str:
    """
        construct a string in MDL rxn format for a ChemicalReaction
    
    """
@typing.overload
def ReactionToSmarts(reaction: ChemicalReaction) -> str:
    """
        construct a reaction SMARTS string for a ChemicalReaction
    
    """
@typing.overload
def ReactionToSmarts(reaction: ChemicalReaction, params: SmilesWriteParams) -> str:
    """
        construct a reaction SMARTS string for a ChemicalReaction
    
    """
@typing.overload
def ReactionToSmiles(reaction: ChemicalReaction, canonical: bool = True) -> str:
    """
        construct a reaction SMILES string for a ChemicalReaction
    
    """
@typing.overload
def ReactionToSmiles(reaction: ChemicalReaction, params: SmilesWriteParams) -> str:
    """
        construct a reaction SMILES string for a ChemicalReaction
    
    """
def ReactionToV3KRxnBlock(reaction: ChemicalReaction, separateAgents: bool = False) -> str:
    """
        construct a string in MDL v3000 rxn format for a ChemicalReaction
    
    """
def ReactionsFromCDXMLBlock(rxnblock: typing.Any, sanitize: bool = False, removeHs: bool = False) -> typing.Any:
    """
        construct a tuple of ChemicalReactions from a string in CDXML format
    
    """
def ReactionsFromCDXMLFile(filename: str, sanitize: bool = False, removeHs: bool = False) -> typing.Any:
    """
        construct a tuple of ChemicalReactions from a CDXML rxn file
    
    """
def ReduceProductToSideChains(product: Mol, addDummyAtoms: bool = True) -> rdkit.Chem.Mol:
    """
        reduce the product of a reaction to the side chains added by the reaction.              The output is a molecule with attached wildcards indicating where the product was attached.              The dummy atom has the same reaction-map number as the product atom (if available).
    
    """
def RemoveMappingNumbersFromReactions(reaction: ChemicalReaction) -> None:
    """
        Removes the mapping numbers from the molecules of a reaction
    
    """
def SanitizeRxn(rxn: ChemicalReaction, sanitizeOps: int = 4294967295, params: AdjustQueryParameters = ..., catchErrors: bool = False) -> SanitizeFlags:
    """
        Does some sanitization of the reactant and product templates of a reaction.
        
            - The reaction is modified in place.
            - If sanitization fails, an exception will be thrown unless catchErrors is set
        
          ARGUMENTS:
        
            - rxn: the reaction to be modified
            - sanitizeOps: (optional) reaction sanitization operations to be carried out
              these should be constructed by or'ing together the
              operations in rdkit.Chem.rdChemReactions.SanitizeFlags
            - optional adjustment parameters for changing the meaning of the substructure
              matching done in the templates.  The default is 
              rdkit.Chem.rdChemReactions.DefaultRxnAdjustParams which aromatizes
              kekule structures if possible.
            - catchErrors: (optional) if provided, instead of raising an exception
              when sanitization fails (the default behavior), the 
              first operation that failed (as defined in rdkit.Chem.rdChemReactions.SanitizeFlags)
              is returned. Zero is returned on success.
        
          The operations carried out by default are:
            1) fixRGroups(): sets R group labels on mapped dummy atoms when possible
            2) fixAtomMaps(): attempts to set atom maps on unmapped R groups
            3) adjustTemplate(): calls adjustQueryProperties() on all reactant templates
            4) fixHs(): merges explicit Hs in the reactant templates that don't map to heavy atoms
        
    
    """
def SanitizeRxnAsMols(rxn: ChemicalReaction, sanitizeOps: int = 268435455) -> None:
    """
        Does the usual molecular sanitization on each reactant, agent, and product of the reaction
    
    """
def UpdateProductsStereochemistry(reaction: ChemicalReaction) -> None:
    """
        Caution: This is an expert-user function which will change a property (molInversionFlag) of your products.          This function is called by default using the RXN or SMARTS parser for reactions and should really only be called if reactions have been constructed some other way.          The function updates the stereochemistry of the product by considering 4 different cases: inversion, retention, removal, and introduction
    
    """
SANITIZE_ADJUST_REACTANTS: SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS
SANITIZE_ALL: SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL
SANITIZE_ATOM_MAPS: SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS
SANITIZE_MERGEHS: SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS
SANITIZE_NONE: SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE
SANITIZE_RGROUP_NAMES: SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES
