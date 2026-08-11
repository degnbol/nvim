"""
Module containing a C++ implementation of 2D molecule drawing
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['ALL', 'ANNOTATIONS', 'ATOMLABELS', 'BONDS', 'Bottom', 'CircleAndLine', 'ContourAndDrawGaussians', 'ContourAndDrawGrid', 'ContourParams', 'DrawElement', 'DrawMoleculeACS1996', 'HIGHLIGHTS', 'IntStringMap', 'Lasso', 'Left', 'LegendPosition', 'MeanBondLength', 'MolDraw2D', 'MolDraw2DCairo', 'MolDraw2DSVG', 'MolDrawOptions', 'MolToACS1996SVG', 'MolToSVG', 'MultiColourHighlightStyle', 'NONE', 'POSTSHAPES', 'PRESHAPES', 'PrepareAndDrawMolecule', 'PrepareMolForDrawing', 'RADICALS', 'Right', 'SetACS1996Mode', 'SetDarkMode', 'SetMonochromeMode', 'Top', 'UpdateDrawerParamsFromJSON', 'UpdateMolDrawOptionsFromJSON', 'map_indexing_suite_IntStringMap_entry']
class ContourParams(Boost.Python.instance):
    """
    Parameters for drawing contours
    """
    __instance_size__: typing.ClassVar[int] = 160
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    def setColourMap(self, colours: typing.Any) -> None:
        """
        """
    def setContourColour(self, colour: tuple) -> None:
        """
        """
    @property
    def colourMap(*args, **kwargs):
        """
        the color map to use when filling the grid
        """
    @colourMap.setter
    def colourMap(*args, **kwargs):
        ...
    @property
    def contourColour(*args, **kwargs):
        """
        the color to use for drawing the contours
        """
    @contourColour.setter
    def contourColour(*args, **kwargs):
        ...
    @property
    def contourWidth(self) -> float:
        """line width of the contours (default: 1.0)"""
    @contourWidth.setter
    def contourWidth(self, value: float) -> None: ...
    @property
    def coordScaleForQuantization(self) -> float:
        """scaling factor used to convert coordinates to ints when forming the continuous lines (default: 1000.0)"""
    @coordScaleForQuantization.setter
    def coordScaleForQuantization(self, value: float) -> None: ...
    @property
    def dashNegative(self) -> bool:
        """use a dashed line for negative contours (default: True)"""
    @dashNegative.setter
    def dashNegative(self, value: bool) -> None: ...
    @property
    def drawAsLines(self) -> bool:
        """draw the contours as continuous lines isntead of line segments (default: True)"""
    @drawAsLines.setter
    def drawAsLines(self, value: bool) -> None: ...
    @property
    def extraGridPadding(self) -> float:
        """extra space (in molecule coords) around the grid (default: 0.0)"""
    @extraGridPadding.setter
    def extraGridPadding(self, value: float) -> None: ...
    @property
    def fillGrid(self) -> bool:
        """colors the grid in addition to drawing contours (default: False)"""
    @fillGrid.setter
    def fillGrid(self, value: bool) -> None: ...
    @property
    def fillThreshold(self) -> float:
        """magnitude threshold to determine if a grid point is filled (default: 0.01)"""
    @fillThreshold.setter
    def fillThreshold(self, value: float) -> None: ...
    @property
    def fillThresholdIsFraction(self) -> bool:
        """if true, fillThreshold is a fraction of the range of the data (default: True)"""
    @fillThresholdIsFraction.setter
    def fillThresholdIsFraction(self, value: bool) -> None: ...
    @property
    def gridResolution(self) -> float:
        """set the resolution of the grid (default: 0.15)"""
    @gridResolution.setter
    def gridResolution(self, value: float) -> None: ...
    @property
    def isovalScaleForQuantization(self) -> float:
        """scaling factor used to convert isovalues to ints when forming the continuous lines (default: 1000000.0)"""
    @isovalScaleForQuantization.setter
    def isovalScaleForQuantization(self, value: float) -> None: ...
    @property
    def setScale(self) -> bool:
        """set the scale of the drawing object (useful if you draw the grid/contours first) (default: True)"""
    @setScale.setter
    def setScale(self, value: bool) -> None: ...
    @property
    def useFillThreshold(self) -> bool:
        """use a magnitude threshold to determine if a grid point is filled (default: False)"""
    @useFillThreshold.setter
    def useFillThreshold(self, value: bool) -> None: ...
class DrawElement(Boost.Python.enum):
    ALL: typing.ClassVar[DrawElement]  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ALL
    ANNOTATIONS: typing.ClassVar[DrawElement]  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ANNOTATIONS
    ATOMLABELS: typing.ClassVar[DrawElement]  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ATOMLABELS
    BONDS: typing.ClassVar[DrawElement]  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.BONDS
    HIGHLIGHTS: typing.ClassVar[DrawElement]  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.HIGHLIGHTS
    NONE: typing.ClassVar[DrawElement]  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.NONE
    POSTSHAPES: typing.ClassVar[DrawElement]  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.POSTSHAPES
    PRESHAPES: typing.ClassVar[DrawElement]  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.PRESHAPES
    RADICALS: typing.ClassVar[DrawElement]  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.RADICALS
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'NONE': rdkit.Chem.Draw.rdMolDraw2D.DrawElement.NONE, 'PRESHAPES': rdkit.Chem.Draw.rdMolDraw2D.DrawElement.PRESHAPES, 'BONDS': rdkit.Chem.Draw.rdMolDraw2D.DrawElement.BONDS, 'ATOMLABELS': rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ATOMLABELS, 'HIGHLIGHTS': rdkit.Chem.Draw.rdMolDraw2D.DrawElement.HIGHLIGHTS, 'ANNOTATIONS': rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ANNOTATIONS, 'RADICALS': rdkit.Chem.Draw.rdMolDraw2D.DrawElement.RADICALS, 'POSTSHAPES': rdkit.Chem.Draw.rdMolDraw2D.DrawElement.POSTSHAPES, 'ALL': rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ALL}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.Draw.rdMolDraw2D.DrawElement.NONE, 1: rdkit.Chem.Draw.rdMolDraw2D.DrawElement.PRESHAPES, 2: rdkit.Chem.Draw.rdMolDraw2D.DrawElement.BONDS, 4: rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ATOMLABELS, 8: rdkit.Chem.Draw.rdMolDraw2D.DrawElement.HIGHLIGHTS, 16: rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ANNOTATIONS, 32: rdkit.Chem.Draw.rdMolDraw2D.DrawElement.RADICALS, 64: rdkit.Chem.Draw.rdMolDraw2D.DrawElement.POSTSHAPES, 2147483647: rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ALL}
class IntStringMap(Boost.Python.instance):
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
class LegendPosition(Boost.Python.enum):
    Bottom: typing.ClassVar[LegendPosition]  # value = rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Bottom
    Left: typing.ClassVar[LegendPosition]  # value = rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Left
    Right: typing.ClassVar[LegendPosition]  # value = rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Right
    Top: typing.ClassVar[LegendPosition]  # value = rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Top
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'Bottom': rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Bottom, 'Top': rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Top, 'Left': rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Left, 'Right': rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Right}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Bottom, 1: rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Top, 2: rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Left, 3: rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Right}
class MolDraw2D(Boost.Python.instance):
    """
    Drawer abstract base class
    """
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def ClearDrawing(self) -> None:
        """
            clears the drawing by filling it with the background color
        
        """
    def DrawArc(self, center: Point2D, radius: float, angle1: float, angle2: float, rawCoords: bool = False) -> None:
        """
            draws an arc with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.The angles are in degrees; angle2 should be > angle1.
        
        """
    def DrawArrow(self, cds1: Point2D, cds2: Point2D, asPolygon: bool = False, frac: float = 0.05, angle: float = 0.5235987755982988, color: typing.Any = None, rawCoords: bool = False) -> None:
        """
            draws an arrow with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels. If asPolygon is true the head of the arrow will be drawn as a triangle, otherwise two lines are used. The fraction of the arrow length to use for the head is given by frac. The angle of the arrowhead (the angle between the main line and each arrowhead line) is given by angle. The color is a tuple of 3 floats (0-1) in red, green, blue (RGB) order.
        
        """
    def DrawAttachmentLine(self, cds1: Point2D, cds2: Point2D, color: tuple, len: float = 1.0, nSegments: int = 16, rawCoords: bool = False) -> None:
        """
            draw a line indicating the presence of an attachment point (normally a squiggle line perpendicular to a bond). The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        
        """
    def DrawEllipse(self, cds1: Point2D, cds2: Point2D, rawCoords: bool = False) -> None:
        """
            draws a triangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        
        """
    def DrawLine(self, cds1: Point2D, cds2: Point2D, rawCoords: bool = False) -> None:
        """
            draws a line with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        
        """
    @typing.overload
    def DrawMolecule(self, mol: Mol, highlightAtoms: typing.Any = None, highlightAtomColors: typing.Any = None, highlightAtomRadii: typing.Any = None, confId: int = -1, legend: str = '') -> None:
        """
            renders a molecule
            
        
        """
    @typing.overload
    def DrawMolecule(self, mol: Mol, highlightAtoms: typing.Any, highlightBonds: typing.Any, highlightAtomColors: typing.Any = None, highlightBondColors: typing.Any = None, highlightAtomRadii: typing.Any = None, confId: int = -1, legend: str = '') -> None:
        """
            renders a molecule
            
        
        """
    def DrawMoleculeWithHighlights(self, mol: Mol, legend: str, highlight_atom_map: typing.Any, highlight_bond_map: typing.Any, highlight_radii: typing.Any, highlight_linewidth_multipliers: typing.Any, confId: int = -1) -> None:
        """
            renders a molecule with multiple highlight colours
            
        
        """
    def DrawMolecules(self, mols: typing.Any, highlightAtoms: typing.Any = None, highlightBonds: typing.Any = None, highlightAtomColors: typing.Any = None, highlightBondColors: typing.Any = None, highlightAtomRadii: typing.Any = None, confIds: typing.Any = None, legends: typing.Any = None) -> None:
        """
            renders multiple molecules
            
        
        """
    def DrawPolygon(self, cds: typing.Any, rawCoords: bool = False) -> None:
        """
            draws a polygon with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        
        """
    def DrawReaction(self, rxn: typing.Any, highlightByReactant: bool = False, highlightColorsReactants: typing.Any = None, confIds: typing.Any = None) -> None:
        """
            renders a reaction
            
        
        """
    def DrawRect(self, cds1: Point2D, cds2: Point2D, rawCoords: bool = False) -> None:
        """
            draws a rectangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        
        """
    @typing.overload
    def DrawString(self, string: str, pos: Point2D, rawCoords: bool = False) -> None:
        """
            add text to the canvas. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        
        """
    @typing.overload
    def DrawString(self, string: str, pos: Point2D, align: int, rawCoords: bool = False) -> None:
        """
            add aligned text to the canvas. The align argument can be 0 (=MIDDLE), 1 (=START), or 2 (=END).The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        
        """
    def DrawTriangle(self, cds1: Point2D, cds2: Point2D, cds3: Point2D, rawCoords: bool = False) -> None:
        """
            draws a triangle with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        
        """
    def DrawWavyLine(self, cds1: Point2D, cds2: Point2D, color1: tuple, color2: tuple, nSegments: int = 16, vertOffset: float = 0.05, rawCoords: bool = False) -> None:
        """
            draw a line indicating the presence of an attachment point (normally a squiggle line perpendicular to a bond). The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        
        """
    def FillPolys(self) -> bool:
        """
            returns whether or not polygons are being filled
        
        """
    def FlexiMode(self) -> bool:
        """
            returns whether or not FlexiMode is being used
        
        """
    def FontSize(self) -> float:
        """
            get the default font size. The units are, roughly, pixels.
        
        """
    @typing.overload
    def GetDrawCoords(self, point: Point2D) -> Point2D:
        """
            get the coordinates in drawing space for a particular point in molecule space
        
        """
    @typing.overload
    def GetDrawCoords(self, atomIndex: int) -> Point2D:
        """
            get the coordinates in drawing space for a particular atom
        
        """
    def GetMolSize(self, mol: Mol, highlightAtoms: typing.Any = None, highlightBonds: typing.Any = None, highlightAtomColors: typing.Any = None, highlightBondColors: typing.Any = None, highlightAtomRadii: typing.Any = None, confId: int = -1, legend: str = '') -> tuple:
        """
            returns the width and height required to draw a molecule at the current size
        
        """
    def Height(self) -> int:
        """
            get the height of the drawing canvas
        
        """
    def LineWidth(self) -> float:
        """
            returns the line width being used
        
        """
    def Offset(self) -> Point2D:
        """
            returns the offset (in drawing coordinates) for the drawing
        
        """
    def SetColour(self, tpl: tuple) -> None:
        """
            set the color being used fr drawing and filling
        
        """
    def SetDrawOptions(self, opts: MolDrawOptions) -> None:
        """
            Copies the drawing options passed in over our drawing options
        
        """
    def SetFillPolys(self, val: bool) -> None:
        """
            sets whether or not polygons are filled
        
        """
    def SetFlexiMode(self, mode: bool) -> None:
        """
            when FlexiMode is set, molecules will always been drawn with the default values for bond length, font size, etc.
        
        """
    def SetFontSize(self, new_size: float) -> None:
        """
            change the default font size. The units are, roughly, pixels.
        
        """
    def SetLineWidth(self, width: float) -> None:
        """
            set the line width being used
        
        """
    def SetOffset(self, x: int, y: int) -> None:
        """
            set the offset (in drawing coordinates) for the drawing
        
        """
    def SetScale(self, width: int, height: int, minv: Point2D, maxv: Point2D, mol: typing.Any = None) -> None:
        """
            uses the values provided to set the drawing scaling
        
        """
    def Width(self) -> int:
        """
            get the width of the drawing canvas
        
        """
    def drawOptions(self) -> MolDrawOptions:
        """
            Returns a modifiable version of the current drawing options
        
        """
class MolDraw2DCairo(MolDraw2D):
    """
    Cairo molecule drawer
    """
    __instance_size__: typing.ClassVar[int] = 1008
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def FinishDrawing(self) -> None:
        """
            add the last bits to finish the drawing
        
        """
    def GetDrawingText(self) -> typing.Any:
        """
            return the PNG data as a string
        
        """
    def WriteDrawingText(self, fName: str) -> None:
        """
            write the PNG data to the named file
        
        """
    def __init__(self, width: int, height: int, panelWidth: int = -1, panelHeight: int = -1, noFreetype: bool = False) -> None:
        ...
class MolDraw2DSVG(MolDraw2D):
    """
    SVG molecule drawer
    """
    __instance_size__: typing.ClassVar[int] = 1280
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddMoleculeMetadata(self, mol: Mol, confId: int = -1) -> None:
        """
            add RDKit-specific information to the bottom of the drawing
        
        """
    def FinishDrawing(self) -> None:
        """
            add the last bits of SVG to finish the drawing
        
        """
    def GetDrawingText(self) -> str:
        """
            return the SVG
        
        """
    def TagAtoms(self, mol: Mol, radius: float = 0.2, events: typing.Any = None) -> None:
        """
            allow atom selection in the SVG
        
        """
    def __init__(self, width: int, height: int, panelWidth: int = -1, panelHeight: int = -1, noFreetype: bool = False) -> None:
        ...
class MolDrawOptions(Boost.Python.instance):
    """
    Drawing options
    """
    __instance_size__: typing.ClassVar[int] = 752
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    def getAnnotationColour(self) -> typing.Any:
        """
            method returning the annotation colour
        
        """
    def getAtomNoteColour(self) -> typing.Any:
        """
            method returning the atom note colour
        
        """
    def getAtomPalette(self) -> dict:
        """
            returns the current atom palette as a dictionary mapping ints to 4-tuples
        
        """
    def getBackgroundColour(self) -> typing.Any:
        """
            method returning the background colour
        
        """
    def getBondNoteColour(self) -> typing.Any:
        """
            method returning the bond note colour
        
        """
    def getHighlightColour(self) -> typing.Any:
        """
            method returning the highlight colour
        
        """
    def getLegendColour(self) -> typing.Any:
        """
            method returning the legend colour
        
        """
    def getQueryColour(self) -> typing.Any:
        """
            method returning the query colour
        
        """
    def getSymbolColour(self) -> typing.Any:
        """
            method returning the symbol colour
        
        """
    def getVariableAttachmentColour(self) -> typing.Any:
        """
            method for getting the colour of variable attachment points
        
        """
    def setAnnotationColour(self, tpl: tuple) -> None:
        """
            method for setting the annotation colour
        
        """
    def setAtomNoteColour(self, tpl: tuple) -> None:
        """
            method for setting the atom note colour
        
        """
    def setAtomPalette(self, cmap: typing.Any) -> None:
        """
            sets the palette for atoms and bonds from a dictionary mapping ints to 3-tuples
        
        """
    def setBackgroundColour(self, tpl: tuple) -> None:
        """
            method for setting the background colour
        
        """
    def setBondNoteColour(self, tpl: tuple) -> None:
        """
            method for setting the bond note colour
        
        """
    def setHighlightColour(self, tpl: tuple) -> None:
        """
            method for setting the highlight colour
        
        """
    def setLegendColour(self, tpl: tuple) -> None:
        """
            method for setting the legend colour
        
        """
    def setQueryColour(self, tpl: tuple) -> None:
        """
            method for setting the query colour
        
        """
    def setSymbolColour(self, tpl: tuple) -> None:
        """
            method for setting the symbol colour
        
        """
    def setVariableAttachmentColour(self, tpl: tuple) -> None:
        """
            method for setting the colour of variable attachment points
        
        """
    def updateAtomPalette(self, cmap: typing.Any) -> None:
        """
            updates the palette for atoms and bonds from a dictionary mapping ints to 3-tuples
        
        """
    def useAvalonAtomPalette(self) -> None:
        """
            use the Avalon renderer palette for atoms and bonds
        
        """
    def useBWAtomPalette(self) -> None:
        """
            use a black and white palette for atoms and bonds
        
        """
    def useCDKAtomPalette(self) -> None:
        """
            use the CDK palette for atoms and bonds
        
        """
    def useDefaultAtomPalette(self) -> None:
        """
            use the default colour palette for atoms and bonds
        
        """
    @property
    def addAtomIndices(self) -> bool:
        """adds atom indices to drawings. Default False. (default: False)"""
    @addAtomIndices.setter
    def addAtomIndices(self, value: bool) -> None: ...
    @property
    def addBondIndices(self) -> bool:
        """adds bond indices to drawings. Default False. (default: False)"""
    @addBondIndices.setter
    def addBondIndices(self, value: bool) -> None: ...
    @property
    def addStereoAnnotation(self) -> bool:
        """adds R/S and E/Z to drawings. Default False. (default: False)"""
    @addStereoAnnotation.setter
    def addStereoAnnotation(self, value: bool) -> None: ...
    @property
    def addStereoGroupAnnotation(self) -> bool:
        """Whether to add the enhanced stereo labels.  Default is True. (default: True)"""
    @addStereoGroupAnnotation.setter
    def addStereoGroupAnnotation(self, value: bool) -> None: ...
    @property
    def additionalAtomLabelPadding(self) -> float:
        """additional padding to leave around atom labels. Expressed as a fraction of the font size. (default: 0.0)"""
    @additionalAtomLabelPadding.setter
    def additionalAtomLabelPadding(self, value: float) -> None: ...
    @property
    def annotationColour(*args, **kwargs):
        """
        the annotation colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """
    @annotationColour.setter
    def annotationColour(*args, **kwargs):
        ...
    @property
    def annotationFontScale(self) -> float:
        """Scale of font for atom and bond annotation relative to atomlabel font.  Default=0.75. (default: 0.5)"""
    @annotationFontScale.setter
    def annotationFontScale(self, value: float) -> None: ...
    @property
    def atomHighlightsAreCircles(self) -> bool:
        """forces atom highlights always to be circles.Default (false) is to put ellipses roundlonger labels. (default: False)"""
    @atomHighlightsAreCircles.setter
    def atomHighlightsAreCircles(self, value: bool) -> None: ...
    @property
    def atomLabelDeuteriumTritium(self) -> bool:
        """labels deuterium as D and tritium as T (default: False)"""
    @atomLabelDeuteriumTritium.setter
    def atomLabelDeuteriumTritium(self, value: bool) -> None: ...
    @property
    def atomLabels(*args, **kwargs):
        """
        maps indices to atom labels
        """
    @atomLabels.setter
    def atomLabels(*args, **kwargs):
        ...
    @property
    def atomNoteColour(*args, **kwargs):
        """
        the atom note colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """
    @atomNoteColour.setter
    def atomNoteColour(*args, **kwargs):
        ...
    @property
    def atomRegions(*args, **kwargs):
        """
        regions to outline
        """
    @atomRegions.setter
    def atomRegions(*args, **kwargs):
        ...
    @property
    def backgroundColour(*args, **kwargs):
        """
        the background colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """
    @backgroundColour.setter
    def backgroundColour(*args, **kwargs):
        ...
    @property
    def baseFontSize(self) -> float:
        """relative size of font.  Defaults to 0.6.  -1 means use default. (default: -1.0)"""
    @baseFontSize.setter
    def baseFontSize(self, value: float) -> None: ...
    @property
    def bondLineWidth(self) -> float:
        """if positive, this overrides the default line width for bonds (default: 2.0)"""
    @bondLineWidth.setter
    def bondLineWidth(self, value: float) -> None: ...
    @property
    def bondNoteColour(*args, **kwargs):
        """
        the bond note colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """
    @bondNoteColour.setter
    def bondNoteColour(*args, **kwargs):
        ...
    @property
    def bracketsAroundAtomLists(self) -> bool:
        """Whether to put brackets round atom lists in query atoms.  Default is true. (default: True)"""
    @bracketsAroundAtomLists.setter
    def bracketsAroundAtomLists(self, value: bool) -> None: ...
    @property
    def centreMoleculesBeforeDrawing(self) -> bool:
        """Moves the centre of the drawn molecule to (0,0).Default False. (default: False)"""
    @centreMoleculesBeforeDrawing.setter
    def centreMoleculesBeforeDrawing(self, value: bool) -> None: ...
    @property
    def circleAtoms(self) -> bool:
        """default: True"""
    @circleAtoms.setter
    def circleAtoms(self, value: bool) -> None: ...
    @property
    def clearBackground(self) -> bool:
        """clear the background before drawing a molecule (default: True)"""
    @clearBackground.setter
    def clearBackground(self, value: bool) -> None: ...
    @property
    def comicMode(self) -> bool:
        """simulate hand-drawn lines for bonds. When combined with a font like Comic-Sans or Comic-Neue, this gives xkcd-like drawings. Default is false. (default: False)"""
    @comicMode.setter
    def comicMode(self, value: bool) -> None: ...
    @property
    def continuousHighlight(self) -> bool:
        """default: True"""
    @continuousHighlight.setter
    def continuousHighlight(self, value: bool) -> None: ...
    @property
    def drawMolsSameScale(self) -> bool:
        """when drawing multiple molecules with DrawMolecules, forces them to use the same scale.  Default is true. (default: True)"""
    @drawMolsSameScale.setter
    def drawMolsSameScale(self, value: bool) -> None: ...
    @property
    def drawingExtentsInclude(self) -> int:
        """Drawing extents are computed taking into account only selected DrawElement items.  Default=DrawElement.ALL (default: 2147483647)"""
    @drawingExtentsInclude.setter
    def drawingExtentsInclude(self, value: int) -> None: ...
    @property
    def dummiesAreAttachments(self) -> bool:
        """default: False"""
    @dummiesAreAttachments.setter
    def dummiesAreAttachments(self, value: bool) -> None: ...
    @property
    def dummyIsotopeLabels(self) -> bool:
        """adds isotope labels on dummy atoms. Default True. (default: True)"""
    @dummyIsotopeLabels.setter
    def dummyIsotopeLabels(self, value: bool) -> None: ...
    @property
    def explicitMethyl(self) -> bool:
        """Draw terminal methyls explictly.  Default is false. (default: False)"""
    @explicitMethyl.setter
    def explicitMethyl(self, value: bool) -> None: ...
    @property
    def fillHighlights(self) -> bool:
        """default: True"""
    @fillHighlights.setter
    def fillHighlights(self, value: bool) -> None: ...
    @property
    def fixedBondLength(self) -> float:
        """If > 0.0, fixes bond length to this number of pixelsunless that would make it too big.  Default -1.0 meansno fix.  If both set, fixedScale takes precedence. (default: -1.0)"""
    @fixedBondLength.setter
    def fixedBondLength(self, value: float) -> None: ...
    @property
    def fixedFontSize(self) -> int:
        """font size in pixels. default=-1 means not fixed.  If set, always used irrespective of scale, minFontSize and maxFontSize. (default: -1)"""
    @fixedFontSize.setter
    def fixedFontSize(self, value: int) -> None: ...
    @property
    def fixedScale(self) -> float:
        """If > 0.0, fixes scale to that fraction of width ofdraw window.  Default -1.0 means adjust scale to fit. (default: -1.0)"""
    @fixedScale.setter
    def fixedScale(self, value: float) -> None: ...
    @property
    def flagCloseContactsDist(self) -> int:
        """default: 3"""
    @flagCloseContactsDist.setter
    def flagCloseContactsDist(self, value: int) -> None: ...
    @property
    def fontFile(self) -> str:
        """Font file for use with FreeType text drawer.  Can also be BuiltinTelexRegular (the default) or BuiltinRobotoRegular. (default: '')"""
    @fontFile.setter
    def fontFile(self, value: str) -> None: ...
    @property
    def highlightBondWidthMultiplier(self) -> int:
        """What to multiply default bond width by for highlighting bonds. Default-8. (default: 8)"""
    @highlightBondWidthMultiplier.setter
    def highlightBondWidthMultiplier(self, value: int) -> None: ...
    @property
    def highlightColour(*args, **kwargs):
        """
        the highlight colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """
    @highlightColour.setter
    def highlightColour(*args, **kwargs):
        ...
    @property
    def highlightRadius(self) -> float:
        """Default radius for highlight circles. (default: 0.3)"""
    @highlightRadius.setter
    def highlightRadius(self, value: float) -> None: ...
    @property
    def includeAtomTags(self) -> bool:
        """include atom tags in output (default: False)"""
    @includeAtomTags.setter
    def includeAtomTags(self, value: bool) -> None: ...
    @property
    def includeChiralFlagLabel(self) -> bool:
        """add a molecule annotation with "ABS" if the chiral flag is set. Default is false. (default: False)"""
    @includeChiralFlagLabel.setter
    def includeChiralFlagLabel(self, value: bool) -> None: ...
    @property
    def includeMetadata(self) -> bool:
        """When possible, include metadata about molecules and reactions to allow them to be reconstructed. Default is true. (default: True)"""
    @includeMetadata.setter
    def includeMetadata(self, value: bool) -> None: ...
    @property
    def includeRadicals(self) -> bool:
        """include radicals in the drawing (it can be useful to turn this off for reactions and queries). Default is true. (default: True)"""
    @includeRadicals.setter
    def includeRadicals(self, value: bool) -> None: ...
    @property
    def isotopeLabels(self) -> bool:
        """adds isotope labels on non-dummy atoms. Default True. (default: True)"""
    @isotopeLabels.setter
    def isotopeLabels(self, value: bool) -> None: ...
    @property
    def legendColour(*args, **kwargs):
        """
        the legend colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """
    @legendColour.setter
    def legendColour(*args, **kwargs):
        ...
    @property
    def legendFontSize(self) -> int:
        """font size in pixels of the legend (if drawn) (default: 16)"""
    @legendFontSize.setter
    def legendFontSize(self, value: int) -> None: ...
    @property
    def legendFraction(self) -> float:
        """fraction of the draw panel to be used for the legend if present (default: 0.1)"""
    @legendFraction.setter
    def legendFraction(self, value: float) -> None: ...
    @property
    def legendPosition(*args, **kwargs):
        """
        legend position enum. Default=Bottom. Values: LegendPosition.Bottom, LegendPosition.Top, LegendPosition.Left, LegendPosition.Right.
        """
    @legendPosition.setter
    def legendPosition(*args, **kwargs):
        ...
    @property
    def legendVerticalText(self) -> bool:
        """when legend is Left or Right, draw text vertically (one char per line) (default: True)"""
    @legendVerticalText.setter
    def legendVerticalText(self, value: bool) -> None: ...
    @property
    def maxFontSize(self) -> int:
        """maximum font size in pixels. default=40, -1 means no maximum. (default: 40)"""
    @maxFontSize.setter
    def maxFontSize(self, value: int) -> None: ...
    @property
    def minFontSize(self) -> int:
        """minimum font size in pixels. default=6, -1 means no minimum. (default: 6)"""
    @minFontSize.setter
    def minFontSize(self, value: int) -> None: ...
    @property
    def multiColourHighlightStyle(*args, **kwargs):
        """
        Either 'CircleAndLine' or 'Lasso', to control style ofmulti-coloured highlighting in DrawMoleculeWithHighlights.Default is CircleAndLine.
        """
    @multiColourHighlightStyle.setter
    def multiColourHighlightStyle(*args, **kwargs):
        ...
    @property
    def multipleBondOffset(self) -> float:
        """offset for the extra lines in a multiple bond as a fraction of mean bond length (default: 0.15)"""
    @multipleBondOffset.setter
    def multipleBondOffset(self, value: float) -> None: ...
    @property
    def noAtomLabels(self) -> bool:
        """disables inclusion of atom labels in the rendering (default: False)"""
    @noAtomLabels.setter
    def noAtomLabels(self, value: bool) -> None: ...
    @property
    def padding(self) -> float:
        """Fraction of empty space to leave around molecule.  Default=0.05. (default: 0.05)"""
    @padding.setter
    def padding(self, value: float) -> None: ...
    @property
    def prepareMolsBeforeDrawing(self) -> bool:
        """call prepareMolForDrawing() on each molecule passed to DrawMolecules() (default: True)"""
    @prepareMolsBeforeDrawing.setter
    def prepareMolsBeforeDrawing(self, value: bool) -> None: ...
    @property
    def queryColour(*args, **kwargs):
        """
        the query colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """
    @queryColour.setter
    def queryColour(*args, **kwargs):
        ...
    @property
    def reagentPadding(self) -> float:
        """Fraction of empty space to leave around each component of a reaction drawing.  Default=0.0. (default: 0.0)"""
    @reagentPadding.setter
    def reagentPadding(self, value: float) -> None: ...
    @property
    def rotate(self) -> float:
        """Rotates molecule about centre by this number of degrees, (default: 0.0)"""
    @rotate.setter
    def rotate(self, value: float) -> None: ...
    @property
    def scaleBondWidth(self) -> bool:
        """Scales the width of drawn bonds using image scaling. (default: False)"""
    @scaleBondWidth.setter
    def scaleBondWidth(self, value: bool) -> None: ...
    @property
    def scaleHighlightBondWidth(self) -> bool:
        """Scales the width of drawn highlighted bonds using image scaling. (default: True)"""
    @scaleHighlightBondWidth.setter
    def scaleHighlightBondWidth(self, value: bool) -> None: ...
    @property
    def scalingFactor(self) -> float:
        """scaling factor for pixels->angstrom when auto scalingbeing used.  Default is 20. (default: 20.0)"""
    @scalingFactor.setter
    def scalingFactor(self, value: float) -> None: ...
    @property
    def showAllCIPCodes(self) -> bool:
        """show all defined CIP codes (no hiding!). Default False. (default: False)"""
    @showAllCIPCodes.setter
    def showAllCIPCodes(self, value: bool) -> None: ...
    @property
    def simplifiedStereoGroupLabel(self) -> bool:
        """if all specified stereocenters are in a single StereoGroup, show a molecule-level annotation instead of the individual labels. Default is false. (default: False)"""
    @simplifiedStereoGroupLabel.setter
    def simplifiedStereoGroupLabel(self, value: bool) -> None: ...
    @property
    def singleColourBonds(self) -> bool:
        """if true all bonds are drawn using symbolColour rather than inheriting their colour from the atoms. Default is false. (default: False)"""
    @singleColourBonds.setter
    def singleColourBonds(self, value: bool) -> None: ...
    @property
    def singleColourWedgeBonds(self) -> bool:
        """if true wedged and dashed bonds are drawn using symbolColour rather than inheriting their colour from the atoms. Default is false. (default: False)"""
    @singleColourWedgeBonds.setter
    def singleColourWedgeBonds(self, value: bool) -> None: ...
    @property
    def splitBonds(self) -> bool:
        """default: False"""
    @splitBonds.setter
    def splitBonds(self, value: bool) -> None: ...
    @property
    def standardColoursForHighlightedAtoms(self) -> bool:
        """If true, highlighted hetero atoms are drawn in standard colours rather than black.  Default=False (default: False)"""
    @standardColoursForHighlightedAtoms.setter
    def standardColoursForHighlightedAtoms(self, value: bool) -> None: ...
    @property
    def stereoGroupAbsLabel(self) -> str:
        """String to use for enhanced stereo 'ABS' groups.  Default='abs'. (default: 'abs')"""
    @stereoGroupAbsLabel.setter
    def stereoGroupAbsLabel(self, value: str) -> None: ...
    @property
    def stereoGroupAndLabel(self) -> str:
        """String to use for enhanced stereo 'AND' groups.  Default='and'. (default: 'and')"""
    @stereoGroupAndLabel.setter
    def stereoGroupAndLabel(self, value: str) -> None: ...
    @property
    def stereoGroupOrLabel(self) -> str:
        """String to use for enhanced stereo 'OR' groups.  Default='or'. (default: 'or')"""
    @stereoGroupOrLabel.setter
    def stereoGroupOrLabel(self, value: str) -> None: ...
    @property
    def symbolColour(*args, **kwargs):
        """
        the symbol colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """
    @symbolColour.setter
    def symbolColour(*args, **kwargs):
        ...
    @property
    def unspecifiedStereoIsUnknown(self) -> bool:
        """if true, double bonds with unspecified stereo are drawn crossed, potential stereocenters with unspecified stereo are drawn with a wavy bond. Default is false. (default: False)"""
    @unspecifiedStereoIsUnknown.setter
    def unspecifiedStereoIsUnknown(self, value: bool) -> None: ...
    @property
    def useComplexQueryAtomSymbols(self) -> bool:
        """replace any atom, any hetero, any halo queries with complex query symbols A, Q, X, M, optionally followed by H if hydrogen is included (except for AH, which stays *). Default is true (default: True)"""
    @useComplexQueryAtomSymbols.setter
    def useComplexQueryAtomSymbols(self, value: bool) -> None: ...
    @property
    def useMolBlockWedging(self) -> bool:
        """If the molecule came from a MolBlock, prefer the wedging information that provides.  If false, use RDKit rules.  Default false (default: False)"""
    @useMolBlockWedging.setter
    def useMolBlockWedging(self, value: bool) -> None: ...
    @property
    def variableAtomRadius(self) -> float:
        """radius value to use for atoms involved in variable attachment points. (default: 0.4)"""
    @variableAtomRadius.setter
    def variableAtomRadius(self, value: float) -> None: ...
    @property
    def variableAttachmentColour(*args, **kwargs):
        """
        the variable attachment colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """
    @variableAttachmentColour.setter
    def variableAttachmentColour(*args, **kwargs):
        ...
    @property
    def variableBondWidthMultiplier(self) -> int:
        """what to multiply standard bond width by for variable attachment points. (default: 16)"""
    @variableBondWidthMultiplier.setter
    def variableBondWidthMultiplier(self, value: int) -> None: ...
class MultiColourHighlightStyle(Boost.Python.enum):
    CircleAndLine: typing.ClassVar[MultiColourHighlightStyle]  # value = rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.CircleAndLine
    Lasso: typing.ClassVar[MultiColourHighlightStyle]  # value = rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.Lasso
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'CircleAndLine': rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.CircleAndLine, 'Lasso': rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.Lasso}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.CircleAndLine, 1: rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.Lasso}
class map_indexing_suite_IntStringMap_entry(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __repr__(arg1: map_indexing_suite_IntStringMap_entry) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def data(self) -> str:
        """
        """
    def key(self) -> int:
        """
        """
def ContourAndDrawGaussians(drawer: MolDraw2D, locs: typing.Any, heights: typing.Any, widths: typing.Any, nContours: int = 10, levels: typing.Any = None, params: ContourParams = ..., mol: typing.Any = None) -> None:
    """
        Generates and draws contours for a set of gaussians
        
          - drawer: the MolDraw2D object to use
          - locs: locations of the gaussians
          - heights: the heights (or weights) of the gaussians
          - widths: the standard deviations of the gaussians
          - nContours: the number of contours to draw
          - levels: the contours to use
          - ps: additional parameters controlling the contouring.
          - mol: molecule used to help set scale.
        
          The values are calculated on a grid with spacing params.gridResolution.
          If params.setScale  is set, the grid size will be calculated based on the
          locations of the gaussians and params.extraGridPadding. Otherwise the current
          size of the viewport will be used.
        
          If the levels argument is empty, the contour levels will be determined
          automatically from the max and min values on the grid and levels will
          be updated to include the contour levels.
        
          If params.fillGrid is set, the data on the grid will also be drawn using
          the color scheme in params.colourMap
        
          If mol is not 0, uses the molecule to help set the scale, assuming that
          it will be drawn over the plot, so needs to fit on it.
        */
    
    """
def ContourAndDrawGrid(drawer: MolDraw2D, data: typing.Any, xcoords: typing.Any, ycoords: typing.Any, nContours: int = 10, levels: typing.Any = None, params: ContourParams = ..., mol: typing.Any = None) -> None:
    """
        Generates and draws contours for data on a grid
        
          - drawer: the MolDraw2D object to use
          - data: numpy array with the data to be contoured
          - xcoords: the x coordinates of the grid
          - ycoords: the y coordinates of the grid
          - nContours: the number of contours to draw
          - levels: the contours to use
          - ps: additional parameters controlling the contouring
          - mol: molecule used to help set scale.
        
          The values are calculated on a grid with spacing params.gridResolution.
          If params.setScale  is set, the grid size will be calculated based on the
          locations of the gaussians and params.extraGridPadding. Otherwise the current
          size of the viewport will be used.
        
          If the levels argument is empty, the contour levels will be determined
          automatically from the max and min values on the grid and levels will
          be updated to include the contour levels.
        
          If params.fillGrid is set, the data on the grid will also be drawn using
          the color scheme in params.colourMap
        
          If mol is not 0, uses the molecule to help set the scale, assuming that
          it will be drawn over the plot, so needs to fit on it.
        */
    
    """
def DrawMoleculeACS1996(drawer: MolDraw2D, mol: Mol, legend: str = '', highlightAtoms: typing.Any = None, highlightBonds: typing.Any = None, highlightAtomColors: typing.Any = None, highlightBondColors: typing.Any = None, highlightAtomRadii: typing.Any = None, confId: int = -1) -> None:
    """
        Draws molecule in ACS 1996 mode.
    
    """
def MeanBondLength(mol: Mol, confId: int = -1) -> float:
    """
        Calculate the mean bond length for the molecule.
    
    """
def MolToACS1996SVG(mol: Mol, legend: str = '', highlightAtoms: typing.Any = None, highlightBonds: typing.Any = None, highlightAtomColors: typing.Any = None, highlightBondColors: typing.Any = None, highlightAtomRadii: typing.Any = None, confId: int = -1) -> str:
    """
        Returns ACS 1996 mode svg for a molecule
    
    """
def MolToSVG(mol: Mol, width: int = 300, height: int = 300, highlightAtoms: typing.Any = None, kekulize: bool = True, lineWidthMult: int = 1, includeAtomCircles: bool = True, confId: int = -1) -> str:
    """
        Returns svg for a molecule
    
    """
def PrepareAndDrawMolecule(drawer: MolDraw2D, mol: Mol, legend: str = '', highlightAtoms: typing.Any = None, highlightBonds: typing.Any = None, highlightAtomColors: typing.Any = None, highlightBondColors: typing.Any = None, highlightAtomRadii: typing.Any = None, confId: int = -1, kekulize: bool = True) -> None:
    """
        Preps a molecule for drawing and actually draws it
        
    
    """
def PrepareMolForDrawing(mol: Mol, kekulize: bool = True, addChiralHs: bool = True, wedgeBonds: bool = True, forceCoords: bool = False, wavyBonds: bool = False) -> rdkit.Chem.Mol:
    """
        Does some cleanup operations on the molecule to prepare it to draw nicely.
        The operations include: kekulization, addition of chiral Hs (so that we can draw
        wedges to them), wedging of bonds at chiral centers, and generation of a 2D
        conformation if the molecule does not already have a conformation
        
        Returns a modified copy of the molecule.
        
    
    """
def SetACS1996Mode(drawOptions: MolDrawOptions, meanBondLength: float) -> None:
    """
        Set the draw options to produce something as close as possible to
        the ACS 1996 guidelines as described at
        https://en.wikipedia.org/wiki/Wikipedia:Manual_of_Style/Chemistry/Structure_drawing
        
         - MolDrawOptions opt - the options what will be changed
         - float meanBondLength - mean bond length of the molecule
        
         Works best if the MolDraw2D object is created with width and height -1 (a
         flexiCanvas).
         The mean bond length may be calculated with MeanBondLength.
         It is used to calculate the offset for the lines in multiple bonds.
        
         Options changed are:
           bondLineWidth = 0.6
           scaleBondWidth = false
           scalingFactor = 14.4 / meanBondLen
           multipleBondOffset = 0.18
           highlightBondWidthMultiplier = 32
           setMonochromeMode - black and white
           fixedFontSize = 10
           additionalAtomLabelPadding = 0.066
           fontFile - if it isn't set already, then if RDBASE is set and the file
                      exists, uses $RDBASE/Data/Fonts/FreeSans.ttf.  Otherwise uses
                      BuiltinRobotoRegular.
         */
        
    
    """
@typing.overload
def SetDarkMode(d2d: MolDrawOptions) -> None:
    """
        set dark mode for a MolDrawOptions object
    
    """
@typing.overload
def SetDarkMode(d2d: MolDraw2D) -> None:
    """
        set dark mode for a MolDraw2D object
    
    """
@typing.overload
def SetMonochromeMode(options: MolDrawOptions, fgColour: tuple, bgColour: tuple) -> None:
    """
        set monochrome mode for a MolDrawOptions object
    
    """
@typing.overload
def SetMonochromeMode(drawer: MolDraw2D, fgColour: tuple, bgColour: tuple) -> None:
    """
        set monochrome mode for a MolDraw2D object
    
    """
def UpdateDrawerParamsFromJSON(drawer: MolDraw2D, json: str) -> None:
    """
    """
def UpdateMolDrawOptionsFromJSON(opts: MolDrawOptions, json: str) -> None:
    """
    """
ALL: DrawElement  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ALL
ANNOTATIONS: DrawElement  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ANNOTATIONS
ATOMLABELS: DrawElement  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.ATOMLABELS
BONDS: DrawElement  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.BONDS
Bottom: LegendPosition  # value = rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Bottom
CircleAndLine: MultiColourHighlightStyle  # value = rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.CircleAndLine
HIGHLIGHTS: DrawElement  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.HIGHLIGHTS
Lasso: MultiColourHighlightStyle  # value = rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.Lasso
Left: LegendPosition  # value = rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Left
NONE: DrawElement  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.NONE
POSTSHAPES: DrawElement  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.POSTSHAPES
PRESHAPES: DrawElement  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.PRESHAPES
RADICALS: DrawElement  # value = rdkit.Chem.Draw.rdMolDraw2D.DrawElement.RADICALS
Right: LegendPosition  # value = rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Right
Top: LegendPosition  # value = rdkit.Chem.Draw.rdMolDraw2D.LegendPosition.Top
