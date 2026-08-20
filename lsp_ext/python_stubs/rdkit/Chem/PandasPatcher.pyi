# fix_pybind_stubs: rdkit 2026.3.5 5beea910
from _typeshed import Incomplete
from rdkit.Chem.rdchem import Mol as Mol

RDK_MOLS_AS_IMAGE_ATTR: str
InteractiveRenderer: None
PrintAsImageString: None
molJustify: None
pandas_frame: None
pandas_formats_name: str
to_html_class_name: str
get_adjustment_module_name: str
get_adjustment_name: str
html_formatter_module_name: str
def is_molecule_image(s): ...

class MolFormatter:
    def __init__(self, orig_formatter: Incomplete | None = ...) -> None: ...
    @staticmethod
    def default_formatter(x): ...
    @staticmethod
    def is_mol(x): ...
    @classmethod
    def get_formatters(cls, df, orig_formatters): ...
    def __call__(self, x): ...
def check_rdk_attr(frame, attr): ...
def set_rdk_attr(frame, attr): ...
def patched_to_html(self, *args, **kwargs): ...
def patched_get_formatter(self, i, *args, **kwargs): ...
def patched_write_cell(self, s, *args, **kwargs): ...
def patched_get_adjustment(): ...
def renderImagesInAllDataFrames(images: bool = ...): ...
def changeMoleculeRendering(frame, renderer: str = ...): ...
def patchPandas(): ...
def unpatchPandas(): ...

# present at runtime, absent from the generated stub:
from xml.parsers.expat import ExpatError as ExpatError
from _io import StringIO as StringIO
def dataframe_applymap(self, func, na_action=..., **kwargs):
    r"""
    Apply a function to a Dataframe elementwise.

    .. versionadded:: 2.1.0

       DataFrame.applymap was deprecated and renamed to DataFrame.map.

    This method applies a function that accepts and returns a scalar
    to every element of a DataFrame.

    Parameters
    ----------
    func : callable
        Python function, returns a single value from a single value.
    na_action : {None, 'ignore'}, default None
        If 'ignore', propagate NaN values, without passing them to func.
    **kwargs
        Additional keyword arguments to pass as keywords arguments to
        `func`.

    Returns
    -------
    DataFrame
        Transformed DataFrame.

    See Also
    --------
    DataFrame.apply : Apply a function along input axis of DataFrame.
    DataFrame.replace: Replace values given in `to_replace` with `value`.
    Series.map : Apply a function elementwise on a Series.

    Examples
    --------
    >>> df = pd.DataFrame([[1, 2.12], [3.356, 4.567]])
    >>> df
           0      1
    0  1.000  2.120
    1  3.356  4.567

    >>> df.map(lambda x: len(str(x)))
       0  1
    0  3  4
    1  5  5

    Like Series.map, NA values can be ignored:

    >>> df_copy = df.copy()
    >>> df_copy.iloc[0, 0] = pd.NA
    >>> df_copy.map(lambda x: len(str(x)), na_action="ignore")
         0  1
    0  NaN  4
    1  5.0  5

    It is also possible to use `map` with functions that are not
    `lambda` functions:

    >>> df.map(round, ndigits=1)
         0    1
    0  1.0  2.1
    1  3.4  4.6

    Note that a vectorized version of `func` often exists, which will
    be much faster. You could square each number elementwise.

    >>> df.map(lambda x: x**2)
               0          1
    0   1.000000   4.494400
    1  11.262736  20.857489

    But it's better to avoid map in that case.

    >>> df**2
               0          1
    0   1.000000   4.494400
    1  11.262736  20.857489
    """
def dataframeformatter_class(frame, columns=..., col_space=..., header=..., index=..., na_rep=..., formatters=..., justify=..., float_format=..., sparsify=..., index_names=..., max_rows=..., min_rows=..., max_cols=..., show_dimensions=..., decimal=..., bold_rows=..., escape=...):
    r"""
    Class for processing dataframe formatting options and data.

    Used by DataFrame.to_string, which backs DataFrame.__repr__.

        Parameters
        ----------
        buf : str, Path or StringIO-like, optional, default None
            Buffer to write to. If None, the output is returned as a string.
        columns : array-like, optional, default None
            The subset of columns to write. Writes all columns by default.
        col_space : %(col_space_type)s, optional
            %(col_space)s
        header : %(header_type)s, optional
            %(header)s.
        index : bool, optional, default True
            Whether to print index (row) labels.
        na_rep : str, optional, default 'NaN'
            String representation of ``NaN`` to use.
        formatters : list, tuple or dict of one-param. functions, optional
            Formatter functions to apply to columns' elements by position or
            name.
            The result of each function must be a unicode string.
            List/tuple must be of length equal to the number of columns.
        float_format : one-parameter function, optional, default None
            Formatter function to apply to columns' elements if they are
            floats. This function must return a unicode string and will be
            applied only to the non-``NaN`` elements, with ``NaN`` being
            handled by ``na_rep``.
        sparsify : bool, optional, default True
            Set to False for a DataFrame with a hierarchical index to print
            every multiindex key at each row.
        index_names : bool, optional, default True
            Prints the names of the indexes.
        justify : str, default None
            How to justify the column labels. If None uses the option from
            the print configuration (controlled by set_option), 'right' out
            of the box. Valid values are

            * left
            * right
            * center
            * justify
            * justify-all
            * start
            * end
            * inherit
            * match-parent
            * initial
            * unset.
        max_rows : int, optional
            Maximum number of rows to display in the console.
        max_cols : int, optional
            Maximum number of columns to display in the console.
        show_dimensions : bool, default False
            Display DataFrame dimensions (number of rows by number of columns).
        decimal : str, default '.'
            Character recognized as decimal separator, e.g. ',' in Europe.

        Returns
        -------
        str or None
            If buf is None, returns the result as a string. Otherwise returns
            None.
    """
from pandas.io.formats import printing as get_adjustment_module
def html_formatter_class(formatter, classes=..., border=..., table_id=..., render_links=...):
    r"""
    Internal class for formatting output data in html.
    This class is intended for shared functionality between
    DataFrame.to_html() and DataFrame._repr_html_().
    Any logic in common with other output formatting methods
    should ideally be inherited from classes in format.py
    and this class responsible for only producing html markup.
    """
from pandas.io.formats import html as html_formatter_module
import importlib as importlib
import logging as logging
from xml.dom import minidom as minidom
def orig_get_adjustment():
    ...
def orig_get_formatter(self, i):
    ...
def orig_to_html(self, buf=..., encoding=..., classes=..., notebook=..., border=..., table_id=..., render_links=...):
    r"""
    Render a DataFrame to an html table.

    Parameters
    ----------
    buf : str, path object, file-like object, or None, default None
        String, path object (implementing ``os.PathLike[str]``), or file-like
        object implementing a string ``write()`` function. If None, the result is
        returned as a string.
    encoding : str, default “utf-8”
        Set character encoding.
    classes : str or list-like
        classes to include in the `class` attribute of the opening
        ``<table>`` tag, in addition to the default "dataframe".
    notebook : {True, False}, optional, default False
        Whether the generated HTML is for IPython Notebook.
    border : int or bool
        When an integer value is provided, it sets the border attribute in
        the opening tag, specifying the thickness of the border.
        If ``False`` or ``0`` is passed, the border attribute will not
        be present in the ``<table>`` tag.
        The default value for this parameter is governed by
        ``pd.options.display.html.border``.
    table_id : str, optional
        A css id is included in the opening `<table>` tag if specified.
    render_links : bool, default False
        Convert URLs to HTML links.
    """
def orig_write_cell(self, s, kind=..., indent=..., tags=...):
    ...
from pandas.io import formats as pandas_formats
import pandas as pd
from pandas.io.formats.printing import pprint_thing as pprint_thing
import re as re
def to_html_class(fmt):
    r"""
    Class for creating dataframe output in multiple formats.

    Called in pandas.core.generic.NDFrame:
        - to_csv
        - to_latex

    Called in pandas.DataFrame:
        - to_html
        - to_string

    Parameters
    ----------
    fmt : DataFrameFormatter
        Formatter with the formatting options.
    """
