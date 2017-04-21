# This program is public domain
# Author: Paul Kienzle
"""
Table plotter
"""

__all__ = ['table_plot']

def table_plot(data, form="line", label=None, title=None):
    """
    Plot periodic table data using element symbol vs. value.

    :Parameters:
        *data* : { Element: float }
            Data values to plot

        *form* = "line" : "line|grid"
            Table layout to use

    :Returns: None
    """
    import pylab
    if form == "line":
        bbox = dict(boxstyle="round", lw=1, ec=(0, 0, 0), fc=(0.85, 0.8, 0.8))
        for el, value in data.items():
            if value is not None:
                pylab.text(el.number, value, el.symbol,
                           bbox=bbox, va='center', ha='center')
        pylab.xlim(0, 100)
        pylab.xlabel('Element number')
        values = [v for v in data.values()]
        minv, maxv = min(values), max(values)
        margin = (maxv - minv)*0.05
        pylab.ylim(minv-margin, maxv+margin)
        if label is not None:
            pylab.ylabel(label)
        if title is not None:
            pylab.title(title)
