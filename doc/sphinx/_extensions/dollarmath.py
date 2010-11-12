# This program is public domain
# Author: Paul Kienzle
r"""
Allow $math$ markup in text and docstrings, ignoring \$.

The $math$ markup should be separated from the surrounding text by spaces.  To
embed markup within a word, place backslash-space before and after.  For
convenience, the final $ can be followed by punctuation (period, comma or
semicolon).
"""

import re

_dollar = re.compile(r"(?:^|(?<=\s))[$]([^\n]*?)(?<![\\])[$](?:$|(?=\s|[.,;\\]))")
_notdollar = re.compile(r"\\[$]")

def replace_dollar(content):
    original = content
    content = _dollar.sub(r":math:`\1`",content)
    content = _notdollar.sub("$", content)
    #if '$' in content:
    #    import sys
    #    sys.stdout.write("\n========> not converted\n")
    #    sys.stdout.write(content)
    #    sys.stdout.write("\n")
    #elif '$' in original:
    #    import sys
    #    sys.stdout.write("\n========> converted\n")
    #    sys.stdout.write(content)
    #    sys.stdout.write("\n")
    return content

def rewrite_rst(app, docname, source):
    source[0] = replace_dollar(source[0])

def rewrite_autodoc(app, what, name, obj, options, lines):
    lines[:] = [replace_dollar(L) for L in lines]

def setup(app):
    app.connect('source-read', rewrite_rst)
    app.connect('autodoc-process-docstring', rewrite_autodoc)


def test_dollar():
    assert replace_dollar(u"no dollar")==u"no dollar"
    assert replace_dollar(u"$only$")==u":math:`only`"
    assert replace_dollar(u"$first$ is good")==u":math:`first` is good"
    assert replace_dollar(u"so is $last$")==u"so is :math:`last`"
    assert replace_dollar(u"and $mid$ too")==u"and :math:`mid` too"
    assert replace_dollar(u"$first$, $mid$, $last$")==u":math:`first`, :math:`mid`, :math:`last`"
    assert replace_dollar(u"dollar\$ escape")==u"dollar$ escape"
    assert replace_dollar(u"dollar \$escape\$ too")==u"dollar $escape$ too"
    assert replace_dollar(u"emb\ $ed$\ ed")==u"emb\ :math:`ed`\ ed"
    assert replace_dollar(u"$first$a")==u"$first$a"
    assert replace_dollar(u"a$last$")==u"a$last$"
    assert replace_dollar(u"a $mid$dle a")==u"a $mid$dle a"

if __name__ == "__main__":
    test_dollar()
