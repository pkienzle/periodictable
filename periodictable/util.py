# This program is in the public domain
# Author: Paul Kienzle
"""
Helper functions
"""

def cell_volume(a=None, b=None, c=None, alpha=None, beta=None, gamma=None):
    r"""
    Compute cell volume from lattice parameters.

    :Parameters:
        *a*, *b*, *c* : float | |Ang|
            Lattice spacings.  *a* is required.
            *b* and *c* default to *a*.
        *alpha*, *beta*, *gamma* : float | |deg|
            Lattice angles.  *alpha* defaults to 90\ |deg|.
            *beta* and *gamma* default to *alpha*.

    :Returns:
        *V* : float | |Ang^3|
            Cell volume

    :Raises:
        *TypeError* : missing or invalid parameters

    The following formula works for all lattice types:

    .. math::

        V = a b c \sqrt{1 - \cos^2 \alpha - \cos^2 \beta - \cos^2 \gamma
                          + 2 \cos \alpha \cos \beta \cos \gamma}
    """
    from math import cos, radians, sqrt
    if a is None:
        raise TypeError('missing lattice parameters')
    if b is None:
        b = a
    if c is None:
        c = a
    calpha = cos(radians(alpha)) if alpha is not None else 0
    cbeta = cos(radians(beta)) if beta is not None else calpha
    cgamma = cos(radians(gamma)) if gamma is not None else calpha
    V = a*b*c*sqrt(1 - calpha**2 - cbeta**2 - cgamma**2 + 2*calpha*cbeta*cgamma)
    return V

def require_keywords(function):
    r"""
    Decorator which forces all keyword arguments to the function to be
    explicitly named.

    For example:

        >>> @require_keywords
        ... def fn(a, b, c=3): pass
        >>> fn(1, 2, 3)
        Traceback (most recent call last):
        ...
        TypeError: name=value required for c
        >>> fn(1, 2, c=6)
        >>> fn(b=1, a=2, c=6)

    Variable arguments are not currently supported:

        >>> @require_keywords
        ... def fn(a, b, c=6, *args, **kw): pass
        Traceback (most recent call last):
        ...
        NotImplementedError: only named arguments for now

    .. Note:: The call signature is not preserved.

    We can't preserve the function signature for the call since the only
    way we can count the number of non-keyword arguments is to
    use the \*args, \*\*kw call style.  Python 3+ provides the '\*' call
    signature element which will force all keywords after '\*' to be named.
    """
    import functools
    try:
        from inspect import signature
        getargspec = _getargspec_from_signature
    except ImportError: # CRUFT: py 2.7 support
        from inspect import getargspec

    args, vararg, varkwd, defaults = getargspec(function)
    if defaults is None:
        defaults = []
    named_args = args[:-len(defaults)]
    named_kwds = args[-len(defaults):]
    # Keep it simple for now
    if vararg or varkwd:
        raise NotImplementedError("only named arguments for now")
    @functools.wraps(function)
    def _require_kwds(*args, **kw):
        if len(args) > len(named_args):
            raise TypeError("name=value required for "+", ".join(named_kwds))
        return function(*args, **kw)
    return _require_kwds

def _getargspec_from_signature(function):
    """
    Reproduce getargspec() interface using newer signature protocol
    """
    from inspect import signature

    args, vararg, varkwd, defaults = [], None, None, []
    sig = signature(function)
    for p in sig.parameters.values():
        args.append(p.name)
        if p.default is not p.empty:
            defaults.append(p.default)
        if p.kind is p.VAR_POSITIONAL:
            vararg = p.name
        elif p.kind is p.VAR_KEYWORD:
            varkwd = p.name
    return args, vararg, varkwd, defaults
