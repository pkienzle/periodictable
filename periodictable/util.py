# This program is in the public domain
# Author: Paul Kienzle
"""
Helper functions
"""
from math import sqrt

def parse_uncertainty(s):
    """
    Given a floating point value plus uncertainty return the pair (val, unc).

    Format is val, val(unc), [nominal] or [low,high].

    The val(unc) form is like 23.0035(12), but also 23(1), 23.0(1.0), or 
    maybe even 23(1.0). This parser does not handle exponential notation
    such as 1.032(4)E10

    The nominal form has zero uncertainty, as does a bare value.

    The [low,high] form is assumed to be a rectangular distribution of 1-sigma
    equivalent width (high-low)/sqrt(12).

    An empty string is returned as None,None rather than 0,inf.
    """
    if s == "": # missing
        # TODO: maybe 0 +/- inf ?
        return None, None

    # Parse [nominal] or [low,high]
    if s.startswith('['):
        s = s[1:-1]
        parts = s.split(',')
        if len(parts) > 1:
            low, high = float(parts[0]), float(parts[1])
            # Use equivalent 1-sigma width for a rectangular distribution
            return (high+low)/2, (high-low)/sqrt(12)
        else:
            return float(parts[0]), 0

    # Parse value(unc) with perhaps '#' at the end
    parts = s.split('(')
    if len(parts) > 1:
        # Split the value and uncertainty.
        value, unc = parts[0], parts[1].split(')')[0]
        # Count digits after the decimal for value and produce
        # 0.00...0{unc} with the right number of zeros.
        # e.g., 23.0035(12) but not 23(1) or 23.0(1.0) or 23(1.0)
        if '.' not in unc and '.' in value:
            zeros = len(value.split('.')[1]) - len(unc)
            unc = "0." + ("0"*zeros) + unc
        return float(value), float(unc)

    # Plain value with no uncertainty
    return float(s), 0

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
