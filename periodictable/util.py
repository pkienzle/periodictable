# This program is in the public domain
# Author: Paul Kienzle
"""
Helper functions
"""

def cell_volume(a=None,b=None,c=None,alpha=None,beta=None,gamma=None):
    """
    Compute cell volume from lattice parameters.

    :Parameters:
        *a*, *b*, *c* : float | |A|
            Lattice spacings.  If *b* or *c* are missing they default to *a*.
        *alpha*, *beta*, *gamma* : float | |deg|
            Lattice angles.  If any are missing they default to 90\ |deg|

    :Returns:
        *V* : float |A^3|
            Cell volume

    :Raises:
        *TypeError* : missing or invalid parameters

    The following formula works for all lattice types:

    .. math::

        V = a b c \sqrt{1 - \cos^2 \alpha - \cos^2 \beta - cos^2 \gamma
                          + 2 \cos \alpha \cos beta \cos gamma}
    """
    from math import cos, radians, sqrt
    if a is None: raise TypeError('missing lattice parameters')
    if b is None: b = a
    if c is None: c = a
    calpha, cbeta, cgamma = [cos(radians(v)) if v is not None else 0
                             for v in alpha, beta, gamma]
    V = a*b*c*sqrt(1 - calpha**2 - cbeta**2 - cgamma**2 + 2*calpha*cbeta*cgamma)
    return V

def require_keywords(function):
    """
    Decorator which forces all keyword arguments to the function to be
    explicitly named.

    For example:

        >>> @require_keywords
        ... def fn(a,b,c=3): pass
        >>> fn(1,2,3)
        Traceback (most recent call last):
        ...
        TypeError: name=value required for c
        >>> fn(1,2,c=6)
        >>> fn(b=1,a=2,c=6)

    Variable arguments are not currently supported:

        >>> @require_keywords
        ... def fn(a,b,c=6,*args,**kw): pass
        Traceback (most recent call last):
        ...
        NotImplementedError: only named arguments for now

    .. Note:: The call signature is not preserved.

    We can't preserve the function signature for the call since the only
    way we can count the number of non-keyword arguments is to
    use the *args, **kw call style.  Python 3+ provides the '*' call
    signature element which will force all keywords after '*' to be named.
    """
    import inspect
    import functools

    args, vararg, varkwd, defaults = inspect.getargspec(function)
    if defaults is None: defaults = []
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


#From werkzeug.utils
#copyright: (c) 2010 by the Werkzeug Team, see AUTHORS for more details.
#license: BSD, see werkzeug LICENSE for more details.
class _Missing(object):
    def __repr__(self):
        return 'no value'
    def __reduce__(self):
        return '_missing'
_missing = _Missing()
class cached_property(object):
    """A decorator that converts a function into a lazy property. The
function wrapped is called the first time to retrieve the result
and then that calculated result is used the next time you access
the value::

class Foo(object):

@cached_property
def foo(self):
# calculate something important here
return 42

The class has to have a `__dict__` in order for this property to
work.
"""

    # implementation detail: this property is implemented as non-data
    # descriptor. non-data descriptors are only invoked if there is
    # no entry with the same name in the instance's __dict__.
    # this allows us to completely get rid of the access function call
    # overhead. If one choses to invoke __get__ by hand the property
    # will still work as expected because the lookup logic is replicated
    # in __get__ for manual invocation.

    def __init__(self, func, name=None, doc=None, writeable=False):
        if writeable:
            from warnings import warn
            warn(DeprecationWarning('the writeable argument to the '
                                    'cached property is a noop since 0.6 '
                                    'because the property is writeable '
                                    'by default for performance reasons'))

        self.__name__ = name or func.__name__
        self.__module__ = func.__module__
        self.__doc__ = doc or func.__doc__
        self.func = func

    def __get__(self, obj, type=None):
        if obj is None:
            return self
        value = obj.__dict__.get(self.__name__, _missing)
        if value is _missing:
            value = self.func(obj)
            obj.__dict__[self.__name__] = value
        return value
