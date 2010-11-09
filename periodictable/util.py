"""
Helper functions
"""
import inspect
import functools

def require_keywords(function):
    """
    Decorator which forces all keyword arguments to the function to be 
    explicitly named.
    
    :Example:

        >>> @require_keywords
        ... def fn(a,b,c=3): pass
        >>> fn(1,2,3)
        Traceback (most recent call last):
        ...
        TypeError: name=value required for c
        >>> fn(1,2,c=6)
        >>> fn(b=1,a=2,c=6)
    
    :Note: The call signature is not preserved.
    
    We can't preserve the function signature for the call since the only 
    way we can count the number of non-keyword arguments is to 
    use the *args, **kw call style.  Python 3+ provides the '*' call 
    signature element which will force all keywords after '*' to be named.
    """
    argspec = inspect.getargspec(function)
    named_args = argspec.args[:-len(argspec.defaults)]
    named_kwds = argspec.args[-len(argspec.defaults):]
    vararg = argspec.varargs
    varkwd = argspec.keywords
    # Keep it simple for now
    if vararg or varkwd:
        raise NotImplementError("only named arguments for now")
    @functools.wraps(function)
    def _require_kwds(*args, **kw):
        if len(args) > len(named_args):
            raise TypeError("name=value required for "+", ".join(named_kwds))
        return function(*args, **kw)
    return _require_kwds
