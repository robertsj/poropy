# Taken from the python documentation and adapted.
class Error(Exception):
    """Base class for exceptions in this package."""
    def __init__(self,msg):
        self.expr = expr
        self.msg = msg

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg

class TransitionError(Error):
    """Raised when an operation attempts a state transition that's not
    allowed.

    Attributes:
        prev -- state at beginning of transition
        next -- attempted new state
        msg  -- explanation of why the specific transition is not allowed
    """

    def __init__(self, prev, nxt, msg):
        self.prev = prev
        self.next = nxt
        self.msg = msg
