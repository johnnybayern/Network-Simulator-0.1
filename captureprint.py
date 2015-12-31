import sys
import cStringIO

__author__ = "Zhuo"


def print_hello():
    print "hello"


def capture_output(func, *args):
    old_stdout = sys.stdout  # Keep track of the previous value.
    stream = cStringIO.StringIO()
    sys.stdout = stream
    # Below you can do whatever you want, import module1, call functions
    func(*args)
    # print_hello()
    # restore the previous stdout.
    sys.stdout = old_stdout
    variable = stream.getvalue()   # This will get the "hello" string inside the variable
    return variable

# print capture_output(print_hello)
