"""Generic object pickler and compressor

This module saves and reloads compressed representations of generic Python
objects to and from the disk.
"""

__author__ = "Bill McNeill <billmcn@speakeasy.net>"
__version__ = "1.0"

import sys
import os.path
import pickle
import gzip

def save(object, filename, protocol = -1):
    """Save an object to a compressed disk file.
       Works well with huge objects.
    """
    with gzip.GzipFile(filename, 'wb') as file:
        pickle.dump(object, file, protocol)

def load(filename):
    """Loads a compressed object from disk
    """
    with gzip.GzipFile(filename, 'rb') as file:
        object = pickle.load(file)
    return object

class Object:
    x = 7
    y = "This is an object."

if __name__ == "__main__":
    filename = sys.argv[1]
    if os.path.isfile(filename):
        o = load(filename)
        print(f"Loaded {o}")
    else:
        o = Object()
        save(o, filename)
        print(f"Saved {o}")
