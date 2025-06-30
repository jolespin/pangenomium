#!/usr/bin/env python
import os
def get_basename(x):
    _, fn = os.path.split(x)
    if fn.endswith(".gz"):
        fn = fn[:-3]
    return ".".join(fn.split(".")[:-1])

