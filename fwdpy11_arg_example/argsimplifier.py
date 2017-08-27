import numpy as np
import msprime


class ArgSimplifier(object):
    __gc_interval = None
    __nodes = msprime.NodeTable()
    __edges = msprime.EdgesetTable()

    def __init__(self, gc_interval):
        self.gc_interval = gc_interval

    def __call__(self, generation, ancestry):
        if generation > 0 and generation % self.gc_interval == 0.0:
            na = np.array(ancestry.nodes, copy=False)
            ea = np.array(ancestry.edges, copy=False)
            samples = np.array(ancestry.samples, copy=False)
            return True

        return False

    @property
    def gc_interval(self):
        return self.__gc_interval

    @gc_interval.setter
    def gc_interval(self, value):
        try:
            int(value)
        except:
            raise ValueError("GC interval must be an integer")
        if value <= 0:
            raise ValueError("GC interval must be and integer > 0")
        self.__gc_interval = int(value)
