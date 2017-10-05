import numpy as np
import msprime
import time


class ArgSimplifier(object):
    """
    Python class to interface between an
    AncestryTracker and msprime
    """

    def __init__(self, gc_interval):
        """
        :param gc_interval: Garbage collection interval
        """
        self.gc_interval = gc_interval
        self.last_gc_time = 0.0
        self.__nodes = msprime.NodeTable()
        self.__edges = msprime.EdgeTable()
        self.__time_sorting = 0.0
        self.__time_appending = 0.0
        self.__time_simplifying = 0.0
        self.__time_prepping = 0.0
        self.__last_edge_start = 0

    def simplify(self, generation, ancestry):
        # update node times:
        if self.__nodes.num_rows > 0:
            tc = self.__nodes.time
            dt = float(generation) - self.last_gc_time
            tc += dt
            self.last_gc_time = generation
            flags = np.empty([self.__nodes.num_rows], dtype=np.uint32)
            flags.fill(1)
            self.__nodes.set_columns(
                flags=flags, population=self.__nodes.population, time=tc)

        start = time.time()
        ancestry.prep_for_gc()
        na = np.array(ancestry.nodes, copy=False)
        ea = np.array(ancestry.edges, copy=False)
        samples = np.array(ancestry.samples, copy=False)
        flags = np.empty([len(na)], dtype=np.uint32)
        flags.fill(1)
        stop = time.time()
        # children_length = np.empty([len(ea)], dtype=np.uint32)
        # children_length.fill(1)
        self.__time_prepping += (stop - start)

        start = time.time()
        self.__nodes.append_columns(flags=flags,
                                    population=na['population'],
                                    time=na['generation'])
        self.__edges.append_columns(left=ea['left'],
                                    right=ea['right'],
                                    parent=ea['parent'],
                                    child=ea['child'])
        # children_length=children_length)
        stop = time.time()
        self.__time_appending += (stop - start)
        start = time.time()
        print("EdgeTable length = ",len(self.__edges))
        print("Sort start is ", self.__last_edge_start)
        msprime.sort_tables(nodes=self.__nodes,
                            edges=self.__edges,
                            edge_start=self.__last_edge_start)
        stop = time.time()
        self.__time_sorting += (stop - start)
        start = time.time()
        msprime.simplify_tables(samples=samples.tolist(), 
                nodes=self.__nodes, edges=self.__edges)
        self.__last_edge_start = len(self.__edges)
        print("EdgeTable length after simplify is ",self.__last_edge_start)
        stop = time.time()
        self.__time_simplifying += (stop - start)
        return (True, self.__nodes.num_rows)

    def __call__(self, generation, ancestry):
        """
        This is called from C++ during a simulation.

        :param generation: Current generation in a simulation.
        :param ancestry: An instance of AncestryTracker

        :rtype: tuple

        :returns: A bool and an int
        """
        if generation > 0 and generation % self.gc_interval == 0.0:
            return self.simplify(generation, ancestry)
        # Keep tuple size constant,
        # for sake of sanity.
        return (False, None)

    @property
    def nodes(self):
        """
        A NumPy record array representing the nodes.
        """
        return self.__nodes

    @property
    def edgesets(self):
        """
        A NumPy record array representing the edge sets.
        """
        return self.__edges

    @property
    def gc_interval(self):
        """
        The GC interval
        """
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

    @property
    def last_gc_time(self):
        """
        The last time GC was performed
        """
        return self.__last_gc_time

    @last_gc_time.setter
    def last_gc_time(self, value):
        self.__last_gc_time = float(value)

    @property
    def times(self):
        """
        A dict representing times spent in various
        steps.
        """
        return {'prepping': self.__time_prepping,
                'sorting': self.__time_sorting,
                'appending': self.__time_appending,
                'simplifying': self.__time_simplifying}
