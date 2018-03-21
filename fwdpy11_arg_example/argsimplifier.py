import numpy as np
import msprime
import time
import itertools
from collections import namedtuple

InitMeta = namedtuple('InitMeta', 'position origin_generation origin')

class ArgSimplifier(object):
    """
    Python class to interface between an
    AncestryTracker and msprime
    """

    from .wfarg import reverse_time

    def __init__(self, gc_interval, trees=None):
        """
        :param gc_interval: Garbage collection interval
        :param trees: An instance of :class:`msprime.TreeSequence`
        """
        self.gc_interval = gc_interval
        self.last_gc_time = 0.0
        self.__nodes = msprime.NodeTable()
        self.__edges = msprime.EdgeTable()
        self.__sites = msprime.SiteTable()
        self.__mutations = msprime.MutationTable()
        self.__process = True
        if trees is not None:
            self.__process = False
            trees.dump_tables(nodes=self.__nodes, edges=self.__edges, sites = self.__sites, mutations = self.__mutations)
            if(len(self.__mutations.metadata) == 0):
            	meta_list = [InitMeta(self.__sites[mut[0]][0], self.__nodes[mut[1]][1], "initial tree") for mut in self.__mutations]
            	encoded, offset = msprime.pack_bytes(list(map(pickle.dumps, meta_list)))
            	self.mutations.set_columns(site=self.__mutations.site, node=self.__mutations.node, derived_state=self.__mutations.derived_state, derived_state_offset=self.__mutations.derived_state_offset, parent=self.__mutations.parent, metadata_offset=offset, metadata=encoded)
        self.__time_sorting = 0.0
        self.__time_appending = 0.0
        self.__time_simplifying = 0.0
        self.__time_prepping = 0.0

    def simplify(self, pop, ancestry):
        # print(type(ancestry))
        generation = pop.generation
        # update node times:
        if self.__nodes.num_rows > 0:
            tc = self.__nodes.time
            dt = float(generation) - self.last_gc_time
            tc += dt
            self.last_gc_time = generation
            flags = np.ones(self.__nodes.num_rows, dtype=np.uint32)
            self.__nodes.set_columns(
                flags=flags, population=self.__nodes.population, time=tc)

        before = time.process_time()
        # Acquire mutex
        ancestry.acquire()
        self.reverse_time(ancestry.nodes)
        ana = np.array(ancestry.nodes, copy=False)
        aea = np.array(ancestry.edges, copy=False)
        ama = np.array(ancestry.mutations, copy=False)
        pma = np.array(pop.mutations) #must be copy
        asa = np.array(ancestry.samples, copy=False)
        flags = np.ones(len(ana), dtype=np.uint32)
        
        self.__time_prepping += time.process_time() - before

        before = time.process_time()
        clen = len(self.__nodes)
        self.__nodes.append_columns(flags=flags,
                                    population=ana['population'],
                                    time=ana['generation'])

        before = time.process_time()
        self.__edges.append_columns(left=aea['left'],
                                    right=aea['right'],
                                    parent=aea['parent'],
                                    child=aea['child'])
        if(len(ama) > 0):
           before = time.process_time()
           self.__sites.append_columns(pma['pos'][ama['mutation_id']],
                                  ancestral_state=np.zeros(len(ama), np.int8) + ord('0'),
                                  ancestral_state_offset=np.arange(len(ama) + 1, dtype=np.uint32))

           before = time.process_time()
           encoded, offset = msprime.pack_bytes(list(map(pickle.dumps,pma[ama['mutation_id']])))
           self.__mutations.append_columns(site=np.arange(len(ama), dtype=np.int32) + self.__mutations.num_rows,
                                      node=ama['node_id'],
                                      derived_state=np.ones(len(ama), np.int8) + ord('0'),
                                      derived_state_offset=np.arange(len(ama) + 1, dtype=np.uint32),
                                      metadata_offset=offset, metadata=encoded)
                                      
        msprime.sort_tables(nodes=self.__nodes, edges=self.__edges, sites=self.__sites, mutations=self.__mutations)
        self.__time_sorting += time.process_time() - before
        before = time.process_time()
        sample_map = msprime.simplify_tables(samples=asa.tolist(),
                                             nodes=self.__nodes, edges=self.__edges, sites=self.__sites, mutations=self.__mutations)
        for i in asa:
            assert(sample_map[i] != -1)
        # Release any locks on the ancestry object
        ancestry.release()
        self.__last_edge_start = len(self.__edges)
        self.__time_simplifying += time.process_time() - before
        self.__process = True
        return (True, self.__nodes.num_rows)

    def __call__(self, pop, ancestry, override):
        """
        This is called from C++ during a simulation.

        :param pop: An instance of SlocusPop
        :param ancestry: An instance of AncestryTracker
        :param override: override the gc interval and forces simplification if there are nodes/edges to simplify

        :rtype: tuple

        :returns: A bool and an int
        """
        if len(ancestry.nodes) > 0 and len(ancestry.edges) > 0:
            if pop.generation > 0 and (pop.generation % self.gc_interval == 0.0 or override):
                return self.simplify(pop, ancestry)
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
    def edges(self):
        """
        A NumPy record array representing the edge sets.
        """
        return self.__edges

    @property
    def sites(self):
        """
        A NumPy record array representing the sites.
        """
        return self.__sites

    @property
    def mutations(self):
        """
        A NumPy record array representing the mutations.
        """
        return self.__mutations

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
