import numpy as np
import msprime
import time
import itertools
import pickle
from collections import namedtuple

InitMeta = namedtuple('InitMeta', 'position origin_generation origin')

class ArgSimplifier(object):
    """
    Python class to interface between an
    AncestryTracker and msprime
    """

    def __init__(self, gc_interval, params, trees=None):
        """
        :param gc_interval: Garbage collection interval
        :param trees: An instance of :class:`msprime.TreeSequence`
        """
        self.__gc_interval = gc_interval
        self.__nodes = msprime.NodeTable()
        self.__edges = msprime.EdgeTable()
        self.__sites = msprime.SiteTable()
        self.__mutations = msprime.MutationTable()
        
        if trees is not None:
            self.__process = False
            trees.dump_tables(nodes=self.__nodes, edges=self.__edges, sites = self.__sites, mutations = self.__mutations)
            
            if self.__nodes.num_rows > 0: #add simulation time to input trees
               total_generations = len(params.demography) 
               tc = self.__nodes.time
               dt = float(total_generations)
               tc += dt
               flags = np.ones(self.__nodes.num_rows, dtype=np.uint32)
               self.__nodes.set_columns(
                   flags=flags, population=self.__nodes.population, time=tc)
                   
            if(self.__mutations.num_rows > 0 and len(self.__mutations.metadata) == 0): #add default mutation metadata if none present
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
        
        before = time.process_time()
        # Acquire mutex
        ancestry.acquire()
        ana = np.array(ancestry.nodes, copy=False)
        aea = np.array(ancestry.edges, copy=False)
        ama = np.array(ancestry.mutations, copy=False)
        pma = np.array(pop.mutations.array()) #must be copy
        node_indexes = ancestry.node_indexes
        flags = np.ones(len(ana), dtype=np.uint32)       
        self.__time_prepping += time.process_time() - before
        
        before = time.process_time()
        self.__nodes.append_columns(flags=flags,
                                    population=ana['population'],
                                    time=ana['generation'])

        self.__edges.append_columns(left=aea['left'],
                                    right=aea['right'],
                                    parent=aea['parent'],
                                    child=aea['child'])
        if(len(ama) > 0):
           self.__sites.append_columns(pma['pos'][ama['mutation_id']],
                                  ancestral_state=np.zeros(len(ama), np.int8) + ord('0'),
                                  ancestral_state_offset=np.arange(len(ama) + 1, dtype=np.uint32))
           ###encodes full mutation info as metadata in mutation table in order of numpy pop.mutations.array dtype 
           ###unpickled and transformed into a tuple can be used to construct a fwdpy11.Mutation
           ###e.g. fwdpy11.Mutation(tuple(pickle.loads(simplifier.mutations[i].metadata))[:-1])
           ###uses everything but the last element
           encoded, offset = msprime.pack_bytes(list(map(pickle.dumps,pma[ama['mutation_id']])))
           self.__mutations.append_columns(site=np.arange(len(ama), dtype=np.int32) + self.__mutations.num_rows,
                                      node=ama['node_id'],
                                      derived_state=np.ones(len(ama), np.int8) + ord('0'),
                                      derived_state_offset=np.arange(len(ama) + 1, dtype=np.uint32),
                                      metadata_offset=offset, metadata=encoded)        
        self.__time_appending += time.process_time() - before
        
        before = time.process_time()                              
        msprime.sort_tables(nodes=self.__nodes, edges=self.__edges, sites=self.__sites, mutations=self.__mutations)
        self.__time_sorting += time.process_time() - before
        before = time.process_time()
        samples = list(range(node_indexes[0],node_indexes[1]))
        sample_map = msprime.simplify_tables(samples= samples,
                                             nodes=self.__nodes, edges=self.__edges, sites=self.__sites, mutations=self.__mutations)
        for i in samples:
            assert(sample_map[i] != -1)
           
        # Release any locks on the ancestry object
        ancestry.release()
        self.__time_simplifying += time.process_time() - before
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
        if ancestry is not None and len(ancestry.nodes) > 0 and len(ancestry.edges) > 0:
            if pop.generation > 0 and (pop.generation % self.__gc_interval == 0.0 or override):
                return self.simplify(pop, ancestry)
        return (False, self.__nodes.num_rows)

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
    def times(self):
        """
        A dict representing times spent in various
        steps.
        """
        return {'prepping': self.__time_prepping,
                'sorting': self.__time_sorting,
                'appending': self.__time_appending,
                'simplifying': self.__time_simplifying}
