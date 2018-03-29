import numpy as np
import msprime
import time
import itertools
import pickle
from collections import namedtuple
from .wfarg import AncestryTracker

InitMeta = namedtuple('InitMeta', 'position origin_generation origin')

class ArgSimplifier(object):
    """
    Python class to interface between an
    AncestryTracker and msprime
    """

    def __init__(self, rng, gc_interval, pop, params, trees=None):
        """
        :param rng: Random Number Generator
        :param gc_interval: Garbage collection interval
        :param pop: An instance of :class:`fwdpy11:SLpop`
        :param params: An instance of :class:`fwdpy11:SLpop`
        :param trees: An instance of :class:`fwdpy11.model_params.SlocusParams`
        """
        self.__gc_interval = gc_interval
        self.__nodes = msprime.NodeTable()
        self.__edges = msprime.EdgeTable()
        self.__sites = msprime.SiteTable()
        self.__mutations = msprime.MutationTable()
        self.__rng = rng
        self.__pop = pop
        self.__params = params
        total_generations = len(params.demography) 
        
        if trees is not None:
            self.__process = False
            trees.dump_tables(nodes=self.__nodes, edges=self.__edges, sites = self.__sites, mutations = self.__mutations)
            
            if self.__nodes.num_rows > 0: #add simulation time to input trees
               
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
        
        
        self._anc_tracker = AncestryTracker(pop.N, self.__nodes.num_rows,
                                            total_generations)
            
        self.__time_sorting = 0.0
        self.__time_appending = 0.0
        self.__time_simplifying = 0.0
        self.__time_prepping = 0.0
        self.__time_simulating = 0.0
        
        from .wfarg import evolve_singlepop_regions_track_ancestry
        from fwdpy11.internal import makeMutationRegions, makeRecombinationRegions
        mm = makeMutationRegions(self.__params.nregions, self.__params.sregions)
        rm = makeRecombinationRegions(self.__params.recregions)
        
        self.__time_simulating = evolve_singlepop_regions_track_ancestry(self.__rng, self.__pop, self._anc_tracker, self,  
                                                       self.__params.demography,
                                                       self.__params.mutrate_s,
                                                       self.__params.recrate, mm, rm,
                                                       self.__params.gvalue, self.__params.pself)


    def simplify(self):
        # print(type(ancestry))
        generation = self.__pop.generation
        
        before = time.process_time()
        # Acquire mutex
        self._anc_tracker.acquire()
        ana = np.array(self._anc_tracker.nodes, copy=False)
        aea = np.array(self._anc_tracker.edges, copy=False)
        ama = np.array(self._anc_tracker.mutations, copy=False)
        pma = np.array(self.__pop.mutations.array()) #must be copy
        node_indexes = self._anc_tracker.node_indexes
        anc_samples = self._anc_tracker.anc_samples
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
        all_samples = samples + [i for i in anc_samples if i not in samples]
        sample_map = msprime.simplify_tables(samples= all_samples,
                                             nodes=self.__nodes, edges=self.__edges, sites=self.__sites, mutations=self.__mutations)
        for idx, val in enumerate(anc_samples):
            print(anc_samples[idx],sample_map[val],sample_map[self._anc_tracker.anc_samples[0]],self._anc_tracker.anc_samples[idx])
            anc_samples[idx] = sample_map[val]
            print(idx,anc_samples[idx],self._anc_tracker.anc_samples[idx])
        
        #print(all_samples,ancestry.anc_samples[0])
           
        # Release any locks on the ancestry object
        self._anc_tracker.release()
        self.__time_simplifying += time.process_time() - before
        return (True, self.__nodes.num_rows)
        
    def __call__(self, override):
        """
        This is called from C++ during a simulation.
        :param override: override the gc interval and forces simplification if there are nodes/edges to simplify

        :rtype: tuple

        :returns: A bool and an int
        """
        if len(self._anc_tracker.nodes) > 0 and len(self._anc_tracker.edges) > 0:
            if self.__pop.generation > 0 and (self.__pop.generation % self.__gc_interval == 0.0 or override):
                return self.simplify()
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
                'simplifying': self.__time_simplifying,
                'simulating': self.__time_simulating}
