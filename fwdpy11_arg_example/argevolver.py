import numpy as np
import msprime
import time
import itertools
import pickle
from collections import namedtuple
from .wfarg import AncestryTracker

InitMeta = namedtuple('InitMeta', 'position origin_generation origin')

class ArgEvolver(object):
    """
    Python class to interface between
    forward simulation and msprime
    """

    def __init__(self, rng, gc_interval, pop, params, anc_sampler=None, trees=None):
        """
        :param rng: Random Number Generator
        :param gc_interval: Garbage collection interval
        :param pop: An instance of :class:`fwdpy11:SLpop`
        :param params: An instance of :class:`fwdpy11:SLpop`
        :param trees: An instance of :class:`fwdpy11.model_params.SlocusParams`
        
        runs forward simulation defined by above params
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
        self._sampler = anc_sampler
        self.__anc_samples = []
        self._new_anc_samples = []
            
        self.__time_sorting = 0.0
        self.__time_appending = 0.0
        self.__time_simplifying = 0.0
        self.__time_prepping = 0.0
        self.__time_simulating = 0.0
        
        from .wfarg import evolve_singlepop_regions_track_ancestry
        from fwdpy11.internal import makeMutationRegions, makeRecombinationRegions
        mm = makeMutationRegions(self.__params.nregions, self.__params.sregions)
        rm = makeRecombinationRegions(self.__params.recregions)
        
        self.__time_simulating = evolve_singlepop_regions_track_ancestry(self.__rng, self.__pop, self._anc_tracker, self, self._anc_sampler,  
                                                       self.__params.demography,
                                                       self.__params.mutrate_s,
                                                       self.__params.recrate, mm, rm,
                                                       self.__params.gvalue, self.__params.pself)


    def _anc_sampler(self):
        temp = []
        
        if(self._sampler):
           new_indiv_samples = self._sampler(self.__pop,self.__params)
           
           if(new_indiv_samples.size > 0):
           
              if(not np.issubdtype(new_indiv_samples.dtype, np.integer)):
                  raise RuntimeError("sample dtype must be an integral type")
                  
              sorted_new_indiv_samples = np.sort(new_indiv_samples)
              if(sorted_new_indiv_samples[0] < 0 or sorted_new_indiv_samples[-1] >= self.__pop.N):
    	          raise RuntimeError("ancestral samples out of bounds")
           
              g1 = lambda val: 2*val + self._anc_tracker.node_indexes[0]
              g2 = lambda val: g1(val) + 1
              temp = [g(val) for val in sorted_new_indiv_samples for g in (g1,g2)]
           
              self._new_anc_samples = temp  
                    
              if(self.__pop.generation == 0):
                  self.__anc_samples = self._new_anc_samples
                  self._new_anc_samples = []  
                       
              return True
           return False
        return False
    
    def _simplify(self):
        generation = self.__pop.generation
        
        before = time.process_time()
        # Acquire mutex
        #self._anc_tracker.acquire()
        ana = np.array(self._anc_tracker.nodes, copy=False)
        aea = np.array(self._anc_tracker.edges, copy=False)
        ama = np.array(self._anc_tracker.mutations, copy=False)
        pma = np.array(self.__pop.mutations.array()) #must be copy
        node_indexes = self._anc_tracker.node_indexes
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
        all_samples = samples + self.__anc_samples #since I force GC during an ancestral sample, I can wait to add them and anc_samples won't overlap with samples
        sample_map = msprime.simplify_tables(samples= all_samples,
                                             nodes=self.__nodes, edges=self.__edges, sites=self.__sites, mutations=self.__mutations)
        
        self.__anc_samples = self._new_anc_samples + self.__anc_samples #doesn't need to be in there before because these ancestral samples will be in the current samples list
        self._new_anc_samples = []
        self.__anc_samples = [sample_map[val] for val in self.__anc_samples]    
           
        # Release any locks on the ancestry object
        #self._anc_tracker.release()
        self._anc_tracker.post_process_gc(self.__nodes.num_rows)
        self.__time_simplifying += time.process_time() - before
        
    def __call__(self, override):
        """
        This is called from C++ during a simulation.
        :param override: override the gc interval and forces simplification if there are nodes/edges to simplify

        :rtype: tuple

        :returns: A bool and an int
        """
        if len(self._anc_tracker.nodes) > 0 and len(self._anc_tracker.edges) > 0:
            if self.__pop.generation > 0 and (self.__pop.generation % self.__gc_interval == 0.0 or override):
                self._simplify()

    @property
    def nodes(self):
        """
        An msprime Node Table.
        """
        return self.__nodes

    @property
    def edges(self):
        """
        An msprime Edge Table.
        """
        return self.__edges

    @property
    def sites(self):
        """
        An msprime Site Table.
        """
        return self.__sites

    @property
    def mutations(self):
        """
        An msprime Mutation Table.
        """
        return self.__mutations
        
    @property
    def anc_samples(self):
        """
        A list of ancestral samples taken during the simulation.
        """
        return self.__anc_samples

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
