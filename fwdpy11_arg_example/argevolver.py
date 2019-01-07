import numpy as np
import msprime
import time
import itertools
import struct
from .wfarg import AncestryTracker

class ArgEvolver(object):
    """
    Python class to interface between
    forward simulation and msprime
    """
        
    def __init__(self, rng, parsed_args, pop, params, anc_sampler=None, trees=None):
        """
        :param rng: Random Number Generator
        :param gc_interval: Garbage collection interval
        :param pop: An instance of :class:`fwdpy11:SLpop`
        :param params: An instance of :class:`fwdpy11:SLpop`
        :param trees: An instance of :class:`fwdpy11.model_params.SlocusParams`
        
        runs forward simulation defined by above params
        """
        self.__gc_interval = parsed_args.gc
        self.__nodes = msprime.NodeTable()
        self.__edges = msprime.EdgeTable()
        self.__sites = msprime.SiteTable()
        self.__mutations = msprime.MutationTable()
        self.__rng = rng
        self.__pop = pop
        self.__params = params
        self.__pop2array = np.array(parsed_args.pop2,dtype=np.uint32)
        self.__migarray = np.array(parsed_args.migration,dtype=np.float32)
        self.__total_generations = len(params.demography) 
        
        if trees is not None:
            self.__process = False
            trees.dump_tables(nodes=self.__nodes, edges=self.__edges, sites = self.__sites, mutations = self.__mutations)
            
            if self.__nodes.num_rows > 0: #add simulation time to input trees
               
               tc = self.__nodes.time
               dt = float(self.__total_generations)
               tc += dt
               flags = np.ones(self.__nodes.num_rows, dtype=np.uint32)
               self.__nodes.set_columns(
                   flags=flags, population=self.__nodes.population, time=tc)
                   
            if(self.__mutations.num_rows > 0 and len(self.__mutations.metadata) == 0): #add default mutation metadata if none present (differentiates these mutations from those generated in simulation)
            	meta_list = np.full(len(self.__mutations),-1,dtype=np.int32)
            	encoded = meta_list.view(np.int8)
            	offset = np.arange(0,4*(len(meta_list)+1),4,dtype=np.uint32)
            	self.mutations.set_columns(site=self.__mutations.site, node=self.__mutations.node, derived_state=self.__mutations.derived_state, derived_state_offset=self.__mutations.derived_state_offset, parent=self.__mutations.parent, metadata_offset=offset, metadata=encoded)
        
        
        self._anc_tracker = AncestryTracker(pop.N, self.__nodes.num_rows, self.__total_generations)
        self._sampler = anc_sampler
        self.__anc_samples = []
        self._gc_anc_samples = []
        self._anc_sampler(False) #sample of generation 0
            
        self.__time_sorting = 0.0
        self.__time_appending = 0.0
        self.__time_simplifying = 0.0
        self.__time_prepping = 0.0
        self.__time_simulating = 0.0
        from .wfarg import evolve_track_ancestry
        
        import fwdpy11.SlocusPop
        from fwdpy11.internal import makeMutationRegions, makeRecombinationRegions
        pneutral = 0
		
        mm = makeMutationRegions(self.__rng, self.__pop, self.__params.nregions,
                                 self.__params.sregions, pneutral)
        rm = makeRecombinationRegions(self.__rng, self.__params.recrate, 
        						 self.__params.recregions)
      	
        self.__time_simulating = evolve_track_ancestry(self.__rng, self.__pop, 
                                                       self._anc_tracker, self,  
                                                       self.__params.demography,
                                                       self.__pop2array,
                                                       self.__migarray, 
                                                       self.__params.mutrate_s, 
                                                       mm, rm)

    def _anc_sampler(self, simplify_generation):
        temp = []
        if(self._sampler and self.__pop.generation <= self.__total_generations):
           pop_size1 = self.__pop.N
           pop_size2 = 0
           if(self.__pop.generation >= self.__pop2array[1] and self.__pop.generation < self.__pop2array[2]):
              pop_size2 = self.__pop2array[0]
           new_indiv_samples = self._sampler(self.__pop.generation, pop_size1, pop_size2, self.__params, self.__total_generations)
           if(new_indiv_samples.size > 0):
              if(not np.issubdtype(new_indiv_samples.dtype, np.integer)):
                  raise RuntimeError("sample dtype must be an integral type")
                  
              sorted_new_indiv_samples = np.sort(new_indiv_samples)
              max = pop_size1 + pop_size2
              if(sorted_new_indiv_samples[0] < 0 or sorted_new_indiv_samples[-1] >= max):
    	          raise RuntimeError("ancestral samples out of bounds")
              
              #note: cannot prevent the recycling of positions from previous simulations
              if(self.__pop.generation > 0):
                  self._anc_tracker.preserve_mutations_sample(sorted_new_indiv_samples,self.__pop)  
                  
              g1 = lambda val: 2*val + self._anc_tracker.node_indexes[0]
              g2 = lambda val: g1(val) + 1
              temp = [g(val) for val in sorted_new_indiv_samples for g in (g1,g2)]
                           
              if(simplify_generation): self._gc_anc_samples = temp
              else: self.__anc_samples = temp + self.__anc_samples
    
    def _simplify(self):
        before = time.process_time()
        generation = self.__pop.generation
        self._anc_tracker.pre_process_gc(self.__pop)
        # Acquire mutex
        #self._anc_tracker.acquire()
        ana = np.array(self._anc_tracker.nodes, copy=False)
        aea = np.array(self._anc_tracker.edges, copy=False)
        ama = np.array(self._anc_tracker.mutations, copy=False)
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
        
        if(ama.size > 0):
            self.__sites.append_columns(ama['pos'],
                                  ancestral_state=np.zeros(len(ama), np.int8) + ord('0'),
                                  ancestral_state_offset=np.arange(len(ama) + 1, dtype=np.uint32))
            self.__mutations.append_columns(site=np.arange(len(ama), dtype=np.int32) + self.__mutations.num_rows,
                                      node=ama['node_id'],
                                      derived_state=np.ones(len(ama), np.int8) + ord('0'),
                                      derived_state_offset=np.arange(len(ama) + 1, dtype=np.uint32))        
        self.__time_appending += time.process_time() - before
        
        before = time.process_time()                              
        msprime.sort_tables(nodes=self.__nodes, edges=self.__edges, sites=self.__sites, mutations=self.__mutations)
        self.__time_sorting += time.process_time() - before
        before = time.process_time()
        
        if(self.__pop.generation < self.__total_generations or len(self._gc_anc_samples) == 0): 
            samples = list(range(node_indexes[0],node_indexes[1]))
        else: 
            samples = self._gc_anc_samples
        
        all_samples = samples + self.__anc_samples #due to sampler behavior, anc_samples won't overlap with samples
        sample_map = msprime.simplify_tables(samples= all_samples,
                                             nodes=self.__nodes, edges=self.__edges, sites=self.__sites, mutations=self.__mutations)
        
        self.__anc_samples = self._gc_anc_samples + self.__anc_samples #doesn't need to be in there before because these ancestral samples will be in the current samples list
        self._gc_anc_samples = []
        self.__anc_samples = [sample_map[val] for val in self.__anc_samples]    
           
        # Release any locks on the ancestry object
        #self._anc_tracker.release()
        self._anc_tracker.post_process_gc(self.__nodes.num_rows)
        self.__time_simplifying += time.process_time() - before
        
    def __call__(self):
        """
        This is called from C++ during a simulation.
        """
        simplification = False
        
        if self.__pop.generation > 0 and (self.__pop.generation % self.__gc_interval == 0.0 or self.__pop.generation == self.__total_generations):
           self._anc_sampler(True)
           simplification = True
           self._simplify()
           
        if(not(simplification)): self._anc_sampler(False) 
        
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
        
    @property
    def total_generations(self):
        """
        Total Number of Generations in Simulation
        """
        return self.__total_generations
        
    @property
    def pop(self):
        """
        ArgEvolver's pop structure
        """
        return self.__pop

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
