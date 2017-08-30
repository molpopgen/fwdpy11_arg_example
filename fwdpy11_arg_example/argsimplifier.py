import numpy as np
import msprime


class ArgSimplifier(object):
    __gc_interval = None
    __nodes = msprime.NodeTable()
    __edges = msprime.EdgesetTable()

    def __init__(self, gc_interval):
        self.gc_interval = gc_interval

    def simplify(self,ancestry):
        # update node times:
        if self.__nodes.num_rows > 0:
            tc = self.__nodes.time
            dt = ancestry.offspring_generation + 1
            tc += self.gc_interval
            # print("current node time adjustment: ")
            # print(self.__nodes.time)
            # print(tc)
            flags = np.empty([self.__nodes.num_rows], dtype = np.uint32)
            flags.fill(1)
            self.__nodes.set_columns(flags=flags,population=self.__nodes.population,time=tc)
            # print("current edge info:")
            # print(self.__edges.parent)
            # print(self.__edges.children)
            # print(self.__nodes)

        ancestry.prep_for_gc()
        na = np.array(ancestry.nodes, copy=False)
        ea = np.array(ancestry.edges, copy=False)
        # print(na[:10])
        # print(na[-10:])
        # print(ea[:10])
        # print(ea[-10:])
        # print("min time: ",na['generation'].min(),na['generation'].max())
        samples = np.array(ancestry.samples, copy=False)
        flags=np.empty([len(na)], dtype=np.uint32)
        flags.fill(1)
        # is_sample=np.empty([len(samples)], dtype = flags.dtype)
        # is_sample.fill(1)
        # flags[-len(samples):]=is_sample
        # X=False
        # if self.__nodes.num_rows>0:
        #     X=True
        #     print("before append:")
        #     print(self.__nodes)
        #     print(self.__edges)
        self.__nodes.append_columns(flags=flags,
                population=na['population'],
                time=na['generation'])
        self.__edges.append_columns(left=ea['left'],
                right=ea['right'],
                parent=ea['parent'],
                children=ea['child'],
                children_length=[1]*len(ea))
        # if X is True:
        #     print(self.__nodes)
        #     print(self.__edges)
        #     print("done before append")
        # print(ancestry.offspring_generation)
        msprime.sort_tables(nodes=self.__nodes, edgesets=self.__edges)
        # print(len(self.__nodes.time))
        # print("here")
        # sample_map = {j:i for i,j in enumerate(samples)}
        #print(sample_map)
        # The rv tuple is : (we did GC, the next index to use for nodes,
        x=msprime.load_tables(nodes=self.__nodes, edgesets=self.__edges)
        # if self.__nodes.num_rows > 0:
        #     print("what are we putting in:")
        #     print(self.__nodes)
        #     print(self.__edges)
        #     print("input done")
        x=x.simplify(samples=samples.tolist())
        x.dump_tables(nodes=self.__nodes, edgesets=self.__edges)
        #print(self.__nodes)
        #print("simplified edges:")
        #print(self.__edges)
        #print("simplified nodes:")
        #print(self.__nodes)
        #print(self.__nodes.num_rows)
        #Do we really need this, or is min/max ok:
        # a map of input to output nodes for the last generation
        return (True,self.__nodes.num_rows)# ,sample_map)

    def __call__(self, generation, ancestry):
        if generation > 0 and generation % self.gc_interval == 0.0:
            return self.simplify(ancestry)
        # Keep tuple size constant,
        # for sake of sanity.
        return (False,None,None)
        
    @property
    def nodes(self):
        return self.__nodes

    @property
    def edgesets(self):
        return self.__edges

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
        self.__gc_interval=int(value)
