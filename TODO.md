# To get msprime integration working:

The C++ side can/should do more to help

* Get a list of sample indexes
* Convert forward time to backward time
* Sort nodes and edges
* Prune edges that where a child ID is never used later as a parent id
* Prune extinct nodes

During edge pruning, we can build a list of non-extinct nodes that we use
to actually prune the node list.  We can used unorderd_set for a constant-time
lookup table.

We will need to book-keep things like:

* The last generation where GC happened
* The last population size, as that relates to "len(samples"), and is the next ID value to use

We can probably get rid of:

* Tracking vectors of parent/offspring IDs. As I've moved to a more "naive" method of node insertion, all I really need is the first ID of the parental generation and its population size, etc.  Likewise for offspring.

We may want to consider:

* Some rudimentary pruning each generation.

Safety stuff:

* We are now using signed integers for some fields.  For "big" sims (large N/rho), we can overflow things pretty quick.  I need to add a check that throws an exception when we overflow.  I'll be strict about this and simply throw the exception.  However, a safer/saner implementation for production code would be to check if the next generation will overflow and do a simplification if so.   The interesting edge case is if the population size and/or 4Nr is so large that we should overflow mid-generation.  That is a nasty thing to ponder...
