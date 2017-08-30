import fwdpy11
import fwdpy11.model_params
import fwdpy11.wright_fisher as wf
import fwdpy11.sampling
import fwdpy11.fitness
import numpy as np

import sys
popsize = int(sys.argv[1])
rho = float(sys.argv[2])
theta = float(sys.argv[3])
seed = int(sys.argv[4])
nreps = int(sys.argv[5])


np.random.seed(seed)
seeds = np.random.choice(3000000,nreps,replace=False)

for i,seedi in enumerate(seeds):
    pop = fwdpy11.SlocusPop(popsize)
    recrate = float(rho) / (4.0 * float(popsize))

    mu = 0.0
    mu_n = theta/float(4*popsize)
    pdict = {'rates': (mu_n, mu, recrate),
             'nregions': [fwdpy11.Region(0,1,1)],
             'sregions': [],
             'recregions': [fwdpy11.Region(0, 1, 1)],
             'gvalue': fwdpy11.fitness.SlocusMult(2.0),
             'demography': np.array([popsize] * 10 * popsize, dtype=np.uint32)
             }

    params = fwdpy11.model_params.SlocusParams(**pdict)
    rng = fwdpy11.GSLrng(seedi)

    wf.evolve(rng,pop,params)
    s = fwdpy11.sampling.sample_separate(rng,pop,10)
    print(len(s[0]))
