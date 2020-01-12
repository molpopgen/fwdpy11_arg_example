import msprime as msp
import scipy.stats as st
import numpy as np

def ancient_sample_test(num_modern=1000, anc_pop = 0, anc_num = 1, anc_time = 200, split_time_anc = 400, Ne0 = 10000, Ne1 = 10000, mu = 1.25e-8, length = 1000, num_rep = 1000, error = None, coverage = False, seed = None):
	if error is None:
		error = np.zeros(anc_num)
	samples = [msp.Sample(population = 0, time = 0)]*num_modern
	samples.extend([msp.Sample(population = anc_pop, time = anc_time)]*(2*anc_num))
	pop_config = [msp.PopulationConfiguration(initial_size = Ne0), msp.PopulationConfiguration(initial_size = Ne1)]
	divergence = [msp.MassMigration(time = split_time_anc, source = 1, destination = 0, proportion = 1.0)]
	sims = msp.simulate(samples=samples,Ne=Ne0,population_configurations=pop_config,demographic_events=divergence,mutation_rate=mu,length=length,num_replicates=num_rep, random_seed = seed)
	freq = []
	reads = []
	GT = []
	sim_num = 0
	poly_count = 0
	for sim in sims:
		poly_count += sim.num_mutations
		for variant in sim.variants():
			var_array = variant.genotypes
			cur_freq = sum(var_array[:-(2*anc_num)])/float(num_modern)
			if cur_freq == 0 or cur_freq == 1: continue
			freq.append(cur_freq)
			reads.append([])
			GT.append([])
			for i in range(anc_num):
				if i == 0: cur_GT = var_array[-2:]
				else: cur_GT = var_array[-(2*(i+1)):-(2*i)]
				cur_GT = sum(cur_GT)
				GT[-1].append(cur_GT)
				reads[-1].append([None,None])
				if coverage:
					num_reads = st.poisson.rvs(coverage)
					p_der = cur_GT/2.*(1-error[i])+(1-cur_GT/2.)*error[i]
					derived_reads = st.binom.rvs(num_reads, p_der)
					reads[-1][-1] = (num_reads-derived_reads,derived_reads)
	print(poly_count)
	return np.array(freq), GT, reads


if __name__ == "__main__":
	np.random.seed(42)
	num_ind = 10
	i = 6803976
	freq_sim, GT_sim, reads_sim = ancient_sample_test(num_modern=100,anc_pop=0,anc_num=num_ind,Ne0=3000,Ne1=3000,anc_time=5419,split_time_anc=5919,length=500,num_rep=13333,coverage=1,error=st.expon.rvs(size=num_ind,scale=.05,random_state=i),seed = i)