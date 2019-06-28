from ancient_genotypes_simulation import *
from ancient_genotypes import *
from numpy import *
import pandas
import pickle
from joblib import Parallel, delayed

def sim_ghost_admixture(anc_pop, i):
	coverage=1
	num_ind=10
	print([coverage, num_ind])
	freq_sim, GT_sim, reads_sim = ancient_sample_test(num_modern=100,anc_pop=0,anc_num=num_ind,Ne0=3000,Ne1=3000,anc_time=5419,split_time_anc=5919,length=500,num_rep=10,coverage=coverage,error=st.expon.rvs(size=num_ind,scale=.05,random_state=i),seed = i)
	freqs_sim, read_list_sim = get_read_dict(freq_sim,reads_sim) 
	params_pop_sim_free = optimize_pop_params_error_parallel(freqs_sim,read_list_sim,num_core=1,detail=0,continuity=False)
	params_pop_sim_continuity = optimize_pop_params_error_parallel(freqs_sim,read_list_sim,num_core=1,detail=0,continuity=True)
	return [params_pop_sim_free,params_pop_sim_continuity]

cov = [1]
results = []
random.seed(42)
for anc_pop in [0,1]:
	print(anc_pop)
	results.append([])
	results[-1].append([])
	results[-1][-1].append(Parallel(n_jobs=50)(delayed(sim_ghost_admixture)(anc_pop,i) for i in random.randint(1000000,size=20)))
	pickle.dump(results,open("results_ghost_admixture_cont_half_4_coverage_1_5_ind_f_0_half_rep.pickle","wb"))
