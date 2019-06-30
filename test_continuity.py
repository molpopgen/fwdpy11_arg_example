from ancient_genotypes_simulation import *
from ancient_genotypes import *
from numpy import *
import concurrent.futures

def sim(tuple):
	anc_pop_id = tuple[0] 
	i = tuple[1]
	coverage=1
	num_ind=10
	print([coverage, num_ind])
	freq_sim, GT_sim, reads_sim = ancient_sample_test(num_modern=100,anc_pop=anc_pop_id,anc_num=num_ind,Ne0=3000,Ne1=3000,anc_time=5419,split_time_anc=5919,length=500,num_rep=13333,coverage=coverage,error=st.expon.rvs(size=num_ind,scale=.05,random_state=i),seed = i)
	freqs_sim, read_list_sim = get_read_dict(freq_sim,reads_sim) 
	params_pop_sim_free = optimize_pop_params_error_parallel(freqs_sim,read_list_sim,num_core=1,detail=0,continuity=False)
	params_pop_sim_continuity = optimize_pop_params_error_parallel(freqs_sim,read_list_sim,num_core=1,detail=0,continuity=True)
	return (params_pop_sim_continuity[0][1], params_pop_sim_free[0][1], params_pop_sim_continuity[0][0][0], params_pop_sim_free[0][0][0], params_pop_sim_continuity[0][0][1], params_pop_sim_free[0][0][1])

results = []
random.seed(42)
num_replicates = 80
list = [i for i in random.choice(range(1000000),size=num_replicates*2,replace=False)]

for anc_pop in [0,1]:
	print(anc_pop)
	list2 = [list[i] for i in range(num_replicates*anc_pop,num_replicates*(anc_pop+1))]
	with concurrent.futures.ProcessPoolExecutor() as pool:
		futures = {pool.submit(sim, (anc_pop,i)) for i in list2}
		for fut in concurrent.futures.as_completed(futures):
			results.append(fut.result())

file = open("results.txt","w")
file.write("(0,1)\t(0,1)\t(0,1)\t(0,1)\t(0,1)\t(0,1)\t(0,2)\t(0,2)\t(0,2)\t(0,2)\t(0,2)\t(0,2)\n")
file.write("Continuity_L\tFree_L\tContinuity_t1\tFree_t1\tContinuity_t2\tFree_t2\tContinuity_L\tFree_L\tContinuity_t1\tFree_t1\tContinuity_t2\tFree_t2\n")
for i in range(num_replicates):
	file.write(str(results[i][0]) + "\t" + str(results[i][1]) + "\t" + str(results[i][2]) + "\t" + str(results[i][3]) + "\t" + str(results[i][4]) + "\t" + str(results[i][5]) + "\t" + str(results[i+num_replicates][0]) + "\t" + str(results[i+num_replicates][1]) + "\t" + str(results[i+num_replicates][2]) + "\t" + str(results[i+num_replicates][3]) + "\t" + str(results[i+num_replicates][4]) + "\t" + str(results[i+num_replicates][5]) + "\n")

file.close()