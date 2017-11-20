import sys
import pandas as pd

dfs=[]
for i in sys.argv[1:]:
    x=pd.read_csv(i,compression='gzip',sep="\t")
    dfs.append(x)

df=pd.concat(dfs)
df.to_csv("detailed_times.txt",sep="\t",index=False)
