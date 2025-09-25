

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import glob


#usage
#python kallistoanalysis.py kallistofilepath pathtooutputdirectory(with/) basenameforfilesaving
#3 outputs - csv of all kallisto ordered with cumulativeercents, csv with top20 , graph with two y axis
#
#define file path
Output = sys.argv[1] # first argument is the output
basename = sys.argv[2]  # second argument is the basenmae
File = sys.argv[3]



#read in the raw kallisto output tsv file
Raw = pd.read_csv(File, sep='\t')

#sort from highest to lowest tpm
Sorted = Raw.sort_values("tpm", ascending = False)
#tpm sum lol i know this should technically be a million but just to check
tpm_sum = Sorted['tpm'].sum()


#adding percentage colum and cumulative frequency colum
Sorted["percent"] = (Sorted["tpm"]/tpm_sum)*100
Sorted["cumulativepercent"] = Sorted["percent"].cumsum()

#Export top top 20 expressed transcripts into a new csv and whole sorted file
top20 = Sorted.head(20)
Sorted.to_csv(f'{Output}{basename}_all.csv') #save to intermediate outputs
top20.to_csv(f'{Output}{basename}_top20.csv')


#visualise all kallisto sorted data
#still working on this, it ends up being too big maybe just do top500?
top500 = Sorted.head(500)
plt.figure(figsize=(15,8))
ax = sns.barplot(top500, x="target_id", y="tpm")
plt.tick_params(
    axis='x',
    which='both',
    bottom=False,
    top=False,
    labelbottom=False)
plt.savefig(f'{Output}{basename}_top500graph.png')
#visualise kallisto data top 20
plt.figure(figsize=(15,8))
ax = sns.barplot(top20, x="target_id", y="tpm")
plt.tick_params(
    axis='x',
    which='both',
    bottom=False,
    top=False,
    labelbottom=False)
ax2 = ax.twinx()
sns.lineplot(top20, x="target_id", y="cumulativepercent", color='red', lw=3)
plt.savefig(f'{Output}{basename}_top20graph.png')

