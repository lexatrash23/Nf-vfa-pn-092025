#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import glob


#usage
#python kallistoanalysis.py kallistofilepath pathtooutputdirectory(with/) basenameforfilesaving
#3 outputs - csv of all kallisto ordered with cumulativeercents, csv with top20 , graph with two y axis
#
#define files from command line arguments
basename = sys.argv[1]  # second argument is the basenmae of output file
File = sys.argv[2] #Input abundance tsv file



#read in the raw kallisto output tsv file
Raw = pd.read_csv(File, sep='\t')

#sort from highest to lowest tpm
Sorted = Raw.sort_values("tpm", ascending = False)
#tpm sum should add up to a million

tpm_sum = Sorted['tpm'].sum()


#adding percentage colum and cumulative frequency column
Sorted["percent"] = (Sorted["tpm"]/tpm_sum)*100
Sorted["cumulativepercent"] = Sorted["percent"].cumsum()

#Export top top 20 expressed transcripts into a new csv and whole sorted file
top20 = Sorted.head(20)
Sorted.to_csv(f'{basename}_Kallisto_Trinity_all.csv')
top20.to_csv(f'{basename}_Kallisto_Trinity_top20.csv')


#visualise all kallisto sorted data
top500 = Sorted.head(500)
plt.figure(figsize=(15,8))
ax = sns.barplot(top500, x="target_id", y="tpm")
plt.tick_params(
    axis='x',
    which='both',
    bottom=False,
    top=False,
    labelbottom=False)
plt.savefig(f'{basename}_Kallisto_Trinity_top500graph.png')
#visualise kallisto data top 20 with cumulative percent
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
plt.savefig(f'{basename}_Kallisto_Trinity_top20graph.png')


