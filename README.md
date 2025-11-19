
**Pipeline and Documentation still in development**

# Overview
___

## Quick Start [Following from previous lexatrash23-vf-092025 run] 
**1.** Download nextflow.config file from lexatrash-vfa-092025  
**2.** Edit config file by providing appropriate input for each described parameter. config file must be in the same directory that the script is run from  
**3.** (Optional)Edit config file by providing appropriate memory allocation for each described parameters  
**4.** run on command line:  
    ``` nextflow pull lexatrash23/Nf-vfa-pn-092025 ```  
**5.** run on command line (or preferably in an sbatch script) :  
   ``` nextflow run lexatrash23/Nf-vfa-pn-092025 -with-conda -with-dag -with-timeline -with-report -with-trace  ```     
**6.** Copy multiqc_report.html from lexatrash-vf-092025 results folder into lexatrash-vfa-092025 results/htmls to complete html folder
___
## Additional notes 
**_Suggested Dir Structure:_**   
Sample_Name>Analysis>.config  

**_config params info_:**  
sampleURL - parameter should link to page from which the user can search through and download available. This could be a custom dropbox/drive link etc.. For example of a Dataframe and Alignment search Rshiny app click [here](https://github.com/lexatrash23)  
genome_id - Genbank genome ID that will be used to create hyperlinks for Blastn searches  
data - path to results folder from lexatrash-vf-092025 nextflow run.   
input_list - Interproscan entry.list  
input_panther - PANTHER entry names csv   
input_toxindomains - TSV file of proteins with IP domains of interest.     
Settings -  settings parameters can be removed as necessary.  
**_When complete:_**     
Check slurm output file (if using sbatch script) to ensure all tasks were run successfully. If completed sucessfully, work directory can be deleted to clear space.  

**_Rerunning and Rerunning after cancellation:_**  
-If desired -resume flag can be used to resume nextflow script when troubleshooting failed steps to avoid repeating successful steps    
-If slurm job running nextflow pipeline is cancelled prior to completion, and subsequent run fails, work directory may need to be deleted prior to rerunning to ensure proper conda environment installation   
___

Pipeline image : 
![Pipelineimagesimple](Pipeline_figures/venomflowanalysis.png)
![Pipelineimagedetailed](Pipeline_figures/venomflowanalysiswithoutputs.png)
___

### Required Inputs
This Nextflow pipeline requires lexatrash23/Nf-vf-pn-092025 to be successfully run first. The files found in the results folder of  lexatrash23/Nf-vf-pn-092025 will serve as inputs for lexatrash23/Nf-vfa-pn-092025  
Apart from this files, the following input files are required:   
1. De novo Trinity Assembly (fasta)  
2. Interproscan entry.list  
3. PANTHER protein names csv  
4. TSV file of Proteins with Interproscan domains of interest    

### Optional Inputs  
5. Massspec analysis file   

### Provided test files
Test run can be down with the following provided test files:    
1. trinity_test: A subset(200 sequences) of the trinity assembly from a Doryteuthis pealeii Posterior Salivary Gland tissue  
2. Results: Results folder from lexatrash23/Nf-vf-pn-092025 run using test files  
3. entry.list, panther_names.csv, toxindomains.tsv  
4. massspec_analysis file    
Download Test_files folder and specify respective file paths in local config
 file.  
 Due to size limitations, these are only sample fasta and hence the results will be limited, for an example of a complete html report, refer to the example below. 
### Output files
Several Output files are produced by this pipeline. Main output files of interest are listed below:   
HTMLs: 
A.html -> starting menu page, all other html pages can be navigated to from this. For full html webpage map, please refer to the pipeline image above. 

Excel: 
if genome was available: transdf_distinct_blastn and overview_csv_filtered
if no genome was available: transdf_distinct and overview_csv


Document with more detailed information on input/output/scripts can be found [here](https://github.com/lexatrash23)

## Example output 
Example html output can be seen [here](https://lexatrash23.github.io/CephTranscriptomicsGUI/Samples/Doryteuthis_pealeii/DP3_PD/A.html)

## Next Step: Compilation of HTMLs and R shiny search apps [OPTIONAL]
For suggested R shiny script to create Apps for parsing and downloading dataframe files and alignment outputs click [here](https://github.com/lexatrash23)
___
## Citing
This pipeline can be cited as follows:

[in preparation]

Achrak E, Ferd J, Schulman J, Dang T, Krampis K, Holford M. VenomFlow: An Automated Bioinformatic Pipeline for Identification of Disulfide-Rich Peptides from Venom Arsenals. Methods Mol Biol. 2022;2498:89-97. doi: 10.1007/978-1-0716-2313-8_6. PMID: 35727542.

Please also include citations for the individual bioinformatic tools utilized in this pipeline including: [fastqc](https://github.com/s-andrews/FastQC), [multiqc](https://github.com/MultiQC/MultiQC), [seqkit](https://bioinf.shenwei.me/seqkit/), [Kallisto](https://github.com/pachterlab/kallisto), [BLAST](https://support.nlm.nih.gov/kbArticle/?pn=KA-03391), [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki), [BUSCO](https://busco.ezlab.org/busco_userguide.html), [SignalP](https://services.healthtech.dtu.dk/services/SignalP-5.0/)and [Interproscan](https://www.ebi.ac.uk/interpro/about/interproscan/) 
___
## Pipeline development and Contact information

This Nextflow pipeline was developed by [the Holford Lab](https://holfordlab.com/) by [Praveena Naidu](https://github.com/lexatrash23) and edited/reviewed by [Add Name/hyperlinks of individuals or bioinformatics teams that have reviewed this code]. The Holford Lab can be contacted at holfordlab@gmail.com.
___


