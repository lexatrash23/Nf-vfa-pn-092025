#!/usr/bin/env nextflow

// Process 1: For kallistoanalysistrinity.py python,pandas,seaborn,matplotlib
process kallistoAnalysisTrinity {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4



    label 'process_low'



    conda 'python=3.8 pandas seaborn matplotlib'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Kallisto/Transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(kallisto_file_trinity)

    output:
    tuple val(sample), path("*_all.csv"), emit: trin_all_csv
    path ("*_top20.csv"), emit: trin_top20_csv
    tuple val(sample), path("*_top500graph.png"), emit: trin_top500_png
    tuple val(sample), path("*_top20graph.png"), emit: trin_top20_png

    script:
    """
    python3 ${workflow.projectDir}/bin/Intermediate_Scripts/kallistoanalysistrinity.py ${sample} ${kallisto_file_trinity}
    """
}

// Process 2: For kallistoanalysistrans.py dependencies:python,pandas,seaborn,matplotlib
process kallistoAnalysisTrans {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'



    conda 'python=3.8 pandas seaborn matplotlib'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Kallisto/ORFs/", mode: 'copy'

    input:
    tuple val(sample), path(kallisto_file_transdecoder)

    output:
    tuple val(sample), path("*_all.csv"), emit: trans_all_csv
    path ("*_top20.csv"), emit: trans_top20_csv
    tuple val(sample), path("*_top500graph.png"), emit: trans_top500_png
    tuple val(sample), path("*_top20graph.png"), emit: trans_top20_png

    script:
    """
    python3 ${workflow.projectDir}/bin/Intermediate_Scripts/kallistoanalysistrans.py "${sample}" ${kallisto_file_transdecoder}
    """
}

// Process 3: Extract Signal Sequences dependencies python biopython
process ExtractSignalSequences {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'



    conda 'python=3.8 biopython'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/SignalSequences/", mode: 'copy'

    input:
    tuple val(sample), path(secreted_pep), path(mature_fasta)

    output:
    tuple val(sample), path("*signalsequences.fasta"), emit: signalsequences

    script:
    """
    python3 ${workflow.projectDir}/bin/Intermediate_Scripts/IS2.py ${secreted_pep} ${mature_fasta} "${sample}_signalsequences.fasta"
    """
}

// Process 4: Create Trinity Dataframe dependecies : R, biocmanager 
process CreateTrinityDataframe {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'



    conda 'r-base bioconductor-biostrings r-tidyr r-dplyr'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Dataframes/Transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), path(blastx_file), path(kallisto_csv)

    output:
    tuple val(sample), path("*TBK.csv"), emit: TBK
    path ("*TBK_distinct.csv.gz"), emit: TBK_distinct

    script:
    """
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS1.R ${trinity_fasta} ${blastx_file} ${kallisto_csv} ${sample}
    """
}

// Process 5: Create Interproscan Dataframe dependecies : R, biocmanager 
process CreateInterproscanDataframe {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base bioconductor-biostrings r-dplyr bioconductor-biostrings bioconductor-biomart r-tidyr'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Dataframes/ORFs/", mode: 'copy'

    input:
    tuple val(sample), path(Interproscan), path(ListFile), path(PantherFile)

    output:
    tuple val(sample), path("*Final_interproscan_dataframe.csv"), emit: Interproscan_dataframe

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS4.R" ${Interproscan} ${ListFile} ${PantherFile} ${sample}
    """
}

// Process 6: Create Transdecoder Dataframe dependecies : R, biocmanager 
process CreateTransdecoderDataframe {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-dplyr r-tidyr bioconductor-biostrings r-stringr'


    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Dataframes/ORFs/", pattern: "*.csv", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Fastas/ORFs/", pattern: "*.fasta", mode: 'copy'

    input:
    tuple val(sample), path(transdecoder_pep), path(transdecoder_cds), path(blastp_file), path(mature_fasta), path(signalsequences), path(Interproscan_dataframe), path(kallistotrans)

    output:
    tuple val(sample), path("*transdf.csv"), emit: transdf
    tuple val(sample), path("*transdf_distinct.csv"), emit: transdf_distinct
    tuple val(sample), path("*_secreted_proteins.pep"), emit: secretedpep
    tuple val(sample), path("*_secreted_proteins.cds"), emit: secretedcds

    script:
    """
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS5.R ${transdecoder_pep} ${transdecoder_cds} ${blastp_file} ${mature_fasta} ${signalsequences} ${Interproscan_dataframe} ${kallistotrans} ${sample}
    """
}



// Process 8: Create BUSCOgraphtranscriptome  
process BUSCOtranscriptome {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/busco/transcriptome/${count}/", mode: 'copy'

    input:
    tuple val(sample), path(buscodirtr), val(count)

    output:
    tuple val(sample), path("*.png"), emit: busco_transcriptome

    script:
    """
    python3 "${workflow.projectDir}/bin/Intermediate_Scripts/IS6.py" -wd "${buscodirtr}" -o ${count}.buscofiguretranscriptome.png
    """
}

// Process 9: Create BUSCOgraphtranslatome 
process BUSCOtranslatome {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/busco/translatome/${count}/", mode: 'copy'

    input:
    tuple val(sample), path(buscodirtl), val(count)

    output:
    tuple val(sample), path("*.png"), emit: busco_translatome

    script:
    """
    python3 "${workflow.projectDir}/bin/Intermediate_Scripts/IS6.py" -wd "${buscodirtl}" -o ${count}.buscofigure.translatome.png
    """
}

// Process 10: Create FigureGenerationTrinity
process FigureGenerationTrinity {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-dplyr r-tidyr r-ggplot2 r-ggalluvial bioconductor-biostrings r-ggrepel r-cowplot r-stringr r-forcats'
    container 'community.wave.seqera.io/library/bioconductor-biostrings_r-base_r-dplyr_r-ggalluvial_pruned:77dba7ba8dae5174'


    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/HTMLFigures/Transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(TBK)

    output:
    tuple val(sample), path("pie1.png"), emit: pie1
    tuple val(sample), path("pie2.png"), emit: pie2
    tuple val(sample), path("pie3.png"), emit: pie3
    tuple val(sample), path("pie4_plot.png"), emit: pie4Plot
    tuple val(sample), path("pie4_legend.png"), emit: pie4Legend
    tuple val(sample), path("alluvial1_plot.png"), emit: alluvial1Plot
    tuple val(sample), path("alluvial1_legend.png"), emit: alluvial1Legend
    tuple val(sample), path("alluvial2_plot.png"), emit: alluvial2Plot
    tuple val(sample), path("alluvial2_legend.png"), emit: alluvial2Legend

    tuple val(sample), path("Table13.csv"), emit: Table13

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Figure_generation_Trinity.R" ${TBK} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"
    """
}

// Process 11: Create FigureGenerationTransdecoder
process FigureGenerationTransdecoder {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-dplyr r-tidyr r-ggplot2 r-ggalluvial bioconductor-biostrings r-ggrepel r-cowplot r-stringr r-forcats'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/HTMLFigures/ORFs/", mode: 'copy'

    input:
    tuple val(sample), path(transdf)

    output:
    tuple val(sample), path("pie5.png"), emit: pie5
    tuple val(sample), path("pie6.png"), emit: pie6
    tuple val(sample), path("pie7.png"), emit: pie7
    tuple val(sample), path("pie8_plot.png"), emit: pie8Plot
    tuple val(sample), path("pie8_legend.png"), emit: pie8Legend
    tuple val(sample), path("alluvial3_plot.png"), emit: alluvial3Plot
    tuple val(sample), path("alluvial3_legend.png"), emit: alluvial3Legend
    tuple val(sample), path("alluvial4_plot.png"), emit: alluvial4Plot
    tuple val(sample), path("alluvial4_legend.png"), emit: alluvial4Legend

    tuple val(sample), path("Table14.csv"), emit: Table14

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Figure_generation_Transdecoder.R" ${transdf} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"
    """
}



// Process 13: Create FigureGenerationSignalp
process FigureGenerationSignalp {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-dplyr r-tidyr r-ggplot2 r-ggalluvial bioconductor-biostrings r-ggrepel r-cowplot r-stringr r-forcats'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/HTMLFigures/Secreted/", mode: 'copy'

    input:
    tuple val(sample), path(transdf)

    output:
    tuple val(sample), path("pie9.png"), emit: pie9
    tuple val(sample), path("pie10.png"), emit: pie10
    tuple val(sample), path("pie11.png"), emit: pie11
    tuple val(sample), path("pie12_plot.png"), emit: pie12Plot
    tuple val(sample), path("pie12_legend.png"), emit: pie12Legend
    tuple val(sample), path("alluvial5_plot.png"), emit: alluvial5Plot
    tuple val(sample), path("alluvial5_legend.png"), emit: alluvial5Legend
    tuple val(sample), path("alluvial6_plot.png"), emit: alluvial6Plot
    tuple val(sample), path("alluvial6_legend.png"), emit: alluvial6Legend
    tuple val(sample), path("Table15.csv"), emit: Table15

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Figure_generation_SignalP.R" ${transdf} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"
    """
}

// Process 12: Create TableGenerationTrinity
process TableGenerationTrinity {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-dplyr r-DT'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/HTMLFigures/Transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(TBK), val(genome_id), val(species)

    output:
    tuple val(sample), path("Table1.csv"), emit: Table1
    tuple val(sample), path("Table2.csv"), emit: Table2
    tuple val(sample), path("Table3.csv"), emit: Table3
    tuple val(sample), path("Table4.csv"), emit: Table4

    script:
    """
	
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Generating_TopTables_Trinity.R" ${TBK} ${genome_id} ${species}
    """
}

// Process 11: Create TableGenerationTransdecoder  
process TableGenerationTransdecoder {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-dplyr r-DT r-stringr'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/HTMLFigures/ORFs/", pattern: "*{5,6,7,8}.csv", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/HTMLFigures/Secreted/", pattern: "*{9,10,11,12}.csv", mode: 'copy'

    input:
    tuple val(sample), path(transdf), val(genome_id), val(species)

    output:
    tuple val(sample), path("Table5.csv"), emit: Table5
    tuple val(sample), path("Table6.csv"), emit: Table6
    tuple val(sample), path("Table7.csv"), emit: Table7
    tuple val(sample), path("Table8.csv"), emit: Table8
    tuple val(sample), path("Table9.csv"), emit: Table9
    tuple val(sample), path("Table10.csv"), emit: Table10
    tuple val(sample), path("Table11.csv"), emit: Table11
    tuple val(sample), path("Table12.csv"), emit: Table12

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Generating_Tables_Transdecoder_SignalP.R" ${transdf} ${genome_id} ${species} 
    """
}


process ToxinVsNonToxin {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'



    conda 'r-base=4.3 bioconductor-go.db r-dplyr r-tidyr r-stringr bioconductor-annotationdbi conda-forge::r-archive r-readr'

    publishDir "CommonIntermediateFiles/Pipelines/Analysis/ToxinVsNonToxinMetaData/", mode: 'copy'

    input:
    tuple path(toxprotblastmetadata), path(nontoxprotmetadata), path(IP_metadata)

    output:
    path "*IP.csv", emit: toxvsnontoxIP
    path "*MF.csv", emit: toxvsnontoxMF
    path "*BP.csv", emit: toxvsnontoxBP

    script:

    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS8.R" ${toxprotblastmetadata} ${nontoxprotmetadata} ${IP_metadata} 
    """
}

//Process Minimap2: For all complete ORFs if genome is available, just to get alignment statistics
process Minimap1 {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_high'




    conda 'minimap2 bioconda::samtools bioconda::stringtie bioconda::bedtools'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Minimap/AllCompleteORFs", mode: 'copy'

    input:
    tuple val(sample), path(genome), path(transdecoder_cds)

    output:
    tuple val(sample), path("*.bam")
    tuple val(sample), path("*.paf"), emit: completecdspaf

    script:
    """
    minimap2 -cx splice:hq -uf ${genome} ${transdecoder_cds} > completecdsalignment.paf
    minimap2 -ax splice:hq -uf ${genome} ${transdecoder_cds} > completecdsminimap.sam
    samtools view -bS completecdsminimap.sam > completecdsminimap.bam
    rm -r *.sam

    """
}

//Process 12: Add Massspec and Blastn6 results where available

process AddMSGenomeIfAvailableAndCreateOverview {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_medium'


    conda 'r-base=4.3 r-dplyr r-ggplot2 r-ggalluvial r-gridbase r-ggvenn bioconductor-genomicranges r-igraph bioconductor-biostrings r-pafr'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/Overview/Dataframes/Unannotated/", pattern: "*.csv", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/Overview/VennDiagrams", pattern: "*.png", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/Overview/Fastas", pattern: "*.pep", mode: 'copy'

    input:
    tuple val(sample), val(species), path(massspec), path(blastn6), path(transdf), path(paf), path(toxvsnontoxIP)

    output:
    tuple val(sample), path("*_Venn_strict.png"), emit: VennPngStrict
    tuple val(sample), path("*filtered_strict.csv"), emit: VennCsvStrict
    // BS > 250;2 KE > 1%;2 Coverage > 50%;2 CR>5%;2 TD OE ;2
    tuple val(sample), path("*_Venn_lax.png"), emit: VennPngLax
    tuple val(sample), path("*filtered_lax.csv"), emit: VennCsvLax
    // BS:50;1 KE > 0%;1 Coverage > 0 ;1 CR: 1%;1  TD Any;1
    path "*.pep"
    //unfiltered
    path "*final_unfiltered.csv"

    script:
    def ms_arg = massspec.name != 'NO_FILE' ? "${massspec}" : "NULL"
    def bn6_arg = blastn6.name != 'NO_FILE' ? "${blastn6}" : "NULL"
    def paf_arg = paf.name != 'NO_FILE' ? "${paf}" : "NULL"
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS15.R" \
        "${sample}" "${species}" ${transdf} ${toxvsnontoxIP} \
        "${bn6_arg}" "${ms_arg}" "${paf_arg}"
    """
}

// Process 7: Create CreateInterproscanFigures
process CreateInterproscanFigures {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-dplyr r-gridbase  r-ggplot2 r-ggrepel r-cowplot r-stringr r-forcats'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/HTMLFigures/Interproscan/", pattern: "*.png", mode: 'copy'

    input:
    tuple val(sample), path(transdf_distinct_csv), path(toxvsnontoxIP), path(toxvsnontoxMF), path(toxvsnontoxBP)

    output:
    tuple val(sample), path("Plot_1_IP.png"), emit: IP_1
    tuple val(sample), path("Legend_1_IP.png"), emit: IP_1_Legend
    tuple val(sample), path("Plot_2_IP.png"), emit: IP_2
    tuple val(sample), path("Legend_2_IP.png"), emit: IP_2_Legend
    tuple val(sample), path("Plot_3_IP.png"), emit: IP_3
    tuple val(sample), path("Legend_3_IP.png"), emit: IP_3_Legend
    tuple val(sample), path("Plot_1_MF.png"), emit: MF_1
    tuple val(sample), path("Legend_1_MF.png"), emit: MF_1_Legend
    tuple val(sample), path("Plot_2_MF.png"), emit: MF_2
    tuple val(sample), path("Legend_2_MF.png"), emit: MF_2_Legend
    tuple val(sample), path("Plot_3_MF.png"), emit: MF_3
    tuple val(sample), path("Legend_3_MF.png"), emit: MF_3_Legend
    tuple val(sample), path("Plot_1_BP.png"), emit: BP_1
    tuple val(sample), path("Legend_1_BP.png"), emit: BP_1_Legend
    tuple val(sample), path("Plot_2_BP.png"), emit: BP_2
    tuple val(sample), path("Legend_2_BP.png"), emit: BP_2_Legend
    tuple val(sample), path("Plot_3_BP.png"), emit: BP_3
    tuple val(sample), path("Legend_3_BP.png"), emit: BP_3_Legend
    tuple val(sample), path("*.csv")

    script:
    """
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS7.R ${transdf_distinct_csv} ${sample} ${toxvsnontoxIP} ${toxvsnontoxMF} ${toxvsnontoxBP} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"
    """
}

// Process 21: RmarkdownB
process RmarkdownB {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-readr r-tidyr r-knitr r-gridExtra r-kableExtra r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(samplesheet)

    output:
    tuple val(sample), path("*.html")

    script:
    """
    samplesheet_abs=\$(readlink -f "${samplesheet}")
   Rscript -e "rmarkdown::render(
  '${workflow.projectDir}/bin/Rmarkdown_scripts/B.Rmd',
  output_dir='.',
  knit_root_dir='.',
  params=list(
    author='${author}',
    sample='${sample}',
    samplesheet='\${samplesheet_abs}'
  )
)"
  

    """
}

process CreateSampleSheet {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'

    input:
    tuple val(sample_name), val(author_name), val(csv_content)

    output:
    tuple val(sample_name), val(author_name), path("${sample_name}_metadata.csv"), emit: samplesheet

    script:
    """
    cat << 'EOF' > ${sample_name}_metadata.csv
    ${csv_content}
    EOF
    """
}



// Process 20: RmarkdownA
process RmarkdownA {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), val(sampleURL)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/A.Rmd', output_dir = '.')" '${workflow.projectDir}/bin/Rmarkdown_scripts/' ${sampleURL} ${sample} ${author}

    """
}



// Process 22:
process RmarkdownCDEGIK {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'




    conda 'r-base=4.3 r-knitr r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/C.Rmd', output_dir = '.')" '${workflow.projectDir}/bin/Rmarkdown_scripts/' ${sample} ${author}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/D.Rmd', output_dir = '.')" '${workflow.projectDir}/bin/Rmarkdown_scripts/' ${sample} ${author}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/E.Rmd', output_dir = '.')" '${workflow.projectDir}/bin/Rmarkdown_scripts/' ${sample} ${author}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/G.Rmd', output_dir = '.')" '${workflow.projectDir}/bin/Rmarkdown_scripts/' ${sample} ${author}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/I.Rmd', output_dir = '.')" '${workflow.projectDir}/bin/Rmarkdown_scripts/' ${sample} ${author}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/K.Rmd', output_dir = '.')" '${workflow.projectDir}/bin/Rmarkdown_scripts/' ${sample} ${author}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/P.Rmd', output_dir = '.')" '${workflow.projectDir}/bin/Rmarkdown_scripts/' ${sample} ${author}


    """
}


// Process 23:
process RmarkdownH {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'




    conda 'r-base=4.3 r-readr r-tidyr r-knitr r-gridExtra r-kableExtra r-png r-gridbase r-DT r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(kallistotop20graphtrinity), path(kallistotop500graphtrinity), path(alluvial1plot), path(alluvial1legend), path(alluvial2plot), path(alluvial2legend), path(pie1), path(pie2), path(pie3), path(pie4plot), path(pie4legend), path(topkallisto), path(busco_figure)

    output:
    tuple val(sample), path("*.html")

    script:
    """
    kallistotop20graphtrinity_abs=\$(readlink -f "${kallistotop20graphtrinity}")
    kallistotop500graphtrinity_abs=\$(readlink -f "${kallistotop500graphtrinity}")
    alluvial1P_abs=\$(readlink -f "${alluvial1plot}")
    alluvial2P_abs=\$(readlink -f "${alluvial2plot}")
    alluvial1L_abs=\$(readlink -f "${alluvial1legend}")
    alluvial2L_abs=\$(readlink -f "${alluvial2legend}")
    pie1_abs=\$(readlink -f "${pie1}")
    pie2_abs=\$(readlink -f "${pie2}")
    pie3_abs=\$(readlink -f "${pie3}")
    pie4P_abs=\$(readlink -f "${pie4plot}")
    pie4L_abs=\$(readlink -f "${pie4legend}")

    topkallisto_abs=\$(readlink -f "${topkallisto}")

    busco_files=(${busco_figure})
    busco_figure1_abs=\${busco_files[0]:-""}
    busco_figure2_abs=\${busco_files[1]:-""}
    busco_figure3_abs=\${busco_files[2]:-""}
    
    # Get absolute paths only for files that exist
    [ -n "\$busco_figure1_abs" ] && busco_figure1_abs=\$(readlink -f "\$busco_figure1_abs")
    [ -n "\$busco_figure2_abs" ] && busco_figure2_abs=\$(readlink -f "\$busco_figure2_abs")
    [ -n "\$busco_figure3_abs" ] && busco_figure3_abs=\$(readlink -f "\$busco_figure3_abs")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/H.Rmd',
      output_dir='.',
      params=list(
        kallistotop20graphtrinity='\$kallistotop20graphtrinity_abs',
        kallistotop500graphtrinity='\$kallistotop500graphtrinity_abs',
        busco_figure1='\$busco_figure1_abs',
        busco_figure2='\$busco_figure2_abs',
        busco_figure3='\$busco_figure3_abs',
        alluvial1P='\$alluvial1P_abs',
        alluvial1L='\$alluvial1L_abs',
        alluvial2P='\$alluvial2P_abs',
        alluvial2L='\$alluvial2L_abs',
        pie1='\$pie1_abs',
        pie2='\$pie2_abs',
        pie3='\$pie3_abs',
        pie4P='\$pie4P_abs',
        pie4L='\$pie4L_abs',
        topkallisto='\$topkallisto_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}


// Process 24:
process RmarkdownJ {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'




    conda 'r-base=4.3 r-readr r-tidyr r-knitr r-gridExtra r-kableExtra r-png r-gridbase r-DT r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(kallistotop20graphtransdecoder), path(kallistotop500graphtransdecoder), path(alluvial3plot), path(alluvial3legend), path(alluvial4plot), path(alluvial4legend), path(pie5), path(pie6), path(pie7), path(pie8plot), path(pie8legend), path(topkallisto_transdecoder), path(busco_figure_transdecoder)

    output:
    tuple val(sample), path("*.html")

    script:
    """
    kallistotop20graphtransdecoder_abs=\$(readlink -f "${kallistotop20graphtransdecoder}")
    kallistotop500graphtransdecoder_abs=\$(readlink -f "${kallistotop500graphtransdecoder}")
    alluvial3P_abs=\$(readlink -f "${alluvial3plot}")
    alluvial4P_abs=\$(readlink -f "${alluvial4plot}")
    alluvial3L_abs=\$(readlink -f "${alluvial3legend}")
    alluvial4L_abs=\$(readlink -f "${alluvial4legend}")
    pie5_abs=\$(readlink -f "${pie5}")
    pie6_abs=\$(readlink -f "${pie6}")
    pie7_abs=\$(readlink -f "${pie7}")
    pie8P_abs=\$(readlink -f "${pie8plot}")
    pie8L_abs=\$(readlink -f "${pie8legend}")
    topkallisto_transdecoder_abs=\$(readlink -f "${topkallisto_transdecoder}")

    busco_files=(${busco_figure_transdecoder})
    busco_figure1_abs=\${busco_files[0]:-""}
    busco_figure2_abs=\${busco_files[1]:-""}
    busco_figure3_abs=\${busco_files[2]:-""}
    
    # Get absolute paths only for files that exist
    [ -n "\$busco_figure1_abs" ] && busco_figure1_abs=\$(readlink -f "\$busco_figure1_abs")
    [ -n "\$busco_figure2_abs" ] && busco_figure2_abs=\$(readlink -f "\$busco_figure2_abs")
    [ -n "\$busco_figure3_abs" ] && busco_figure3_abs=\$(readlink -f "\$busco_figure3_abs")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/J.Rmd',
      output_dir='.',
      params=list(
        kallistotop20graphtransdecoder='\$kallistotop20graphtransdecoder_abs',
        kallistotop500graphtransdecoder='\$kallistotop500graphtransdecoder_abs',
        busco_figure1='\$busco_figure1_abs',
        busco_figure2='\$busco_figure2_abs',
        busco_figure3='\$busco_figure3_abs',
        alluvial3P='\$alluvial3P_abs',
        alluvial3L='\$alluvial3L_abs',
        alluvial4P='\$alluvial4P_abs',
        alluvial4L='\$alluvial4L_abs',
        pie5='\$pie5_abs',
        pie6='\$pie6_abs',
        pie7='\$pie7_abs',
        pie8P='\$pie8P_abs',
        pie8L='\$pie8L_abs',
        topkallisto_transdecoder='\$topkallisto_transdecoder_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 25:
process RmarkdownL {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-readr r-tidyr r-knitr r-gridExtra r-kableExtra r-png r-gridbase r-DT r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(alluvial5plot), path(alluvial5legend), path(alluvial6plot), path(alluvial6legend), path(pie9), path(pie10), path(pie11), path(pie12plot), path(pie12legend), path(topkallisto_signalp)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    alluvial5P_abs=\$(readlink -f "${alluvial5plot}")
    alluvial6P_abs=\$(readlink -f "${alluvial6plot}")
    alluvial5L_abs=\$(readlink -f "${alluvial5legend}")
    alluvial6L_abs=\$(readlink -f "${alluvial6legend}")
    pie9_abs=\$(readlink -f "${pie9}")
    pie10_abs=\$(readlink -f "${pie10}")
    pie11_abs=\$(readlink -f "${pie11}")
    pie12P_abs=\$(readlink -f "${pie12plot}")
    pie12L_abs=\$(readlink -f "${pie12legend}")
    topkallisto_signalp_abs=\$(readlink -f "${topkallisto_signalp}")

    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/L.Rmd',
      output_dir='.',
      params=list(
        alluvial5P='\$alluvial5P_abs',
        alluvial5L='\$alluvial5L_abs',
        alluvial6P='\$alluvial6P_abs',
        alluvial6L='\$alluvial6L_abs',
        pie9='\$pie9_abs',
        pie10='\$pie10_abs',
        pie11='\$pie11_abs',
        pie12P='\$pie12P_abs',
        pie12L='\$pie12L_abs',
        topkallisto_signalp='\$topkallisto_signalp_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 26:
process RmarkdownM {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Table1)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Table1_abs=\$(readlink -f "${Table1}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/M.Rmd',
      output_dir='.',
      params=list(
        Table1='\$Table1_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 27:
process RmarkdownN {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Table2), path(Table3)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Table2_abs=\$(readlink -f "${Table2}")
    Table3_abs=\$(readlink -f "${Table3}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/N.Rmd',
      output_dir='.',
      params=list(
        Table2='\$Table2_abs',
        Table3='\$Table3_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 28:
process RmarkdownO {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Table4)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Table4_abs=\$(readlink -f "${Table4}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/O.Rmd',
      output_dir='.',
      params=list(
        Table4='\$Table4_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 29:
process RmarkdownQ {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Table5)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Table5_abs=\$(readlink -f "${Table5}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/Q.Rmd',
      output_dir='.',
      params=list(
        Table5='\$Table5_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 30:
process RmarkdownR {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Table6), path(Table7)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Table6_abs=\$(readlink -f "${Table6}")
    Table7_abs=\$(readlink -f "${Table7}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/R.Rmd',
      output_dir='.',
      params=list(
        Table6='\$Table6_abs',
        Table7='\$Table7_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 31:
process RmarkdownS {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Table8)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Table8_abs=\$(readlink -f "${Table8}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/S.Rmd',
      output_dir='.',
      params=list(
        Table8='\$Table8_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}



// Process 32:
process RmarkdownV {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Table9)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Table9_abs=\$(readlink -f "${Table9}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/V.Rmd',
      output_dir='.',
      params=list(
        Table9='\$Table9_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 33:
process RmarkdownW {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Table10), path(Table11)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Table10_abs=\$(readlink -f "${Table10}")
    Table11_abs=\$(readlink -f "${Table11}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/W.Rmd',
      output_dir='.',
      params=list(
        Table10='\$Table10_abs',
        Table11='\$Table11_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 34:
process RmarkdownX {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr r-downloadthis r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Table12)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    Table12_abs=\$(readlink -f "${Table12}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/X.Rmd',
      output_dir='.',
      params=list(
        Table12='\$Table12_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

// Process 35 STRICT:
process RmarkdownZ {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4



    label 'process_low'




    conda 'conda-forge::r-base=4.3 conda-forge::r-rmarkdown r-DT R-dplyr r-knitr r-png r-gridbase r-downloadthis r-gridExtra'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(Vennlax), path(VennStrict), path(table), val(protspace)

    output:
    tuple val(sample), path("*.html")

    script:
    """
    Venn_abs1=\$(readlink -f "${Vennlax}")
    Venn_abs2=\$(readlink -f "${VennStrict}")
    table_abs=\$(readlink -f "${table}")

   version=\$([[ "${protspace}" == "TRUE" || "${protspace}" == "True" || "${protspace}" == "true" ]] && echo V1 || echo V2)
    Rmarkdown="${workflow.projectDir}/bin/Rmarkdown_scripts/\${version}/Z.Rmd"


    Rscript -e "rmarkdown::render(
      '\$Rmarkdown',
      output_dir='.',
      params=list(
        VENN1='\$Venn_abs1',
        VENN2='\$Venn_abs2',
        TABLE='\$table_abs',
        AuthorName='${author}',
        SampleName='${sample}',
        Path= '${workflow.projectDir}/bin/Rmarkdown_scripts'
      )
    )"
    """
}

process RmarkdownF {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'




    conda 'r-base=4.3 r-readr r-tidyr r-knitr r-gridExtra r-png r-gridbase r-rmarkdown'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(IP_1), path(IP_2), path(IP_3), path(IP_1_Legend), path(IP_2_Legend), path(IP_3_Legend), path(MF_1), path(MF_2), path(MF_3), path(MF_1_Legend), path(MF_2_Legend), path(MF_3_Legend), path(BP_1), path(BP_2), path(BP_3), path(BP_1_Legend), path(BP_2_Legend), path(BP_3_Legend)

    output:
    tuple val(sample), path("*.html")

    script:
    """
    IP_1_abs=\$(readlink -f "${IP_1}")
    IP_2_abs=\$(readlink -f "${IP_2}")
    IP_3_abs=\$(readlink -f "${IP_3}")
    IP_1_Legend_abs=\$(readlink -f "${IP_1_Legend}")
    IP_2_Legend_abs=\$(readlink -f "${IP_2_Legend}")
    IP_3_Legend_abs=\$(readlink -f "${IP_3_Legend}")
    MF_1_abs=\$(readlink -f "${MF_1}")
    MF_2_abs=\$(readlink -f "${MF_2}")
    MF_3_abs=\$(readlink -f "${MF_3}")
    MF_1_Legend_abs=\$(readlink -f "${MF_1_Legend}")
    MF_2_Legend_abs=\$(readlink -f "${MF_2_Legend}")
    MF_3_Legend_abs=\$(readlink -f "${MF_3_Legend}")
    BP_1_abs=\$(readlink -f "${BP_1}")
    BP_2_abs=\$(readlink -f "${BP_2}")
    BP_3_abs=\$(readlink -f "${BP_3}")
    BP_1_Legend_abs=\$(readlink -f "${BP_1_Legend}")
    BP_2_Legend_abs=\$(readlink -f "${BP_2_Legend}")
    BP_3_Legend_abs=\$(readlink -f "${BP_3_Legend}")

    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/F.Rmd',
      output_dir='.',
      params=list(
        IP_2='\$IP_2_abs',
        IP_2_Legend='\$IP_2_Legend_abs',
        MF_2='\$MF_2_abs',
        MF_2_Legend='\$MF_2_Legend_abs',
        BP_2='\$BP_2_abs',
        BP_2_Legend='\$BP_2_Legend_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
        Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/T.Rmd',
      output_dir='.',
      params=list(
        IP_1='\$IP_1_abs',
        IP_1_Legend='\$IP_1_Legend_abs',
        MF_1='\$MF_1_abs',
        MF_1_Legend='\$MF_1_Legend_abs',
        BP_1='\$BP_1_abs',
        BP_1_Legend='\$BP_1_Legend_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
        Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/Y.Rmd',
      output_dir='.',
      params=list(
        IP_3='\$IP_3_abs',
        IP_3_Legend='\$IP_3_Legend_abs',
        MF_3='\$MF_3_abs',
        MF_3_Legend='\$MF_3_Legend_abs',
        BP_3='\$BP_3_abs',
        BP_3_Legend='\$BP_3_Legend_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}
// Process 36:
process Blast0Chunks {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Blast0/", mode: 'copy'

    input:
    tuple val(sample), path(blastx0), path(blastp0)

    output:
    tuple val(sample), path("*"), emit: blast0chunks

    script:

    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS11.R" ${sample} ${blastx0} ${blastp0}
    """
}



process Blast0Chunksn {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-kableExtra r-DT r-dplyr'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Blast0/", mode: 'copy'

    input:
    tuple val(sample), path(blastx0), path(blastp0), path(blastn0)

    output:
    tuple val(sample), path("*"), emit: blast0chunks

    script:

    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS11.R" ${sample} ${blastx0} ${blastp0} ${blastn0}
    """
}



process Annotate {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'


    conda 'r-base=4.3 r-knitr r-dplyr bioconductor-biostrings r-rentrez r-stringr r-xml2 r-tidyr conda-forge::r-archive r-readr r-purrr'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/Overview/Dataframes/Annotated/", pattern: "*Annotated_df.csv", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/ProtSpace/Input/", pattern: "*ProtSpaceAnnotation.csv", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/ProtSpace/Input/", pattern: "*.fasta", mode: 'copy'
    publishDir "${params.outdir}/${sample}/FinalOutputs/", pattern: "*Select_Annotated_df.csv", mode: 'copy'

    input:
    tuple val(sample), path(final_filtered_lax), path(toxprotblast6), path(nontoxprotblast6), path(Diamondblast6), path(toxprotblastmetadata), path(nontoxprotmetadata), path(toxvsnontoxIP), path(toxvsnontoxMF), path(toxvsnontoxBP)

    output:
    tuple val(sample), path("*all_Annotated_df.csv"), emit: Annotated_df
    tuple val(sample), path("*Select_Annotated_df.csv")
    tuple val(sample), path("*ProtSpaceAnnotation.csv"), emit: ProtSpaceAnnotatedCSV
    tuple val(sample), path("*_ProtSpacePEP.fasta"), emit: FilteredLaxPep
    tuple val(sample), path("*_ProtSpaceCDS.fasta"), emit: FilteredLaxcds

    script:

    """
    # Check if file exists and is not empty placeholder
    if [ -s ${Diamondblast6} ] && [ "\$(cat ${Diamondblast6})" != "NULL" ]; then
        BLASTN_ARG="${Diamondblast6}"
    else
        BLASTN_ARG="NULL"
    fi


    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS16.R" ${final_filtered_lax} ${toxprotblast6} ${toxprotblastmetadata} ${nontoxprotblast6} ${nontoxprotmetadata} ${toxvsnontoxIP} ${toxvsnontoxMF} ${toxvsnontoxBP} ${sample} ${Diamondblast6}
    """
}

process ProtSpace {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'





    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/ProtSpace/", mode: 'copy'
    publishDir "${params.outdir}/${sample}/FinalOutputs/", pattern: "*.parquetbundle", mode: 'copy'

    input:
    tuple val(sample), path(filteredlaxfasta), path(ProtSpaceAnnotatedCSV)

    output:
    tuple val(sample), path("projections/projections_data.parquet"), emit: ProtSpaceParquet
    tuple val(sample), path("*.parquetbundle")

    script:
    // esm2_650m is used instead of prot_t5 as there are permission errors when prot_t5 is used on the test cluster
    """
    ProtSpaceAnnotatedCSV_abs=\$(readlink -f "${ProtSpaceAnnotatedCSV}")

    protspace embed -i ${filteredlaxfasta} -e esm2_3b -o embeddings/
    protspace project -i embeddings/esm2_3b.h5 -m umap2 -o projections/
    protspace prepare -i embeddings/esm2_3b.h5 -a ${ProtSpaceAnnotatedCSV}



    """
}

process ProtSpaceToxin {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4


    label 'process_low'



    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/ProtSpace/", mode: 'copy'
    publishDir "${params.outdir}/${sample}/FinalOutputs/", pattern: "*.parquetbundle", mode: 'copy'

    input:
    tuple val(sample), path(final_filtered_lax), path(table), path(ToxinFasta), path(ToxinMetadata)

    output:
    tuple val(sample), path("*withtoxins.parquetbundle")

    script:
    // esm2_650m is used instead of prot_t5 as there are permission errors when prot_t5 is used on the test cluster
    """
    awk -F'|' '/^>/ { print ">" \$3 } 1' ${ToxinFasta} > ToxProtFasta.int.fasta 
    awk '/^>/ {print \$1; next} {print}' ToxProtFasta.int.fasta  > ToxProtFasta.cleaned.fasta
    cat ${final_filtered_lax} ToxProtFasta.cleaned.fasta > ProtspaceInput.fasta
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS17.R ${table} ${ToxinMetadata}
    protspace embed -i ProtspaceInput.fasta -e esm2_3b -o embeddings/
    protspace project -i embeddings/esm2_3b.h5 -m umap2 -o projections/
    protspace prepare -i embeddings/esm2_3b.h5 -a ProtSpaceToxin.csv
    mv data.parquetbundle data.withtoxins.parquetbundle



    """
}

// Process 35 ProtSpace Page // Optional:
process RmarkdownU {
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_low'

    conda 'r-base=4.3 r-rmarkdown r-plotly  r-dplyr r-arrow r-umap r-forcats'

    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

    input:
    tuple val(sample), val(author), path(parquet), path(metadata)

    output:
    tuple val(sample), path("*.html")

    script:
    """

    parquet_abs=\$(readlink -f "${parquet}")
    metadata_abs=\$(readlink -f "${metadata}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/U.Rmd',
      output_dir='.',
      params=list(
        parquet_path='\$parquet_abs',
        metadata_path='\$metadata_abs',
        AuthorName = '${author}',
        SampleName = '${sample}'
      )
    )"
    """
}

//Minimap for just our candidates that fall in category to visualise against genome 
process Minimap {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }

    maxRetries 4




    label 'process_high'


    conda 'minimap2 bioconda::samtools bioconda::stringtie bioconda::bedtools'

    publishDir "${params.outdir}/${sample}/Pipelines/Analysis/IntermediateFiles/Minimap/CandidateToxins/", mode: 'copy'

    input:
    tuple val(sample), path(genome), path(FilteredLaxcds)

    output:
    tuple val(sample), path("minimaptrial.bam"), emit: bam
    tuple val(sample), path("minimaptrial.sorted.bam"), emit: selectbam
    tuple val(sample), path("minimaptrial.sorted.bam.bai"), emit: selectbambai
    tuple val(sample), path("expandedgenomeregions.fa"), emit: selectfasta
    tuple val(sample), path("expandedgenomeregions.fa.fai"), emit: selectfastafai

    script:
    """
    minimap2 -ax splice:hq -uf ${genome} ${FilteredLaxcds} > minimaptrial.sam
    samtools view -bS minimaptrial.sam > minimaptrial.bam
    samtools sort minimaptrial.bam -o minimaptrial.sorted.bam
    stringtie -o minimaptrial.gtf minimaptrial.sorted.bam
    bedtools bamtobed -split -i minimaptrial.bam > minimaptrial.bed
    rm -r minimaptrial.sam
    samtools faidx ${genome}
    cut -f 1,2 ${genome}.fai > chrom.sizes
    bedtools slop -i minimaptrial.bed -g chrom.sizes -b 500 > minimaptrial.slop.bed
    bedtools getfasta -fi ${genome} -bed minimaptrial.slop.bed -s -name -split -fo expandedgenomeregions.fa
    samtools faidx expandedgenomeregions.fa
    samtools index minimaptrial.sorted.bam

    """
}
// Define input file patterns via parameters



workflow {

    //Define CSV channel.. channel factory creates the channel from a csv file. can be defined in the config or command line 
    csv_channel = channel.fromPath(params.input_csv).splitCsv(header: true, sep: ',').map { row -> row.collectEntries { key, value -> [key.replaceAll('"', ''), value?.toString()?.replaceAll('"', '')] } }
    //Define results folder from the Venomflow 
    venomflowfiles = csv_channel.map { row ->
        def samplename = row.Sample_name
        //0
        def results_path = file(row.Venomflowresultsfolder).toAbsolutePath()

        def kallisto_trinity = results_path ? file("${results_path}/kallisto/trinity/output/abundance.tsv") : ''
        //1
        def kallisto_trans = results_path ? file("${results_path}/kallisto/transdecoder/output/abundance.tsv") : ''
        //2
        def combined_pep = results_path ? file("${results_path}/ORFprediction/Combined/All/*cd95.pep") : ''
        //3
        def combined_cds = results_path ? file("${results_path}/ORFprediction/Combined/All/*cd95pep.cds") : ''
        //4
        def mature_fasta = results_path ? file("${results_path}/Secreted/Mature/Signalp/*_mature.fasta") : ''
        //5
        def blastx_files = results_path ? file("${results_path}/Blast/Blastx/*.6.txt") : ''
        //6
        def blastp_files = results_path ? file("${results_path}/Blast/Blastp_Toxin/*.6.txt") : ''
        //7
        def Interproscan_file = results_path ? file("${results_path}/Interproscan/*.tsv") : ''
        //8
        def signalp_summary = results_path ? file("${results_path}/Secreted/Mature/Signalp/*_summary.signalp5") : ''
        //9
        def Blastn6 = results_path ? file("${results_path}/Blast/Blastn/*.blastn.db.6.txt") : ''
        //10
        def blastn0txt = results_path ? file("${results_path}/Blast/Blastn/*.blastn.db.0.txt") : ''
        //11
        def blastx0txt = results_path ? file("${results_path}/Blast/Blastx/*.blastx.db.0.txt") : ''
        //12
        def blastp0txt = results_path ? file("${results_path}/Blast/Blastp_Toxin/*.blastp.db.0.txt") : ''
        //13
        def basename = row.basename
        //14
        def busco_transcriptome_dir = results_path ? file("${results_path}/BUSCO/transcriptome/Transcriptome1") : ''
        //15
        def busco_translatome_dir = results_path ? file("${results_path}/BUSCO/translatome/Transdecoder/") : ''
        //16
        def genomeid = row.NCBI_Genome_id ?: ''
        //17
        def species = row.Species ?: ''
        //18 
        def MS = row.massspec_csv ? file(row.massspec_csv) : ''
        //19
        def Gavailability = row.isgenomeavailble ?: ''
        //20
        def MSavailability = row.ismassspecavailable ?: ''
        //21
        def ToxinDomains = row.Toxin_domains ? file(row.Toxin_domains) : ''
        //22
        def busco_transcriptome_dir2 = results_path ? file("${results_path}/BUSCO/transcriptome/Transcriptome2") : ''
        //23
        def busco_transcriptome_dir3 = results_path ? file("${results_path}/BUSCO/transcriptome/Combined") : ''
        //24
        def busco_translatome_dir2 = results_path ? file("${results_path}/BUSCO/translatome/TD2/") : ''
        //25
        def busco_translatome_dir3 = results_path ? file("${results_path}/BUSCO/translatome/Combined/") : ''
        //26
        def Transcriptome1 = row.Transcriptome1 ? file(row.Transcriptome1) : ''
        //27
        def Transcriptome2 = row.Transcriptome2 ? file(row.Transcriptome2) : ''
        //28
        def TranscriptomeC = results_path ? file("${results_path}/Combined_Transcriptome/*.cdhit95.fasta") : ''
        //29
        def complete_pep = results_path ? file("${results_path}/ORFprediction/Combined/Complete/*.pep") : ''
        //30
        def complete_cds = results_path ? file("${results_path}/ORFprediction/Combined/Complete/*.cds") : ''
        //31
        def combined_mature = results_path ? file("${results_path}/Secreted/Mature/Combined/*.combined.mature.deduplicated.pep.fasta") : ''
        //32
        def SampleURL = row.SearchAndDownloadURL ?: ''
        //33
        def TD_pep = results_path ? file("${results_path}/ORFprediction/Transdecoder/*.pep") : ''
        //34
        def TD2_pep = results_path ? file("${results_path}/ORFprediction/TD2/*.pep") : ''
        //35
        def TD_cds = results_path ? file("${results_path}/ORFprediction/Transdecoder/*.cds") : ''
        //36
        def TD2_cds = results_path ? file("${results_path}/ORFprediction/TD2/*.cds") : ''
        //37
        def Diamondblast6 = row.Diamondblast6Path ? file(row.Diamondblast6Path) : ''
        //38
        def BlastpNonToxin = results_path ? file("${results_path}/Blast/Blastp_NonToxin/*.blastp.db.6.txt") : ''
        //39
        def AuthorName = row.AnalysisdataAuth ?: ''
        //40
        def genome = row.Genome_fasta_path ? file(row.Genome_fasta_path) : ''
        return [samplename, kallisto_trinity, kallisto_trans, combined_pep, combined_cds, mature_fasta, blastx_files, blastp_files, Interproscan_file, signalp_summary, Blastn6, blastn0txt, blastx0txt, blastp0txt, basename, busco_transcriptome_dir, busco_translatome_dir, genomeid, species, MS, Gavailability, MSavailability, ToxinDomains, busco_transcriptome_dir2, busco_transcriptome_dir3, busco_translatome_dir2, busco_translatome_dir3, Transcriptome1, Transcriptome2, TranscriptomeC, complete_pep, complete_cds, combined_mature, SampleURL, TD_pep, TD2_pep, TD_cds, TD2_cds, Diamondblast6, BlastpNonToxin, AuthorName, genome]
    }

    // Metadafiles 
    // IP entry list
    def interproscan_metadata_file = file(params.input_interproscan, checkIfExists: false).exists()
        ? file(params.input_interproscan)
        : file(params.Fallback_interproscan)
    InterproscanMetadata = Channel.fromPath(interproscan_metadata_file)
    // Panther names entry list
    def panther_metadata_file = file(params.input_panther, checkIfExists: false).exists()
        ? file(params.input_panther)
        : file(params.Fallback_panther)
    PantherMetadata = Channel.fromPath(panther_metadata_file)

    // Toxin metadata TSV
    def toxin_metadata_file = file(params.input_toxin_metadata, checkIfExists: false).exists()
        ? file(params.input_toxin_metadata)
        : file(params.Fallback_toxin_metadata)
    ToxinMetadata = Channel.fromPath(toxin_metadata_file)

    // NonToxin metadata TSV
    def nontoxin_metadata_file = file(params.input_nontoxin_metadata, checkIfExists: false).exists()
        ? file(params.input_nontoxin_metadata)
        : file(params.Fallback_nontoxin_metadata)
    NonToxinMetadata = Channel.fromPath(nontoxin_metadata_file)


    // Process 1: Define Inputs for kallistoanalysistrinity  sample id + kallisto_file trinity tuple 
    kallistoanalysistrinityinput = venomflowfiles.map { [it[0], it[1]] }
    //Run process kallistoAnalysisTrinity
    kallistoanalysistrinityinput | kallistoAnalysisTrinity


    //Process 2: Define Inputs for kallistoanalysistrans sample id + kallisto_file trans tuple 
    kallistoanalysistransinput = venomflowfiles.map { [it[0], it[2]] }
    //Run process kallistoanalysistrans
    kallistoanalysistransinput | kallistoAnalysisTrans

    //Process 3: Define Inputs for ExtractSignalSequences sample id + transdecoderpep + maturefasta tuple 

    def ExtractSignalSequencesinput = venomflowfiles.map { item ->
        def maturefasta = item[32] && item[32] != '' ? item[32] : item[5]
        [item[0], item[30], maturefasta]
    }
    //Run process ExtractSignalSequences
    ExtractSignalSequencesinput | ExtractSignalSequences


    //Process 4: Define Inputs for ExtractCreateTrinityDataframe  sample + trinity_fasta + blastx_file + kallisto_csv  tuple 

    CreateTrinityDataframeinput_single = venomflowfiles
        .filter { !it[28] }
        .map { [it[0], it[27], it[6]] }
        .join(kallistoAnalysisTrinity.out.trin_all_csv)

    CreateTrinityDataframeinput_combined = venomflowfiles
        .filter { it[28] }
        .map { [it[0], it[29], it[6]] }
        .join(kallistoAnalysisTrinity.out.trin_all_csv)

    CreateTrinityDataframeinput = CreateTrinityDataframeinput_single.mix(CreateTrinityDataframeinput_combined)
    //Run process ExtractCreateTrinityDataframe
    CreateTrinityDataframeinput | CreateTrinityDataframe


    //Process 5: Define Inputs for CreateInterproscanDataframe  sample + Interproscan + ListFile + PantherFile  tuple 
    CreateInterproscanDataframeinput = venomflowfiles
        .map {
            return [it[0], it[8]]
        }
        .combine(InterproscanMetadata)
        .combine(PantherMetadata)
    //Run process CreateInterproscanDataframe
    CreateInterproscanDataframeinput | CreateInterproscanDataframe


    //Process 6: Define Inputs for CreateTransdecoderDataframe  sample + transdecoder_pep + transdecoder_cds + blastp6_file  mature_fasta, Signalp_summary, signalsequences, Interproscan_dataframe, kallistotrans
    CreateTransdecoderDataframeinput = venomflowfiles
        .map { item ->
            def maturefasta = item[32] && item[32] != '' ? item[32] : item[5]
            def combinedpep2 = item[3] && item[3] != ''
                ? item[3]
                : (item[34] && item[34] != '' ? item[34] : item[35])
            def combinedcds2 = item[4] && item[4] != ''
                ? item[4]
                : (item[36] && item[36] != '' ? item[36] : item[37])
            [item[0], combinedpep2, combinedcds2, item[7], maturefasta]
        }
        .join(ExtractSignalSequences.out.signalsequences)
        .join(CreateInterproscanDataframe.out.Interproscan_dataframe)
        .join(kallistoAnalysisTrans.out.trans_all_csv)
    //Run process ExtractCreateTrinityDataframe
    CreateTransdecoderDataframeinput | CreateTransdecoderDataframe


    //Process 8: Define Input BUSCOtranscriptome
    BUSCOtranscriptomeinput_1 = venomflowfiles.map { [it[0], it[15], "1"] }
    BUSCOtranscriptomeinput_2 = venomflowfiles.filter { it[28] }.map { [it[0], it[23], "2"] }
    BUSCOtranscriptomeinput_C = venomflowfiles.filter { it[28] }.map { [it[0], it[24], "C"] }
    BUSCOtranscriptomeinput = BUSCOtranscriptomeinput_1.mix(BUSCOtranscriptomeinput_2).mix(BUSCOtranscriptomeinput_C)
    //Run Process BUSCOtranscriptome
    BUSCOtranscriptomeinput | BUSCOtranscriptome

    //Process 9: Define Input BUSCOtranslatome
    BUSCOtranslatomeinput_1 = venomflowfiles.map { [it[0], it[25], "TD2"] }
    BUSCOtranslatomeinput_2 = venomflowfiles.filter { it[16] }.map { [it[0], it[16], "TD"] }
    BUSCOtranslatomeinput_C = venomflowfiles.filter { it[26] }.map { [it[0], it[26], "C"] }
    BUSCOtranslatomeinput = BUSCOtranslatomeinput_1.mix(BUSCOtranslatomeinput_2).mix(BUSCOtranslatomeinput_C)

    //Run Process BUSCOtranslatome
    BUSCOtranslatomeinput | BUSCOtranslatome


    //Process 10: Define Input for TableGenerationTrinity sample+TBK + genomeid + species 
    def genome = venomflowfiles.map {
        return [it[0], it[17]]
    }
    def species = venomflowfiles.map {
        return [
            it[0],
            it[18],
        ]
    }
    def TableGenerationTrinityinput = CreateTrinityDataframe.out.TBK.join(genome).join(species)
    //Run Process TableGenerationTrinity
    TableGenerationTrinityinput | TableGenerationTrinity


    //Process 12: Define Input for FigureGenerationTrinity 
    def FigureGenerationTrinityinput = CreateTrinityDataframe.out.TBK
    // Run process FigureGenerationTrinity
    FigureGenerationTrinityinput | FigureGenerationTrinity


    //Process 13: Define Input FigureGenerationTransdecoder
    def FigureGenerationTransdecoderinput = CreateTransdecoderDataframe.out.transdf
    // Run Process FigureGenerationTransdecoder
    FigureGenerationTransdecoderinput | FigureGenerationTransdecoder

    // Run Process14:  FigureGenerationSignalp
    FigureGenerationTransdecoderinput | FigureGenerationSignalp


    //Process 11: Define 4 possible sample types for table generation transdecoder sample+transdf + genomeid + species 
    SampleGenomeSpecies = venomflowfiles.map {
        return [it[0], it[17], it[18]]
    }


    TableGenerationTransdecoderInput = CreateTransdecoderDataframe.out.transdf.join(SampleGenomeSpecies)
    //Run Process TableGenerationTrinity
    TableGenerationTransdecoderInput | TableGenerationTransdecoder

    //Input 

    //Define Input for ToxinVsNonToxin
    UniProtMetadata = ToxinMetadata.combine(NonToxinMetadata).combine(InterproscanMetadata)

    //Run ToxinvsNonToxin
    UniProtMetadata | ToxinVsNonToxin

    //Define Input CreateInterproscanFigures
    def CreateInterproscanFiguresInput = CreateTransdecoderDataframe.out.transdf_distinct
        .combine(ToxinVsNonToxin.out.toxvsnontoxIP)
        .combine(ToxinVsNonToxin.out.toxvsnontoxMF)
        .combine(ToxinVsNonToxin.out.toxvsnontoxBP)

    // Run CreateInterproscanFigures
    CreateInterproscanFiguresInput | CreateInterproscanFigures

    //Minimap1 for all complete ORFs 

    Minimap1Input = venomflowfiles
        .filter { it[41] }
        .map { item ->
            def combinedcds2 = item[4] && item[4] != ''
                ? item[4]
                : (item[36] && item[36] != '' ? item[36] : item[37])
            [item[0], item[41], combinedcds2]
        }
    //Minimap run 
    Minimap1Input | Minimap1
    // Overview 

    AddMSGenomeIfAvailableAndCreateOverviewInput = venomflowfiles
        .map { it -> [it[0], it[18], it[19] ?: file('NO_FILE'), it[10] ?: file('NO_FILE')] }
        .join(CreateTransdecoderDataframe.out.transdf_distinct)
        .join(Minimap1.out.completecdspaf, remainder: true)
        .map { it -> it[0..-2] + [it[-1] ?: file('NO_FILE')] }
        .combine(ToxinVsNonToxin.out.toxvsnontoxIP)

    AddMSGenomeIfAvailableAndCreateOverviewInput | AddMSGenomeIfAvailableAndCreateOverview
    //Define Input for Annotate
    WithDiamond = venomflowfiles.filter { it[38] }.map { [it[0], it[7], it[39], it[38]] }
    //sample, blastptox, blastpnontox, diamond

    WithoutDiamond = venomflowfiles.filter { !it[38] }.map { [it[0], it[7], it[39], []] }
    //sample, blastptox, blastpnontox, null

    SpeciesSpecificBlastInputs = WithDiamond.mix(WithoutDiamond)

    VennCsvLax = AddMSGenomeIfAvailableAndCreateOverview.out.VennCsvLax
    SpeciesSpecificAnnotateInputs = VennCsvLax.join(SpeciesSpecificBlastInputs)
    CommonMetadata = ToxinMetadata
        .combine(NonToxinMetadata)
        .combine(ToxinVsNonToxin.out.toxvsnontoxIP)
        .combine(ToxinVsNonToxin.out.toxvsnontoxMF)
        .combine(ToxinVsNonToxin.out.toxvsnontoxBP)

    Annotate_Input = SpeciesSpecificAnnotateInputs.combine(CommonMetadata)
    // Run Annotate
    Annotate_Input | Annotate

    Toxin_fasta_file = file(params.toxprot_fasta, checkIfExists: false).exists()
        ? file(params.toxprot_fasta)
        : file(params.Fallback_toxin_fasta)
    ToxinFasta = Channel.fromPath(Toxin_fasta_file)
    ToxinFastaAll = ToxinFasta.join(ToxinMetadata)

    //Define Input for ProtSpace
    if (params.protspace) {
        ProtSpace_input = Annotate.out.FilteredLaxPep.join(Annotate.out.ProtSpaceAnnotatedCSV)
        //Run Process for ProtSpace
        ProtSpace_input | ProtSpace
        // Full protspace alongside toxprot
        ProtSpaceToxinInput = VennCsvLax.join(Annotate.out.Annotated_df).combine(ToxinFastaAll)
        ProtSpaceToxinInput | ProtSpaceToxin
    }

    // Define Sample Metadata
    csv_file = file(params.input_csv)
    header_line = csv_file.text.readLines()[0]

    MetadataInput = csv_channel.map { row ->
        def sample_name = row.Sample_name
        def author_name = row.AnalysisdataAuth
        def csv_content = header_line + '\n' + row.collect { _k, v -> v }.join(',')
        tuple(sample_name, author_name, csv_content)
    }


    samplesheets = CreateSampleSheet(MetadataInput)
    // Then pass to RmarkdownB
    samplesheets | RmarkdownB

    //Minimap define 
    Minimap_input = venomflowfiles
        .filter { it[41] }
        .map { [it[0], it[41]] }
        .join(Annotate.out.FilteredLaxcds)
    //Minimap run 
    Minimap_input | Minimap
    //RmarkdownA 
    RmarkdownAInput = venomflowfiles.map {
        return [it[0], it[40], it[33]]
    }
    RmarkdownAInput | RmarkdownA
    //IntermediateHTMLS
    RmarkdownCDEGIKInput = (venomflowfiles.map {
        return [it[0], it[40]]
    })
    RmarkdownCDEGIKInput | RmarkdownCDEGIK
    groupedbuscotranscriptome = BUSCOtranscriptome.out.busco_transcriptome.groupTuple()
    groupedbuscotranslatome = BUSCOtranslatome.out.busco_translatome.groupTuple()
    //groupedbuscotranscriptome = BUSCOtranscriptome.out.busco_transcriptome
    //groupedbuscotranslatome = BUSCOtranslatome.out.busco_translatome
    // RmarkdownH
    RmarkdownHInput = RmarkdownCDEGIKInput
        .join(kallistoAnalysisTrinity.out.trin_top20_png)
        .join(kallistoAnalysisTrinity.out.trin_top500_png)
        .join(FigureGenerationTrinity.out.alluvial1Plot)
        .join(FigureGenerationTrinity.out.alluvial1Legend)
        .join(FigureGenerationTrinity.out.alluvial2Plot)
        .join(FigureGenerationTrinity.out.alluvial2Legend)
        .join(FigureGenerationTrinity.out.pie1)
        .join(FigureGenerationTrinity.out.pie2)
        .join(FigureGenerationTrinity.out.pie3)
        .join(FigureGenerationTrinity.out.pie4Plot)
        .join(FigureGenerationTrinity.out.pie4Legend)
        .join(FigureGenerationTrinity.out.Table13)
        .join(groupedbuscotranscriptome)

    RmarkdownHInput | RmarkdownH
    // RmarkdownJ
    RmarkdownJInput = RmarkdownCDEGIKInput
        .join(kallistoAnalysisTrans.out.trans_top20_png)
        .join(kallistoAnalysisTrans.out.trans_top500_png)
        .join(FigureGenerationTransdecoder.out.alluvial3Plot)
        .join(FigureGenerationTransdecoder.out.alluvial3Legend)
        .join(FigureGenerationTransdecoder.out.alluvial4Plot)
        .join(FigureGenerationTransdecoder.out.alluvial4Legend)
        .join(FigureGenerationTransdecoder.out.pie5)
        .join(FigureGenerationTransdecoder.out.pie6)
        .join(FigureGenerationTransdecoder.out.pie7)
        .join(FigureGenerationTransdecoder.out.pie8Plot)
        .join(FigureGenerationTransdecoder.out.pie8Legend)
        .join(FigureGenerationTransdecoder.out.Table14)
        .join(groupedbuscotranslatome)

    RmarkdownJInput | RmarkdownJ
    // RmarkdownL
    RmarkdownLInput = RmarkdownCDEGIKInput
        .join(FigureGenerationSignalp.out.alluvial5Plot)
        .join(FigureGenerationSignalp.out.alluvial5Legend)
        .join(FigureGenerationSignalp.out.alluvial6Plot)
        .join(FigureGenerationSignalp.out.alluvial6Legend)
        .join(FigureGenerationSignalp.out.pie9)
        .join(FigureGenerationSignalp.out.pie10)
        .join(FigureGenerationSignalp.out.pie11)
        .join(FigureGenerationSignalp.out.pie12Plot)
        .join(FigureGenerationSignalp.out.pie12Legend)
        .join(FigureGenerationSignalp.out.Table15)

    RmarkdownLInput | RmarkdownL

    // RmarkdownM 
    RmarkdownMInput = RmarkdownCDEGIKInput.join(TableGenerationTrinity.out.Table1)
    RmarkdownMInput | RmarkdownM
    // RmarkdownN
    RmarkdownNInput = RmarkdownCDEGIKInput
        .join(TableGenerationTrinity.out.Table2)
        .join(TableGenerationTrinity.out.Table3)

    RmarkdownNInput | RmarkdownN

    // RmarkdownO
    RmarkdownOInput = RmarkdownCDEGIKInput.join(TableGenerationTrinity.out.Table4)
    RmarkdownOInput | RmarkdownO

    // RmarkdownQ 
    RmarkdownQInput = RmarkdownCDEGIKInput.join(TableGenerationTransdecoder.out.Table5)
    RmarkdownQInput | RmarkdownQ

    // RmarkdownR
    RmarkdownRInput = RmarkdownCDEGIKInput
        .join(TableGenerationTransdecoder.out.Table6)
        .join(TableGenerationTransdecoder.out.Table7)
    RmarkdownRInput | RmarkdownR
    // RmarkdownS
    RmarkdownSInput = RmarkdownCDEGIKInput.join(TableGenerationTransdecoder.out.Table8)

    RmarkdownSInput | RmarkdownS
    // RmarkdownV
    RmarkdownVInput = RmarkdownCDEGIKInput.join(TableGenerationTransdecoder.out.Table9)

    RmarkdownVInput | RmarkdownV
    // RmarkdownW
    RmarkdownWInput = RmarkdownCDEGIKInput
        .join(TableGenerationTransdecoder.out.Table10)
        .join(TableGenerationTransdecoder.out.Table11)
    RmarkdownWInput | RmarkdownW
    // RmarkdownX
    RmarkdownXInput = RmarkdownCDEGIKInput.join(TableGenerationTransdecoder.out.Table12)
    RmarkdownXInput | RmarkdownX


    RmarkdownUInput = RmarkdownCDEGIKInput.join(ProtSpace.out.ProtSpaceParquet).join(Annotate.out.ProtSpaceAnnotatedCSV)
    RmarkdownUInput | RmarkdownU

    /*
    Rmarkdown


*/
    // Input RmarkdownF
    def RmarkdownF_input = RmarkdownCDEGIKInput
        .join(CreateInterproscanFigures.out.IP_1)
        .join(CreateInterproscanFigures.out.IP_2)
        .join(CreateInterproscanFigures.out.IP_3)
        .join(CreateInterproscanFigures.out.IP_1_Legend)
        .join(CreateInterproscanFigures.out.IP_2_Legend)
        .join(CreateInterproscanFigures.out.IP_3_Legend)
        .join(CreateInterproscanFigures.out.MF_1)
        .join(CreateInterproscanFigures.out.MF_2)
        .join(CreateInterproscanFigures.out.MF_3)
        .join(CreateInterproscanFigures.out.MF_1_Legend)
        .join(CreateInterproscanFigures.out.MF_2_Legend)
        .join(CreateInterproscanFigures.out.MF_3_Legend)
        .join(CreateInterproscanFigures.out.BP_1)
        .join(CreateInterproscanFigures.out.BP_2)
        .join(CreateInterproscanFigures.out.BP_3)
        .join(CreateInterproscanFigures.out.BP_1_Legend)
        .join(CreateInterproscanFigures.out.BP_2_Legend)
        .join(CreateInterproscanFigures.out.BP_3_Legend)

    // Rmarkdown F
    RmarkdownF_input | RmarkdownF
    // Input RmarkdownZ 
    def Protspace = channel.value(params.protspace ?: "FALSE")

    def RmarkdownZ_input = RmarkdownCDEGIKInput
        .join(AddMSGenomeIfAvailableAndCreateOverview.out.VennPngLax)
        .join(AddMSGenomeIfAvailableAndCreateOverview.out.VennPngStrict)
        .join(Annotate.out.Annotated_df)
        .combine(Protspace)

    // Rmarkdown Z
    RmarkdownZ_input | RmarkdownZ
    BlastChunksInput = venomflowfiles.map { [it[0], it[12], it[13], it[11] ? it[11] : []] }

    BlastChunksInput | Blast0Chunksn
}
