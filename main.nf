#!/usr/env nextflow


// Process 1: For kallistoanalysistrinity.py python,pandas,seaborn,matplotlib
process kallistoAnalysisTrinity {
    errorStrategy 'ignore'

    conda 'python=3.8 pandas seaborn matplotlib'

    publishDir "${sample}/VenomFlowAnalysis/results/KallistoAnalysis/Trinity/", mode: 'copy'

    input:
    tuple val(sample), path (kallisto_file_trinity)

    output:
    tuple val(sample), path ("*_all.csv"), emit: trin_all_csv
    path "*_top20.csv", emit: trin_top20_csv
    tuple val(sample), path ("*_top500graph.png"), emit: trin_top500_png
    tuple val(sample), path ("*_top20graph.png"), emit: trin_top20_png

    script:

    """
    
    python3 ${workflow.projectDir}/bin/Intermediate_Scripts/kallistoanalysistrinity.py ./ ${params.basename} ${kallisto_file_trinity}
    """
}

// Process 2: For kallistoanalysistrans.py dependencies:python,pandas,seaborn,matplotlib
process kallistoAnalysisTrans {
    errorStrategy 'ignore'

    conda 'python=3.8 pandas seaborn matplotlib'

    publishDir "${sample}/VenomFlowAnalysis/results/KallistoAnalysis/Transdecoder/", mode: 'copy'

    input:
    tuple val(sample) , path (kallisto_file_transdecoder)

    output:
    tuple val(sample), path ("*_all.csv"), emit: trans_all_csv
    path "*_top20.csv", emit: trans_top20_csv
    tuple val(sample), path ("*_top500graph.png"), emit: trans_top500_png
    tuple val(sample), path ("*_top20graph.png"), emit: trans_top20_png

    script:
    """
    
    python3 ${workflow.projectDir}/bin/Intermediate_Scripts/kallistoanalysistrans.py ./ "${params.basename}trans" ${kallisto_file_transdecoder}
    """
}

// Process 3: Extract Signal Sequences dependencies python biopython
process ExtractSignalSequences {
    errorStrategy 'ignore'

    conda 'python=3.8 biopython'

    publishDir "${sample}/VenomFlowAnalysis/results/Signal_sequences/", mode: 'copy'

    input:
    tuple val(sample), path(transdecoder_pep), path(mature_fasta)

    output:
    tuple val(sample) path ("signalsequences.fasta") , emit: signalsequences

    script:
    """
    
    python3 ${workflow.projectDir}/bin/Intermediate_Scripts/IS2.py ${transdecoder_pep} ${mature_fasta} "${sample}_signalsequences.fasta"
    """
}

// Process 4: Create Trinity Dataframe dependecies : R, biocmanager 
process CreateTrinityDataframe {
    errorStrategy 'ignore'

    conda 'r-base bioconductor-biostrings r-tidyr r-dplyr'

    publishDir "${sample}/VenomFlowAnalysis/results/TrinityDataframe/", mode: 'copy'
    publishDir "${sample}/VenomFlowAnalysis/results/RappData/TrinityDataframe/", pattern: "*.gz", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), path(blastx_file), path(kallisto_csv)

    output:
    tuple val(sample), path ("*TBK.csv"), emit: TBK
    path "*TBK_distinct.csv.gz", emit: TBK_distinct

    script:
    """
    
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS1.R ${trinity_fasta} ${blastx_file} ${kallisto_csv} ${sample}
    """
}

// Process 5: Create Interproscan Dataframe dependecies : R, biocmanager 
process CreateInterproscanDataframe {
    errorStrategy 'ignore'

    conda 'r-base bioconductor-biostrings r-dplyr bioconductor-go.db bioconductor-biomart r-tidyr'

    publishDir "${sample}/VenomFlowAnalysis/results/InterproscanDataframe/", mode: 'copy'

    input:
    tuple val(sample), path(Interproscan), path(ListFile), path(PantherFile)

    output:
    tuple val(sample), path ("*Final_interproscan_dataframe.csv"), emit: Interproscan_dataframe

    script:
    """
    
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS4.R" ${Interproscan} ${ListFile} ${PantherFile} ${sample}
    """
}

// Process 6: Create Transdecoder Dataframe dependecies : R, biocmanager 
process CreateTransdecoderDataframe {
    errorStrategy 'ignore'

    conda 'r-base=4.3 r-dplyr r-tidyr bioconductor-biostrings'


    publishDir "${sample}/VenomFlowAnalysis/results/TransdecoderDataframe/", pattern: "*.csv", mode: 'copy'
    publishDir "${sample}/VenomFlowAnalysis/results/Secreted_proteins_fasta/", pattern: "*.fasta", mode: 'copy'

    input:
    tuple val(sample), path(transdecoder_pep), path(transdecoder_cds), path(blastp_file), path(mature_fasta), path(Signalp_summary), path(signalsequences), path(Interproscan_dataframe), path(kallistotrans)

    output:
    tuple val(sample), path ("*transdf.csv"), emit: transdf
    tuple val(sample), path ("*transdf_distinct.csv"), emit: transdf_distinct
    tuple val(sample), path ("*secreted_proteins.pep.fasta"), emit: toxin_fastapep
    tuple val(sample), path ("*secreted_proteins.cds.fasta"), emit: toxin_fastacds


    script:
    """
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS5.R ${transdecoder_pep} ${transdecoder_cds} ${blastp_file} ${mature_fasta} ${Signalp_summary} ${signalsequences} ${Interproscan_dataframe} ${kallistotrans} ${sample}
    """
}

// Process 7: Create Interproscantoxinplotly
process CreateInterproscanToxinPlotly {
    errorStrategy 'ignore'

    conda 'bioconductor-biostrings r-dplyr r-tidyr r-htmlwidgets r-plotly'


    publishDir "${sample}/VenomFlowAnalysis/results/htmls/", pattern: "*.html", mode: 'copy'
    publishDir "${sample}/VenomFlowAnalysis/results/Secreted_proteins_fasta/", pattern: "*.fasta", mode: 'copy'
    publishDir "${sample}/VenomFlowAnalysis/results/Toxin_Domain_Filtering/", pattern: "*.csv", mode: 'copy'


    input:
    tuple val(sample), path(transdf_distinct_csv), path(transtoxinfastacds), path(transtoxinfastapep) path(Toxin_domains)

    output:
    tuple val(sample), path ("secreted_sequences_with_toxin_related_domains.pep.fasta"), emit: filtered_sequencespep
    tuple val(sample), path ("secreted_sequences_with_toxin_related_domains.cds.fasta"), emit: filtered_sequences
    tuple val(sample), path ("plotly_graph.html"), emit: plotly_graph
    tuple val(sample), path ("*.csv"), emit: annotationcsvs


    script:
    """
    
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS7.R ${transdf_distinct_csv} ${Toxin_domains} ${transtoxinfastacds} ${transtoxinfastapep} ${sample}
    """
}

// Process 8: Create BUSCOgraphtranscriptome  
process BUSCOtranscriptome {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco2.yaml"

    publishDir "${sample}/VenomFlowAnalysis/results/busco/transcriptome/${count}/", mode: 'copy'

    input:
    tuple val (sample), path (buscodirtr), val(count)

    output:
    uple val(sample), path ("*.png"), emit: busco_transcriptome

    script:
    """
     
    python3 "${workflow.projectDir}/bin/Intermediate_Scripts/IS6.py" -wd "${buscodirtr}"
    """
}

// Process 9: Create BUSCOgraphtranslatome 
process BUSCOtranslatome {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco2.yaml"

    publishDir "${sample}/VenomFlowAnalysis/results/busco/translatome/${count}/", mode: 'copy'

    input:
    val (sample), path (buscodirtl), val(count)

    output:
    tuple val(sample), path ("*.png"), emit: busco_translatome

    script:
    """
     
    python3 "${workflow.projectDir}/bin/Intermediate_Scripts/IS6.py" -wd "${buscodirtl}"

    """
}

// Process 10: Create FigureGenerationTrinity
process FigureGenerationTrinity {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/Figures/Images/Trinity", pattern: ".png", mode: 'copy'
    publishDir "${sample}/results/Figures/Tables/Trinity", pattern: "*.csv", mode: 'copy'


    input:
    tuple val(sample), path (TBK)

    output:
    tuple val(sample),path "pie1.png", emit: pie1
    tuple val(sample),path "pie2.png", emit: pie2
    tuple val(sample),path "pie3.png", emit: pie3
    tuple val(sample),path "pie4.png", emit: pie4
    tuple val(sample),path "alluvial1.png", emit: alluvial1
    tuple val(sample),path "alluvial2.png", emit: alluvial2
    tuple val(sample),path "Table13.csv", emit: Table13

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Figure_generation_Trinity.R" ${TBK} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"

    """
}

// Process 11: Create FigureGenerationTransdecoder
process FigureGenerationTransdecoder {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/Figures/Images/Transdecoder/",pattern: ".png",  mode: 'copy'
    publishDir "${sample}/results/Figures/Tables/Transdecoder", pattern: "*.csv", mode: 'copy'


    input:
    tuple val(sample),path (transdf)

    output:
    tuple val(sample), path ("pie5.png"), emit: pie5
    tuple val(sample), path ("pie6.png"), emit: pie6
    tuple val(sample), path ("pie7.png"), emit: pie7
    tuple val(sample), path ("pie8.png"), emit: pie8
    tuple val(sample), path ("alluvial3.png"), emit: alluvial3
    tuple val(sample), path ("alluvial4.png"), emit: alluvial4
    tuple val(sample), path ("Table14.csv"), emit: Table14

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Figure_generation_Transdecoder.R" ${transdf} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"

    """
}



// Process 13: Create FigureGenerationSignalp
process FigureGenerationSignalp {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/Figures/Images/SignalP/", pattern: ".png", mode: 'copy'
    publishDir "${sample}/results/Figures/Tables/SignalP/", pattern: "*.csv", mode: 'copy'

    input:
    tuple val(sample), path (transdf)

    output:
    tuple val(sample), path ("pie9.png"), emit: pie9
    tuple val(sample), path ("pie10.png"), emit: pie10
    tuple val(sample), path ("pie11.png"), emit: pie11
    tuple val(sample), path ("pie12.png"), emit: pie12
    tuple val(sample), path ("alluvial5.png"), emit: alluvial5
    tuple val(sample), path ("alluvial6.png"), emit: alluvial6
    tuple val(sample), path ("Table15.csv"), emit: Table15

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Figure_generation_SignalP.R" ${transdf} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"

    """
}

// Process 12: Create TableGenerationTrinity
process TableGenerationTrinity {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/Figures/Tables/Trinity", mode: 'copy'

    input:
    tuple val(sample), path(TBK), val(genome_id), val(species)

    output:
    tuple val(sample), path ("Table1.csv"), emit: Table1
    tuple val(sample), path ("Table2.csv"), emit: Table2
    tuple val(sample), path ("Table3.csv"), emit: Table3
    tuple val(sample), path ("Table4.csv"), emit: Table4

    script:
    """
	
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Generating_TopTables_Trinity.R" ${TBK} ${genome_id} ${species}
    """
}

// Process 11: Create TableGenerationTransdecoder  
process TableGenerationTransdecoder {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/Figures/Tables/Transdecoder", pattern: "*{5,6,7,8}.csv", mode: 'copy'
    publishDir "${sample}/results/Figures/Tables/SignalP", pattern: "*{9,10,11,12}.csv", mode: 'copy'


    input:
    tuple val(sample), path(transdf), val(genome_id), val(species)

    output:
    tuple val(sample),path ("Table5.csv"), emit: Table5
    tuple val(sample),path ("Table6.csv"), emit: Table6
    tuple val(sample),path ("Table7.csv"), emit: Table7
    tuple val(sample),path ("Table8.csv"), emit: Table8
    tuple val(sample),path ("Table9.csv"), emit: Table9
    tuple val(sample),path ("Table10.csv"), emit: Table10
    tuple val(sample),path ("Table11.csv"), emit: Table11
    tuple val(sample),path ("Table12.csv"), emit: Table12

    script:
    """
	
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Generating_Tables_Transdecoder_SignalP.R" ${transdf} ${genome_id} ${species} 
    """
}

//Process 12: Add Massspec and Blastn6 results where available

process AddMSGenomeIfAvailableAndCreateOverview {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/FinalDataFrames/", pattern: "*.csv", mode: 'copy'
    publishDir "${sample}/results/Overview/", pattern: "*.csv", mode: 'copy'


    input:
    tuple val(sample), val(species), path(massspec), , path(blastn6), path(transdf) path(Toxindomains)

    output:
    path("*.csv")
    tuple val(sample), path(*_Venn_strict.png), emit VennPngStrict
    tuple val(sample), path(*_Venn_Diagram_union_strict.csv), emit VennCsvStrict
    tuple val(sample), path(*_Venn_lax.png), emit VennPngLax
    tuple val(sample), path(*_Venn_Diagram_union_lax.csv), emit VennCsvLax

    script:
    """
	    # Check if file exists and is not empty placeholder
    if [ -s ${blastn6} ] && [ "\$(cat ${blastn6})" != "NULL" ]; then
        BLASTN_ARG="${blastn6}"
    else
        BLASTN_ARG="NULL"
    fi

     if [ -s ${massspec} ] && [ "\$(cat ${massspec})" != "NULL" ]; then
        MS_ARG="${massspec}"
    else
        MS_ARG="NULL"
    fi

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS14.R" ${sample} ${species} ${transdf} ${Toxindomains}  \$BLASTN_ARG \$MS_ARG
    """
}



// Process 20: RmarkdownA
process RmarkdownA {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), val (sampleURL)

    output:
    path "*.html"

    script:
    """

    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/A.Rmd', output_dir = '.')" "${workflow.projectDir}/bin/Rmarkdown_scripts/" "${sampleURL}" "${params.name}"
    """
}

// Process 21: RmarkdownB
process RmarkdownB {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), val(Name), path(samplesheet)

    output:
    path "*.html"

    script:
    """

    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/B.Rmd',
      output_dir='.',
      knit_root_dir='.',
      args=c('${Name}', '${samplesheet}')
    )"
  

    """
}

// Process 22:
process RmarkdownCDEGIK {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    output:
    path "*.html"

    script:
    """

    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/C.Rmd', output_dir = '.')" "${workflow.projectDir}/bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/D.Rmd', output_dir = '.')" "${workflow.projectDir}/bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/E.Rmd', output_dir = '.')" "${workflow.projectDir}/bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/G.Rmd', output_dir = '.')" "${workflow.projectDir}/bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/I.Rmd', output_dir = '.')" "${workflow.projectDir}/bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('${workflow.projectDir}/bin/Rmarkdown_scripts/K.Rmd', output_dir = '.')" "${workflow.projectDir}/bin/Rmarkdown_scripts/" ${params.name}


    """
}

// Process 23:
process RmarkdownH {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (kallistotop20graphtrinity),
    path (kallistotop500graphtrinity),
    path (busco_figure),
    path (alluvial1),
    path (alluvial2),
    path (pie1),
    path (pie2),
    path (pie3),
    path (pie4),
    path (topkallisto),

    output:
    path "*.html"

    script:
    """
    kallistotop20graphtrinity_abs=\$(readlink -f "${kallistotop20graphtrinity}")
    kallistotop500graphtrinity_abs=\$(readlink -f "${kallistotop500graphtrinity}")
    busco_figure_abs=\$(readlink -f "${busco_figure}")
    alluvial1_abs=\$(readlink -f "${alluvial1}")
    alluvial2_abs=\$(readlink -f "${alluvial2}")
    pie1_abs=\$(readlink -f "${pie1}")
    pie2_abs=\$(readlink -f "${pie2}")
    pie3_abs=\$(readlink -f "${pie3}")
    pie4_abs=\$(readlink -f "${pie4}")
    topkallisto_abs=\$(readlink -f "${topkallisto}")

    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/H.Rmd',
      output_dir='.',
      params=list(
        kallistotop20graphtrinity='\$kallistotop20graphtrinity_abs',
        kallistotop500graphtrinity='\$kallistotop500graphtrinity_abs',
        busco_figure='\$busco_figure_abs',
        alluvial1='\$alluvial1_abs',
        alluvial2='\$alluvial2_abs',
        pie1='\$pie1_abs',
        pie2='\$pie2_abs',
        pie3='\$pie3_abs',
        pie4='\$pie4_abs',
        topkallisto='\$topkallisto_abs',
        name = '${params.name}'
      )
    )"
    """
}


// Process 24:
process RmarkdownJ {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (kallistotop20graphtransdecoder),
    path (kallistotop500graphtransdecoder),
    path (busco_figure_transdecoder),
    path (alluvial3),
    path (alluvial4),
    path (pie5),
    path (pie6),
    path (pie7),
    path (pie8),
    path (topkallisto_transdecoder)

    output:
    path "*.html"

    script:
    """
    kallistotop20graphtransdecoder_abs=\$(readlink -f "${kallistotop20graphtransdecoder}")
    kallistotop500graphtransdecoder_abs=\$(readlink -f "${kallistotop500graphtransdecoder}")
    busco_figure_transdecoder_abs=\$(readlink -f "${busco_figure_transdecoder}")
    alluvial3_abs=\$(readlink -f "${alluvial3}")
    alluvial4_abs=\$(readlink -f "${alluvial4}")
    pie5_abs=\$(readlink -f "${pie5}")
    pie6_abs=\$(readlink -f "${pie6}")
    pie7_abs=\$(readlink -f "${pie7}")
    pie8_abs=\$(readlink -f "${pie8}")
    topkallisto_transdecoder_abs=\$(readlink -f "${topkallisto_transdecoder}")

    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/J.Rmd',
      output_dir='.',
      params=list(
        kallistotop20graphtransdecoder='\$kallistotop20graphtransdecoder_abs',
        kallistotop500graphtransdecoder='\$kallistotop500graphtransdecoder_abs',
        busco_figure_transdecoder='\$busco_figure_transdecoder_abs',
        alluvial3='\$alluvial3_abs',
        alluvial4='\$alluvial4_abs',
        pie5='\$pie5_abs',
        pie6='\$pie6_abs',
        pie7='\$pie7_abs',
        pie8='\$pie8_abs',
        topkallisto_transdecoder='\$topkallisto_transdecoder_abs',
        name = '${params.name}'
      )
    )"
    """
}

// Process 25:
process RmarkdownL {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (alluvial5),
    path (alluvial6),
    path (pie9),
    path (pie10),
    path (pie11),
    path (pie12),
    path (topkallisto_signalp),

    output:
    path "*.html"

    script:
    """

    alluvial5_abs=\$(readlink -f "${alluvial5}")
    alluvial6_abs=\$(readlink -f "${alluvial6}")
    pie9_abs=\$(readlink -f "${pie9}")
    pie10_abs=\$(readlink -f "${pie10}")
    pie11_abs=\$(readlink -f "${pie11}")
    pie12_abs=\$(readlink -f "${pie12}")
    topkallisto_signalp_abs=\$(readlink -f "${topkallisto_signalp}")

    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/L.Rmd',
      output_dir='.',
      params=list(
        alluvial5='\$alluvial5_abs',
        alluvial6='\$alluvial6_abs',
        pie9='\$pie9_abs',
        pie10='\$pie10_abs',
        pie11='\$pie11_abs',
        pie12='\$pie12_abs',
        topkallisto_signalp='\$topkallisto_signalp_abs',
        name = '${params.name}'
      )
    )"
    """
}

// Process 26:
process RmarkdownM {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Table1)

    output:
    path "*.html"

    script:
    """

    Table1_abs=\$(readlink -f "${Table1}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/M.Rmd',
      output_dir='.',
      params=list(
        Table1='\$Table1_abs',
        name = '${params.name}'
      )
    )"
    """
}

// Process 27:
process RmarkdownN {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Table2), path (Table3)

    output:
    path "*.html"

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
        name = '${params.name}'
      )
    )"
    """
}

// Process 28:
process RmarkdownO {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Table4)

    output:
    path "*.html"

    script:
    """

    Table4_abs=\$(readlink -f "${Table4}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/O.Rmd',
      output_dir='.',
      params=list(
        Table4='\$Table4_abs',
        name = '${params.name}'
      )
    )"
    """
}

// Process 29:
process RmarkdownQ {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Table5)

    output:
    path "*.html"

    script:
    """

    Table5_abs=\$(readlink -f "${Table5}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/Q.Rmd',
      output_dir='.',
      params=list(
        Table5='\$Table5_abs',
        name = '${params.name}'
      )
    )"
    """
}

// Process 30:
process RmarkdownR {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Table6),
    path (Table7)

    output:
    path "*.html"

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
        name = '${params.name}'
      )
    )"
    """
}

// Process 31:
process RmarkdownS {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Table8)

    output:
    path "*.html"

    script:
    """

    Table8_abs=\$(readlink -f "${Table8}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/S.Rmd',
      output_dir='.',
      params=list(
        Table8='\$Table8_abs',
        name = '${params.name}'
      )
    )"
    """
}



// Process 32:
process RmarkdownV {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Table9)

    output:
    path "*.html"

    script:
    """

    Table9_abs=\$(readlink -f "${Table9}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/V.Rmd',
      output_dir='.',
      params=list(
        Table9='\$Table9_abs',
        name = '${params.name}'
      )
    )"
    """
}

// Process 33:
process RmarkdownW {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Table10),
    path (Table11)

    output:
    path "*.html"

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
        name = '${params.name}'
      )
    )"
    """
}

// Process 34:
process RmarkdownX {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Table12)

    output:
    path "*.html"

    script:
    """

    Table12_abs=\$(readlink -f "${Table12}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/X.Rmd',
      output_dir='.',
      params=list(
        Table12='\$Table12_abs',
        name = '${params.name}'
      )
    )"
    """
}

// Process 35:
process RmarkdownZ {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/htmls", mode: 'copy'

    input:
    tuple val(sample), path (Venn),
    path (table)

    output:
    path "*.html"

    script:
    """

    Venn_abs=\$(readlink -f "${Venn}")
    table_abs=\$(readlink -f "${table}")


    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/Z.Rmd',
      output_dir='.',
      params=list(
        VENN='\$Venn_abs',
        TABLE='\$table_abs',
        name = '${params.name}'
      )
    )"
    """
}

// Process 36:
process Blast0Chunks {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/RappData/Alignmentapp", mode: 'copy'

    input:
    tuple val(sample), path (blastx0),
    path (blastp0)

    output:
    path "*", emit: blast0chunks

    script:

    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS11.R" ${params.basename} ${blastx0} ${blastp0}
    """
}



process Blast0Chunksn {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "${sample}/results/RappData/Alignmentapp", mode: 'copy'

    input:
    tuple val(sample), path (blastx0), path (blastp0), path (blastn0)

    output:
    path "*", emit: blast0chunks

    script:

    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS11.R" ${params.basename} ${blastx0} ${blastp0} ${blastn0}
    """
}



// Define input file patterns via parameters



workflow {

    //Define CSV channel.. Channel factory creates the channel from a csv file. can be defined in the config or command line 
    csv_channel = Channel.fromPath(params.input_csv).splitCsv(header: true, sep: ',').map { row -> row.collectEntries { key, value -> [key.replaceAll('"', ''), value?.toString()?.replaceAll('"', '')]}}
    //Define results folder from the Venomflow 
    venomflowfiles = csv_channel.map { row ->
        def samplename = row.Sample_name //0
        def results_path = file(row.Venomflowresultsfolder).toAbsolutePath() 
        def kallisto_trinity = file("${results_path}/kallisto/trinity/output/abundance.tsv") //1
        def kallisto_trans = file("${results_path}/kallisto/transdecoder/output/abundance.tsv") //2
        def combined_pep = file("${results_path}/ORFprediction/Combined/*combined.deduplicated.pep").first() //3
        def combined_cds = file("${results_path}/ORFprediction/Combined/*combined.deduplicated.cds").first() //4
        def mature_fasta = file("${results_path}/Signalp/*_mature.fasta").first() //5
        def blastx_files = file("${results_path}/Blast/Blastx/*.blastx.db.6.txt").first() //6
        def blastp_files = file("${results_path}/Blast/Blastp/*.blastp.db.6.txt").first() //7
        def Interproscan_file = file("${results_path}/Interproscan/*.cleaned.pep.tsv").first() //8
        def signalp_summary = file("${results_path}/Signalp/*_summary.signalp5").first() //9
        def Blastn6 = file("${results_path}/Blast/Blastn/*.blastn.db.6.txt").first() //10
        def blastn0txt = file("${results_path}/Blast/Blastn/*.blastn.db.0.txt").first() //11
        def blastx0txt = file("${results_path}/Blast/Blastx/*.blastx.db.0.txt").first() //12
        def blastp0txt = file("${results_path}/Blast/Blastp/*.blastp.db.0.txt").first() //13
        def basename = row.basename //14
        def busco_transcriptome_dir = file("${results_path}/BUSCO/transcriptome/Transcriptome1") //15
        def busco_translatome_dir = file("${results_path}/BUSCO/translatome/Transdecoder/") //16
        def genomeid = row.NCBI_Genome_id //17
        def species = row.Species //18 
        def MS = file(row.massspec_csv) //19
        def Gavailability = row.isgenomeavailble //20
        def MSavailability = row.ismassspecavailable //21
        def ToxinDomains = file(row.Toxin_domains) //22
        def busco_transcriptome_dir2 = file("${results_path}/BUSCO/transcriptome/Transcriptome2") //23
        def busco_transcriptome_dir3 = file("${results_path}/BUSCO/transcriptome/Combined") //24
        def busco_translatome_dir2 = file("${results_path}/BUSCO/translatome/TD2/") //25
        def busco_translatome_dir3 = file("${results_path}/BUSCO/translatome/Combined/") //26
        def Transcriptome1 = file(row.Transcriptome1) //27
        def Transcriptome2 = file(row.Transcriptome2) //28
        def TranscriptomeC = file("${results_path}/Transcriptome/*transcriptome_combined.deduplicated.fasta").first() //29
        def complete_pep = file("${results_path}/ORFprediction/Combined/*combine.complete.pep").first() //30
        def complete_cds = file("${results_path}/ORFprediction/Combined/*combine.complete.cds").first() //31
       
        return [row.sample_id, kallisto_trinity, kallisto_trans, combined_pep, combined_cds, mature_fasta, blastx_files, blastp_files, Interproscan_file, signalp_summary, Blastn6, blastn0txt, blastx0txt, blastp0txt, basename, busco_transcriptome_dir, busco_translatome_dir,genomeid,species, MS, Gavailability, MSavailability,ToxinDomains,busco_transcriptome_dir2, busco_transcriptome_dir3, busco_translatome_dir2, busco_translatome_dir3,Transcriptome1,Transcriptome2,TranscriptomeC,complete_pep,complete_cds]
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

    // Metadata HTML (B)
    MetadataInput | RmarkdownB

    def isEmpty(value) {
    return !value || 
           value.toString().trim() == "" || 
           value.toString().trim().toLowerCase() == "null"
}


    // Process 1: Define Inputs for kallistoanalysistrinity  sample id + kallisto_file trinity tuple 
    kallistoanalysistrinityinput = venomflowfiles.map{ [it[0], it[1]]}
    //Run process kallistoAnalysisTrinity
    kallistoanalysistrinityinput | kallistoAnalysisTrinity


    //Process 2: Define Inputs for kallistoanalysistrans sample id + kallisto_file trans tuple 
    kallistoanalysistransinput = venomflowfiles.map{ [it[0], it[2]]}
    //Run process kallistoanalysistrans
    kallistoanalysistransinput | kallistoAnalysisTrans



    //Process 3: Define Inputs for ExtractSignalSequences sample id + transdecoderpep + maturefasta tuple 
    ExtractSignalSequencesinput = venomflowfiles.map{ [it[0], it[3], it[5]]}
    //Run process ExtractSignalSequences
    ExtractSignalSequencesinput | ExtractSignalSequences


    //Process 4: Define Inputs for ExtractCreateTrinityDataframe  sample + trinity_fasta + blastx_file + kallisto_csv  tuple 
    
    CreateTrinityDataframeinput_single = venomflowfiles.filter{isEmpty(it[28])}
    .map{[it[0], it[27], it[6]]}.join(kallistoAnalysisTrinity.out.trin_all_csv)
    
    CreateTrinityDataframeinput_combined = venomflowfiles.filter{!isEmpty(it[28])}
    .map{[it[0], it[28], it[6]]}.join(kallistoAnalysisTrinity.out.trin_all_csv)

  CreateTrinityDataframeinput = CreateTrinityDataframeinput_single.mix(CreateTrinityDataframeinput_combined)
    //Run process ExtractCreateTrinityDataframe
    CreateTrinityDataframeinput | CreateTrinityDataframe


    //Process 5: Define Inputs for CreateInterproscanDataframe  sample + Interproscan + ListFile + PantherFile  tuple 
    def ListFile = Channel.fromPath(params.input_interproscan)
    def PantherFile = Channel.fromPath(params.input_panther)
    CreateInterproscanDataframeinput = venomflowfiles.map{return [it[0], it[8]]}.combine([ListFile, PantherFile])
    //Run process CreateInterproscanDataframe
    CreateInterproscanDataframeinput | CreateInterproscanDataframe


    //Process 6: Define Inputs for CreateTransdecoderDataframe  sample + transdecoder_pep + transdecoder_cds + blastp6_file  mature_fasta, Signalp_summary, signalsequences, Interproscan_dataframe, kallistotrans
    CreateTransdecoderDataframeinput = venomflowfiles.map {return [it[0], it[3], it[4], t[7], it[5], it[9]]}
                                      .join(ExtractSignalSequences.out.signalsequences)
                                       .join(CreateInterproscanDataframe.out.Interproscan_dataframe)
                                       .join(kallistoAnalysisTrans.out.trans_all_csv)
    //Run process ExtractCreateTrinityDataframe
    CreateTransdecoderDataframeinput | CreateTransdecoderDataframe



    //Process 7: Define Inputs for CreateInterproscanToxinPlotly sample + transdf_distinct_csv  + transtoxinfasta + toxindomains
    def toxindomains = Channel.fromPath(params.input_toxin_domains)
    CreateInterproscanToxinPlotlyinput = CreateTransdecoderDataframe.out.transdf_distinct.join(CreateTransdecoderDataframe.out.toxin_fastacds).join(CreateTransdecoderDataframe.out.toxin_fastapep).combine([toxindomains])
    // Run process CreateInterproscanToxinPlotly
    CreateInterproscanToxinPlotlyinput | CreateInterproscanToxinPlotly

    //Process 8: Define Input BUSCOtranscriptome
    BUSCOtranscriptomeinput_1 = venomflowfiles.map {[it[0], it[15], "1"]}
    BUSCOtranscriptomeinput_2 = venomflowfiles.filter{!isEmpty(it[28])}.map {[it[0], it[23], "2"]}
    BUSCOtranscriptomeinput_C = venomflowfiles.filter{!isEmpty(it[28])}.map {[it[0], it[24], "C"]}
    BUSCOtranscriptomeinput = BUSCOtranscriptomeinput_1.mix(BUSCOtranscriptomeinput_2).mix(BUSCOtranscriptomeinput_C)
    //Run Process BUSCOtranscriptome
    BUSCOtranscriptomeinput | BUSCOtranscriptome

    //Process 9: Define Input BUSCOtranslatome
    BUSCOtranslatomeinput_1 = venomflowfiles.map {[it[0], it[16], "1"]}
    BUSCOtranslatomeinput_2 = venomflowfiles.filter{!isEmpty(it[28])}.map {[it[0], it[25], "2"]}
    BUSCOtranslatomeinput_C = venomflowfiles.filter{!isEmpty(it[28])}.map {[it[0], it[26], "C"]}
    BUSCOtranslatomeinput = BUSCOtranscriptomeinput_1.mix(BUSCOtranscriptomeinput_2).mix(BUSCOtranscriptomeinput_C)

    //Run Process BUSCOtranslatome
    BUSCOtranslatomeinput | BUSCOtranslatome

    groupedbuscotranscriptome = BUSCOtranscriptome.out.busco_transcriptome.groupTuple()
    groupedbuscotranslatome = BUSCOtranslatomeinput.out.busco_transcriptome.groupTuple()

    //Process 10: Define Input for TableGenerationTrinity sample+TBK + genomeid + species 
    def genome = venomflowfiles.map{return [it[0], it[17]]}
    def species = venomflowfiles.map{return [it[0], it[18], ]}
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
    SampleGenomeSpecies = venomflowfiles.map { return [it[0], it[17], it[18]] }
  

    TableGenerationTransdecoderInput = CreateTransdecoderDataframe.out.transdf.join(SampleGenomeSpecies)
    //Run Process TableGenerationTrinity
    TableGenerationTransdecoderInput | TableGenerationTransdecoder

      //Input 

      AddMSGenomeIfAvailableAndCreateOverviewInput_MG = venomflowfiles.filter { it[20] == 'Y' && it[21] == 'Y' }
                                                        .map([it[0], it[18], it[19],it[10] ])
                                                        .join(CreateTransdecoderDataframe.out.transdf_distinct)
                                                        .combine(toxindomains)
      AddMSGenomeIfAvailableAndCreateOverviewInput_M = venomflowfiles.filter { it[20] == 'N' && it[21] == 'Y' }
                                                        .map([it[0], it[18], it[19], "NULL" ])
                                                        .join(CreateTransdecoderDataframe.out.transdf_distinct)
                                                        .combine(toxindomains)
      AddMSGenomeIfAvailableAndCreateOverviewInput_G = venomflowfiles.filter { it[20] == 'Y' && it[21] == 'N' }
                                                         .map([it[0], it[18], "NULL" ,it[10] ])
                                                        .join(CreateTransdecoderDataframe.out.transdf_distinct)
                                                        .combine(toxindomains)
      AddMSGenomeIfAvailableAndCreateOverviewInput_N = venomflowfiles.filter { it[20] == 'N' && it[21] == 'N' }
                                                        .map([it[0], it[18], "NULL", "NULL" ])
                                                        .join(CreateTransdecoderDataframe.out.transdf_distinct)
                                                        .combine(toxindomains)

    
      AddMSGenomeIfAvailableAndCreateOverviewInput_MG = venomflowfiles.filter { it[20] == 'Y' && it[21] == 'Y' }
    
    AddMSGenomeIfAvailableAndCreateOverviewInput.mix(AddMSGenomeIfAvailableAndCreateOverviewInput_M).mix(AddMSGenomeIfAvailableAndCreateOverviewInput_G).mix(AddMSGenomeIfAvailableAndCreateOverviewInput_N)
    // Run process AddMSGenomeIfAvailableAndCreateOverview
    AddMSGenomeIfAvailableAndCreateOverviewInputAddMSGenomeIfAvailableAndCreateOverview

  /*
    //RmarkdownCDEGIK()


    //RmarkdownH(
        kallistoAnalysisTrinity.out.trin_top20_png,
        kallistoAnalysisTrinity.out.trin_top500_png,
        BUSCOtranscriptome.out.busco_transcriptome,
        FigureGenerationTrinity.out.alluvial1,
        FigureGenerationTrinity.out.alluvial2,
        FigureGenerationTrinity.out.pie1,
        FigureGenerationTrinity.out.pie2,
        FigureGenerationTrinity.out.pie3,
        FigureGenerationTrinity.out.pie4,
        FigureGenerationTrinity.out.Table13,
    //)

    //RmarkdownJ(
        kallistoAnalysisTrans.out.trans_top20_png,
        kallistoAnalysisTrans.out.trans_top500_png,
        BUSCOtranslatome.out.busco_translatome,
        FigureGenerationTransdecoder.out.alluvial3,
        FigureGenerationTransdecoder.out.alluvial4,
        FigureGenerationTransdecoder.out.pie5,
        FigureGenerationTransdecoder.out.pie6,
        FigureGenerationTransdecoder.out.pie7,
        FigureGenerationTransdecoder.out.pie8,
        FigureGenerationTransdecoder.out.Table14,
    //)

   // RmarkdownL(
        FigureGenerationSignalp.out.alluvial5,
        FigureGenerationSignalp.out.alluvial6,
        FigureGenerationSignalp.out.pie9,
        FigureGenerationSignalp.out.pie10,
        FigureGenerationSignalp.out.pie11,
        FigureGenerationSignalp.out.pie12,
        FigureGenerationSignalp.out.Table15,
    //)

    //RmarkdownM(
        TableGenerationTrinity.out.Table1
    //)

    //RmarkdownN(
        TableGenerationTrinity.out.Table2,
        TableGenerationTrinity.out.Table3,
    //)

    //RmarkdownO(
        TableGenerationTrinity.out.Table4
   // )

    //RmarkdownQ(
        TableGenerationTransdecoder.out.Table5
    //)

    //RmarkdownR(
        TableGenerationTransdecoder.out.Table6,
        TableGenerationTransdecoder.out.Table7,
   //)

   // RmarkdownS(
        TableGenerationTransdecoder.out.Table8
    //)

    RmarkdownV(
        TableGenerationTransdecoder.out.Table9
    )

    RmarkdownW(
        TableGenerationTransdecoder.out.Table10,
        TableGenerationTransdecoder.out.Table11,
    )

    RmarkdownX(
        TableGenerationTransdecoder.out.Table12
    )

    if (params.isgenomeavailable == 'Y') {


        def blastn0txt = Channel.fromPath(params.input_blastn0_files)
        def blastx0txt = Channel.fromPath(params.input_blastx0_files)
        def blastp0txt = Channel.fromPath(params.input_blastp0_files)
        Blast0Chunksn(blastx0txt, blastp0txt, blastn0txt)
    }
    else {
        def blastx0txt = Channel.fromPath(params.input_blastx0_files)
        def blastp0txt = Channel.fromPath(params.input_blastp0_files)
        Blast0Chunks(blastx0txt, blastp0txt)
    }
}
row.Sample_name	row.Trinity_fasta			BUSCO_lin1	BUSCO_lin2	Protein_fasta_path_for_Blast	Protein_fasta_name	Strandedness	Genome_fasta_path	Genome_fasta_name	Author_Name	NCBI_Genome_id	Species	isgenomeavailble	ismassspecavailable	basename	SearchAndDownloadURL	massspec_csv	SRA	Access	Published	Age	Animal Length/Size (cm)	Tissue	Tissue mass (mg)	Sex	Source	Bred/Wildcaught	Date Collected	Collected by	Storage conditions  before extraction 	Homogenization method	RNA extraction Method 	Concentration(ng/ul)	Extraction date	Extracted by	Submission for QC/Seq Date	QC/Seq by 	Submission Order Number	RIN	Library Prep	Selection	Strandedness	Library Layout	Platform	ReadLength	Results Received Date	Preprocessing  by 	Trinity Assembly  by 	Venomflow by 	VenomflowAnalysis by 	Related Publications	Lab 	Contact email 	Misc Notes

*//