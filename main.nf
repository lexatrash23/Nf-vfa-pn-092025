#!/usr/env nextflow

// Process 1: Creating the settings, notes and triminfo files
process MetadataCreation {
    errorStrategy 'ignore'

    publishDir "results/metadata", mode: 'copy'

    input:
    tuple val(notes), val(samplename), val(SRA), val(Tissue), val(Age), val(Sex), val(RNAExtractionMethod), val(LibraryPrep), val(Strandedness), val(Selection), val(Readlength)

    output:
    path "notes.txt", emit: notes
    path "settings.csv", emit: settings

    script:
    """
    #notes
    cat <<EOF > "notes.txt"
${notes}
EOF
   
    
    #settings

    # Write CSV header
    echo "Field,Value" > "settings.csv"

    # Write key-value rows
    echo "Sample Name,${samplename}" >> "settings.csv"
    echo "SRA,${SRA}" >> "settings.csv"
    echo "Tissue,${Tissue}" >> "settings.csv"
    echo "Age,${Age}" >> "settings.csv"
    echo "Sex,${Sex}" >> "settings.csv"
    echo "RNA Extraction Method,${RNAExtractionMethod}" >> "settings.csv"
    echo "Library Prep,${LibraryPrep}" >> "settings.csv"
    echo "Strandedness,${Strandedness}" >> "settings.csv"
    echo "Selection,${Selection}" >> "settings.csv"
    echo "Read Length,${Readlength}" >> "settings.csv"
    
    
    """
}

// Process 2: For kallistoanalysistrinity.py python,pandas,seaborn,matplotlib
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

// Process 3: For kallistoanalysistrans.py dependencies:python,pandas,seaborn,matplotlib
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

// Process 4: Extract Signal Sequences dependencies python biopython
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
    
    python3 ${workflow.projectDir}/bin/Intermediate_Scripts/IS2.py ${transdecoder_pep} ${mature_fasta} signalsequences.fasta
    """
}

// Process 5: Create Trinity Dataframe dependecies : R, biocmanager 
process CreateTrinityDataframe {
    errorStrategy 'ignore'

    conda 'r-base bioconductor-biostrings r-tidyr r-dplyr'

    publishDir "${sample}/VenomFlowAnalysis/results/TrinityDataframe/", mode: 'copy'
    publishDir "${sample}/VenomFlowAnalysis/results/RappData/TrinityDataframe/", pattern: "*.gz", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), path(blastx_file), path(kallisto_csv)

    output:
    tuple val(sample), path ("_TBK.csv"), emit: TBK
    path "_TBK_distinct.csv.gz", emit: TBK_distinct

    script:
    """
    
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS1.R ${trinity_fasta} ${blastx_file} ${kallisto_csv}
    """
}

// Process 6: Create Interproscan Dataframe dependecies : R, biocmanager 
process CreateInterproscanDataframe {
    errorStrategy 'ignore'

    conda 'r-base bioconductor-biostrings r-dplyr bioconductor-go.db bioconductor-biomart r-tidyr'

    publishDir "${sample}/VenomFlowAnalysis/results/InterproscanDataframe/", mode: 'copy'

    input:
    tuple val(sample), path(Interproscan), path(ListFile), path(PantherFile)

    output:
    tuple val(sample), path ("Final_interproscan_dataframe.csv"), emit: Interproscan_dataframe

    script:
    """
    
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS4.R" ${Interproscan} ${ListFile} ${PantherFile}
    """
}

// Process 7: Create Transdecoder Dataframe dependecies : R, biocmanager 
process CreateTransdecoderDataframe {
    errorStrategy 'ignore'

    conda 'r-base=4.3 r-dplyr r-tidyr bioconductor-biostrings'


    publishDir "${sample}/VenomFlowAnalysis/results/TransdecoderDataframe/", mode: 'copy'
    publishDir "${sample}/VenomFlowAnalysis/results/Secreted_proteins_fasta/", pattern: "*.fasta", mode: 'copy'

    input:
    tuple val(sample), path(transdecoder_pep), path(transdecoder_cds), path(blastp_file), path(mature_fasta), path(Signalp_summary), path(signalsequences), val(basename), path(Interproscan_dataframe), path(kallistotrans)

    output:
    tuple val(sample), path ("_transdf.csv"), emit: transdf
    tuple val(sample), path ("_transdf_distinct.csv"), emit: transdf_distinct
    tuple val(sample), path ("secreted_proteins.fasta"), emit: toxin_fasta

    script:
    """
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS5.R ${transdecoder_pep} ${transdecoder_cds} ${blastp_file} ${mature_fasta} ${Signalp_summary} ${signalsequences} ${Interproscan_dataframe} ${kallistotrans} ${basename}
    """
}

// Process 8: Create Interproscantoxinplotly
process CreateInterproscanToxinPlotly {
    errorStrategy 'ignore'

    conda 'bioconductor-biostrings r-dplyr r-tidyr r-htmlwidgets r-plotly'


    publishDir "${sample}/VenomFlowAnalysis/results/htmls/", pattern: "*.html", mode: 'copy'
    publishDir "${sample}/VenomFlowAnalysis/results/Secreted_proteins_fasta/", pattern: "*.fasta", mode: 'copy'

    input:
    tuple val(sample), path(transdf_distinct_csv), path(transtoxinfasta), path(Toxin_domains)

    output:
    tuple val(sample), path ("filtered_sequences.fasta"), emit: filtered_sequences
    tuple val(sample), path ("plotly_graph.html"), emit: plotly_graph

    script:
    """
    
    Rscript ${workflow.projectDir}/bin/Intermediate_Scripts/IS7.R ${transdf_distinct_csv} ${Toxin_domains} ${transtoxinfasta}
    """
}

// Process 9: Create BUSCOgraphtranscriptome  
process BUSCOtranscriptome {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco2.yaml"

    publishDir "${sample}/VenomFlowAnalysis/results/busco/transcriptome/", mode: 'copy'

    input:
    tuple val (sample), path (buscodirtr)

    output:
    uple val(sample), (path "*.png"), emit: busco_transcriptome

    script:
    """
     
    python3 "${workflow.projectDir}/bin/Intermediate_Scripts/IS6.py" -wd "${buscodirtr}"
    """
}

// Process 10: Create BUSCOgraphtranslatome 
process BUSCOtranslatome {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco2.yaml"

    publishDir "${sample}/VenomFlowAnalysis/results/busco/translatome", mode: 'copy'

    input:
    val (sample), path (buscodirtl)

    output:
    tuple val(sample), path ("*.png"), emit: busco_translatome

    script:
    """
     
    python3 "${workflow.projectDir}/bin/Intermediate_Scripts/IS6.py" -wd "${buscodirtl}"

    """
}

// Process 11: Create TableGenerationTrinity
process TableGenerationTrinity {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts2_outputs", mode: 'copy'

    input:
    tuple path(TBK), val(genome_id), val(species)

    output:
    path "Table1.csv", emit: Table1
    path "Table2.csv", emit: Table2
    path "Table3.csv", emit: Table3
    path "Table4.csv", emit: Table4

    script:
    """
	
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Generating_TopTables_Trinity.R" ${TBK} ${genome_id} ${species}
    """
}

// Process 12: Create TableGenerationTransdecoder  
process TableGenerationTransdecoder {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts2_outputs", mode: 'copy'

    input:
    tuple path(transdf), val(genome_id), val(species)

    output:
    path "Table5.csv", emit: Table5
    path "Table6.csv", emit: Table6
    path "Table7.csv", emit: Table7
    path "Table8.csv", emit: Table8
    path "Table9.csv", emit: Table9
    path "Table10.csv", emit: Table10
    path "Table11.csv", emit: Table11
    path "Table12.csv", emit: Table12

    script:
    """
	
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Generating_Tables_Transdecoder_SignalP.R" ${transdf} ${genome_id} ${species}
    """
}


// Process 13: Create FigureGenerationTrinity
process FigureGenerationTrinity {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts2_outputs", mode: 'copy'

    input:
    path TBK

    output:
    path "pie1.png", emit: pie1
    path "pie2.png", emit: pie2
    path "pie3.png", emit: pie3
    path "pie4.png", emit: pie4
    path "alluvial1.png", emit: alluvial1
    path "alluvial2.png", emit: alluvial2
    path "Table13.csv", emit: Table13

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Figure_generation_Trinity.R" ${TBK} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"

    """
}

// Process 14: Create FigureGenerationTransdecoder
process FigureGenerationTransdecoder {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts2_outputs", mode: 'copy'

    input:
    path transdf

    output:
    path "pie5.png", emit: pie5
    path "pie6.png", emit: pie6
    path "pie7.png", emit: pie7
    path "pie8.png", emit: pie8
    path "alluvial3.png", emit: alluvial3
    path "alluvial4.png", emit: alluvial4
    path "Table14.csv", emit: Table14

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Figure_generation_Transdecoder.R" ${transdf} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"

    """
}

// Process 15: Create FigureGenerationSignalp
process FigureGenerationSignalp {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts2_outputs", mode: 'copy'

    input:
    path transdf

    output:
    path "pie9.png", emit: pie9
    path "pie10.png", emit: pie10
    path "pie11.png", emit: pie11
    path "pie12.png", emit: pie12
    path "alluvial5.png", emit: alluvial5
    path "alluvial6.png", emit: alluvial6
    path "Table15.csv", emit: Table15

    script:
    """
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts2/Figure_generation_SignalP.R" ${transdf} "${workflow.projectDir}/bin/Intermediate_Scripts2/color_palette.rds"

    """
}

// Process 16: Create AddMassSpec
process AddMassSpec {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/RappData/Single", pattern: "*.gz", mode: 'copy'
    publishDir "results/RappData/Combined", pattern: "*.csv", mode: 'copy'

    input:
    tuple path(transdf), path(massspecdata), val(species), val(basename)

    output:
    path "*_filtered_masspec_select.csv", emit: filtered_massspec
    path "*_distinct_masspec.csv.gz", emit: distinct_massspec

    script:
    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS8.R" ${transdf} ${massspecdata} "${species}" ${basename}

    """
}

// Process 17: Create SkipMassSpec
process SkipMassSpec {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/RappData/Single", pattern: "*.gz", mode: 'copy'
    publishDir "results/RappData/Combined", pattern: "*.csv", mode: 'copy'

    input:
    tuple path(transdf), val(species), val(basename)

    output:
    path "*_filtered_nomasspec_.csv", emit: filtered_nomasspec
    path "*_distinct_nomasspec.csv.gz", emit: distinct_nomasspec

    script:
    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS9.R" ${transdf} "${species}" ${basename}

    """
}

// Process 18: Create VennOverview
process VennOverviewMS {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(filtered_massspec), path(Toxin_domains)

    output:
    path "*.csv", emit: overviewcsv
    path "*.png", emit: overviewpng

    script:
    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS10.R" ${filtered_massspec} "${Toxin_domains}"

    """
}

// Process 19: Create VennOverview
process VennOverviewNoMS {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(filtered_nomasspec), path(Toxin_domains)

    output:
    path "*.csv", emit: overviewcsv
    path "*.png", emit: overviewpng

    script:
    """
    
    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS10.R" ${filtered_nomasspec} "${Toxin_domains}"
    
    """
}

// Process 20: RmarkdownA
process RmarkdownA {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    val sampleURL

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path settings
    path notes

    output:
    path "*.html"

    script:
    """

    settings_abs=\$(readlink -f "${settings}")
    notes_abs=\$(readlink -f "${notes}")

    Rscript -e "rmarkdown::render(
      '${workflow.projectDir}/bin/Rmarkdown_scripts/B.Rmd',
      output_dir='.',
      params=list(
        settings='\$settings_abs',
        notes='\$notes_abs',
        name = '${params.name}'
      )
      )"
    """
}

// Process 22:
process RmarkdownCDEGIK {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path kallistotop20graphtrinity
    path kallistotop500graphtrinity
    path busco_figure
    path alluvial1
    path alluvial2
    path pie1
    path pie2
    path pie3
    path pie4
    path topkallisto

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path kallistotop20graphtransdecoder
    path kallistotop500graphtransdecoder
    path busco_figure_transdecoder
    path alluvial3
    path alluvial4
    path pie5
    path pie6
    path pie7
    path pie8
    path topkallisto_transdecoder

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path alluvial5
    path alluvial6
    path pie9
    path pie10
    path pie11
    path pie12
    path topkallisto_signalp

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Table1

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Table2
    path Table3

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Table4

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Table5

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Table6
    path Table7

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Table8

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Table9

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Table10
    path Table11

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Table12

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

    publishDir "results/htmls", mode: 'copy'

    input:
    path Venn
    path table

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

    publishDir "results/RappData/Alignmentapp", mode: 'copy'

    input:
    path blastx0
    path blastp0

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

    publishDir "results/RappData/Alignmentapp", mode: 'copy'

    input:
    path blastx0
    path blastp0
    path blastn0

    output:
    path "*", emit: blast0chunks

    script:

    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS11.R" ${params.basename} ${blastx0} ${blastp0} ${blastn0}
    """
}

process BlastnIntegration {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    path venndiagram
    path trandfdistinctmass
    path blastn6

    output:
    path "transdf_distinct_blastn.csv", emit: transdf_distinct_blastn
    path "venn_overview_blastn.csv", emit: venn_overview_blastn
    path "venn_overview_blastn_filtered.csv", emit: venn_overview_blastn_filtered

    script:

    """

    Rscript "${workflow.projectDir}/bin/Intermediate_Scripts/IS12.R" ${venndiagram} ${trandfdistinctmass} ${blastn6}
    """
}

// Define input file patterns via parameters



workflow {

    //Define CSV channel.. Channel factory creates the channel from a csv file. can be defined in the config or command line 
    csv_channel = Channel.fromPath(params.input_csv).splitCsv(header: true, sep: ',')
    //Define results folder from the Venomflow 
    venomflowfiles = csv_channel.map { row ->
        def samplename = row.Sample_name //0
        def results_path = file(row.Venomflowresultsfolder).toAbsolutePath() 
        def kallisto_trinity = file("${results_path}/kallisto/trinity/output/abundance.tsv") //1
        def kallisto_trans = file("${results_path}/kallisto/transdecoder/output/abundance.tsv") //2
        def transdecoder_pep = file("${results_path}/Transdecoder/*.transdecoder.pep").first() //3
        def transdecoder_cds = file("${results_path}/Transdecoder/*.transdecoder.cds").first() //4
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
        def busco_transcriptome_dir = file("${results_path}/BUSCO/transcriptome/") //15
        def busco_translatome_dir = file("${results_path}/BUSCO/translatome/") //16
        return [row.sample_id, kallisto_trinity, kallisto_trans, transdecoder_pep, transdecoder_cds, mature_fasta, blastx_files, blastp_files, Interproscan, signalp_summary, Blastn6, blastn0txt, blastx0txt, blastp0txt, basename, busco_transcriptome_dir, busco_translatome_dir]
    }



    //Define Inputs for kallistoanalysistrinity  sample id + kallisto_file trinity tuple 
    kallistoanalysistrinityinput = venomflowfiles.map{return [it[0], it[1]]}
    //Run process kallistoAnalysisTrinity
    kallistoanalysistrinityinput | kallistoAnalysisTrinity


    //Define Inputs for kallistoanalysistrans sample id + kallisto_file trans tuple 
    kallistoanalysistransinput = venomflowfiles.map{return [it[0], it[2]]}
    //Run process kallistoanalysistrans
    kallistoanalysistransinput | kallistoAnalysisTrans



    //Define Inputs for ExtractSignalSequences sample id + transdecoderpep + maturefasta tuple 
    ExtractSignalSequencesinput = venomflowfiles.map{return [it[0], it[3], it[5]]}
    //Run process ExtractSignalSequences
    ExtractSignalSequencesinput | ExtractSignalSequences


    //Define Inputs for ExtractCreateTrinityDataframe  sample + trinity_fasta + blastx_file + kallisto_csv  tuple 
    ExtractSignalSequencesinput = venomflowfiles.map{return [it[0], it[3], it[5]]}.join(kallistoAnalysisTrinity.out.trin_all_csv)
    //Run process ExtractCreateTrinityDataframe
    ExtractSignalSequencesinput | ExtractSignalSequences


    //Define Inputs for CreateInterproscanDataframe  sample + Interproscan + ListFile + PantherFile  tuple 
    def ListFile = Channel.fromPath(params.input_list)
    def PantherFile = Channel.fromPath(params.input_panther)
    CreateInterproscanDataframeinput = venomflowfiles.map{return [it[0], it[8]]}.combine([ListFile, PantherFile])
    //Run process CreateInterproscanDataframe
    CreateInterproscanDataframeinput | CreateInterproscanDataframe

    //Define Inputs for CreateTransdecoderDataframe  sample + transdecoder_pep + transdecoder_cds + blastp6_file  mature_fasta, Signalp_summary, signalsequences, Interproscan_dataframe, kallistotrans
    CreateTransdecoderDataframeinput = venomflowfiles.map {return [it[0], it[3], it[4], t[7], it[6], it[9], it[14]]}
                                       .join(CreateInterproscanDataframe.out.Interproscan_dataframe)
                                       .join(kallistoAnalysisTrans.out.trans_all_csv)
    //Run process ExtractCreateTrinityDataframe
    CreateTransdecoderDataframeinput | CreateTransdecoderDataframe



    //Define Inputs for CreateInterproscanToxinPlotly sample + transdf_distinct_csv  + transtoxinfasta + toxindomains
    def toxindomains = Channel.fromPath(params.input_toxindomains)
    CreateInterproscanToxinPlotlyinput = CreateTransdecoderDataframe.out.transdf_distinct.join(CreateTransdecoderDataframe.out.toxin_fasta).combine([toxindomains])
    // Run process CreateInterproscanToxinPlotly
    CreateInterproscanToxinPlotlyinput | CreateInterproscanToxinPlotly

    //Define Input BUSCOtranscriptome
    BUSCOtranscriptomeinput = venomflowfiles.map {return [it[0], it[15]]}
    //Run Process BUSCOtranscriptome
    BUSCOtranscriptomeinput | BUSCOtranscriptome

    //Define Input BUSCOtranslatome
    BUSCOtranslatomeinput = venomflowfiles.map {return [it[0], it[16]]}
    //Run Process BUSCOtranslatome
    BUSCOtranslatomeinput | BUSCOtranslatome



    InterproscanToxinData | CreateInterproscanToxinPlotly
    BUSCOtranscriptome()
    BUSCOtranslatome()

    def genome_id = Channel.value(params.genome_id)
    def species = Channel.value(params.species)
    def TBK = CreateTrinityDataframe.out.TBK

    def TBK_table_data = TBK.combine(genome_id).combine(species).view()

    TBK_table_data | TableGenerationTrinity

    def transdf = CreateTransdecoderDataframe.out.transdf

    def Trans_table_data = transdf
        .combine(genome_id)
        .combine(species)

    Trans_table_data | TableGenerationTransdecoder
    def TBK_only = CreateTrinityDataframe.out.TBK

    TBK_only | FigureGenerationTrinity
    transdf | FigureGenerationTransdecoder
    transdf | FigureGenerationSignalp
    def basename = Channel.value(params.basename)

    if (params.ismassspecavailable == 'Y') {

        def massspecdata = Channel.fromPath(params.massspecdata)

        def massspecdatatable = transdf
            .combine(massspecdata)
            .combine(species)
            .combine(basename)

        massspecdatatable | AddMassSpec
        def filtered_massspec = AddMassSpec.out.filtered_massspec
        def overviewms = filtered_massspec.combine(Toxin_domains)
        overviewms | VennOverviewMS
        RmarkdownZ(VennOverviewMS.out.overviewpng, VennOverviewMS.out.overviewcsv)
        if (params.isgenomeavailable == 'Y') {
            def Blastn6 = Channel.fromPath(params.input_blastn_files)
            BlastnIntegration(VennOverviewMS.out.overviewcsv, AddMassSpec.out.distinct_massspec, Blastn6)
        }
    }
    else {

        def nomassspecdatatable = transdf
            .combine(species)
            .combine(basename)

        nomassspecdatatable | SkipMassSpec
        def filtered_nomasspec = SkipMassSpec.out.filtered_nomasspec
        def overviewnoms = filtered_nomasspec.combine(Toxin_domains)
        overviewnoms | VennOverviewNoMS
        RmarkdownZ(VennOverviewNoMS.out.overviewpng, VennOverviewNoMS.out.overviewcsv)
        if (params.isgenomeavailable == 'Y') {
            def Blastn6 = Channel.fromPath(params.input_blastn_files)
            BlastnIntegration(VennOverviewNoMS.out.overviewcsv, SkipMassSpec.out.distinct_nomasspec, Blastn6)
        }
    }
    def sampleurl = Channel.value(params.sampleURL)
    sampleurl | RmarkdownA

    def N = Channel.value(params.notes)
    def SN = Channel.value(params.samplename)
    def R = Channel.value(params.SRA)
    def Ti = Channel.value(params.Tissue)
    def A = Channel.value(params.Age)
    def SX = Channel.value(params.Sex)
    def REM = Channel.value(params.RNAExtractionMethod)
    def LP = Channel.value(params.LibraryPrep)
    def SD = Channel.value(params.Strandedness)
    def SL = Channel.value(params.Selection)
    def RL = Channel.value(params.Readlength)
    def Metadata = N.combine(SN).combine(R).combine(Ti).combine(A).combine(SX).combine(REM).combine(LP).combine(SD).combine(SL).combine(RL)

    Metadata | MetadataCreation

    RmarkdownB(
        MetadataCreation.out.settings,
        MetadataCreation.out.notes,
    )

    RmarkdownCDEGIK()


    RmarkdownH(
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
    )

    RmarkdownJ(
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
    )

    RmarkdownL(
        FigureGenerationSignalp.out.alluvial5,
        FigureGenerationSignalp.out.alluvial6,
        FigureGenerationSignalp.out.pie9,
        FigureGenerationSignalp.out.pie10,
        FigureGenerationSignalp.out.pie11,
        FigureGenerationSignalp.out.pie12,
        FigureGenerationSignalp.out.Table15,
    )

    RmarkdownM(
        TableGenerationTrinity.out.Table1
    )

    RmarkdownN(
        TableGenerationTrinity.out.Table2,
        TableGenerationTrinity.out.Table3,
    )

    RmarkdownO(
        TableGenerationTrinity.out.Table4
    )

    RmarkdownQ(
        TableGenerationTransdecoder.out.Table5
    )

    RmarkdownR(
        TableGenerationTransdecoder.out.Table6,
        TableGenerationTransdecoder.out.Table7,
    )

    RmarkdownS(
        TableGenerationTransdecoder.out.Table8
    )

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

