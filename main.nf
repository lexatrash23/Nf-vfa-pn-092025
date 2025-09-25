#!/usr/bin/env nextflow

def launchDir = System.getProperty('user.dir')

// Process 0: Creating the settings, notes and triminfo files
process MetadataCreation {

    publishDir "results/metadata", mode: 'copy'

    input:
    tuple val(notes), val(triminfo), val(samplename), val(SRA), val(Tissue), val(Age), val(Sex), val(RNAExtractionMethod), val(LibraryPrep), val(Strandedness), val(Selection), val(Readlength)

    output:
    path "notes.txt", emit: notes
    path "triminfo.txt", emit: triminfo
    path "settings.csv", emit: settings

    script:
    """
    #notes
    cat <<EOF > "notes.txt"
${notes}
EOF
    
    #triminfo

    cat <<EOF > "triminfo.txt"
${triminfo}
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

// Process 1: For kallistoanalysistrinity.py python,pandas,seaborn,matplotlib
process kallistoAnalysisTrinity {

    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    path (kallisto_file_trinity)

    output:
    path "*_all.csv", emit: trin_all_csv
    path "*_top20.csv", emit: trin_top20_csv
    path "*_top500graph.png", emit: trin_top500_png
    path "*_top20graph.png", emit: trin_top20_png

    script:

    """
    
    python3 bin/Intermediate_Scripts/kallistoanalysistrinity.py ./ ${params.basename} ${kallisto_file_trinity}
    """
}

// Process 2: For kallistoanalysistrans.py dependencies:python,pandas,seaborn,matplotlib
process kallistoAnalysisTrans {

    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    path (kallisto_file_transdecoder)

    output:
    path "*_all.csv", emit: trans_all_csv
    path "*_top20.csv", emit: trans_top20_csv
    path "*_top500graph.png", emit: trans_top500_png
    path "*_top20graph.png", emit: trans_top20_png


    script:
    """
    
    python3 bin/Intermediate_Scripts/kallistoanalysistrans.py ./ ${params.basename} ${kallisto_file_transdecoder}
    """
}

// Process 3: Extract Signal Sequences dependencies python biopython
process ExtractSignalSequences {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(transdecoder_pep), path(mature_fasta)

    output:
    path "signalsequences.fasta",  emit: signalsequences

    script: 
    """
    
    python3 bin/Intermediate_Scripts/IS2.py ${transdecoder_pep} ${mature_fasta} signalsequences.fasta
    """

}

// Process 4: Create Trinity Dataframe dependecies : R, biocmanager 
process CreateTrinityDataframe {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(trinity_fasta), path(blastx_file), path(kallisto_csv)

    output:
    path "_TBK.csv", emit: TBK
    path "_TBK_distinct.csv.gz", emit: TBK_distinct

    script: 
    """
    
    Rscript bin/Intermediate_Scripts/IS1.R ${trinity_fasta} ${blastx_file} ${kallisto_csv}
    """

}

// Process 4: Create Interproscan Dataframe dependecies : R, biocmanager 
process CreateInterproscanDataframe {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(Interproscan), path(ListFile), path(PantherFile)

    output:
    path "Final_interproscan_dataframe.csv",  emit: Interproscan_dataframe

    script: 
    """
    
    Rscript "bin/Intermediate_Scripts/IS4.R" ${Interproscan} ${ListFile} ${PantherFile}
    """

}

// Process 5: Create Transdecoder Dataframe dependecies : R, biocmanager 
process CreateTransdecoderDataframe {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(transdecoder_pep), path(transdecoder_cds), path(blastp_file), path(mature_fasta), path(Signalp_summary), path(signalsequences), path(Interproscan_dataframe), path(kallistotrans)

    output:
    path "_transdf.csv", emit: transdf
    path "_transdf_distinct.csv", emit: transdf_distinct   
    path "FINAL_CSV_distinct_filtered_putative_toxins.fasta", emit: toxin_fasta

    script: 
    def basename = transdecoder_pep.getSimpleName()
    """
    
    Rscript bin/Intermediate_Scripts/IS5.R ${transdecoder_pep} ${transdecoder_cds} ${blastp_file} ${mature_fasta} ${Signalp_summary} ${signalsequences} ${Interproscan_dataframe} ${kallistotrans} ${basename}
    """

}

// Process 6: Create Interproscantoxinplotly
process CreateInterproscanToxinPlotly {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    tuple path(transdf_distinct_csv), path(Toxin_domains), path(transtoxinfasta)

    output:
    path "filtered_sequences.fasta", emit: filtered_sequences
    path "plotly_graph.html", emit: plotly_graph
 

    script: 
    """
    
    Rscript bin/Intermediate_Scripts/IS7.R ${transdf_distinct_csv} ${Toxin_domains} ${transtoxinfasta}
    """

}

// Process 7: Create BUSCOgraphtranscriptome  Literally the busco provided script for visualisation
process BUSCOtranscriptome {
    conda "bin/Setup/busco2.yaml"

    publishDir "results/busco/transcriptome", mode: 'copy'

    input:

    output:
    path "*.png", emit: busco_transcriptome
 

    script: 
    """
     
    python3 "bin/Intermediate_Scripts/IS6.py" -wd "${params.data}/busco_transcriptome/"
    """

}

// Process 8: Create BUSCOgraphtranslatome  Literally the busco provided script for visualisation
process BUSCOtranslatome {
    conda "bin/Setup/busco2.yaml"

    publishDir "results/busco/translatome", mode: 'copy'

    input:

    output:
    path "*.png", emit: busco_translatome
 

    script: 
    """
     
    python3 "bin/Intermediate_Scripts/IS6.py" -wd "${params.data}/busco_translatome/"

    """

}

// Process 9: Create TableGenerationTrinity
process TableGenerationTrinity {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

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
	
    Rscript "bin/Intermediate_Scripts2/Generating_TopTables_Trinity.R" ${TBK} ${genome_id} ${species}
    """

}

// Process 10: Create TableGenerationTransdecoder   visualisation
process TableGenerationTransdecoder {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

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
	
    Rscript "bin/Intermediate_Scripts2/Generating_Tables_Transdecoder_SignalP.R" ${transdf} ${genome_id} ${species}
    """

}


// Process 11: Create FigureGenerationTrinity
process FigureGenerationTrinity {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts2_outputs", mode: 'copy'

    input:
    tuple path(TBK), path(colour)

    output:
    path "pie1.png", emit: pie1
    path "pie2.png", emit: pie2
    path "pie3.png", emit: pie3
    path "pie4.png", emit: pie4
    path "alluvial1.png", emit: alluvial1
    path "alluvial2.png", emit: alluvial2
    path "Table13.csv", emit:Table13

 
    script: 
    """
	
    Rscript "bin/Intermediate_Scripts2/Figure_generation_Trinity.R" ${TBK} ${colour}

    """

}

// Process 12: Create FigureGenerationTransdecoder
process FigureGenerationTransdecoder {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts2_outputs", mode: 'copy'

    input:
    tuple path(transdf), path(colour)

    output:
    path "pie5.png", emit: pie5
    path "pie6.png", emit: pie6
    path "pie7.png", emit: pie7
    path "pie8.png", emit: pie8
    path "alluvial3.png", emit: alluvial3
    path "alluvial4.png", emit: alluvial4
    path "Table14.csv", emit:Table14
 

    script: 
    """
	
    Rscript "bin/Intermediate_Scripts2/Figure_generation_Transdecoder.R" ${transdf} ${colour}

    """

}

// Process 13: Create FigureGenerationSignalp
process FigureGenerationSignalp {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts2_outputs", mode: 'copy'

    input:
    tuple path(transdf), path(colour)

    output:
    path "pie9.png", emit: pie9
    path "pie10.png", emit: pie10
    path "pie11.png", emit: pie11
    path "pie12.png", emit: pie12
    path "alluvial5.png", emit: alluvial5
    path "alluvial6.png", emit: alluvial6
    path "Table15.csv", emit:Table15
 

    script: 
    """

    Rscript "bin/Intermediate_Scripts2/Figure_generation_SignalP.R" ${transdf} ${colour}

    """

}

// Process 14: Create AddMassSpec
process AddMassSpec {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(transdf), path(massspecdata), val(species), val(basename)

    output:
    path "*_filtered_masspec_select.csv", emit: filtered_massspec
    path "*_distinct_masspec.csv.gz" , emit: distinct_massspec
 

    script:
    """

    Rscript "bin/Intermediate_Scripts/IS8.R" ${transdf} ${massspecdata} "${species}" ${basename}

    """

}

// Process 15: Create SkipMassSpec
process SkipMassSpec {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(transdf), val(species), val(basename)

    output:
    path "*_filtered_nomasspec_.csv", emit: filtered_nomasspec
    path "*_distinct_nomasspec.csv.gz" , emit: distinct_nomasspec
 

    script:
    """

    Rscript "bin/Intermediate_Scripts/IS9.R" ${transdf} "${species}" ${basename}

    """

}

// Process 16: Create VennOverview
process VennOverviewMS {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(filtered_massspec), path(Toxin_domains)

    output:
    path "*.csv", emit: overviewcsv
    path "*.png" , emit: overviewpng
 

    script:
    """

    Rscript "bin/Intermediate_Scripts/IS10.R" ${filtered_massspec} "${Toxin_domains}"

    """

}

// Process 16b: Create VennOverview
process VennOverviewNoMS {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    tuple path(filtered_nomasspec), path(Toxin_domains)

    output:
    path "*.csv", emit: overviewcsv
    path "*.png" , emit: overviewpng
 

    script:
    """
    
    Rscript "bin/Intermediate_Scripts/IS10.R" ${filtered_nomasspec} "${Toxin_domains}"
    
    """

}

// Process 17: RmarkdownA
process RmarkdownA {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    val(sampleURL)

    output:
    path "*.html"
 

    script:
    """

    Rscript -e "rmarkdown::render('bin/Rmarkdown_scripts/A.Rmd', output_dir = '.')" "bin/Rmarkdown_scripts/" "${sampleURL}" "${params.name}"
    """

}

// Process 18: RmarkdownB
process RmarkdownB {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path settings
    path notes
    path triminfo
    
    
    output:
    path "*.html"

    script:
    """

    settings_abs=\$(readlink -f "${settings}")
    notes_abs=\$(readlink -f "${notes}")
    triminfo_abs=\$(readlink -f "${triminfo}")

    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/B.Rmd',
      output_dir='.',
      params=list(
        settings='\$settings_abs',
        notes='\$notes_abs',
        triminfo='\$triminfo_abs',
        name = '${params.name}'
      )
      )"
    """
}

process RmarkdownCDEGIK {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:

    
    output:
    path "*.html"

    script:
    """

    Rscript -e "rmarkdown::render('bin/Rmarkdown_scripts/C.Rmd', output_dir = '.')" "bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('bin/Rmarkdown_scripts/D.Rmd', output_dir = '.')" "bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('bin/Rmarkdown_scripts/E.Rmd', output_dir = '.')" "bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('bin/Rmarkdown_scripts/G.Rmd', output_dir = '.')" "bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('bin/Rmarkdown_scripts/I.Rmd', output_dir = '.')" "bin/Rmarkdown_scripts/" ${params.name}
    Rscript -e "rmarkdown::render('bin/Rmarkdown_scripts/K.Rmd', output_dir = '.')" "bin/Rmarkdown_scripts/" ${params.name}


    """
}

process RmarkdownH {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (kallistotop20graphtrinity)
    path (kallistotop500graphtrinity)
    path (busco_figure)
    path (alluvial1)
    path (alluvial2)
    path (pie1)
    path (pie2)
    path (pie3)
    path (pie4)
    path (topkallisto)

    
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
      'bin/Rmarkdown_scripts/H.Rmd',
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


process RmarkdownJ {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (kallistotop20graphtransdecoder)
    path (kallistotop500graphtransdecoder)
    path (busco_figure_transdecoder)
    path (alluvial3)
    path (alluvial4)
    path (pie5)
    path (pie6)
    path (pie7)
    path (pie8)
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
      'bin/Rmarkdown_scripts/J.Rmd',
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

process RmarkdownL {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (alluvial5)
    path (alluvial6)
    path (pie9)
    path (pie10)
    path (pie11)
    path (pie12)
    path (topkallisto_signalp)

    
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
      'bin/Rmarkdown_scripts/L.Rmd',
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

process RmarkdownM {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Table1)
    
    output:
    path "*.html"

    script:
    """

    Table1_abs=\$(readlink -f "${Table1}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/M.Rmd',
      output_dir='.',
      params=list(
        Table1='\$Table1_abs',
        name = '${params.name}'
      )
    )"
    """
}

process RmarkdownN {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Table2)
    path (Table3)
    
    output:
    path "*.html"

    script:
    """

    Table2_abs=\$(readlink -f "${Table2}")
    Table3_abs=\$(readlink -f "${Table3}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/N.Rmd',
      output_dir='.',
      params=list(
        Table2='\$Table2_abs',
        Table3='\$Table3_abs',
        name = '${params.name}'
      )
    )"
    """
}

process RmarkdownO {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Table4)
    
    output:
    path "*.html"

    script:
    """

    Table4_abs=\$(readlink -f "${Table4}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/O.Rmd',
      output_dir='.',
      params=list(
        Table4='\$Table4_abs',
        name = '${params.name}'
      )
    )"
    """
}

process RmarkdownQ {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Table5)
    
    output:
    path "*.html"

    script:
    """

    Table5_abs=\$(readlink -f "${Table5}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/Q.Rmd',
      output_dir='.',
      params=list(
        Table5='\$Table5_abs',
        name = '${params.name}'
      )
    )"
    """
}

process RmarkdownR {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Table6)
    path (Table7)
    
    output:
    path "*.html"

    script:
    """

    Table6_abs=\$(readlink -f "${Table6}")
    Table7_abs=\$(readlink -f "${Table7}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/R.Rmd',
      output_dir='.',
      params=list(
        Table6='\$Table6_abs',
        Table7='\$Table7_abs',
        name = '${params.name}'
      )
    )"
    """
}

process RmarkdownS {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Table8)
    
    output:
    path "*.html"

    script:
    """

    Table8_abs=\$(readlink -f "${Table8}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/S.Rmd',
      output_dir='.',
      params=list(
        Table8='\$Table8_abs',
        name = '${params.name}'
      )
    )"
    """
}



process RmarkdownV {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Table9)
    
    output:
    path "*.html"

    script:
    """

    Table9_abs=\$(readlink -f "${Table9}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/V.Rmd',
      output_dir='.',
      params=list(
        Table9='\$Table9_abs',
        name = '${params.name}'
      )
    )"
    """
}

process RmarkdownW {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Table10)
    path (Table11)
    
    output:
    path "*.html"

    script:
    """

    Table10_abs=\$(readlink -f "${Table10}")
    Table11_abs=\$(readlink -f "${Table11}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/W.Rmd',
      output_dir='.',
      params=list(
        Table10='\$Table10_abs',
        Table11='\$Table11_abs',
        name = '${params.name}'
      )
    )"
    """
}

process RmarkdownX {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Table12)
    
    output:
    path "*.html"

    script:
    """

    Table12_abs=\$(readlink -f "${Table12}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/X.Rmd',
      output_dir='.',
      params=list(
        Table12='\$Table12_abs',
        name = '${params.name}'
      )
    )"
    """
}

process Movefilteredseq {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/Intermediate_Scripts1_outputs", mode: 'copy'

    input:
    path (filtered_sequences)
    
    output:

    script:
    """
    mv "results/htmls/filtered_sequences.fasta" "results/Intermediate_Scripts1_outputs/filtered_sequences.fasta"
    """
}


process RmarkdownZ {
    conda "bin/Setup/VenomFlowAnalysis2.yaml"

    publishDir "results/htmls", mode: 'copy'

    input:
    path (Venn)
    path (table)
    
    output:
    path "*.html"

    script:
    """

    Venn_abs=\$(readlink -f "${Venn}")
    table_abs=\$(readlink -f "${table}")


    Rscript -e "rmarkdown::render(
      'bin/Rmarkdown_scripts/Z.Rmd',
      output_dir='.',
      params=list(
        VENN='\$Venn_abs',
        TABLE='\$table_abs',
        name = '${params.name}'
      )
    )"
    """
}

// Define input file patterns via parameters
params.input_kallisto_trinity = "${params.data}/*_kalltrin.tsv"
params.input_kallisto_trans   = "${params.data}/*_kalltrans.tsv"
params.input_transdecoder_pep = "${params.data}/*_transpep.pep"
params.input_transdecoder_cds = "${params.data}/*_transcds.cds"
params.input_mature_fasta   = "${params.data}/*_mature.fasta"
params.input_trinity_fasta = "${params.data}/*_trinity.fasta"
params.input_blastx_files   = "${params.data}/*_blastxunitox6.txt"
params.input_blastp_files   = "${params.data}/*_blastpunitox6.txt"
params.input_interproscan = "${params.data}/*.cleaned.pep.tsv"
params.input_signalp_summary   = "${params.data}/*_summary.signalp5"
params.genome_id = null
params.species = null 
params.colour = file('bin/Intermediate_Scripts2/color_palette.rds')
params.ismassspecavailable = 'N'
params.basename = null
params.sampleURL = 'NULL'
workflow {

    def kallisto_file_trinity = Channel.fromPath(params.input_kallisto_trinity)

    def kallisto_file_trans = Channel.fromPath(params.input_kallisto_trans)

    def transdecoder_pep = Channel.fromPath(params.input_transdecoder_pep)

    def transdecoder_cds = Channel.fromPath(params.input_transdecoder_cds)
 

    def mature_fasta = Channel.fromPath(params.input_mature_fasta)


    def signal_inputs = transdecoder_pep.combine(mature_fasta)


    kallisto_file_trinity | kallistoAnalysisTrinity
    kallisto_file_trans   | kallistoAnalysisTrans
    signal_inputs | ExtractSignalSequences

    def kallisto_trin_csv = kallistoAnalysisTrinity.out.trin_all_csv


    def trinity_fasta = Channel 
        .fromPath(params.input_trinity_fasta)
    
    def blastx_files = Channel 
        .fromPath(params.input_blastx_files)

    def blastp_files = Channel 
        .fromPath(params.input_blastp_files)

    def signalp_summary = Channel 
        .fromPath(params.input_signalp_summary)

    def trinity_data = trinity_fasta.combine(blastx_files).combine(kallisto_trin_csv)

    trinity_data | CreateTrinityDataframe

    def Interproscan = Channel.fromPath(params.input_interproscan)
    def ListFile = Channel.fromPath(params.input_list)
    def PantherFile = Channel.fromPath(params.input_panther)

    def Interproscan_data = Interproscan
                            .combine(ListFile)
                            .combine(PantherFile)


    Interproscan_data | CreateInterproscanDataframe

    def signalpsequences = ExtractSignalSequences.out.signalsequences

    def Interproscan_dataframe = CreateInterproscanDataframe.out.Interproscan_dataframe

    def kallisto_trans_csv = kallistoAnalysisTrans.out.trans_all_csv

    def Transdecoderdf_data = transdecoder_pep
                              .combine(transdecoder_cds)
                              .combine(blastp_files)
                              .combine(mature_fasta)
                              .combine(signalp_summary)
                              .combine(signalpsequences)
                              .combine(Interproscan_dataframe)
                              .combine(kallisto_trans_csv)

    
    Transdecoderdf_data | CreateTransdecoderDataframe

    def transdf_distinct_csv = CreateTransdecoderDataframe.out.transdf_distinct
    def Toxin_domains = Channel.fromPath(params.input_toxindomains)
    def transtoxinfasta = CreateTransdecoderDataframe.out.toxin_fasta
    def InterproscanToxinData = transdf_distinct_csv
                                .combine(Toxin_domains)
                                .combine(transtoxinfasta)

    InterproscanToxinData | CreateInterproscanToxinPlotly
   BUSCOtranscriptome()
   BUSCOtranslatome()

   def genome_id = Channel.value(params.genome_id)
   def species = Channel.value(params.species)
   def info = genome_id.combine(species)
   def TBK = CreateTrinityDataframe.out.TBK

    def TBK_table_data  = TBK.combine(genome_id).combine(species).view()

    TBK_table_data | TableGenerationTrinity

    def transdf = CreateTransdecoderDataframe.out.transdf
    
    def Trans_table_data  = transdf    
                         .combine(genome_id)
                         .combine(species)

    Trans_table_data | TableGenerationTransdecoder
    def colours = Channel.value(params.colour)
    def TBK_only = CreateTrinityDataframe.out.TBK

    def Tbkfigures = TBK_only
                     .combine(colours)
    def Transdffigures = transdf
                     .combine(colours)
    Tbkfigures | FigureGenerationTrinity
    Transdffigures | FigureGenerationTransdecoder
    Transdffigures | FigureGenerationSignalp
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
            RmarkdownZ (VennOverviewMS.out.overviewpng,VennOverviewMS.out.overviewcsv)

        } else {

            def nomassspecdatatable = transdf
                                    .combine(species)
                                    .combine(basename)

            nomassspecdatatable | SkipMassSpec
            def filtered_nomasspec = SkipMassSpec.out.filtered_nomasspec
            def overviewnoms = filtered_nomasspec.combine(Toxin_domains)
            overviewnoms | VennOverviewNoMS
            RmarkdownZ (VennOverviewNoMS.out.overviewpng, VennOverviewNoMS.out.overviewcsv)
            
        }
    def sampleurl = Channel.value(params.sampleURL)
    sampleurl | RmarkdownA
    
    def N = Channel.value(params.notes)
    def T = Channel.value(params.triminfo)
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
    def Metadata = N.combine(T).combine(SN).combine(R).combine(Ti).combine(A).combine(SX).combine(REM).combine(LP).combine(SD).combine(SL).combine(RL)
    
    Metadata | MetadataCreation
    
    RmarkdownB(
        MetadataCreation.out.settings,
        MetadataCreation.out.notes,
        MetadataCreation.out.triminfo
    )
 
    RmarkdownCDEGIK()

    
    RmarkdownH(
        kallistoAnalysisTrinity.out.trin_top20_png,
        kallistoAnalysisTrinity.out.trin_top500_png,
        BUSCOtranscriptome.out.busco_transcriptome,
        FigureGenerationTrinity.out.pie1,
        FigureGenerationTrinity.out.pie2,
        FigureGenerationTrinity.out.pie3,
        FigureGenerationTrinity.out.pie4,
        FigureGenerationTrinity.out.alluvial1,
        FigureGenerationTrinity.out.alluvial2,
        FigureGenerationTrinity.out.Table13
        
    )
    
    RmarkdownJ(
        kallistoAnalysisTrans.out.trans_top20_png,
        kallistoAnalysisTrans.out.trans_top500_png,
        BUSCOtranslatome.out.busco_translatome,
        FigureGenerationTransdecoder.out.pie5,
        FigureGenerationTransdecoder.out.pie6,
        FigureGenerationTransdecoder.out.pie7,
        FigureGenerationTransdecoder.out.pie8,
        FigureGenerationTransdecoder.out.alluvial3,
        FigureGenerationTransdecoder.out.alluvial4,
        FigureGenerationTransdecoder.out.Table14
        
    )
    
    RmarkdownL(
        FigureGenerationSignalp.out.pie9,
        FigureGenerationSignalp.out.pie10,
        FigureGenerationSignalp.out.pie11,
        FigureGenerationSignalp.out.pie12,
        FigureGenerationSignalp.out.alluvial5,
        FigureGenerationSignalp.out.alluvial6,
        FigureGenerationSignalp.out.Table15
        
    )
    
    RmarkdownM(
    TableGenerationTrinity.out.Table1
)

    RmarkdownN(
    TableGenerationTrinity.out.Table2,
    TableGenerationTrinity.out.Table3
)

    RmarkdownO(
    TableGenerationTrinity.out.Table4
)

    RmarkdownQ(
    TableGenerationTransdecoder.out.Table5
)

    RmarkdownR(
    TableGenerationTransdecoder.out.Table6,
    TableGenerationTransdecoder.out.Table7
)

    RmarkdownS(
    TableGenerationTransdecoder.out.Table8
)

    RmarkdownV(
    TableGenerationTransdecoder.out.Table9
)

    RmarkdownW(
    TableGenerationTransdecoder.out.Table10,
    TableGenerationTransdecoder.out.Table11
)

    RmarkdownX(
    TableGenerationTransdecoder.out.Table12
)

 Movefilteredseq (CreateInterproscanToxinPlotly.out.filtered_sequences)
 
 
}







