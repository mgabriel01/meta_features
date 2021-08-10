# meta_features
get meta-features/profiles from bigwig files and gff annotations

## Prerequisites

- awk (GNU awk), version >= 4

- R version >= 3.4, with packages : 
            
            ggplot2
            GenomicFeatures
            GenomicAlignments
            GenomicRanges
            foreach
            Biostrings
            doParallel
            data.table
            rtracklayer
            ggpubr
            plyr
            gridExtra
            optparse
            compiler
            gtools
            
- UCSC tools for bigwig/tdf files (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) :
           
          bedGraphToBigWig

          wigToBigWig

          bigWigToWig

          bigWigToBedGraph 

- igvtools (https://software.broadinstitute.org/software/igv/download)
- bedtools version >= 2.29.0 (https://bedtools.readthedocs.io/en/latest/)

## Steps

- 1-) Use the script **get_tdf_bw_from_bam.sh** [https://github.com/mgabriel01/meta_features/blob/main/get_tdf_bw_from_bam.sh] to have bigwig/tdf files raw & normalized in rpm (the bigwig files will be used for the meta-profile). Modify the section of the `input data` according to your situation.


- 2-) Use the script **getMergedExonsPerGenes.sh** [https://github.com/mgabriel01/meta_features/blob/main/getMergedExonsPerGenes.sh] to have meta-transcripts (for a given gene, exons from all transcripts are merged, to have a single line in genome browsers like IGV). Run the script without any options to see the inputs. Dont forget to redirect the stdout to a file.[**skip this step if you already have your defined features in bed format**]

              Example on GAPDH :
   ![](https://github.com/mgabriel01/meta_features/blob/main/igv_snapshot_gapdh_metatranscript.png)

- 3-) Use the script **getTranscriptByExonsByStrand.sh** [https://github.com/mgabriel01/meta_features/blob/main/getTranscriptByExonsByStrand.sh] on the result of the script **getMergedExonsPerGenes.sh** (empty run to see the inputs), it will give you the higher level (by default it's "transcript") ; **run it with the option `-u "yes"`**, in order to avoid to have many genes with the same gene_id, that could create issues with the next script (paralogous genes), then concatenate both files (it's in order to create the intronic features).[**skip this step if you already have your defined features in bed format**]


- 4-) Use the script **getIntronsByTranscripts.R** [https://github.com/mgabriel01/meta_features/blob/main/getIntronsByTranscripts.R] to create the "intron" feature. Run the script without any options to see the inputs. Dont forget to redirect the stdout to a file.[**skip this step if you already have your defined features in bed format**]

- 5-) Use this script **getExonsIntronsNumbered.sh**[https://github.com/mgabriel01/meta_features/blob/main/getExonsIntronsNumbered.sh] to give position number to exons and/or introns (depending on what you want, you can have the last intron or exon, by using the pattern "position=Last_intron" or "position=Last_exon" on the output file. Modify the section of the input data according to your situation.[**skip this step if you already have your defined features in bed format**]
              
             Remark : Use the unix command "comm -12 IDs_file1.txt IDs_file2.txt" to have common IDs between exons and introns between two files after selecting the IDs for both, if you want to have features at the same positions (like first introns and first exons, last introns and last exons, etc). 
             Here IDs_file1.txt & IDs_file2.txt are text files with IDs you would like to intersect.

- 6-) Convert the gff files in bed6 format (don't forget, bed format is 0-based), to supply to the script that will do the profiles.[**skip this step if you already have your defined features in bed format**]


- 7-) Fill the tables & associative tables in the script **insert script here**,  in the section `input data`, then run it.


      The actual design in the script will create a plot with 3 panels (in the example, we have 3 fractions : chromatin, cytoplasm, total) ; for each of them, we want the exact same parts : Exon1, Intron1, Exon2, Intron2, Last_intron, Last_exon
      So we supply for each of them : 

                         - the same name of features -> variable "all_features"
                         - the same annotations in bed format -> variable "all_bed"
                         - the same labels for the start and the end of the features -> variable "all_delim"
                         - their own file containing their bigwig files (= coverage files ; supply only the plus strand, the minus one will be found automatically in the same directory) -> variable "all_design"
                         - the associated prefix for each panel, and for each feature in the panel -> variable "all_prefix"
                         
      At the end, the 3 variables, should have the same length (in the example, it's 18)


The design file(s), from the variable `all_design` in the script, should be like this (field separator : tabulation) :

            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_1_1_plus_strand_normalization_RPM.bw	clone_2_1	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_1_2_plus_strand_normalization_RPM.bw	clone_2_1	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_1_3_plus_strand_normalization_RPM.bw	clone_2_1	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_1_ad_1_plus_strand_normalization_RPM.bw	clone_2_1_ad	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_1_ad_2_plus_strand_normalization_RPM.bw	clone_2_1_ad	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_1_ad_3_plus_strand_normalization_RPM.bw	clone_2_1_ad	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_7_1_plus_strand_normalization_RPM.bw	clone_2_7	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_7_2_plus_strand_normalization_RPM.bw	clone_2_7	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_7_3_plus_strand_normalization_RPM.bw	clone_2_7	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_7_ad_1_plus_strand_normalization_RPM.bw	clone_2_7_ad	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_7_ad_2_plus_strand_normalization_RPM.bw	clone_2_7_ad	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_2_7_ad_3_plus_strand_normalization_RPM.bw	clone_2_7_ad	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_wt_1_plus_strand_normalization_RPM.bw	wt	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_wt_2_plus_strand_normalization_RPM.bw	wt	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_wt_3_plus_strand_normalization_RPM.bw	wt	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_wt_ad_1_plus_strand_normalization_RPM.bw	wt_ad	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_wt_ad_2_plus_strand_normalization_RPM.bw	wt_ad	1
            /media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chr_wt_ad_3_plus_strand_normalization_RPM.bw	wt_ad	1
      
    
   1 st column : full path to the bigwig plus strand (**the minus strand will be automatically loaded : it should have the pattern "_minus_", and be in the same directory as the plus strand**)
  
   2 nd column : name of the condition
  
   last column : the normalization factor (if it's already normalized, or you don't want to, just put 1 as in the example)
   

- 8-) Example of plot from the design in the script **insert script here** :


![](https://github.com/mgabriel01/meta_features/blob/main/metagene_plot_multipanel_PCG_chroma_PCG_cyto_PCG_tot_Exon1_Exon2_Intron1_Intron2_Last_exon_Last_intron_model2.png)

## To Do list :
- create a config file for **insert script here**
  


             
              
     
