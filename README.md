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
            
            
            




1) Use the script **insert script here** to have bigwig/tdf files raw & normalized in rpm (the bigwig files will be used for the meta-profile). Modify the section of the input data according to your situation.

2)* Use the script **insert script here** to have meta-transcripts (for a given gene, exons from all transcripts are merged, to have a single line in genome browsers like IGV). Run the script without any options to see the inputs. Dont forget to redirect the stdout to a file.

              Example on GAPDH :
   ![](https://github.com/mgabriel01/meta_features/blob/main/igv_snapshot_gapdh_metatranscript.png)

3)* Use the script **insert script here** to create the "intron" feature. Run the script without any options to see the inputs. Dont forget to redirect the stdout to a file.

4)* Use this script **getExonsIntronsNumbered.sh**[https://github.com/mgabriel01/meta_features/blob/main/getExonsIntronsNumbered.sh] to give position number to exons and/or introns (depending on what you want, you can have the last intron or exon, by using the pattern "position=Last_intron" or "position=Last_exon" . Modify the section of the input data according to your situation.
              
             Remark : Use the unix command "comm -12 IDs_file1.txt IDs_file2.txt" to have common IDs between exons and introns between two files after selecting the IDs for both, if you want to have features at the same positions (like first introns and first exons, last introns and last exons, etc). 
             Here IDs_file1.txt & IDs_file2.txt are text files with IDs you would like to intersect.

5) convert the gff files in bed6 format (don't forget, bed format is 0-based), to supply to the script that will do the profiles.


7) fill the tables & associative tables in the script **insert script here** , then run the script

The design file(s) should be like this :

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
      
  1 st column : full path to the bigwig plus strand (the minus strand will be automatically loaded : it should have the pattern "_minus_" at the same place of "_plus_" in the file in the design)
  2 nd column : name of the condition
  last column : the normalization factor (if it's already normalized, or you don't want to, just put 1 as in the example)



**\* Those steps are optional if you already have your defined features**


             
              
     
