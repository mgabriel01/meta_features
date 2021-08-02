# meta_features
get meta-features/profiles from bigwig files and gff annotations


1) Use the script **insert script here** to have bigwig/tdf files raw & normalized in rpm (the bigwig files will be used for the meta-profile). Modify the section of the input data according to your situation.

2) Use the script **insert script here** to have meta-transcripts (for a given gene, exons from all transcripts are merged, to have a single line in genome browsers like IGV). Run the script without any options to see the inputs. Dont forget to redirect the stdout to a file.

3) Use the script **insert script here** to create the "intron" feature. Run the script without any options to see the inputs. Dont forget to redirect the stdout to a file.

4) Use this script **getExonsIntronsNumbered.sh**[https://github.com/mgabriel01/meta_features/blob/main/getExonsIntronsNumbered.sh] to give position number to exons and/or introns (depending on what you want, you can have the last intron or exon, by using the pattern "position=Last_intron" or "position=Last_exon" . Modify the section of the input data according to your situation.
              
              -> extra : Use the unix command "comm -12" to have common IDs between exons and introns between two files after selecting the IDs for both, if you take features at the same positions (like first introns and first exons, last introns and last exons, etc)

5) convert the gff files in bed6 format (don't forget bed format is 0-based), to supply to the script that will do the profiles.


             
              
     
