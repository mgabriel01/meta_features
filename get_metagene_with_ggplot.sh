#!/bin/bash

################    Explanation ###############
###############################################

#This design will create 3 panels (in the example, we have 3 fractions : chromatin, cytoplasm, total) ; for each of them, we want the exact same parts : Exon1, Intron1, Exon2, Intron2, Last_intron, Last_exon
#So we supply for each of them : 
#
#                         - the same name of features -> variable "all_features"
#                         - the same annotations in bed format -> variable "all_bed"
#                         - the same labels for the start and the end of the features -> variable "all_delim"
#                         - their own file containing their bigwig files (= coverage files ; supply only the plus strand, the minus one will be found automatically in the same directory) -> variable "all_design"
#                         - the associated prefix for each panel, and for each feature in the panel -> variable "all_prefix"
# At the end, the 3 variables, should have the same length (in the example, it's 18)
#
##############################################
##############################################



########## input data ################


#check if the R scripts are in the same directory as the main script
metagene_rscript=$(dirname "$0")/get_metagene_with_ggplot.R
if [ ! -f $metagene_rscript ];then echo -e "get_metagene_with_ggplot.R script is missing ! ";exit;fi
chmod 755 $metagene_rscript


metaplot_metagene=$(dirname "$0")/getMultipanelMetagenePlot.R
if [ ! -f $metaplot_metagene ];then echo -e "getMultipanelMetagenePlot.R is missing ! ";exit;fi
chmod 755 $metaplot_metagene


#output_dir="/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/PCG_plots/"
output_dir="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/PCG_plots/"

#it's parallelized on the length of the tables (they all have the same length, here in the example 18)
process_number_limit=4

                                 
#name of the features, for each panel
all_features=("Exon1" "Intron1" "Exon2" "Intron2" "Last_intron" "Last_exon"
              "Exon1" "Intron1" "Exon2" "Intron2" "Last_intron" "Last_exon"
              "Exon1" "Intron1" "Exon2" "Intron2" "Last_intron" "Last_exon"
                  )


#bed file of the features, for each panel (should be the same length as the other tables)
all_bed=(

"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/chroma_fraction_upregulated_1Exons.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/chroma_fraction_upregulated_1Introns.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/chroma_fraction_upregulated_2Exons.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/chroma_fraction_upregulated_2Introns.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/chroma_fraction_upregulated_LastIntrons.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/chroma_fraction_upregulated_LastExons.bed"

"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/cyto_fraction_upregulated_1Exons.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/cyto_fraction_upregulated_1Introns.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/cyto_fraction_upregulated_2Exons.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/cyto_fraction_upregulated_2Introns.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/cyto_fraction_upregulated_LastIntrons.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/cyto_fraction_upregulated_LastExons.bed"

"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/tot_fraction_upregulated_1Exons.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/tot_fraction_upregulated_1Introns.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/tot_fraction_upregulated_2Exons.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/tot_fraction_upregulated_2Introns.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/tot_fraction_upregulated_LastIntrons.bed"
"/media/marcgabriel/saylar2/Dominika_XRN1_25_03_2020/bed_for_metagene_PCG/tot_fraction_upregulated_LastExons.bed"

)

#name of the start and end of the features, for each panel (should be the same length as the other tables)
all_delim=("start,end" "start,end" "start,end" "start,end" "start,end" "start,end"
		   "start,end" "start,end" "start,end" "start,end" "start,end" "start,end"
		   "start,end" "start,end" "start,end" "start,end" "start,end" "start,end"
		  )
		  
#link to the file containing the bigwig files, for each panel (should be the same length as the other tables)
#for of the design file :
#  </full/path/to/bigwig> <tabulation> <condition> <tabulation> <normalization factor>
all_design=("/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chroma_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chroma_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chroma_design.tsv" 
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chroma_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chroma_design.tsv" 
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/chroma_design.tsv"
			
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/cyto_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/cyto_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/cyto_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/cyto_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/cyto_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/cyto_design.tsv"
		   
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/total_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/total_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/total_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/total_design.tsv"
			"/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/total_design.tsv"
		   "/media/marcgabriel/homeborn/dominika_hxrn1_cyto_chroma_tot_bigwig_files/total_design.tsv"
		   )		  

   
#name of the panels (should be the same length as the other tables)     
all_prefix=("PCG_chroma"
            "PCG_chroma"
            "PCG_chroma"
            "PCG_chroma"
            "PCG_chroma"
            "PCG_chroma"
            
            "PCG_cyto"
            "PCG_cyto"
            "PCG_cyto"
            "PCG_cyto"
            "PCG_cyto"
            "PCG_cyto"
            
            "PCG_tot"
            "PCG_tot"
            "PCG_tot"
            "PCG_tot"
            "PCG_tot"
            "PCG_tot"
                    )
                    
        
######### end of input data #############



length_of_tables=$(echo -e ${#all_features[*]} ${#all_bed[*]}|sed 's/ /\n/g'|sort -u|wc -l|awk '{print $1}') 

if [[ $length_of_tables -gt 1 ]];then


  echo -e "\n\nAt least one of the table doesn't have the same length as the others, check again ! \n"
  
  
  exit

fi
             
output_dir="${output_dir}/"
output_dir=$(echo $output_dir |sed 's/\/\//\//g')

if [ ! -d $output_dir ];then mkdir $output_dir;fi


subscripts_dir="${output_dir}subscripts/"

if [ -d $subscripts_dir ];then rm -r $subscripts_dir;fi

mkdir ${subscripts_dir}


all_RDATA_files=()

all_parts_Mean=()

#loop across the number of annotations (bed files)
for i in $(seq 0 $((${#all_bed[*]}-1)));do

  #take one bed file
  one_bed=${all_bed[$i]}
  
  #take one design
  one_design=${all_design[$i]}
  
  #take one prefix
  one_prefix=${all_prefix[$i]}
  
  #take one feature
  one_feat=${all_features[$i]}
  
  #one delim
  one_delim=${all_delim[$i]}

        #if no RDATA file (ggplot data)) for this feature, with this prefix, do it
        if [ ! -f "${output_dir}${i}_metapart_metagene_data_${one_prefix}_${one_feat}.RDATA" ];then
     
 
			echo -e "#!/bin/bash\n\n" >${subscripts_dir}${i}_metagene_all_fractions.sh
			


			echo -e "$metagene_rscript -b $one_bed -f $one_design -p $one_prefix -m $one_feat -d \"$one_delim\" -i $i -o \"${output_dir}\"" >>${subscripts_dir}${i}_metagene_all_fractions.sh
			
			chmod 775 ${subscripts_dir}${i}_metagene_all_fractions.sh
			
			all_RDATA_files+=("${output_dir}${i}_metapart_metagene_data_${one_prefix}_${one_feat}.RDATA")
			
			all_parts_Mean+=("${output_dir}${i}_${one_prefix}_allPartsMean_${one_feat}.tsv")

         else
         
           all_RDATA_files+=("${output_dir}${i}_metapart_metagene_data_${one_prefix}_${one_feat}.RDATA")
         
           echo -e "${output_dir}${i}_metapart_metagene_data_${one_prefix}_${one_feat}.RDATA already exists..."
           
           all_parts_Mean+=("${output_dir}${i}_${one_prefix}_allPartsMean_${one_feat}.tsv")
         
         
         fi

done

find "${subscripts_dir}" -name "*_metagene_all_fractions.sh*" | xargs -n 1 -P $process_number_limit bash

echo -e "\nall metagenes in RDATA file are : \n${all_RDATA_files[*]}\n"


all_RDATA_files=$(echo -e "${all_RDATA_files[*]}"|sed 's/ /,/g')

all_prefix=$(echo -e "${all_prefix[*]}"|sed 's/ /\n/g'|sort -u|tr '\n' '_'|sed 's/__/_/g'|sed 's/_$//g')

all_features=$(echo -e "${all_features[*]}"|sed 's/ /\n/g'|sort -u|tr '\n' '_'|sed 's/__/_/g'|sed 's/_$//g')

all_parts_Mean=$(echo -e "${all_parts_Mean[*]}"|sed 's/ /,/g')


$metaplot_metagene $all_RDATA_files $output_dir $all_prefix $all_features $all_parts_Mean




