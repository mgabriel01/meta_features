#!/bin/bash

#this script will print for each gene, the merged exons (across all its transcripts) in gff3 format

getTranscriptByExonsByStrand="/home/marcgabriel/Desktop/scripts/getTranscriptByExonsByStrand.sh"


usage() { echo -e "Usage: $0 <Required arguments>\n\n
 \tRequired arguments :\n
                  -a < annotation (gff3 format) >\n                 
                  -o < output dir >\n
 \tOptional arguments :\n
                  -f <feature (default : gene_name ; other choice : gene_id)>\n
                  -l <level (default : exon)>\n
                  -p < process_number_limit (default=0 ; it means it will use all the available cores to run the subscripts) >\n" 1>&2; exit 1;}

[[ $# -eq 0 ]] && usage

while getopts ":a:p:f:l:o:" opt; do
  case $opt in
  
      p)
      
	      process_number_limit=$OPTARG
	      
	      #echo -e "\n#process_number_limit is : $process_number_limit\n" >&2
	      
              ;;
      a)
      
	      annotation=$OPTARG
	      
	      #echo -e "\n#annotation is : $annotation\n" >&2
	      
              ;;
 
      f)
      
	      feature_to_keep=$OPTARG
	      
	      #echo -e "\n#feature to check is : $feature_to_keep\n" >&2
	      
              ;;   
              
     l)
      
	      level_to_keep=$OPTARG
	      
	      #echo -e "\n#feature to check is : $feature_to_keep\n" >&2
	      
              ;;                 
                       
      
      o)
      
          
	      output_dir=$OPTARG
	      
      
      	;;
     

      #invalid options (options not in the list)
      ######################
      
      
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    
  esac
done


if [ "$process_number_limit" == "" ];then

  process_number_limit=0

fi

if [ "$feature_to_keep" == "" ];then

  feature_to_keep="gene_name"

fi

if [ "$level_to_keep" == "" ];then

  level_to_keep="exon"

fi



if [ "$annotation" == "" ] || [ "$output_dir" == "" ]; then
      
echo -e "\none required argument is missing or is wrong !!\n"
      usage
      exit 1
fi


#echo -e "annotation is : ${annotation}\n
         #output dir is : ${output_dir}\n
         #feature to check is : ${feature_to_keep}\n***\n"

#### inputs #########

#bedtools program
bedtools="bedtools"

#script that will allow to construct introns for each transcript of each gene
getIntronsByTranscripts="/home/marcgabriel/Desktop/scripts/getIntronsByTranscripts.R"

#gff file
#annotation="/home/marcgabriel/Documents/gencode27lift37/gencode.v27lift37.annotation.gff3"

#output_dir="/home/marcgabriel/Documents/gencode27lift37/"

#number of chromosomes to process in parallel
#if 0, that means it will process all the subscripts (genes) in parallel : let it as it is, the script doesn't use a big memory
#process_number_limit=0


#####################


####### process inputs #######

output_dir="${output_dir}/"
output_dir=$(echo $output_dir |sed 's/\/\//\//g')

#if no output dir, then create it
if [ ! -d $output_dir ];then mkdir $output_dir;fi


#if temp fir exists, remove it, and re-create it
temp_dir="${output_dir}tmp_process/"
if [ -d $temp_dir ];then rm -rf $temp_dir;fi
mkdir $temp_dir

#directory of subscripts (to be laucnhed in parallel)
sample_scripts_dir="${temp_dir}subscripts_dir/"
if [ -d $sample_scripts_dir ]; then rm -rf $sample_scripts_dir; fi 
mkdir $sample_scripts_dir

#create introns
gff_file="${temp_dir}modified_annot.gff3"

if [ ! -f $gff_file ];then


  if [ "$include_introns" == "yes" ];then

    sort -T ${temp_dir} -S "4G" --parallel=8 -k1,1 -k4,4n <(grep -v "^#" $annotation) <($getIntronsByTranscripts $annotation $temp_dir "yes") >$gff_file
  
  else
  
    grep -v "^#" $annotation | sort -T ${temp_dir} -S "4G" --parallel=8 -k1,1 -k4,4n >$gff_file
  
  fi

fi



############################


#process exons

#gff of only exons
grep -P "\t${level_to_keep}\t" $gff_file|sort -T ${temp_dir} -k1,1 -k4,4n|awk 'OFS="\t"{$3="exon";print}' >${temp_dir}only_exons.gff



#be carefull ! bedtools 2.27-2.29 don't put the strand when using merge anymore !!
#build merged exons (keep the gene_id info to know if we have a mix of genes); order by gene_id, then by chr, then by start
paste <(awk 'OFS="\t"{print $1,$4-1,$5,$7}' ${temp_dir}only_exons.gff) <(cut -f9 ${temp_dir}only_exons.gff|sed 's/;/\n/g'|grep "${feature_to_keep}"|sed "s/${feature_to_keep}=//g"|sed "s/${feature_to_keep} \"//g"|sed "s/\"//g"|sed "s/ //g")|awk 'OFS="\t"{print $1,$2,$3,$5,".",$4}' |sort -T ${temp_dir} -k1,1 -k2,2n | bedtools merge -s -i stdin -c 4,6 -o distinct,distinct -delim "---" 2>/dev/null|awk 'OFS="\t"{print $1,$2,$3,$4,".",$5}' |sort -T ${temp_dir} -k4,4 -k1,1 -k2,2n >${temp_dir}merged_exons.bed

#store the ones that have overlaps with other genes
#warning ! we have to ensure the research of the exact IDs (to avoid to lose some of them, e.g, with just a "grep -v" on NF1, we may lose all the IDs with this pattern, like ZNF169) !
#|awk '{print $1"\t"}'
grep -E "\-\-\-" ${temp_dir}merged_exons.bed|cut -f4|sed 's/---/\n/g'|while read line;do echo -e "\t${line}---\n---${line}\t\n---${line}---\n\t${line}\t" ;done|sort -T ${temp_dir} -u >${temp_dir}overlapped_genes.txt

sort -u -T ${temp_dir} ${temp_dir}overlapped_genes.txt >${temp_dir}overlapped_genes.tmp && mv ${temp_dir}overlapped_genes.tmp ${temp_dir}overlapped_genes.txt

sed 's/---/\n/g' ${temp_dir}overlapped_genes.txt|grep -v "^$"|sort -u -T ${temp_dir}|awk '{print "\t"$1"\t"}' |sort -u -T ${temp_dir} >${temp_dir}overlapped_genes2.txt

cat ${temp_dir}overlapped_genes2.txt ${temp_dir}overlapped_genes.txt >${temp_dir}overlapped_genes.tmp && mv ${temp_dir}overlapped_genes.tmp ${temp_dir}overlapped_genes.txt && rm ${temp_dir}overlapped_genes2.txt


#select the ones that have no overlaps with other genes
LC_ALL=C grep -F -v -f ${temp_dir}overlapped_genes.txt ${temp_dir}merged_exons.bed >${temp_dir}merged_exons.tmp && mv ${temp_dir}merged_exons.tmp ${temp_dir}merged_exons.bed

#loop over IDs that have overlaps (they will be processed separately)
for one_ID in $(grep -v -P "\t" ${temp_dir}overlapped_genes.txt |sed 's/---/\n/g'|sort -u -T ${temp_dir}|grep -v "^$");do 


      echo -e "#!/bin/bash\n\n" >${sample_scripts_dir}subscript_${one_ID}.sh
   
   
      #print the fasta header (">" + gene_id) ; print the sequences (strand specific, & respect order of exons of minus genes : sort -T ${temp_dir} -k1,1 -k3,3nr) of all the merged exons of all the transcripts for the gene (200 characters max per line)
      echo -e "grep -E \"=$one_ID;|=$one_ID$\" ${temp_dir}only_exons.gff|sort -T ${temp_dir} -k1,1 -k4,4n|$bedtools merge -s -i stdin -c 7 -o distinct -delim \"---\" 2>/dev/null|awk -v ID=$one_ID 'OFS=\"\\\t\"{print \$1,\$2,\$3,ID,\".\",\$4}' >${temp_dir}${one_ID}_sub_chr.bed\n" >>${sample_scripts_dir}subscript_${one_ID}.sh
   
   chmod 755 ${sample_scripts_dir}subscript_${one_ID}.sh

done

#run the subscripts in parallel
find ${sample_scripts_dir} -name "*subscript_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "subscript failure !" 1>&2; exit; }

#when it's done remove the subscripts
find ${sample_scripts_dir} -name "*subscript_*\.sh" | while read F; do rm ${F} ; done

>${temp_dir}merged_exons.tmp
find ${temp_dir} -name "*_sub_chr.bed" | while read F; do cat ${F} >>${temp_dir}merged_exons.tmp; done

cat ${temp_dir}merged_exons.bed ${temp_dir}merged_exons.tmp >${temp_dir}merged_exons.tmp2 && mv ${temp_dir}merged_exons.tmp2 ${temp_dir}merged_exons.bed && rm ${temp_dir}merged_exons.tmp


find ${temp_dir} -name "*_sub_chr.bed" | while read F; do rm ${F} ; done


awk 'OFS="\t"{print $1,".","exon",$2+1,$3,".",$6,".","ID="$4";gene_id="$4";gene_name="$4";Parent="$4}' ${temp_dir}merged_exons.bed >${temp_dir}merged_exons.gff3



sort -T ${temp_dir} -k1,1 -k4,4n ${temp_dir}merged_exons.gff3


#exit
rm -rf ${temp_dir}


