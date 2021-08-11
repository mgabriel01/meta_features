#!/bin/bash


###### input data ###########

gff="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/gencode26_DIS3_scallop_work_on_introns/gencode26_DIS3_scallop_introns.gff3"


output_dir="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/gencode26_DIS3_scallop_work_on_introns/"

#feat_to_number="exon"
feat_to_number="intron"

#do you want to add and additional attribute to the gff to tag first and last features ?
#it will have the pattern "position=Last_*"
tag_last="yes"

#number of subscripts to run in parallel
process_number_limit=30

##########  end of input data #######


output_dir="${output_dir}/"
output_dir=$(echo $output_dir |sed 's/\/\//\//g')

if [ ! -d $output_dir ];then mkdir $output_dir;fi


#directory with all final subscripts to run
final_subscript="${output_dir}final_subscripts/"
if [ -d $final_subscript ]; then rm -rf $final_subscript ;fi

mkdir $final_subscript

parent_list=($(grep -P "\t${feat_to_number}" $gff |cut -f9|sed 's/;/\n/g'|grep "Parent"|grep -v -E -i "rRNA|tRNA|snoRNA|snRNA|scarRNA|7SK|Y_RNA|misc|paralo|PAR_Y|PAR_X" |awk -v awk_end=$awk_end '{print $1awk_end}'|sort -u))		


for i in ${parent_list[*]};do 

	  strand=($(LC_ALL=C grep -E "${i}$|${i};" $gff|cut -f7|sort -u))
	  
	  #check the uniqueness of the IDs, if they are not, don't process them
	  if [[ ${#strand[*]} -eq 1 ]];then
	  
	         clean_id=$(echo "${i}"|sed -E "s/;//g"|sed -E 's/Parent=//g')
	         
	  
			  echo -e "#!/bin/bash\n" >${final_subscript}final_subscript_${clean_id}.sh
			  
			  if [[ "${strand}" == "+" ]];then
			  
			         #if it's a transcript on plus strand, order by ascending start, then descending end
					 echo -e "readarray all_lines <<< \$(LC_ALL=C grep \"$i\" $gff|grep -P \"\\\t${feat_to_number}\\\t\"|sort -k4,4n -k5,5nr)\n" >>${final_subscript}final_subscript_${clean_id}.sh
							   
			  else
			  
			    #if it's a transcript on minus strand, order by descending end, then ascending start
				echo -e "readarray all_lines <<< \$(LC_ALL=C grep \"$i\" $gff|grep -P \"\\\t${feat_to_number}\\\t\"|sort -k5,5nr -k4,4n)\n" >>${final_subscript}final_subscript_${clean_id}.sh
			  
			  fi
			  
			  #give numbers according to the number of rows
			  #if tag_last is set to yes, we give as attribute o the last feature "position=Last_*" (e.g : position=Last_intron or position=Last_exon)
			  echo -e "nb_rows=\$(echo \"\${all_lines[*]}\"|sed 's/^ //g'|grep -v \"^\$\"|wc -l|awk '{print \$1}');echo \"\${all_lines[*]}\"|sed 's/^ //g'|grep -v \"^\$\"|awk -v nb_rows=\$nb_rows -v strand=$strand -v tag_last=$tag_last 'OFS=\"\\\t\"{if(tag_last!=\"yes\"){print \$0\";${feat_to_number}_number=\"NR}else{if(NR==nb_rows){print \$0\";${feat_to_number}_number=\"NR\";position=Last_\"${feat_to_number}}else{print \$0\";${feat_to_number}_number=\"NR}}}'|sed 's/;;/;/g' |sort -k4,4n >${output_dir}${clean_id}_with_${feat_to_number}_num.gff\n" >>${final_subscript}final_subscript_${clean_id}.sh	
			   
	  
	  else
	  
	  
		 echo -e "$i is on 2 strands !\n"
	  
	  
	  fi
 
 

done

find ${final_subscript} -name "final_subscript_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "parallel script failure !" 1>&2; exit; }

find ${output_dir} -name "*with_*.gff"|sort -k1,1 -k4,4n | xargs cat >${output_dir}${feat_to_number}_numbered.gff

find ${output_dir} -name "*with_*.gff"|while read f;do rm $f;done

echo -e "\n\n-> check file : ${output_dir}${feat_to_number}_numbered.gff\n\n"
	



