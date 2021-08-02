#!/bin/bash

#Purpose : this script will build bigwig & tdf files, stranded raw and normalized, and all kind of coverage files in bedgraph, thanks to flags

#this version is faster and takes less space on the disk


#############  Input data ################################

#path to bedtools scripts
bedtools="bedtools"

#path to samtools
samtools="samtools"

#path to tophat outputs (input directory)
input_dir="/media/marcgabriel/saylar4/urines_ffpe_files/bam_files/ /media/marcgabriel/9d38aaa1-ce5d-4422-950f-0f083b797d13/exosomes_dekupl_result/bam_files/"


#library type from the mapper (don't forget if R1 is in the same direction as the mRNA and R2 is reverse, it's fr-secondstrand ; if R2 is in the same direction as the mRNA and R1 is reverse, it's fr-firststrand)
library_type="fr-firststrand"
#library_type="fr-secondstrand"
#library_type="unstranded"


#if yes, we will process the bam files in order 
#already_unique="no"
already_unique="yes"

single="yes"
#single="no"

#make median of raw files
do_summarized_raw="F"

#refine the previous bigwig files if we have to do tdf files
refine_bw="F"

#path to output directory
output_dir="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/visu_files/"

#you can find it on UCSC
bedGraphToBigWig="/home/marcgabriel/Desktop/scripts/bedGraphToBigWig"

#you can find it on UCSC
wigToBigWig="/home/marcgabriel/Desktop/scripts/wigToBigWig"


#number of subscripts to run in parallel
process_number_limit=6


#number of threads to give to samtools (but also "sort" command)
samtools_threads=8


#we increase the size of the buffer for the sort (in multithreaded mode), otherwise it will be single threaded
sort_buffer="8G"

igvtools="/home/marcgabriel/Downloads/IGV_2.8.2/igvtools"

bigWigToWig="/home/marcgabriel/Desktop/scripts/bigWigToWig"

bigWigToBedGraph="/home/marcgabriel/Desktop/scripts/bigWigToBedGraph"


#actually, it's a script that allows
getRowMedians="/home/marcgabriel/Desktop/scripts/getRowMedians.R"

genome="/home/marcgabriel/Documents/gencode32/GRCh38.primary_assembly_official_chromosomes.fa"



#declare -A condition_list=(

	#[D268T24]="FFPE_normal_1"
	#[D268T40]="FFPE_normal_2"
	#[D268T32]="FFPE_normal_3"
	#[D268T48]="FFPE_normal_4"
	#[D268T56]="FFPE_normal_5"
	#[D268T95]="FFPE_normal_6"

	

#)

#declare -A condition_list2=(

	#[FFPE_normal_1]="FFPE_normal"
	#[FFPE_normal_2]="FFPE_normal"
	#[FFPE_normal_3]="FFPE_normal"
	#[FFPE_normal_4]="FFPE_normal"
	#[FFPE_normal_5]="FFPE_normal"
	#[FFPE_normal_6]="FFPE_normal"
	

#)

################## end of input data #######################################



output_dir="${output_dir}/"
output_dir=$(echo $output_dir |sed 's/\/\//\//g')
if [ ! -d $output_dir ]; then mkdir $output_dir; fi

declare -p condition_list2 >${output_dir}saved.sh

temp="${output_dir}temp/"

if [ -d $temp ]; then rm -rf $temp;fi;mkdir $temp


#directory of scripts for each sample (to be run in parallel)
sample_scripts_dir="${output_dir}/sample_scripts/"
if [ -d $sample_scripts_dir ]; then rm -rf $sample_scripts_dir; fi 
mkdir $sample_scripts_dir

all_report=${output_dir}global_report.txt

bedGraphToBigWig_output=${output_dir}BigWigs/

if [ ! -d $bedGraphToBigWig_output ]; then mkdir $bedGraphToBigWig_output ; fi




for one_condition in ${!condition_list[*]};do

  new_name=${condition_list[$one_condition]}
  
  
  cond_name="${condition_list2[$new_name]}"
  
  
  
  echo -e "condition $one_condition &  new name $new_name\n "
  

  
  report=${bedGraphToBigWig_output}${one_condition}_report.txt
  
  if [ ! -f $report ];then

    >$report
    
  fi
  
  
  #list of bam
  file=$(find $input_dir -name "*bam" |grep -v -i "unmapped"| grep "${one_condition}\.unique\.bam$")
  
  if [[ ! -f $file ]];then
  
    file=$(find $input_dir -name "*bam" |grep -v -i "unmapped"| grep "${one_condition}"|grep "bam$")|| { echo "no bam file !" 1>&2; exit; }
    
  fi
  
  echo -e "file is : $file\n"
  
  
  echo -e "#!/bin/bash\n\n" >${sample_scripts_dir}subscript_${one_condition}.sh
  
   echo -e "\n$samtools view -H $file |grep -E -v \"^@PG|^@HD|^@RG\"|grep -v \"user\" | awk 'OFS=\"\\\t\"{print \$2,\$3}' | sed 's/SN\://;s/LN\://g' | LANG=en_EN sort -k1,1 >${output_dir}${one_condition}.genome_size.txt\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
  echo -e "\nchrs_length=${output_dir}${one_condition}.genome_size.txt\n" >>${sample_scripts_dir}subscript_${one_condition}.sh

  
  #select only reads not mapped anywhere else (2 = ID of mate 1 & 2, if it's less, that means it's a singleton)
  #|cut -f1|sort --parallel=8 -T $temp |uniq -c|awk '{if($1<=2){print}}'|awk '{a=a+$1}''END{print a}')
  
  echo -e "\n	if [ ! -f ${bedGraphToBigWig_output}${new_name}_plus_strand_normalization_RPM.tdf ] || [ ! -f ${bedGraphToBigWig_output}${new_name}_minus_strand_normalization_RPM.tdf ] ; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
   echo -e "\necho -e \"$file \n  -----\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
   
    #initialization of the signal to compute the scaling factor
    doit="no"
    
    
    #if the report file is there, look for the scaling factor
    if [[ -f "$report" ]];then
    
  
       my_factor=$(grep "scaling factor is" $report|sed 's/scaling factor is //g')
       
       if [ "$my_factor" == "" ];then
       
       
       
         if [[ -f "$all_report" ]];then
       
          #if the result is empty, check the global file
          my_factor=$(grep -P "$file\t" $all_report|head -n1|cut -f2)
          
         else
         
           my_factor=""
         
         
         fi
          
       fi
       
       #if the scaling factor is there, sore it, otherwise, change the signal
       if [ "$my_factor" != "" ];then
      
         echo -e "\nscaling_factor=$my_factor\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
         
         doit="no"
       
       else
       
         doit="yes"
       
       
       fi
       
    else
    
       doit="yes"
    
    fi
   
    #if the scaling factor isn't there, compute it  
    if [ "$doit" == "yes" ];then
    
      echo -e "\necho \"Analysis date :\" >>$report\n" >>${sample_scripts_dir}subscript_${one_condition}.sh

      echo -e "\ndate \"+%A %d %B %Y %H:%M:%S\" | xargs echo >>$report\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
      echo -e "\necho -e \"\ncondition is $cond_name \n \" >>$report\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
      echo -e "\necho -e \"\nfile is $file \n \" >>$report\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
               
            if [[ "$already_unique" == "yes" ]];then
            
            
            if [[ $single == "no" ]];then
            
            
				
				echo -e "\nmapped_paired=\$($samtools view -@ $samtools_threads -F 0x8 -F 0x4 -f 0x1 -c $file|awk '{print \$1}')\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
				echo -e "\nmapped_singletons=\$($samtools view -@ $samtools_threads -F 0x4 -f 0x8 -c $file |awk '{print \$1}')\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
				
				
		  else
		  
		  
		  		echo -e "\nmapped_paired=\$($samtools view -@ $samtools_threads -F 0x4 -c $file|awk '{print \$1}')\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
		  		
		  		echo -e "\nmapped_singletons=0\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
		  
		  
		  
		  fi
            
            
            else
            
             if [[ $single == "no" ]];then
      
        
				echo -e "\nmapped_paired=\$($samtools view -@ $samtools_threads -F 0x8 -F 0x4 -f 0x1 $file|cut -f1|sort -S $sort_buffer --parallel=$samtools_threads -T $temp -u |wc -l|awk '{print \$1*2}')\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
				echo -e "\nmapped_singletons=\$($samtools view -@ $samtools_threads -F 0x4 -f 0x8 $file|cut -f1|sort -S $sort_buffer --parallel=$samtools_threads -T $temp -u |wc -l|awk '{print \$1}')\n" >>${sample_scripts_dir}subscript_${one_condition}.sh


             else
             
             
             		  		echo -e "\nmapped_paired=\$($samtools view -@ $samtools_threads -F 0x4 -c $file|awk '{print \$1}')\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
             		  		
             		  		echo -e "\nmapped_singletons=0\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
             
             
             
             fi
             
             
  
			fi
			
			echo -e "\ntotal_reads=\$((mapped_paired+mapped_singletons))\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
     
      echo -e "\necho -e \"\ntotal mapped reads is \${total_reads}\n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
      echo -e "\nscaling_factor=\$(echo \"scale=8;x=1000000/\$total_reads ;if(x<1) print 0;x\" | bc)\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
      echo -e "\necho -e \"\nscaling factor is \${scaling_factor}\n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
      echo -e "\necho -e \"\ntotal mapped reads is \${total_reads}\n\" >>$report\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
      echo -e "\necho -e \"\nscaling factor is \${scaling_factor}\n\" >>$report\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
      echo -e "\necho -e \"\n=================================\n\" >>$report\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
   fi
  
 
      
        

  


    echo -e "\nsign_list=(plus minus)\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
    
    echo -e "\nfor sign in \${sign_list[*]};do\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
    
    
      echo -e "\n	if [ ! -f ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_normalization_RPM.tdf ] ; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
    
    
    
       #if no bam file sorted properly, create it
       echo -e "\n	if [ ! -f ${output_dir}${one_condition}_XS_\${sign}_sorted.bam.bai ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
  
				echo -e "\n	  echo -e \"\n no sorted bam file for XS+, we're going to make it\n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
			   
				echo -e "\n	if [ \"\$sign\" == \"plus\" ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
				
				
				if [[ "$library_type" == "fr-firststrand" ]];then
				
				
					if [[ $single == "no" ]];then 
				
						#XS+
						echo -e "cat <($samtools view -h -@ $samtools_threads -f128 -F16 $file) <($samtools view -@ $samtools_threads -f64 -f16 $file) |$samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
						
					
					else
					
					 #XS+
					echo -e "$samtools view -h -@ $samtools_threads -f16 $file|$samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam|| { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
					
					
					fi
				
				elif [[ "$library_type" == "fr-secondstrand" ]];then
				
				
					if [[ $single == "no" ]];then 
					
						#XS+
						echo -e "cat <($samtools view -h -@ $samtools_threads -f128 -f16 $file ) <($samtools view -@ $samtools_threads -f64 -F16 $file) |$samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam|| { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
						
					
					else
					
					 #XS+
					echo -e "$samtools view -h -@ $samtools_threads -F16 $file | $samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
						
					
					fi           
				
				else
				
				 
					echo -e "$samtools view -h -@ $samtools_threads $file  | $samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
					
				
				fi
			   
			  
				echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
				
				echo -e "\n	if [ \"\$sign\" == \"minus\" ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
				
				
				if [[ "$library_type" == "fr-firststrand" ]];then
				
				
					  if [[ $single == "no" ]];then  
					
						#XS-
						echo -e "cat <($samtools view -h -f128 -f16 $file) <($samtools view -@ $samtools_threads -f64 -F16 $file)  | $samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
						
					 else
					 
						#XS-
						echo -e "$samtools view -h -F16 $file  |$samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
							
					 
					 fi

			   
				elif [[ "$library_type" == "fr-secondstrand" ]];then
				
				
				
					 if [[ $single == "no" ]];then 
					
						#XS-
						echo -e "cat <($samtools view -h -f128 -F16 $file) <($samtools view -@ $samtools_threads -f64 -f16 $file)|$samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
						
					 else
					 
						#XS-
						echo -e "$samtools view -h -f16 $file | $samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; } \n" >>${sample_scripts_dir}subscript_${one_condition}.sh
									 
					 fi          
				
				else
				
				
				 
					#XS-
					echo -e "$samtools view -h $file | $samtools sort -@ $samtools_threads -T ${temp}${one_condition}_XS_\${sign}_sorted_tmp - >${output_dir}${one_condition}_XS_\${sign}_sorted.bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
					
							   
				
				fi
			  
				echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh


				echo -e "\n	echo -e \"${output_dir}${one_condition}_XS_\${sign}_sorted.bam done \n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh

		echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
		
		
		echo -e "\n	if [ ! -f ${output_dir}${one_condition}_XS_\${sign}_sorted.bam.bai ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
		
		
		   echo -e "\n	  $samtools index ${output_dir}${one_condition}_XS_\${sign}_sorted.bam\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
		
		echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
	
		echo -e "\n	if [ ! -f ${output_dir}${one_condition}_coverage_raw_\${sign}_strand.bed ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
		
			echo -e "\n	  echo -e \"\n no coverage on \${sign} strand from genomeCoverageBed, we're going to make it\n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
			
			#-g ${chrs_length}
		   echo -e "\n	  ${bedtools} genomecov -bg -split -ibam ${output_dir}${one_condition}_XS_\${sign}_sorted.bam|grep -v -E \"^chr[0-9]+_|^chr[A-Z]+_\" >${output_dir}${one_condition}_coverage_raw_\${sign}_strand.bed ||{ echo \"genomecov 1 failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
	  
		   echo -e "\n	  echo -e \"coverage on strand \$sign done\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
	  
		echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
		



		echo -e "\n	LANG=en_EN sort -T $temp -k1,1 -k2,2n ${output_dir}${one_condition}_coverage_raw_\${sign}_strand.bed >${output_dir}${one_condition}_coverage_raw_\${sign}_strand.tmp && mv ${output_dir}${one_condition}_coverage_raw_\${sign}_strand.tmp ${output_dir}${one_condition}_coverage_raw_\${sign}_strand.bed\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
        
        echo -e "\n	if [ ! -f ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_raw.bw ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
           echo -e "\n	  echo -e \"\n no bigwig file for ${one_condition} on \${sign} strand from bedGraphToBigWig, we're going to make it\n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  $bedGraphToBigWig ${output_dir}${one_condition}_coverage_raw_\${sign}_strand.bed \${chrs_length} ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_raw.bw ||{ echo \"bedGraphToBigWig 1 failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  echo -e \"\n bigwig done for ${one_condition} on strand \$sign \n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
        echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
        echo -e "\n	if [ ! -f ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_raw.tdf ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
           echo -e "\n	  echo -e \"\n no tdf file for ${one_condition} on \${sign} strand from totdf, we're going to make it\n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
           
           
           echo -e "\n	  $bigWigToWig ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_raw.bw ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_raw.wig||{ echo \"bigWigToWig 1 failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  $igvtools toTDF --tmpDir ${bedGraphToBigWig_output} ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_raw.wig ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_raw.tdf $genome && rm ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_raw.wig||{ echo \"totdf 1 failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  echo -e \"\n tdf done for ${one_condition} on strand \$sign \n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
        echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
        echo -e "\n	if [ ! -f ${output_dir}${one_condition}_coverage_normalization_RPM_\${sign}_strand.bed ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  echo -e \"\n no coverage normalization_RPM on \${sign} strand from genomeCoverageBed, we're going to make it\n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
      
          #-g ${chrs_length}
          echo -e "\n	  ${bedtools} genomecov -bg -split -ibam ${output_dir}${one_condition}_XS_\${sign}_sorted.bam -scale \$scaling_factor|grep -v -E \"^chr[0-9]+_|^chr[A-Z]+_\"  >${output_dir}${one_condition}_coverage_normalization_RPM_\${sign}_strand.bed && rm ${output_dir}${one_condition}_XS_\${sign}_sorted.bam||{ echo \"genomecov 2 failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
          echo -e "\n	  echo -e \"coverage normalization_RPM on strand \$sign done\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
      
        echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
        echo -e "\n	LANG=en_EN sort -T $temp -k1,1 -k2,2n ${output_dir}${one_condition}_coverage_normalization_RPM_\${sign}_strand.bed >${output_dir}${one_condition}_coverage_normalization_RPM_\${sign}_strand.tmp && mv ${output_dir}${one_condition}_coverage_normalization_RPM_\${sign}_strand.tmp ${output_dir}${one_condition}_coverage_normalization_RPM_\${sign}_strand.bed\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
      
        echo -e "\n	if [ ! -f ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_normalization_RPM.bw ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  echo -e \"\n no bigwig file normalization_RPM for ${one_condition} on ${sign} strand from bedGraphToBigWig, we're going to make it\n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  $bedGraphToBigWig ${output_dir}${one_condition}_coverage_normalization_RPM_\${sign}_strand.bed \${chrs_length} ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_normalization_RPM.bw && rm ${output_dir}${one_condition}_coverage_normalization_RPM_\${sign}_strand.bed ||{ echo \"bedGraphToBigWig 2 failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  echo -e \"\n bigwig normalization_RPM done for ${one_condition} on strand \$sign \n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
        echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
        echo -e "\n	if [ ! -f ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_normalization_RPM.tdf ]; then\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
           echo -e "\n	  echo -e \"\n no tdf file for ${one_condition} on \${sign} strand from totdf, we're going to make it\n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
           
           
           echo -e "\n	  $bigWigToWig ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_normalization_RPM.bw ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_normalization_RPM.wig||{ echo \"bigWigToWig 2 failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  $igvtools toTDF --tmpDir ${bedGraphToBigWig_output} ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_normalization_RPM.wig ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_normalization_RPM.tdf $genome && rm ${bedGraphToBigWig_output}${new_name}_\${sign}_strand_normalization_RPM.wig||{ echo \"totdf 2 failure !\" 1>&2; exit; }\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
           echo -e "\n	  echo -e \"\n tdf done for ${one_condition} on strand \$sign \n\"\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
      
        echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
        
             echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
        
       
    echo -e "\ndone\n" >>${sample_scripts_dir}subscript_${one_condition}.sh
    
    echo -e "\n	fi\n" >>${sample_scripts_dir}subscript_${one_condition}.sh  
    
   
done


#exit

find ${sample_scripts_dir} -name "*subscript_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "executing bash sorting and/or counting subsritpts failure !" 1>&2; exit; }


report_files=$(find $output_dir -name "*_report.txt"|grep -v "global")

for i in $report_files;do

	my_file=$(grep "file is" $i|sort -u)

	my_cond=$(grep "condition is" $i|sort -u)

	my_factor=$(grep "scaling factor is" $i|sort -u)

	if [ "$my_file" == "" ];then

	  continue
	  
	fi

	echo -e "${my_file}\t${my_cond}\t${my_factor}" |sed 's/file is //g'|sed 's/condition is //g'|sed 's/scaling factor is //g' >>${all_report}

done



sort -u -k1,1 ${all_report} >${all_report}.tmp && mv ${all_report}.tmp ${all_report}

files_to_remove=($(find ${output_dir} -name "*.bam"))
files_to_remove+=($(find ${output_dir} -name "*.bam.bai"))
files_to_remove+=($(find ${output_dir} -name "*strand.bed"))

echo -e "\nfiles to remove : \n${files_to_remove[*]}\n##############\n"

for i in ${files_to_remove[*]};do rm $i ;done

rm -rf $temp


######################### make one file (+ & -) per condition ###############

BigWigSummarizedByConditionsDir="${output_dir}BigWigSummarizedByConditions/"

if [ ! -d $BigWigSummarizedByConditionsDir ]; then mkdir $BigWigSummarizedByConditionsDir; fi

#directory of scripts for each sample (to be run in parallel)
sample_scripts_dir="${BigWigSummarizedByConditionsDir}/sample_scripts/"
if [ -d $sample_scripts_dir ]; then rm -rf $sample_scripts_dir; fi 
mkdir $sample_scripts_dir


chrs_length=$(find ${output_dir} -name "*genome_size.txt"|head -n1)


if [[ "$do_summarized_raw" == "T" ]];then

	norm_list=("normalization_RPM" "raw")

else

	norm_list=("normalization_RPM")

fi

all_strands=(plus_strand minus_strand)


for one_norm in ${norm_list[*]};do


	for one_cond in $(echo -e "${condition_list2[*]}"|tr ' ' '\n'|sort -u);do
	
	
					   echo -e "- $one_cond ($one_norm) :\n"
					   
					   
					   echo -e "#!/bin/bash\n\n" >${sample_scripts_dir}subscript_${one_cond}.sh
					   
					   #load the associative array
					   echo -e "source ${output_dir}saved.sh\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
					   
					   
					   echo -e "for one_strand in ${all_strands[*]};do\n\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
					   
								  								   
									echo -e "	if [[ ! -f ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bw ]];then\n\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
								   
								   
											 
											 echo -e "	files_one_strand_one_cond=()\n" >>${sample_scripts_dir}subscript_${one_cond}.sh


											   											 
												echo -e "	for one_sample in $(echo -e "${!condition_list2[*]}");do\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
																							   
												 
												      echo -e "		if [[ \"\${condition_list2[\$one_sample]}\" == \"$one_cond\" ]];then\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												 
												 
												   
												      echo -e "			echo -e \"\\\t- \$one_sample (+)\\\n\"\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												   
												   												   
												   
												   echo -e "			files_one_strand_one_cond+=(\$(find ${bedGraphToBigWig_output} -name \"\${one_sample}_\${one_strand}_${one_norm}.bw\"))\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												   
												  
												 
												   echo -e "		fi\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												 
												 
											   
											   
											   echo -e "	done\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
											   
											   
											   
												
											   echo -e "	echo -e \"\\\t\\\t- all files + : \${files_one_strand_one_cond[*]}\\\n\"\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												
											   
												
											   echo -e "all_per_pos_counts=()\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												
												
											   
											   echo -e "	count=1\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
											   
											   
											   
											   
											   
											   echo -e "	for one_file in \${files_one_strand_one_cond[*]};do\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
											   
											   												 
												 
												echo -e "		$bigWigToBedGraph \$one_file ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_\${count}.tmp\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												 
												 												 												 
												 
												echo -e "		all_per_pos_counts+=(${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_\${count}.tmp)\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												 
												 
												 
												 
												echo -e "		count=\$((count+1))\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												 

											   
											   echo -e "	done\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
											   

											   

											   echo -e "	$bedtools unionbedg -i \${all_per_pos_counts[*]} >${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_union.bed\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
											   
											   	
																		   
											   #remove pos with 0 values to avoid having a big file
											   
											   echo -e "	paste <(cut -f1-3 ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_union.bed) <($getRowMedians ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_union.bed 4) |awk 'OFS=\"\\\t\"{if(\$4!=0){print}}' >${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bed && $bedGraphToBigWig ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bed ${chrs_length} ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bw && rm ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bed ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_union.bed \${all_per_pos_counts[*]}\n" >>${sample_scripts_dir}subscript_${one_cond}.sh


									   
									  echo -e "fi\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
									   
									
									   
									  echo -e "if [[ ! -f ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.tdf ]];then\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
									   
											


											echo -e "if [[ \"\$refine_bw\" == \"T\" ]];then\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
											 

											echo -e "$bigWigToBedGraph ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bw ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bed && cat ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bed |grep -v \"^#\"  >${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.tmp && mv ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.tmp ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bed && $bedGraphToBigWig ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bed $chrs_length ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bw && $bigWigToWig ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bw ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.wig && $igvtools toTDF --tmpDir ${BigWigSummarizedByConditionsDir} ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.wig ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.tdf $genome && rm ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.wig ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bed\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												 
												
											 
												echo -e "else\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
											 
											 
												 
												echo -e "	$bigWigToWig ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.bw ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.wig\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												 
												 
												 
												echo -e "	$igvtools toTDF --tmpDir ${BigWigSummarizedByConditionsDir} ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.wig ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.tdf $genome && rm ${BigWigSummarizedByConditionsDir}${one_cond}_\${one_strand}_${one_norm}_Mean.wig\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
												 
											 
											 
											echo -e "fi\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
						         
										  
									   
										echo -e "fi\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
									   
				   

						   
						echo -e "done\n" >>${sample_scripts_dir}subscript_${one_cond}.sh
						
					   

			   
	done
	
	
	find ${sample_scripts_dir} -name "*subscript_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "executing bash sorting and/or counting subsritpts failure !" 1>&2; exit; }

done


echo -e "\n====\n"

echo -e "\n check the mean of coverage for each condition in : $BigWigSummarizedByConditionsDir\n\n"
echo -e "\n check the individual files for each sample in : $bedGraphToBigWig_output\n\n"



