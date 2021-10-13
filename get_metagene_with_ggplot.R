#!/usr/bin/env Rscript

##!/usr/bin/env /bioinfo/local/build/R/R-3.3.0_centos/bin/Rscript

#cluster : /bioinfo/local/build/R/R-3.3.0_centos/bin/Rscript

#remark
#example of interactive session on cluster : qsub -I -l mem=50gb,nodes=1:ppn=20


suppressPackageStartupMessages(suppressMessages(library(ggplot2)))
suppressPackageStartupMessages(suppressMessages(library(GenomicFeatures)))
suppressPackageStartupMessages(suppressMessages(library(GenomicAlignments)))
suppressPackageStartupMessages(suppressMessages(library(GenomicRanges)))
suppressPackageStartupMessages(suppressMessages(library(foreach)))
suppressPackageStartupMessages(suppressMessages(library(Biostrings)))
suppressPackageStartupMessages(suppressMessages(library(doParallel)))
suppressPackageStartupMessages(suppressMessages(library(data.table)))
suppressPackageStartupMessages(suppressMessages(library(rtracklayer)))
suppressPackageStartupMessages(suppressMessages(library(ggpubr)))
suppressPackageStartupMessages(suppressMessages(library(plyr)))
suppressPackageStartupMessages(suppressMessages(library(gridExtra)))
suppressPackageStartupMessages(suppressMessages(library(optparse)))
suppressPackageStartupMessages(suppressMessages(library(compiler)))
suppressPackageStartupMessages(suppressMessages(library(gtools)))


####### purpose : plot a metagene from a list of genes in bed format (Marc G.) ###################
#this function takes an indexed bam file, a bed file, and give back for each condition a list of numerics (reads), on each strand for the genes reoriented 5'-3', and then for their antisense

option_list<-list(
  make_option(c("-f","--files_descriptor"),type="character",help = "tab delimited file with 3 columns, 1st : the file ; 2nd : condition ; 3rd : norm. factor"),
  make_option(c("-b","--bed"),type="character",help = "annotation in bed format (required)"),
  make_option(c("-o","--output_dir"),type="character",default="",help = "directory to store the results (required)"),
  make_option(c("-p","--prefix"),type="character",default="",help = "prefix of analysis (=title)"),
  make_option(c("-m","--feature"),type="character",default="",help = "name of the features to use in metagene"),
  make_option(c("-d","--delim"),type="character",default="TSS,TES",help = "delimitation to use"),
  make_option(c("-i","--id"),type="numeric",default=0,help = "id of the job")
)

opt_parser = OptionParser(option_list=option_list)

opt <- parse_args(OptionParser(option_list=option_list))

#if one required parameter is missing, print a warning and stop the script
if(length(opt$bed)==0 |length(opt$delim)==0| length(opt$output_dir)==0 | length(opt$files_descriptor)==0){
  
  cat("one required parameter is missing !","\n")
  print_help(opt_parser)
  q("no")
}

design_file<-opt$files_descriptor

bed_file<-opt$bed

one_prefix<-opt$prefix

one_feature<-opt$feature

one_delim<-opt$delim

one_id<-opt$id

nb_threads<-10

########### parameters ###########
##################################

home<-opt$output_dir
#cluster : /data/tmp/mgabriel
#on the cluster /tmp/mgabriel_outputs/
setwd(home)

isbedGraph=T

#nb process to run in parallel ; be careful !! This number will be multiplied by 2 (nested foreach)
nb_process=1

###################

######### functions ################
####################################

#check the script with ggbio, there's a clear (I think) explanation on it, but globally, we are extracting the reads according to their orientation fr-firststrand (R2 same direction as RNA, & R1 is reverse)
plus_parameters<-list(c(T,T),c(F,F))
minus_parameters<-list(c(T,F),c(F,T))

getReads<-function(my_parameters,my_bam,my_chr,my_start,my_end){
  
  my_start<-as.integer(my_start)
  
  my_end<-as.integer(my_end)
  
  mate_1<-""
  mate_2<-""
  for(one_par in 1:length(my_parameters)){
    
    #isSecondaryAlignment=F
    parameters=ScanBamParam(what=scanBamWhat(),which=GRanges(my_chr, IRanges(my_start,my_end)),
                            flag =scanBamFlag(isUnmappedQuery=F,isMinusStrand=my_parameters[[one_par]][1],isFirstMateRead=my_parameters[[one_par]][2]),tag=c("MD", "MT"))
    
    
    R1R2_all<-readGAlignmentsList(my_bam,use.names=T,param=parameters)
    
    if(one_par==1){
      
      mate_1<-R1R2_all
      
      
    }else{
      
      mate_2<-R1R2_all
      
      #both_mates<-GRangesList(mate_1,mate_2)
      
      both_mates<-c(mate_1,mate_2)
      
      #both_mates<-unlist(both_mates)
      
    }
    
  }
  
  #return(unlist(both_mates))
  
  #http://girke.bioinformatics.ucr.edu/CSHL_RNAseq/mydoc/mydoc_Rgraphics_7/ 
  pos<-sapply(coverage(unlist(both_mates))[my_chr], as.numeric)
  
  my_coverage<-data.frame(coverage=pos[my_start:my_end],position=my_start:my_end)
  
  return(my_coverage)
  
}


getReads<-cmpfun(getReads)

#function that will give a list of averaged signal across genes for the sense & antisense strand, in 5'->3' orientation
getDataforSignal<-function(LocInBed,BamFile,normalization){
  
  bam<-BamFile
  
  genes_in_bed<-read.table(LocInBed,header=F,check.names=F)
  
  #add +1 because bed file is 0 based
  genes_in_bed[2]<-genes_in_bed[2]+1
 
  print(paste("file : ",basename(as.character(bam)),", prefix : ",one_prefix,", feature : ",one_feature,"\n",sep=""))


  
  #list for "list of numerics" from genes 
  one_condition_genes<-list()
  
  #list for "list of numerics" from antisense of genes
  one_condition_antisense<-list()
  
  #system(paste("samtools index"," ",bam,sep="" ))
  
  #list of numerics for sense on both strands
  list_of_numerics<-list()
  
  #list of numerics for antisense on both strands
  list_of_numericsAntisense<-list()
  
  registerDoParallel(cores=nb_threads)
  
  #replace this loop by a foreach
  #list_of_numerics_bothsense_parallelized<-foreach(i=1:nrow(genes_in_bed)) %do% {
  #,.verbose=T
  invisible(foreach(i=1:nrow(genes_in_bed),.inorder=F,.errorhandling = "pass") %dopar% {
    
    #print(paste("gene number :",i,sep=""))
    
    #loc to print in bed format for bedtools
    loc_for_bedtools<-genes_in_bed[i,]
    
    #print(paste("processing location :",loc_for_samtools,"strand",loc_for_bedtools[6],bam,sep=" "))
    
    
    #dev
    
    if(loc_for_bedtools[6]=="+"){
      
      if(isbedGraph==F){
        
        senseReads_frame<-getReads(plus_parameters,as.character(bam),my_chr = as.character(unlist(loc_for_bedtools[1])),my_start=loc_for_bedtools[2],my_end=loc_for_bedtools[3])
        
      }else{
        
        
        regions<-GRanges(seqnames=as.character(unlist(loc_for_bedtools[1])), ranges=IRanges(as.numeric(loc_for_bedtools[2]),as.numeric(loc_for_bedtools[3])), strand="*")
        senseReads_frame<-import(as.character(bam),format="BigWig",which=regions)
        senseReads_frame<-as.data.frame(senseReads_frame)
        
        #if the counts are not 0, detail them by nuc
        if(nrow(senseReads_frame)>0){
          
          tmp_frame2<-data.frame()
          for(y in 1:nrow(senseReads_frame)){
            
            tmp_frame2<-rbind(tmp_frame2,data.frame(position=unique(senseReads_frame[y,]$start:senseReads_frame[y,]$end),coverage=senseReads_frame[y,]$score))
            
          }
          
          senseReads_frame<-tmp_frame2
          
        }else{
          
          senseReads_frame<-data.frame(position=numeric(),coverage=integer())
          
        }
        
        #rebuild the original range, and do the correspondance with the counts, put 0 to the ones are are not present for the given position (NA)
        tmp_frame<-data.frame(position=start(regions):end(regions))
        
        senseReads_frame<-merge(tmp_frame,senseReads_frame, 
                                by.x=colnames(tmp_frame[1]),
                                by.y=colnames(senseReads_frame[1]),
                                all.x=TRUE,
                                all.y=FALSE)
        
        senseReads_frame$coverage[is.na(senseReads_frame$coverage)==T]<-0
        
        senseReads_frame<-senseReads_frame[,c("coverage","position")]
        
      }
      
      
    }else{
      
      if(isbedGraph==F){
        
        senseReads_frame<-getReads(minus_parameters,as.character(bam),my_chr = as.character(unlist(loc_for_bedtools[1])),my_start=loc_for_bedtools[2],my_end=loc_for_bedtools[3])
        
      }else{
        
        regions<-GRanges(seqnames=as.character(unlist(loc_for_bedtools[1])), ranges=IRanges(as.numeric(loc_for_bedtools[2]),as.numeric(loc_for_bedtools[3])), strand="*")
        senseReads_frame<-import(gsub("_plus_","_minus_",as.character(bam)),format="BigWig",which=regions)
        senseReads_frame<-as.data.frame(senseReads_frame)
        
        
        if(nrow(senseReads_frame)>0){
          
          tmp_frame2<-data.frame()
          for(y in 1:nrow(senseReads_frame)){
            
            tmp_frame2<-rbind(tmp_frame2,data.frame(position=unique(senseReads_frame[y,]$start:senseReads_frame[y,]$end),coverage=senseReads_frame[y,]$score))
            
          }
          
          senseReads_frame<-tmp_frame2
          
        }else{
          
          senseReads_frame<-data.frame(position=numeric(),coverage=integer())
          
        }
        
        
        tmp_frame<-data.frame(position=start(regions):end(regions))
        
        senseReads_frame<-merge(tmp_frame,senseReads_frame, 
                                by.x=colnames(tmp_frame[1]),
                                by.y=colnames(senseReads_frame[1]),
                                all.x=TRUE,
                                all.y=FALSE)
        
        senseReads_frame$coverage[is.na(senseReads_frame$coverage)==T]<-0
        
        senseReads_frame<-senseReads_frame[,c("coverage","position")]
        
      }
      
    }
    
    
    #normalize the sense coverage
    senseReads_frame$coverage<-senseReads_frame$coverage*normalization
    
    if(loc_for_bedtools[6]=="+"){
      
      if(isbedGraph==F){
        
        antisenseReads_frame<-getReads(minus_parameters,as.character(bam),my_chr =as.character(unlist(loc_for_bedtools[1])),my_start=loc_for_bedtools[2],my_end=loc_for_bedtools[3])
        
      }else{
        
        regions<-GRanges(seqnames=as.character(unlist(loc_for_bedtools[1])), ranges=IRanges(as.numeric(loc_for_bedtools[2]),as.numeric(loc_for_bedtools[3])), strand="*")
        antisenseReads_frame<-import(gsub("_plus_","_minus_",as.character(bam)),format="BigWig",which=regions)
        antisenseReads_frame<-as.data.frame(antisenseReads_frame)
        
        if(nrow(antisenseReads_frame)>0){
          
          tmp_frame2<-data.frame()
          for(y in 1:nrow(antisenseReads_frame)){
            
            tmp_frame2<-rbind(tmp_frame2,data.frame(position=unique(antisenseReads_frame[y,]$start:antisenseReads_frame[y,]$end),coverage=antisenseReads_frame[y,]$score))
            
          }
          
          antisenseReads_frame<-tmp_frame2
          
        }else{
          
          antisenseReads_frame<-data.frame(position=numeric(),coverage=integer())
          
          
        }
        
        tmp_frame<-data.frame(position=start(regions):end(regions))
        
        antisenseReads_frame<-merge(tmp_frame,antisenseReads_frame, 
                                    by.x=colnames(tmp_frame[1]),
                                    by.y=colnames(antisenseReads_frame[1]),
                                    all.x=TRUE,
                                    all.y=FALSE)
        
        antisenseReads_frame$coverage[is.na(antisenseReads_frame$coverage)==T]<-0
        
        antisenseReads_frame<-antisenseReads_frame[,c("coverage","position")]
        
        
      }
      
      
    }else{
      #rebuild the original range, and do the correspondance with the counts, put 0 to the
      if(isbedGraph==F){
        
        antisenseReads_frame<-getReads(plus_parameters,as.character(bam),my_chr = as.character(unlist(loc_for_bedtools[1])),my_start=loc_for_bedtools[2],my_end=loc_for_bedtools[3])
        
      }else{
        
        regions<-GRanges(seqnames=as.character(unlist(loc_for_bedtools[1])), ranges=IRanges(as.numeric(loc_for_bedtools[2]),as.numeric(loc_for_bedtools[3])), strand="*")
        antisenseReads_frame<-import(as.character(bam),format="BigWig",which=regions)
        antisenseReads_frame<-as.data.frame(antisenseReads_frame)
        
        if(nrow(antisenseReads_frame)>0){
          
          tmp_frame2<-data.frame()
          for(y in 1:nrow(antisenseReads_frame)){
            
            tmp_frame2<-rbind(tmp_frame2,data.frame(position=unique(antisenseReads_frame[y,]$start:antisenseReads_frame[y,]$end),coverage=antisenseReads_frame[y,]$score))
            
          }
          
          antisenseReads_frame<-tmp_frame2
          
        }else{
          
          antisenseReads_frame<-data.frame(position=numeric(),coverage=integer())
          
          
        }
        
        
        
        tmp_frame<-data.frame(position=start(regions):end(regions))
        
        antisenseReads_frame<-merge(tmp_frame,antisenseReads_frame, 
                                    by.x=colnames(tmp_frame[1]),
                                    by.y=colnames(antisenseReads_frame[1]), 
                                    all.x=TRUE,
                                    all.y=FALSE)
        
        antisenseReads_frame$coverage[is.na(antisenseReads_frame$coverage)==T]<-0
        
        antisenseReads_frame<-antisenseReads_frame[,c("coverage","position")]
        
      }
      
    }
    
    
    #normalize the antisense coverage
    antisenseReads_frame$coverage<-antisenseReads_frame$coverage*normalization
    
    one_list<-senseReads_frame$coverage    
    #print(paste("counting sense for",loc_for_samtools,"strand",loc_for_bedtools[6],bam,"done",sep=" "))
    
    #as.matrix to convert in numeric
    one_list_antisense<-antisenseReads_frame$coverage        
    #print(paste("counting antisense for",loc_for_samtools,"strand",loc_for_bedtools[6],bam,"done",sep=" "))
    
    #for the counting on gene -, reverse to 5'->3'        
    if (loc_for_bedtools[6]=="-"){one_list<-rev(one_list)}
    
    one_list<-approx(one_list,n=100,method="linear",ties="ordered")$y
    
    #for the counting on antisense gene -, reverse to 3'<-5' (because the antisense of gene - is like this 5->3')
    if (loc_for_bedtools[6]=="-"){one_list_antisense<-rev(one_list_antisense)}
    
    one_list_antisense<-approx(one_list_antisense,n=100,method="linear",ties="ordered")$y
    
    
    table_sense_antisense<-data.frame(sense=one_list,antisense=one_list_antisense)
    
    write.table(table_sense_antisense,file=paste(i,"_",basename(as.character(BamFile)),"_",one_feature,"_",one_prefix,".tsv",sep=""), sep='\t',row.names=F, col.names=T, quote=F)
    
    
    
  })#end of the counting on the whole bed file with genes
  
  
  #retrieve for each gene, for the corresponding file, the table of values
  sense_antisense_all_genes_one_file<-grep(paste("[0-9]+","_",basename(as.character(BamFile)),"_",one_feature,"_",one_prefix,".tsv",sep=""),list.files(),value=T,perl=T)
  
  #store the as list of list of values
  list_of_numerics_bothsense_parallelized<-list()
  for(one_gene_pos in 1:length(sense_antisense_all_genes_one_file)){
    
    one_gene_table<-read.delim(sense_antisense_all_genes_one_file[one_gene_pos])
    
    list_of_numerics_bothsense_parallelized[[length(list_of_numerics_bothsense_parallelized)+1]]<-list(t(one_gene_table)[1,],t(one_gene_table)[2,])
  }
  
  #remove the files
  file.remove(sense_antisense_all_genes_one_file)
  
  list_of_numerics<-list()
  list_of_numericsAntisense<-list()
  for(one_cov_gene in 1:length(list_of_numerics_bothsense_parallelized)){
    
    list_of_numerics[[length(list_of_numerics)+1]]<-list_of_numerics_bothsense_parallelized[[one_cov_gene]][[1]]
    list_of_numericsAntisense[[length(list_of_numericsAntisense)+1]]<-list_of_numerics_bothsense_parallelized[[one_cov_gene]][[2]]
  }
  
  #each row is the genes values, on x positions, then make the mean at each pos, to have just one row
  list_of_numericsScaled<-as.data.frame(t(sapply(list_of_numerics, rbind)))
  list_of_numericsScaled<-colMeans(list_of_numericsScaled)
  
  list_of_numericsAntisenseScaled<-as.data.frame(t(sapply(list_of_numericsAntisense,rbind)))
  list_of_numericsAntisenseScaled<-colMeans(list_of_numericsAntisenseScaled)
  
  MyListOfListsOfNumScaled<-list(list_of_numericsScaled,list_of_numericsAntisenseScaled)
  
  return(MyListOfListsOfNumScaled)
  
  
  
}

getDataforSignal<-cmpfun(getDataforSignal)
#don't use a list, print the result directly in the disk


getAverageCov<-function(mydataframe){
  
  my_result<-data.frame()
  
  for( i in 1:length(unique(mydataframe$pos))){
    
    one_pos=unique(mydataframe$pos)[i]
    
    for(j in 1:length(unique(mydataframe$condition))){
      
      one_cond=unique(mydataframe$condition)[j]
      
      tmp1<-mydataframe[which(mydataframe$pos==one_pos & mydataframe$condition==one_cond),]
      
      
      tmp2<-tmp1[1,]
      
      tmp2$count<-mean(tmp1$count)
      
      tmp2$the_min<-min(tmp1$count)
      
      tmp2$the_max<-max(tmp1$count)
      
      
      my_result<-rbind(my_result,tmp2)
      
      
    }
    
  }
  
  return(my_result)
  
}

getAverageCov<-cmpfun(getAverageCov)

rescaleList<-function(b="",a="",mylist="",mymax="",mymin=""){
  
  mylist<-unlist(mylist)
  scaled_values<-list()
  
  for(x in 1:length(mylist)){
    
    scaled_values[length(scaled_values)+1]<-(((b-a)*(mylist[x] - mymin))/ (mymax - mymin)) +a
    
  }
  
  return(scaled_values)
  
}

rescaleList<-cmpfun(rescaleList)

#########################################################
#########################################################



one_delim<-unlist(strsplit(one_delim,","))

######### process and check the inputs ################
#######################################################
#process the design file
my_design<-read.delim(design_file,header=F)
names(my_design)<-c("file","condition","norm_factors")
my_design$sample<-""
#create replicates per condition
tmp_design<-data.frame()
for(i in 1:length(unique(my_design$condition))){
  
  one_cond<-as.character(unique(my_design$condition)[i])
  
  table_oneCond<-my_design[which(my_design$condition==one_cond),]
  
  table_oneCond$sample<-paste(table_oneCond$condition,1:nrow(table_oneCond),sep="_")
  
  tmp_design<-rbind(tmp_design,table_oneCond)
  
}
my_design<-tmp_design

#############################################################
#############################################################

all_conditions_numerics_gene<-list()

#registerDoParallel(cores=nb_process)

#sometimes random issue there : try to export the function and packages...
#don't parallelize here, the order is used for the conditions
all_conditions_numerics_gene<-foreach(l=1:nrow(my_design)) %do% {
  
  
  getDataforSignal(LocInBed=bed_file,
                   BamFile=my_design$file[l],
                   normalization=my_design$norm_factors[l]
                   )
  
  
}

#we have something like this :
#x=condition(1,...n),y=scaled counts(1 : 5->3 or 2 : 3<-5)
#all_conditions_numerics_gene[[x]][[y]]
#try plot(all_conditions_numerics_gene[[1]][[1]],type="l")
##############################################



sense_gene<-data.frame()
antisense_gene<-data.frame()
for(i in 1:length(all_conditions_numerics_gene)){
  
  
  sense_gene<-rbind(sense_gene,data.frame(count=unlist(all_conditions_numerics_gene[[i]][1]),pos=1:length(unlist(all_conditions_numerics_gene[[i]][1])),sample=my_design$sample[i],condition=my_design$condition[i]))
  antisense_gene<-rbind(antisense_gene,data.frame(count=unlist(all_conditions_numerics_gene[[i]][2]),pos=1:length(unlist(all_conditions_numerics_gene[[i]][2])),sample=my_design$sample[i],condition=my_design$condition[i]))
  
}


sense_gene_mean<-getAverageCov(sense_gene)
antisense_gene_mean<-getAverageCov(antisense_gene)


correspondance<-data.frame(pos=unique(sense_gene_mean$pos),new_pos=unlist(rescaleList(b=100,a=1,mylist=list(1:length(unique(sense_gene_mean$pos))),mymax=length(unique(sense_gene_mean$pos)),mymin=1)))

for(i in 1:nrow(sense_gene_mean)){
  
  sense_gene_mean$pos[i]<-correspondance[which(correspondance$pos==sense_gene_mean$pos[i]),]$new_pos
}



all_colors<-rainbow(length(unique(my_design$condition)))

#here we fixd the color of the conditions, but you could do it automatically
both_conds<-c("cornflowerblue","red","cyan","deeppink","#3eb489","orange","gold")
names(both_conds)<-unique(my_design$condition)

white_background<-ggplot2::theme(axis.line = element_line(colour = "black"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.background = element_blank())

sense_gene_mean$condition <- factor(sense_gene_mean$condition, levels = unique(as.character(sense_gene_mean$condition)))

sense_gene_mean$feature<-one_feature

sense_gene_mean$project<-one_prefix


sense_gene_mean$extra<-as.numeric(system(paste("wc -l ",bed_file," |awk '{print $1}'",sep=""),intern=T))


#write the table to avoid to re-do all these steps
#rbind(sense_five_mean,sense_gene_mean,sense_three_mean)
write.table(sense_gene_mean,file=paste(home,one_id,"_",one_prefix,"_allPartsMean_",one_feature,".tsv",sep=""), sep='\t',row.names=F, col.names=T, quote=F)

#rbind(sense_five_mean,sense_gene_mean,sense_three_mean)       
ggplot_img_to_save<-ggplot(sense_gene_mean, aes(pos, count,color=condition))+
  geom_line(size=2)+
  geom_ribbon(aes(ymin=the_min, ymax=the_max,fill=condition), alpha=0.2)+
  #with _smooth, we go under 0 sometimes ...find a way to avoid this
  #stat_smooth(aes(x = pos, y = count), formula = y ~ s(x, k =20), method = "gam", se = FALSE,size=2)+
  theme(axis.text.x=element_text(color="black",size=14,face="bold"))+
  xlim(c(min(sense_gene_mean$pos),max(sense_gene_mean$pos)))+
  #fix the limite, if you don't want the smoothing make the y axis go under 0
  scale_y_continuous(expand=c(0,0))+
  
  scale_x_continuous(expand=c(0,0),breaks = c(1, 100), 
                     
                     labels = c(one_delim))+
  geom_hline(yintercept = 0,color="black",size=1)+
  #" upregulated "
  ggtitle(paste("metagene on ",as.numeric(system(paste("wc -l ",bed_file," |awk '{print $1}'",sep=""),intern=T))," ",one_feature, "\ncond : ",one_prefix,sep=""))+
  theme(axis.text.x=element_text(color="black",size=16,hjust=1,face="bold",angle = 45),axis.text.y = element_text(color="black",size=16,face="bold"),plot.title = element_text(hjust = 0.5,size=18),axis.title=element_text(size=16,face="bold"),legend.text=element_text(size=20),legend.title=element_text(size=18,face="bold"))+
  scale_color_manual(values=c(both_conds),name="condition")+
  scale_fill_manual(values=c(both_conds),name="condition")+
  xlab(paste(one_feature," (5'->3')",sep=""))+
  ylab("Mean normalized coverage (rpm)")+
  white_background+
  geom_vline(xintercept =c(1,100), linetype="dashed")


#plot for only genebody (or only the feature that you have entered in the bed file, for example exons, or introns)  
png(filename=paste(home,"metagene_plot_upstream_genebody_downstream_",as.numeric(system(paste("wc -l ",bed_file," |awk '{print $1}'",sep=""),intern=T)),"_",one_prefix,"_",one_feature,".png",sep=""),width=1200,height=1000)

print(
  
  ggplot_img_to_save
  
)

dev.off()


save(ggplot_img_to_save, file =paste(home,one_id,"_metapart_metagene_data_",one_prefix,"_",one_feature,".RDATA",sep=""))

cat("find the results in : ",home,"\n")








