#!/usr/bin/env Rscript
##!/usr/bin/env /bioinfo/local/build/R/R-3.3.0_centos/bin/Rscript

args <- commandArgs(TRUE)


library(gridExtra)
library(ggplot2)
library(ggpubr)
library(egg)
library(gtools)

#all the RData files (comma separated)
ggplot_objects<-args[1]
output_dir<-args[2]
prefix<-args[3]
feature<-args[4]
metagene_tables<-args[5]



ggplot_objects<-as.list(unlist(strsplit(ggplot_objects,",")))

metagene_tables<-unlist(strsplit(metagene_tables,","))


setwd(output_dir)

#### this step is not needed anymore
# all_max<-c()
# 
# ggplot_list<-list()
# 
# 
# for(i in 1:length(ggplot_objects)){
#   
#   #one_ggplot<-
#   load(ggplot_objects[[i]])
#   
#   #extract the legend (ggpubr::ggarrange is buggy)
#   legend <- cowplot::get_legend(ggplot_img_to_save+
#                                   guides(color = guide_legend(nrow = 1)) +
#                                   theme(legend.position = "top"))
#   ggplot_list[[length(ggplot_list)+1]]<-ggplot_img_to_save
#   
#   all_max<-c(all_max,max(ggplot_build(ggplot_img_to_save)$data[[1]]$y))
#   
# }
# 
# my_max<-max(all_max)
# 
# cat("my max is : ",my_max,"\n")
# 
# cat ("number of ggplot objects : ",length(ggplot_list),"\n")


######  recreate plot from the dataframes of each run #################################

all_files<-gtools::mixedsort(metagene_tables)

cat(all_files,"\n")

full_dataframe<-data.frame()
for(i in 1:length(all_files)){
  
  full_dataframe<-rbind(full_dataframe,read.delim(all_files[i],header=T,sep="\t"))
  
}

full_dataframe$condition <- factor(full_dataframe$condition, levels = unique(as.character(full_dataframe$condition)))
full_dataframe$feature <- factor(full_dataframe$feature, levels = unique(as.character(full_dataframe$feature)))
full_dataframe$project <- factor(full_dataframe$project, levels = unique(as.character(full_dataframe$project)))

full_dataframe$project<-paste(full_dataframe$project,paste("n = ",full_dataframe$extra,sep=""),sep="\n")

full_dataframe$project <- factor(full_dataframe$project, levels = unique(as.character(full_dataframe$project)))

both_conds<-c("cornflowerblue","red","cyan","deeppink","#3eb489","orange","gold")

names(both_conds)<-unique(full_dataframe$condition)


white_background<-ggplot2::theme(axis.line = element_line(colour = "black"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.background = element_blank())




my_file_name<-paste(output_dir,"metagene_plot_multipanel_",prefix,"_",feature,".png",sep="")

#pdf(file=paste(output_dir,"metagene_plot_multipanel_",prefix,"_",feature,"_model2",".pdf",sep=""),width=15,height=10)

png(filename=my_file_name,width=1800,height=1000)

print(
  
  ggplot(full_dataframe, aes(pos, log10(count+1),color=condition))+
    geom_line(size=1)+
    geom_ribbon(aes(ymin=log10(the_min+1), ymax=log10(the_max+1),fill=condition), alpha=0.2)+
    xlab("")+
    #with _smooth, we go under 0 sometimes ...find a way to avoid this
    #stat_smooth(aes(x = pos, y = count), formula = y ~ s(x, k =100), method = "gam", se = FALSE,size=2)+
    scale_y_continuous(expand=c(0,0))+
    
    scale_x_continuous(expand=c(0,0))+
    geom_hline(yintercept = 0,color="black",size=1)+
   
    theme(axis.text.x=element_text(color="black",size=16,hjust=1,face="bold",angle = 45),axis.text.y = element_text(color="black",size=25,face="bold"),plot.title = element_text(hjust = 0.5,size=18),axis.title=element_text(size=25,face="bold"),legend.text=element_text(size=20),legend.title=element_text(size=18,face="bold"))+
    scale_color_manual(values=c(both_conds),name="condition")+
    scale_fill_manual(values=c(both_conds),name="condition")+

    ylab("Mean normalized coverage\nlog10(rpm+1)")+
    
    white_background+
    geom_vline(xintercept =c(1,100), linetype="dashed",color="grey40")+
    facet_grid(project~feature,scales="free_y")+
    #
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5,
                                          linetype="solid"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_blank(),
          strip.text.y =element_text(size=16)
    )+
    
   
    coord_cartesian( clip="off")+
    
    
  geom_text(data=full_dataframe[which(full_dataframe$project==unique(full_dataframe$project)[length(unique(full_dataframe$project))]),],
            x=50,
            y=-(log10(max(full_dataframe$count)+1)/20),
            label=full_dataframe[which(full_dataframe$project==unique(full_dataframe$project)[length(unique(full_dataframe$project))]),]$feature,color="black",
            size=8)+
    
    #expand the bottom border
    theme(plot.margin = unit(c(1,1,4,1), "lines"))
  
  
)


dev.off()

cat("\ncheck file ",my_file_name,"\n\n",sep="")

