########### TRACERx proccessing############
###### written by Fares Darawshy - fares.darawshy@nyulangone.org

library(fgsea)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(ggalt)
library(forcats)
library(car)
library(magrittr)
library(lubridate)
library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(gplots)
library(ggrepel)
library(pathfindR)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(omu)
library(maptools)
library(phyloseq)
library("vegan")
library(ade4)
library("reshape2")
library(dplyr)	
library(Matrix)
library(scales)
library(cowplot)
library(randomForest)
library(caret)
library("mlbench")
library(RCurl)
library(VennDiagram)
library(ranger)
library(xgboost)
library(rmarkdown)
library(glue)
library(ggtext)
library(ComplexHeatmap)
library(TCGAbiolinks)
library(table1)
library(pROC)
library(survival)
library(survminer)
library(pdftools)
library(grid)
library(gridExtra)
library(fst)

# set your working directory 
setwd("~/NYU Langone Health Dropbox/Fares Darawshy/Fares Darawshy’s files/Home/Projects/TRACERx")

##### varaibles were adjusted and created for this reserach from TRACERx orignial data. 
#### all secondary primary recurrence events were removed from analyiss 
#### only patients with LUAD stage I and without secondary primary were included. 
### 171 samples from tumor overall (>1 sample per patient, that were sum when needed)

################################################################################
################### analysis for stage I LUAD only ####################

##### comparison old vs young 

#read metadata 
tumor_data <- read.csv(file = "TRACER_tumor_LUAD_Stage_I_meta.csv")
rownames(tumor_data) <- tumor_data$X
tumor_data <- tumor_data[,-1]

#read counts 
mycounts_Tumor <- read.delim2(file = "TRACERx_tumor_LUAD_stage_I_counts.txt")
rownames(mycounts_Tumor) <- mycounts_Tumor$X
mycounts_Tumor <- mycounts_Tumor[,-1]

#check difference 
setdiff(rownames(tumor_data), colnames(mycounts_Tumor))
#correct colnames 
colnames(mycounts_Tumor) <- gsub("(T\\d+)[.-](R\\d+)", "\\1-\\2", colnames(mycounts_Tumor))
#check again
setdiff(rownames(tumor_data), colnames(mycounts_Tumor))
#check they align 
table(rownames(tumor_data)==colnames(mycounts_Tumor))

#create dds object 
#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts_Tumor, function(x) as.numeric(as.character(x))), check.names=F, row.names = rownames(mycounts_Tumor))

#Convert Data to Integers to Run DESEq
d1[] <- lapply(d1, as.integer)

dds <- DESeqDataSetFromMatrix(countData = d1, colData = tumor_data, design = ~ age_grp)

#Normalization Step 
dds <- estimateSizeFactors(dds)

#Retrive normalized counts matrix 
normalized_counts <- counts(dds, normalized=TRUE)
#save it 
#write.table(normalized_counts, file="RNA_analysis/Results/normalized_counts.txt", sep="\t", quote=F, col.names=NA)


#Remove genes with normalized read counts below 5 in 20% of the samples 
norm_counts <- counts(dds, normalized = TRUE)
samples_thresh <- round(0.2*ncol(norm_counts))
genes_to_keep <- norm_counts > 5
genes_to_keep <- rownames(genes_to_keep)[rowSums(genes_to_keep) >= samples_thresh]
dds_tumour <- dds[genes_to_keep, ]

#Transform Data
vsd <- varianceStabilizingTransformation(dds_tumour)


dds_tumour <- dds_tumour
vsd_tumour <- vsd

###############table 1 analysis ########

#sub table that include LUAD and stage I only without secondary primary 

Table.1.data <- tumor_data %>%   distinct(patient, .keep_all = TRUE)

table1(~ 
         #set variables you want to display. Categorical as factors, continous as numeric 
         factor(clinical_sex) + factor(Histology) + as.numeric(age) + factor(luad_subtype) 
       + factor(ethnicity) +factor(Histology_per_tumour_id_muttable) + factor(LUAD_pred_subtype_with.IMA_per_tumour) + 
         factor(site_per_lesion) + factor(vascular_invasion_per_lesion)+factor(pleural_invasion_per_lesion)+
         factor(smoking_status_merged)+as.numeric(cigs_perday)+as.numeric(pack_years) +factor(is.family.lung)+factor(ECOG_PS)+
         as.numeric(size_pathology_per_patient)+factor(Surgery_type)+factor(adjuvant_treatment_YN)+factor(adjuvant_treatment_given)+
         factor(num_cycle_na.added)+factor(first_dfs_any_event_rec.or.new.primary)+factor(first_event_during_followup)+as.numeric(dfs_time_any_event)+
         factor(Relapse_cat_new)+factor(stage) +factor(death_y_n)+five_y_time
       # set the variable you want to stratify by 
       | age_grp, 
       #set data 
       data=Table.1.data, 
       #set stat display options for continous variables (options from stat.default)
       render.continuous = c((.="Median [Q1, Q3]")))
#now add stats (same as above except where notes added ) for these two groups only 
table1(~ 
         factor(clinical_sex) + factor(Histology) + as.numeric(age) + factor(luad_subtype) 
       + factor(ethnicity) +factor(Histology_per_tumour_id_muttable) + factor(LUAD_pred_subtype_with.IMA_per_tumour) + 
         factor(site_per_lesion) + factor(vascular_invasion_per_lesion)+factor(pleural_invasion_per_lesion)+
         factor(smoking_status_merged)+as.numeric(cigs_perday)+as.numeric(pack_years) +factor(is.family.lung)+factor(ECOG_PS)+
         as.numeric(size_pathology_per_patient)+factor(Surgery_type)+factor(adjuvant_treatment_YN)+factor(adjuvant_treatment_given)+
         factor(num_cycle_na.added)+factor(first_dfs_any_event_rec.or.new.primary)+factor(first_event_during_followup)+as.numeric(dfs_time_any_event)+
         factor(Relapse_cat_new)+factor(stage)+factor(death_y_n)+five_y_time
       | age_grp, 
       data=Table.1.data, 
       # don't display overall column 
       overall = F,
       render.continuous = c((.="Median [Q1, Q3]")),
       #add p value as extra column (defined in function above)
       extra.col = list("P-value" = pvalue))

###### stats for numerical variables 
chisq.test(Table.1.data$age, Table.1.data$age_grp)
chisq.test(Table.1.data$cigs_perday, Table.1.data$age_grp)
chisq.test(Table.1.data$pack_years, Table.1.data$age_grp)
chisq.test(Table.1.data$size_pathology_per_patient, Table.1.data$age_grp)
chisq.test(Table.1.data$dfs_time_any_event, Table.1.data$age_grp)
chisq.test(Table.1.data$os_time, Table.1.data$age_grp)
chisq.test(Table.1.data$lung_event_time, Table.1.data$age_grp)




#that leaves 171 samples and 83 patients for analysis 
#choose counts
mycounts_Tumor <- assay(dds_tumour)
mycounts_Tumor <- mycounts_Tumor[,tumor_data$region]
ncol(mycounts_Tumor)
#same for vst 
vst_Tumor <- assay(vsd_tumour)


#### compare old vs young in tumor#####

#choose variable 
v= "age_grp"

#Set Reference Level for Comparison (Control Group)
tumor_data[[v]] <- factor(tumor_data[[v]])
tumor_data[[v]] <- relevel(tumor_data[[v]], ref = "young")

#create appropraite metadata for edgeR
metadata_RNA_analysis <- tumor_data

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "blue3"
col2 <- "green4"

######################Plot PCoA

#define variable 
tumor_data[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(vst_Tumor), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = tumor_data, by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$age_grp
newResults$v <- factor(newResults$v, levels = c("young", "old"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)

####save your figure and output it 
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Young ", "Old"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")



##############Differential Analysis 
mycounts.analysis <- mycounts_Tumor

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Results/edgeR.results_", paste0(v, paste0("_tumor_stage_I_no_second_primary.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_tumor_stage_I_no_second_prim.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_tumor_stage_I_no_second_prim.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()

#######comparison of rec vs no rec in old, tumor samples ########

#create a new recurrence variable
tumor_data$Relapse_cat_y_n <- ifelse(tumor_data$Relapse_cat_new=="No rec", "no_rec", "recurrence")

#choose only old 
tumor_data_old <- tumor_data %>% filter(age_grp=="old")

#choose counts
mycounts_Tumor <- assay(dds_tumour)
mycounts_Tumor <- mycounts_Tumor[,tumor_data_old$region]

#same for vst 
vst_Tumor <- assay(vsd_tumour)[,tumor_data_old$region]

#choose variable 

v= "Relapse_cat_y_n"

#Set Reference Level for Comparison (Control Group)
tumor_data_old[[v]] <- factor(tumor_data_old[[v]])
tumor_data_old[[v]] <- relevel(tumor_data_old[[v]], ref = "no_rec")

#create appropraite metadata for edgeR
metadata_RNA_analysis <- tumor_data_old

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "red"
col2 <- "orange"

######################Plot PCoA

#define variable 
tumor_data_old[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(vst_Tumor), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = tumor_data_old, by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$Relapse_cat_y_n
newResults$v <- factor(newResults$v, levels = c("no_rec", "recurrence"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
#write.table(x, file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Results/Beta.Diversity.Bray.", paste0(v, paste0("_tumor_old_stage_I_no_second_prim.txt")))
#            , sep = "\t", row.names = T)


pdf(file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_tumor_old_stage_I_no_second_prim.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence ", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis 
mycounts.analysis <- mycounts_Tumor

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Results/edgeR.results_", paste0(v, paste0("_tumor_old_stage_I_no_second_prim.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_tumor_old_stage_I_no_second_prim.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_tumor_old_stage_I_no_second_prim.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()





#######comparison of rec vs no rec in young #########
tumor_data$Relapse_cat_y_n <- ifelse(tumor_data$Relapse_cat_new=="No rec", "no_rec", "recurrence")

#choose only young 
tumor_data_young <- tumor_data %>% filter(age_grp=="young")

#choose counts
mycounts_Tumor <- assay(dds_tumour)
mycounts_Tumor <- mycounts_Tumor[,tumor_data_young$region]

#same for vst 
vst_Tumor <- assay(vsd_tumour)[,tumor_data_young$region]


#choose variable 

v= "Relapse_cat_y_n"

#Set Reference Level for Comparison (Control Group)
tumor_data_young[[v]] <- factor(tumor_data_young[[v]])
tumor_data_young[[v]] <- relevel(tumor_data_young[[v]], ref = "no_rec")

#create appropraite metadata for edgeR
metadata_RNA_analysis <- tumor_data_young

#set colors, col1 is for log FC > 0, and col2 if log FC < 0
col1 <- "blue"
col2 <- "lightblue"

######################Plot PCoA

#define variable 
tumor_data_young[[v]] 
#Create Distance Matrix
vegdist = vegdist(t(vst_Tumor), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = tumor_data_young, by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

newResults$v <- newResults$Relapse_cat_y_n
newResults$v <- factor(newResults$v, levels = c("no_rec", "recurrence"))

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean) 
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
#write.table(x, file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Results/Beta.Diversity.Bray.", paste0(v, paste0("_tumor_young_stage_I_no_second_prim.txt")))
#            , sep = "\t", row.names = T)


pdf(file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Figures/Beta.Diversity.Bray_", paste0(v, paste0("_tumor_young_stage_I_no_second_prim.pdf"))), 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= v)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c(col2, col1)) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= v), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= v)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence ", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


##############Differential Analysis 
mycounts.analysis <- mycounts_Tumor

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_RNA_analysis, v)

#create edgeR list 
dgeFull <- DGEList(counts = x, 
                   group = group, 
                   remove.zeros = TRUE)

# calculatee normalization method
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
# Estimate dispersions
dgeFull <- estimateCommonDisp(dgeFull)

dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

#create table
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table), adjust.method = "BH", sort.by = "PValue")

res <- resNoFilt$table
#Reverse Directionality if you need to  
#res$logFC <- res$logFC*(-1)

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#create gene symbol column 
res<- res %>% 
  dplyr::mutate(Gene.symbol = rownames(.))

# save results (can be used in IPA also)
write.csv(res, file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Results/edgeR.results_", paste0(v, paste0("_tumor_young_stage_I_no_second_prim.csv"))))

#save for GSEA a ranked file
# create rank file 
x <-
  res %>% 
  dplyr::filter(!is.na(FDR)) 
x$fcSign<-sign(x$logFC)
x$logP   <- -log10(x$FDR)
x$metric <- x$logP/x$fcSign

y<-x %>% 
  dplyr::select(Gene.symbol, metric) %>% 
  dplyr::rename(SYMBOL=Gene.symbol)
baseline.rnk=y
write_tsv(y, file= paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Results/", paste0("GSEA_edgeR_", paste0(v, paste0("_tumor_young_stage_I_no_second_prim.rnk")))), col_names=FALSE)


###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("~/Dropbox (NYU Langone Health)/Fares Darawshy’s files/Home/Projects/TRACERx/Figures/edgeR_", paste0(v, paste0(alpha, paste0("_tumor_young_stage_I_no_second_prim.pdf")))), height = 8 , width = 6)
ggplot(res, aes(x= logFC, y=sig, label= Gene.symbol))+ 
  geom_point(color = cols)+
  geom_text_repel(label=ifelse(res$logFC>0 & res$FDR < alpha , as.character(res$Gene.symbol),
                               ifelse(res$logFC<0 & res$FDR < alpha, as.character(res$Gene.symbol),'')),
                  size= 4, 
                  force = 25,
                  segment.colour="grey",
                  segment.alpha=0.5) + #Label values based on parameters, including pcal and logFC
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
  xlab("Effect size: log2(fold-change)") + #label X Axis
  ylab("-log10(adjusted p-value)")+  #label Y Axis
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
        axis.text.y=element_text(colour = "grey80", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()


################################################################################
###################### IPA heatmaap for tumor : rec vs no rec in old and young, using age_grp no second primary #######
### this code is just an example to produce a heatmap. you should clean IPA resuls first 

#read IPA results as provided in the supplementary table 
IPA_res <- read.csv(file = "your_IPA_results_sheet.csv")
#####plot heat map using complex heat maps 

# Define colors for each levels of qualitative variables
library(circlize)
library(ComplexHeatmap)
max(IPA_res$Old, na.rm = TRUE)
min(IPA_res$Old,na.rm = TRUE)
max(IPA_res$Young,na.rm = TRUE)
min(IPA_res$Young, na.rm = TRUE)

col_fun = colorRamp2(c(-3.725, 0, 11.147), c("blue", "white", "orange"))
col_fun(seq(-3.725, 0, 11.147))

#convert data to matrix 
IPA_res_mat <- IPA_res

#sort 
IPA_res_mat <- IPA_res %>% arrange(desc(Old))

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

library(ComplexHeatmap)
######ploting
#na colors as white 
#add lines between cells 
#bold rownames
#set size 
ComplexHeatmap::Heatmap(IPA_res_mat, 
                        rect_gp = gpar(col = "black", lwd = 1),
                        na_col = "white",
                        row_names_side = "left",
                        column_names_side = "top",
                        column_title = "", row_title = "", 
                        col = col_fun, 
                        row_names_gp = gpar(fontface="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))




####comparison of old vs young in IPA stage I no second priamry ###
IPA_res <- read.csv(file = "IPA_res_file_old_vs_young")

IPA_res$Pathway <- factor(IPA_res$Pathway, levels = unique(IPA_res$Pathway))
IPA_res <- IPA_res %>% arrange(desc(z.score)) %>% dplyr::filter(z.score != 0) %>% filter(z.score!="#NUM!")
IPA_res$neg_log_p_value <- IPA_res$X.log.p.value. 

#leave only pathways with p value < 0.2 (FDR)
IPA_res <- IPA_res %>% filter(neg_log_p_value > -log(0.05))

#set colors of pthways. oragnge is upregulatedm, blue is downregulated 
IPA_res <- IPA_res %>%  mutate(col=ifelse(z.score >0, "blue3", "green4"))

ggplot(IPA_res, aes(x=as.numeric(z.score), y=fct_reorder(Pathway,as.numeric(z.score), .fun = max), fill=col))+
  geom_col()+
  scale_fill_manual(values = c("blue3", "green4"))+
  xlab("Z Score")+ylab("")+
  guides(fill="none")+
  ggtitle("Top Upregulated Pathways in patients with Dead vs Aalive Patients")+
  theme_classic()+
  theme(axis.title=element_text(size=28,face="bold"), 
        axis.text.x = element_text(size = 24, face = "bold"), 
        axis.text.y = element_text(size = 24, face = "bold"))


######## save your priect 