# 01/2026 ; Fares Darawshy ; Segal Lab 
#  Signatures of Early Lung Cancer. RNA sequencing. Analysis by age 

#load libraries 
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
library(cluster)
library(umap)
###### NYU early stage data #########

# set your WD 
setwd("~/NYU Langone Health Dropbox/Fares Darawshy/Fares Darawshy’s files/Home/Projects/Early.Lung.Cancer/Analysis_by_Age/code_for_upload")

###### read data counts and metadata 

RNA.data <- read.csv(file = "NYU_stage_I_LUAD_meta_no_2nd_prim.csv")

mycounts <-read.delim2("NYU_stage_I_LUAD_counts_no_2nd_prim.txt", sep="\t", row.names=1)

#align tables and arrange them 
rownames(RNA.data) <- RNA.data$X
RNA.data <- RNA.data[,-1]
#
colnames(mycounts) <- gsub("X", "", colnames(mycounts))

#check difference 
setdiff(rownames(RNA.data), colnames(mycounts))
table(rownames(RNA.data)==colnames(mycounts))

#set RNA data for tumor and lung speratly 
RNA_data_tumor <- RNA.data %>% filter(Sample_Type_Involved=="Lung.Tissue.In")
RNA_data_lung <- RNA.data %>% filter(Sample_Type_Involved=="Lung.Tissue.UnIn")


#survival plot according to median age (Supp FIgure 1)
# Load the necessary libraries
library(survival)
library(survminer)
surv_data <- RNA_data_tumor %>% filter(ProgType_Lab!="Secondary.Primary")

# Create a survival object
surv_obj <- Surv(time = as.numeric(surv_data$New_TTP), event = as.numeric(surv_data$Progression))  # Make sure to adjust column names accordingly

# Create a Kaplan-Meier survival curve
surv_fit <- survfit(surv_obj ~ age_grp, data = surv_data)

# Plot the Kaplan-Meier curve
survival.plot <- ggsurvplot(
  surv_fit,
  data = surv_data,
  conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
  fun = "pct", ggtheme = theme_pubr(), 
  palette=c("darkblue", "darkgreen"), size=1,
  legend.labs=c("Old", "Young"), legend.title="Age",  
  # tables.col = "strata", #color of numbers inside risk table
  tables.height = 0.15,
  fontsize = 4,
  risk.table.y.text.col = T,
  cumcensor = TRUE, 
  tables.theme = theme_cleantable(),
  # risk.table.y.text = FALSE,
  risk.table.title = "Number at risk (cumulative number of recurrence)",
  cumcensor.title = "Cumulative number of censored subjects",
  #labs
  title="Kaplan-Meier Curve for Recurrence According to Age",
  xlab = "Time (Days)",
  ylab = "Recurrence Probability",
)

survival.plot





#####differential analysis of tumor and lung tissue (Figures 1 and Figure 3, NYU cohort)####

RNA.data$age_grp <- ifelse(RNA.data$Age <= median(RNA.data$Age), "less_equal_70", "greater_70")

meta <- data.frame(RNA.data)

##### create deseq object 
#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(mycounts, function(x) as.numeric(as.character(x))), check.names=F, row.names = rownames(mycounts))

#Convert Data to Integers to Run DESEq
d1[] <- lapply(d1, as.integer)

#Convert Model Variable into Factor
meta$age_grp <- factor(meta$age_grp)
meta$Sample_Type_Involved <- as.factor(meta$Sample_Type_Involved)
meta$Progression_Lab_Inv <- as.factor(meta$Progression_Lab_Inv)


dds <- DESeqDataSetFromMatrix(countData = d1, colData = meta, design = ~ age_grp + Progression_Lab_Inv)

#Normalization Step 
dds <- estimateSizeFactors(dds)

#Retrive normalized counts matrix 
normalized_counts <- counts(dds, normalized=TRUE)
#save it 
#write.table(normalized_counts, file="Results/normalized_counts_age_grp.txt", sep="\t", quote=F, col.names=NA)


#Filtering
#filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
idx <- rowSums( counts(dds, normalized=TRUE) >= 100 ) >= 3
dds <- dds[idx,]

#Transform Data
vsd <- varianceStabilizingTransformation(dds)

#Drop Levels
dds$age_grp   <- droplevels(dds$age_grp)
vsd$age_grp   <- droplevels(vsd$age_grp)



###############################################################################
################Running Differential Analysis old vs. young ##################

#run old vs young 

###Subset for Tumor

dds_tumor <- dds[,dds$Sample_Type_Involved=="Lung.Tissue.In"]
vsd_tumor <- vsd[,vsd$Sample_Type_Involved=="Lung.Tissue.In"]

#remove second primary 
dds_tumor <- dds_tumor[,dds_tumor$ProgType_Lab!="Secondary.Primary"]
vsd_tumor <- vsd_tumor[,vsd_tumor$ProgType_Lab!="Secondary.Primary"]


#### compare old vs young in tumor without second primary according to median age ####


#define variables 
dds_tumor$age_grp
v= "age_grp"

#Set Reference Level for Comparison (Control Group)
dds_tumor[[v]] <- factor(dds_tumor[[v]])
dds_tumor[[v]] <- relevel(dds_tumor[[v]], ref = "less_equal_70")

vsd_tumor[[v]] <- factor(vsd_tumor[[v]])
vsd_tumor[[v]] <- relevel(vsd_tumor[[v]], ref = "less_equal_70")


######################Plot PCoA


#Create Distance Matrix
vegdist = vegdist(t(assay(vsd_tumor)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd_tumor), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ age_grp,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="age_grp",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ age_grp, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray.Tumor_old_vs_young", paste0(v, paste0("_no_second_prim.txt")))
            , sep = "\t", row.names = T)


pdf(file = "Figures/RNA/Beta.Diversity.Bray.Tumor.samples.age_grp_RNA_no_second_prim.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= age_grp)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green4", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= age_grp), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= age_grp)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Young", "Old"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis

#first select Lung samples only from metadata and counts table 
mycounts_Tumor <- assay(dds_tumor)

#get genes table and add 1
x= as(mycounts_Tumor, "matrix")
x= x+1 

#get Lung only table 
RNA_data_analysis <- meta %>% filter(Sample_Type_Involved=="Lung.Tissue.In") %>% filter(ProgType_Lab!="Secondary.Primary")


# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA_data_analysis, "age_grp")

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

############continue edgeR here#############

res <- resNoFilt$table
#Reverse Directionality if you need to  
res$logFC <- res$logFC*(-1)

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
write.csv(res, file = "Results/RNA/edgeR.results_age_grp_tumor_no_second_primary.csv")

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "blue3"
cols[res$logFC < 0 & res$FDR < alpha ] <- "green4"


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR_Old_vs_Young_RNA_tumor_no_second_primary", paste0(alpha, paste0(".pdf"))), height = 8 , width = 6)
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

########## repeat for Lung subset : old vs young, removing second primary, age median
#run old vs young 

#Subset for lung

dds_lung <- dds[,dds$Sample_Type_Involved=="Lung.Tissue.UnIn"]
vsd_lung <- vsd[,vsd$Sample_Type_Involved=="Lung.Tissue.UnIn"]

#remove second primary 
dds_lung <- dds_lung[,dds_lung$ProgType_Lab!="Secondary.Primary"]
vsd_lung <- vsd_lung[,vsd_lung$ProgType_Lab!="Secondary.Primary"]


#### compare old vs young in lung without second primary according to median age ####

#define variables 
dds_lung$age_grp
v= "age_grp"

#Set Reference Level for Comparison (Control Group)
dds_lung[[v]] <- factor(dds_lung[[v]])
dds_lung[[v]] <- relevel(dds_lung[[v]], ref = "less_equal_70")

vsd_lung[[v]] <- factor(vsd_lung[[v]])
vsd_lung[[v]] <- relevel(vsd_lung[[v]], ref = "less_equal_70")


######################Plot PCoA


#Create Distance Matrix
vegdist = vegdist(t(assay(vsd_lung)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd_lung), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ age_grp,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="age_grp",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ age_grp, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray.lung_old_vs_young", paste0(v, paste0("_no_second_prim.txt")))
            , sep = "\t", row.names = T)


pdf(file = "Figures/RNA/Beta.Diversity.Bray.lung.samples.age_grp_RNA_no_second_prim.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= age_grp)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("green4", "blue3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= age_grp), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= age_grp)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Young", "Old"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis

#first select Lung samples only from metadata and counts table 
mycounts_lung <- assay(dds_lung)

#get genes table and add 1
x= as(mycounts_lung, "matrix")
x= x+1 

#get Lung only table 
RNA_data_analysis <- meta %>% filter(Sample_Type_Involved=="Lung.Tissue.UnIn") %>% filter(ProgType_Lab!="Secondary.Primary")


# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA_data_analysis, "age_grp")

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
res$logFC <- res$logFC*(-1)

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
write.csv(res, file = "Results/RNA/edgeR.results_age_grp_lung_no_second_primary.csv")

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "blue3"
cols[res$logFC < 0 & res$FDR < alpha ] <- "green4"


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR_Old_vs_Young_RNA_lung_no_second_primary", paste0(alpha, paste0(".pdf"))), height = 8 , width = 6)
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


#plot IPA heatmap. this is an exmaple of how you can generate heatmap. You should read the IPA data from supp table and re generate them using the following structure 

IPA_res <- read.csv(file = "IPA/IPA_comparison_NYU_early_stage_tumor_lung_old_vs_young_no_second_primary.csv")


#get max and min values 
max(IPA_res$Tumor, na.rm = TRUE)
min(IPA_res$Tumor, na.rm = TRUE)
max(IPA_res$Lung, na.rm = TRUE)
min(IPA_res$Lung, na.rm = TRUE)

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
col_fun = colorRamp2(c(-5, 0, 5.735), c("blue", "white", "orange"))
col_fun(seq(-5, 0, 5.735))

IPA_res <- IPA_res %>% arrange(desc(Lung))

#convert data to matrix 
IPA_res_mat <- IPA_res


#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

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
                        column_names_gp = gpar(fontace="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))



##### analyzing recurrence vs no recurrence without second primary, using median age as cutoff (Figure 1 and Figure 3 for NYU cohort ) ######


###Subset for Tumor

dds_tumor <- dds[,dds$Sample_Type_Involved=="Lung.Tissue.In"]
vsd_tumor <- vsd[,vsd$Sample_Type_Involved=="Lung.Tissue.In"]

#remove second primary 
dds_tumor <- dds_tumor[,dds_tumor$ProgType_Lab!="Secondary.Primary"]
vsd_tumor <- vsd_tumor[,vsd_tumor$ProgType_Lab!="Secondary.Primary"]

#subset old only 
dds_tumor_old <- dds_tumor[,dds_tumor$age_grp=="greater_70"]
vsd_tumor_old <- vsd_tumor[,vsd_tumor$age_grp=="greater_70"]


#### compare recurrence vs no recurrence in old in tumor without second primary according to median age ####

#change desgin 
dds.analysis <- dds_tumor_old
vsd.analysis <- vsd_tumor_old

design(dds.analysis) <- ~ Progression_Lab_Inv

#define variables 
dds.analysis$Progression_Lab_Inv
v= "Progression_Lab_Inv"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "In.No.Recurrence")

vsd.analysis[[v]] <- factor(vsd.analysis[[v]])
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "In.No.Recurrence")


######################Plot PCoA


#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Progression_Lab_Inv,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Progression_Lab_Inv",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ Progression_Lab_Inv, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray.Tumor_old_rec_vs_norec", paste0(v, paste0("_no_second_prim.txt")))
            , sep = "\t", row.names = T)


pdf(file = "Figures/RNA/Beta.Diversity.Bray.Tumor.samples.old_rec_vs_no_rec_RNA_no_second_prim.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Progression_Lab_Inv)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("orange", "red")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression_Lab_Inv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression_Lab_Inv)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis###################

#first select Lung samples only from metadata and counts table 
mycounts_Tumor <- assay(dds.analysis)

ncol(mycounts_Tumor)
#get genes table and add 1
x= as(mycounts_Tumor, "matrix")
x= x+1 

#get Lung only table 
RNA_data_analysis <- meta %>% filter(Sample_Type_Involved=="Lung.Tissue.In") %>% filter(ProgType_Lab!="Secondary.Primary") %>% filter(age_grp=="greater_70")


# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA_data_analysis, "Progression_Lab_Inv")

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

############continue edgeR here#############

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
write.csv(res, file = "Results/RNA/edgeR.results_old_rec_vs_no_rec_tumor_no_second_primary.csv")

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "red"
cols[res$logFC < 0 & res$FDR < alpha ] <- "orange"


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR_Old_rec_vs_no_rec_RNA_tumor_no_second_primary", paste0(alpha, paste0(".pdf"))), height = 8 , width = 6)
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


#repeat for young: rec vs no recu in tumor 

dds_tumor_young <- dds_tumor[,dds_tumor$age_grp=="less_equal_70"]
vsd_tumor_young <- vsd_tumor[,vsd_tumor$age_grp=="less_equal_70"]

dds.analysis <- dds_tumor_young
vsd.analysis <- vsd_tumor_young

design(dds.analysis) <- ~ Progression_Lab_Inv

#define variables 
dds.analysis$Progression_Lab_Inv
v= "Progression_Lab_Inv"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "In.No.Recurrence")

vsd.analysis[[v]] <- factor(vsd.analysis[[v]])
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "In.No.Recurrence")


######################Plot PCoA


#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Progression_Lab_Inv,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Progression_Lab_Inv",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ Progression_Lab_Inv, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray.Tumor_young_rec_vs_norec", paste0(v, paste0("_no_second_prim.txt")))
            , sep = "\t", row.names = T)


pdf(file = "Figures/RNA/Beta.Diversity.Bray.Tumor.samples.young_rec_vs_no_rec_RNA_no_second_prim.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Progression_Lab_Inv)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("skyblue2", "blue2")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression_Lab_Inv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression_Lab_Inv)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis

#first select Lung samples only from metadata and counts table 
mycounts_Tumor <- assay(dds.analysis)

ncol(mycounts_Tumor)
#get genes table and add 1
x= as(mycounts_Tumor, "matrix")
x= x+1 

#get Lung only table 
RNA_data_analysis <- meta %>% filter(Sample_Type_Involved=="Lung.Tissue.In") %>% filter(ProgType_Lab!="Secondary.Primary") %>% filter(age_grp=="less_equal_70")


# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA_data_analysis, "Progression_Lab_Inv")

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

############continue edgeR here#############

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
write.csv(res, file = "Results/RNA/edgeR.results_young_rec_vs_no_rec_tumor_no_second_primary.csv")

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "blue2"
cols[res$logFC < 0 & res$FDR < alpha ] <- "skyblue2"


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR_young_rec_vs_no_rec_RNA_tumor_no_second_primary", paste0(alpha, paste0(".pdf"))), height = 8 , width = 6)
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




###Subset for Lung###

dds_lung <- dds[,dds$Sample_Type_Involved=="Lung.Tissue.UnIn"]
vsd_lung <- vsd[,vsd$Sample_Type_Involved=="Lung.Tissue.UnIn"]

#remove second primary 
dds_lung <- dds_lung[,dds_lung$ProgType_Lab!="Secondary.Primary"]
vsd_lung <- vsd_lung[,vsd_lung$ProgType_Lab!="Secondary.Primary"]

#subset old only 
dds_lung_old <- dds_lung[,dds_lung$age_grp=="greater_70"]
vsd_lung_old <- vsd_lung[,vsd_lung$age_grp=="greater_70"]


#### compare recurrence vs no recurrence in old in lung without second primary according to median age ####

#change desgin 
dds.analysis <- dds_lung_old
vsd.analysis <- vsd_lung_old

design(dds.analysis) <- ~ Progression_Lab_Inv

#define variables 
dds.analysis$Progression_Lab_Inv
v= "Progression_Lab_Inv"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "UnIn.No.Recurrence")

vsd.analysis[[v]] <- factor(vsd.analysis[[v]])
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "UnIn.No.Recurrence")


######################Plot PCoA


#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Progression_Lab_Inv,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Progression_Lab_Inv",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ Progression_Lab_Inv, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray.lung_old_rec_vs_norec", paste0(v, paste0("_no_second_prim.txt")))
            , sep = "\t", row.names = T)


pdf(file = "Figures/RNA/Beta.Diversity.Bray.lung.samples.old_rec_vs_no_rec_RNA_no_second_prim.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Progression_Lab_Inv)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("orangered", "darkred")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression_Lab_Inv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression_Lab_Inv)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis###################

#first select Lung samples only from metadata and counts table 
mycounts_lung <- assay(dds.analysis)

ncol(mycounts_lung)
#get genes table and add 1
x= as(mycounts_lung, "matrix")
x= x+1 

#get Lung only table 
RNA_data_analysis <- meta %>% filter(Sample_Type_Involved=="Lung.Tissue.UnIn") %>% filter(ProgType_Lab!="Secondary.Primary") %>% filter(age_grp=="greater_70")


# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA_data_analysis, "Progression_Lab_Inv")

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

############continue edgeR here

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
write.csv(res, file = "Results/RNA/edgeR.results_old_rec_vs_no_rec_lung_no_second_primary.csv")

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "darkred"
cols[res$logFC < 0 & res$FDR < alpha ] <- "orangered"


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR_Old_rec_vs_no_rec_RNA_lung_no_second_primary", paste0(alpha, paste0(".pdf"))), height = 8 , width = 6)
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


#repeat for young: rec vs no recu in lung 

dds_lung_young <- dds_lung[,dds_lung$age_grp=="less_equal_70"]
vsd_lung_young <- vsd_lung[,vsd_lung$age_grp=="less_equal_70"]

dds.analysis <- dds_lung_young
vsd.analysis <- vsd_lung_young

design(dds.analysis) <- ~ Progression_Lab_Inv

#define variables 
dds.analysis$Progression_Lab_Inv
v= "Progression_Lab_Inv"

#Set Reference Level for Comparison (Control Group)
dds.analysis[[v]] <- factor(dds.analysis[[v]])
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "UnIn.No.Recurrence")

vsd.analysis[[v]] <- factor(vsd.analysis[[v]])
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "UnIn.No.Recurrence")


######################Plot PCoA


#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Progression_Lab_Inv,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Progression_Lab_Inv",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ Progression_Lab_Inv, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray.lung_young_rec_vs_norec", paste0(v, paste0("_no_second_prim.txt")))
            , sep = "\t", row.names = T)


pdf(file = "Figures/RNA/Beta.Diversity.Bray.lung.samples.young_rec_vs_no_rec_RNA_no_second_prim.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Progression_Lab_Inv)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("skyblue", "darkblue")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Progression_Lab_Inv), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Progression_Lab_Inv)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis

#first select Lung samples only from metadata and counts table 
mycounts_lung <- assay(dds.analysis)

ncol(mycounts_lung)
#get genes table and add 1
x= as(mycounts_lung, "matrix")
x= x+1 

#get Lung only table 
RNA_data_analysis <- meta %>% filter(Sample_Type_Involved=="Lung.Tissue.UnIn") %>% filter(ProgType_Lab!="Secondary.Primary") %>% filter(age_grp=="less_equal_70")


# get your group variable (the condition you want to analyze according to it)
group = get_variable(RNA_data_analysis, "Progression_Lab_Inv")

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
write.csv(res, file = "Results/RNA/edgeR.results_young_rec_vs_no_rec_lung_no_second_primary.csv")

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- "darkblue"
cols[res$logFC < 0 & res$FDR < alpha ] <- "skyblue"


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR_young_rec_vs_no_rec_RNA_lung_no_second_primary", paste0(alpha, paste0(".pdf"))), height = 8 , width = 6)
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






#plot IPA heatmap from IPA results 

IPA_res <- read.csv(file = "IPA/IPA_comparison_NYU_early_tumor_lung_rec_vs_norec_age_no_second_prim.csv")


#get max and min values 
max(IPA_res$Old_Tumor, na.rm = TRUE)
min(IPA_res$Old_Tumor, na.rm = TRUE)
max(IPA_res$Young_Tumor, na.rm = TRUE)
min(IPA_res$Young_Tumor, na.rm = TRUE)
max(IPA_res$Old_Lung, na.rm = TRUE)
min(IPA_res$Old_Lung, na.rm = TRUE)
max(IPA_res$Young_Lung, na.rm = TRUE)
min(IPA_res$Young_Lung, na.rm = TRUE)

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
col_fun = colorRamp2(c(-2.887, 0, 5.578), c("blue", "white", "orange"))
col_fun(seq(-2.887, 0, 5.578))

#convert data to matrix 
IPA_res_mat <- IPA_res

IPA_res_mat <- IPA_res_mat %>% arrange(desc(Old_Lung))

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

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
                        column_names_gp = gpar(fontace="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))


#subset for tumor only (FIgure 1)

IPA_res <- read.csv(file = "IPA/IPA_comparison_NYU_early_tumor_rec_vs_norec_age_no_second_prim.csv")

#get max and min values 
max(IPA_res$Old_Tumor, na.rm = TRUE)
min(IPA_res$Old_Tumor, na.rm = TRUE)
max(IPA_res$Young_Tumor, na.rm = TRUE)
min(IPA_res$Young_Tumor, na.rm = TRUE)

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
col_fun = colorRamp2(c(-1.633, 0, 4.914), c("blue", "white", "orange"))
col_fun(seq(-1.633, 0, 4.914))

#convert data to matrix 
IPA_res <- IPA_res %>% arrange(desc(Old_Tumor))

IPA_res_mat <- IPA_res


#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

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
                        column_names_gp = gpar(fontace="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))


#####lung subset (Figure 3)

IPA_res <- read.csv(file = "IPA/IPA_comparison_NYU_early_lung_rec_vs_norec_age_no_second_prim.csv")


#get max and min values 
max(IPA_res$Old_Lung, na.rm = TRUE)
min(IPA_res$Old_Lung, na.rm = TRUE)
max(IPA_res$Young_Lung, na.rm = TRUE)
min(IPA_res$Young_Lung, na.rm = TRUE)

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
col_fun = colorRamp2(c(-3.162, 0, 6.272), c("blue", "white", "orange"))
col_fun(seq(-3.162, 0, 6.272))

#convert data to matrix 
IPA_res_mat <- IPA_res

IPA_res_mat <- IPA_res_mat %>% arrange(desc(Old_Lung))

#set pathways as rownames
rownames(IPA_res_mat) <- IPA_res_mat$Pathway

#get rid of extra columns
IPA_res_mat <- IPA_res_mat[,-1]

#convert to matrix 
IPA_res_mat <- as.matrix(IPA_res_mat)

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
                        column_names_gp = gpar(fontace="bold"),
                        row_title_side = "left", 
                        column_title_side = "top", 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE, 
                        heatmap_legend_param = list(title="Z-Score", title_position = "lefttop-rot"),
                        row_labels =IPA_res$Pathway, 
                        width = unit(2, "cm"))



