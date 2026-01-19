# 05/08/23 ; Fares Darawshy ; Segal Lab 
#  Signatures of Early Lung Cancer. RNA sequencing. Analysis by age 

setwd("____set your wd")
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


############################# TCGA Cohort analysis ########################


####### you can either download TCGA data as detailed in the following, 
####### OR use the following files in the repository and skip directly to downbstreama analysis 



############# download and prepare TCGA data #####

#install TCGA package 
BiocManager::install("TCGAbiolinks")
#load the package 
library(TCGAbiolinks)

#load the data you want from TCGA 
proj <- "TCGA-LUAD"
query <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
#download the files you want 
GDCdownload(query)
data <- GDCprepare(query)

#get counts and metadata from summarized experiment object 
tcga_counts <- data.frame(assay(data))

#get metadata
tcga_metadata <- data.frame(data@colData)

#make counts table rownames genes names and not gene id 
names <- make.unique(data@rowRanges$gene_name)
rownames(tcga_counts) <- names

#replace metadata rownames format to match them to the columns of counts table 
rownames(tcga_metadata) <- gsub("-", ".", rownames(tcga_metadata))

#make sure counts table columns and rows of metadata are the same 
table(colnames(tcga_counts)==rownames(tcga_metadata))

#get information from metadata 
colnames(tcga_metadata)

#add clinical information you need from other data downloaded form www.cbioportal.org 
#read the file downloaded 
tcga_progression_data <- read.delim2("TCGA/KM_Plot__Disease_Free_(months).txt", header = TRUE)
tcga_progression_data<- tcga_progression_data %>% 
  dplyr::rename(patient=Patient_ID)

#add progression data to metadata 
tcga_metadata$patient
tcga_progression_data$patient
#before joining, get rownames of tcga data so you don't lose them and use them later 
tcga_rownames <- rownames(tcga_metadata)

#join and correct rownames 
tcga_metadata <- inner_join(tcga_metadata, tcga_progression_data, by="patient")
rownames(tcga_metadata) <- tcga_rownames

tcga_metadata$definition <- factor(tcga_metadata$definition)
tcga_metadata %>% 
  dplyr::group_by(definition) %>% 
  dplyr::summarise(n=dplyr::count(.))

#we have 539 primary tumors, 2 recurrent tumor and 59 normal tissue 
#will remove 2 recurrent tumor from analsys which wil leave out 598 samples  
tcga_metadata <- tcga_metadata %>% 
  dplyr::filter(definition != "Recurrent Solid Tumor")


#filter tcga counts so thery include only human genes (For these, use TCGA coutns ifle or NYU counts file, and call it mycoounts)
mycounts <- read.delim2("NYU/ TCGA counts file")
tcga_counts_filtered <- tcga_counts %>% 
  dplyr::filter(rownames(.)%in%rownames(mycounts))

nrow(mycounts)
nrow(tcga_counts_filtered)


#remove those from counts table 
tcga_counts_filtered <- tcga_counts_filtered %>% 
  dplyr::select(rownames(tcga_metadata))

#finally validate that count and clinical data table match 
table(rownames(tcga_metadata)==colnames(tcga_counts_filtered))

#create 5 year survival columns 
tcga_metadata$five_year <- (tcga_metadata$year_of_diagnosis +5)
tcga_metadata$year_of_death
tcga_metadata$dead_or_not_five_years <- ifelse(tcga_metadata$year_of_death < tcga_metadata$five_year, "dead", "alive")
tcga_metadata$dead_or_not_five_years[is.na(tcga_metadata$dead_or_not_five_years)] <- "alive"

#correct recurrence column 
tcga_metadata$DFS_STATUS <- gsub("0:DiseaseFree", "No.Recurrence",tcga_metadata$DFS_STATUS )
tcga_metadata$DFS_STATUS <- gsub("1:Recurred/Progressed", "Recurrence",tcga_metadata$DFS_STATUS)

#create smoking variable 
tcga_metadata$smoking_ever <- ifelse(tcga_metadata$cigarettes_per_day =="NA", "never", "ever")






######## you alterntaively can start from here by reading the files 
tcga_metadata <- read.csv("TCGA_meta_stage_I_LUAD.csv")
tcga_counts_filtered <- read.delim2("TCGA_counts_stage_I_LUAD.txt")


#set age cutoff = 70
tcga_metadata_stage_I

#get clinical metadata of early stage (overall n=325 samples)
tcga_metadata_stage_I <- tcga_metadata %>% 
  dplyr::filter(ajcc_pathologic_stage %in% c("Stage_I", "Stage_IA", "Stage_IB"))

#get the count table accordingly 
tcga_counts_stage_I <- tcga_counts_filtered %>% 
  dplyr::select(rownames(tcga_metadata_stage_I))

ncol(tcga_counts_stage_I) == nrow(tcga_metadata_stage_I)

############## Create summarized experiment object #########################

#Convert Model Variable into Factor
tcga_metadata_stage_I$age_grp_new <- ifelse(tcga_metadata_stage_I$age_at_index<=70, "less_equal_70", "greater_70")

tcga_metadata_stage_I$age_grp <- factor(tcga_metadata_stage_I$age_grp)
tcga_metadata_stage_I$DFS_STATUS <- factor(tcga_metadata_stage_I$DFS_STATUS)

#remove any missing data without age or recurrence information (overall 301 for analysis)
tcga_metadata_stage_I_no_NA <- tcga_metadata_stage_I %>% 
  dplyr::filter(age_grp_new != "NA") %>% 
  dplyr::filter(DFS_STATUS != "NA")


#now select the final counts table 
tcga_counts_stage_I_no_NA <- tcga_counts_filtered %>% 
  dplyr::select(rownames(tcga_metadata_stage_I_no_NA))


#Deseq 
#Convert Count Table into a Numeic Data Frame
d1 = data.frame(lapply(tcga_counts_stage_I_no_NA, function(x) as.numeric(as.character(x))), check.names=F, row.names = rownames(tcga_counts_stage_I_no_NA))

#Convert Data to Integers to Run DESEq
d1[] <- lapply(d1, as.integer)


dds <- DESeqDataSetFromMatrix(countData = d1, colData = tcga_metadata_stage_I_no_NA, design = ~ age_grp_new + DFS_STATUS)

#Normalization Step 
dds <- estimateSizeFactors(dds)

#Retrive normalized counts matrix 
normalized_counts <- counts(dds, normalized=TRUE)
#save it 
#write.table(normalized_counts, file="Results/RNA/normalized_counts_tcga_stage_I_median_age_new.txt", sep="\t", quote=F, col.names=NA)


#Filtering
#filter out genes where there are less than 3 samples with normalized counts greater than or equal to 100.
idx <- rowSums( counts(dds, normalized=TRUE) >= 100 ) >= 3
dds <- dds[idx,]

#Transform Data
vsd <- varianceStabilizingTransformation(dds)

#Drop Levels
dds$age_grp_new   <- droplevels(dds$age_grp_new)
vsd$age_grp_new   <- droplevels(vsd$age_grp_new)

dds$DFS_STATUS   <- droplevels(dds$DFS_STATUS)
vsd$DFS_STATUS   <- droplevels(vsd$DFS_STATUS)

#subset tumor and lung objects 
dds_tumor <- dds[,dds$sample_type=="Primary Tumor"]
dds_lung <- dds[, dds$sample_type=="Solid Tissue Normal"]

vsd_tumor <- vsd[,vsd$sample_type=="Primary Tumor"]
vsd_lung <- vsd[, vsd$sample_type=="Solid Tissue Normal"]



######### Table 1 of early lung cancer by median age #######


#get tabel 1 information 

Table.1.data <- tcga_metadata_stage_I_no_NA

#filter for number of subjects only (363 patients)
Table.1.data.subjects <- Table.1.data %>% 
  distinct(patient,.keep_all = TRUE)

#calculate median age = 66 
median(Table.1.data.subjects$age_at_index, na.rm = TRUE)

#create table 1 using table 1 package 
table1(~ 
         #set variables you want to display. Categorical as factors, continous as numeric 
         factor(synchronous_malignancy) + factor(ajcc_pathologic_stage) + as.numeric(age_at_index) + factor(race) + as.numeric(pack_years_smoked) 
       + factor(ajcc_pathologic_stage) +factor(gender) + factor(vital_status) + as.numeric(days_to_death) + factor(smoking_ever) +factor(dead_or_not_five_years)+factor(DFS_STATUS)
       # set the variable you want to stratify by 
       | age_grp_new, 
       #set data 
       data=Table.1.data.subjects, 
       #set stat display options for continous variables (options from stat.default)
       render.continuous = c((.="Median [Q1, Q3]")))

#now add stats (same as above except where notes added )
table1(~ factor(synchronous_malignancy) + factor(ajcc_pathologic_stage) + as.numeric(age_at_index) + factor(race) + as.numeric(pack_years_smoked) 
       + factor(ajcc_pathologic_stage) +factor(gender) + factor(vital_status) + as.numeric(days_to_death) + factor(smoking_ever)+factor(dead_or_not_five_years)+factor(DFS_STATUS)
       | age_grp_new, 
       data=Table.1.data.subjects, 
       # don't display overall column 
       overall = F,
       render.continuous = c((.="Median [Q1, Q3]")),
       #add p value as extra column (defined in function above)
       extra.col = list("P-value" = pvalue))



################ Kaplan Meier for Recurrence early stage, by median age ########


#use patients data 
Table.1.data.subjects

#We need data from our metadata that include time to recurrence, status of recurrence, mortality
metadata <- Table.1.data.subjects

#search for columns to keep 
#DFS status give recurrence or progression status 
metadata$DFS_STATUS
#DFS_MONTHS give time to recurrenc ein months, should convert this to days 
metadata$DFS_MONTHS
metadata$DFS_days <-  as.numeric(metadata$DFS_MONTHS) * 30.4167

###plot for recurrence 
#Column to be used as time (in days) - DFS_days 
# column to be used as status is DFS_STATUS
survival.data <- metadata[, c("DFS_days","DFS_STATUS", "age_grp_new" )]

#rename the new df 
colnames(survival.data) <- c("time", "status", "Age")

survival.data$status <- ifelse(survival.data$status=="Recurrence", 1, 0)

# turn the data into data frame so it can be used in survival function. 
#turn the class of columns so its either factor or numeric to be used in survival functions 
survival.data<- data.frame(survival.data)
survival.data$time <- as.numeric(survival.data$time)
survival.data$status <- as.numeric(survival.data$status)
survival.data$Age <- as.factor(survival.data$Age)

#remove NA 
survival.data <- na.omit(survival.data)

#get median survival in each group 
survival.data %>% 
  group_by(Age) %>% 
  summarise(median=median(time), 
            Q1=quantile(time, 0.25), 
            Q3=quantile(time, 0.75))

#create survival object

surv <- Surv(survival.data$time, survival.data$status)

#create survival df
sfit <- survfit(Surv(time, status)~Age, data=survival.data)
sfit

survival.plot <- ggsurvplot(sfit, data = survival.data,
                            fun = "pct", ggtheme = theme_pubr(), 
                            conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
                            legend.labs=c("Old", "Young"), legend.title="Age",  
                            palette=c("blue3", "green4"), size=1,
                            title="Kaplan-Meier Curve for Recurrence According to Age",
                            ylab="Overall Recurrence Probability",
                            xlab="Time(Days)",
                            # tables.col = "strata", #color of numbers inside risk table
                            tables.height = 0.15,
                            fontsize = 4,
                            risk.table.y.text.col = T,
                            cumcensor = TRUE, 
                            tables.theme = theme_cleantable(),
                            # risk.table.y.text = FALSE,
                            risk.table.title = "Number at risk (cumulative number of recurrence)",
                            cumcensor.title = "Cumulative number of censored subjects"
)


#save the plot (export as PNG)





#################### Tumor Subset Analysis  #run old vs young ##############

dds.analysis <- dds_tumor
vsd.analysis <- vsd_tumor

#define variables 
dds.analysis$age_grp_new
v= "age_grp_new"
tissue <- "tumor"
stage <- "stage_I_TCGA"

#convert to factor 
dds.analysis[[v]] <- factor(dds.analysis[[v]])
vsd.analysis[[v]] <- factor(vsd.analysis[[v]])

#set reference levels 
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "less_equal_70")
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "less_equal_70")

#define metadata 
metadata_analysis <- data.frame(colData(dds.analysis))

#define colors 
# old
col1 <- "blue3"

#young
col2 <- "green4"
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
colnames(newResults)[colnames(newResults)=="Row.names"] <- "names"

newResults$v <- newResults$age_grp_new

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray.", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0(".txt")))))))
            , sep = "\t", row.names = T)

pdf(file = paste0("Figures/RNA/Beta.Diversity.Bray_PCoA", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0(".pdf"))))))),
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
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Young", "Old"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_analysis, v)

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
write.csv(res, file = paste0("Results/RNA/edgeR.results_",paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0(".csv"))))))))

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR.results_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0(alpha, paste0("_.pdf")))))))), 
    height = 8 , width = 6)
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


####################### Repeat for Lung old vs young, early stage, by median age ####


dds.analysis <- dds_lung
vsd.analysis <- vsd_lung

#define variables 
dds.analysis$age_grp_new
v= "age_grp_new"
tissue <- "lung"
stage <- "stage_I_TCGA"

#convert to factor 
dds.analysis[[v]] <- factor(dds.analysis[[v]])
vsd.analysis[[v]] <- factor(vsd.analysis[[v]])

#set reference levels 
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "less_equal_70")
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "less_equal_70")

#define metadata 
metadata_analysis <- data.frame(colData(dds.analysis))

#define colors 
# old
col1 <- "blue3"

#young
col2 <- "green4"
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
colnames(newResults)[colnames(newResults)=="Row.names"] <- "names"

newResults$v <- newResults$age_grp_new

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray.", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0(".txt")))))))
            , sep = "\t", row.names = T)

pdf(file = paste0("Figures/RNA/Beta.Diversity.Bray_PCoA_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0(".pdf"))))))),
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
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Young", "Old"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_analysis, v)

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
write.csv(res, file = paste0("Results/RNA/edgeR.results_",paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0(".csv"))))))))

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR.results_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0(alpha, paste0("_.pdf")))))))), 
    height = 8 , width = 6)
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


###################plot IPA heatmap, old vs young 

##### this is only an example to generate heatmaps, you should use the results from supplementary tables that match TCGA IPA results according to comparison 

IPA_res <- read.csv(file = "IPA/TCGA/IPA_comparison_TCGA_stage_I_age_70_tumor_lung_old_vs_young.csv")


#get max and min values 
max(IPA_res$TCGA_Tumor, na.rm = TRUE)
min(IPA_res$TCGA_Tumor, na.rm = TRUE)
max(IPA_res$TCGA_Lung, na.rm = TRUE)
min(IPA_res$TCGA_Lung, na.rm = TRUE)

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
col_fun = colorRamp2(c(-1.732, 0, 3.9), c("blue", "white", "orange"))
col_fun(seq(-1.732, 0, 3.9))

IPA_res <- IPA_res %>% arrange(desc(TCGA_Tumor))

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
pdf(file = "Figures/RNA/IPA_comparison_TCGA_stage_I_age_70_tumor_lung_old_vs_young.pdf", height = 20, width = 16)
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
dev.off()



################ Early stage Tumor tissue, comparing recurrence in age groups (according to median)#####

### subset for old and young tables 
dds_old <- dds[,dds$age_grp_new=="greater_70"]
dds_young <- dds[,dds$age_grp_new=="less_equal_70"]

vsd_old <- vsd[,vsd$age_grp_new=="greater_70"]
vsd_young <- vsd[,vsd$age_grp_new=="less_equal_70"]


################# Old subset#### 
#old tumor comparing recu to no rec 
dds.analysis <- dds_old[,dds_old$definition=="Primary solid Tumor"]
vsd.analysis <- vsd_old[,vsd_old$definition=="Primary solid Tumor"]

#change design 
design(dds.analysis) <- ~ DFS_STATUS

#define variables 
dds.analysis$DFS_STATUS
v= "DFS_STATUS"
tissue <- "tumor"
stage <- "stage_I_TCGA"
age_group <- "old"

#convert to factor 
dds.analysis[[v]] <- factor(dds.analysis[[v]])
vsd.analysis[[v]] <- factor(vsd.analysis[[v]])

#set reference levels 
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "No.Recurrence")
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "No.Recurrence")

#define metadata 
metadata_analysis <- data.frame(colData(dds.analysis))

#define colors 
# rec
col1 <- "red2"

#no rec
col2 <- "orange2"
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
colnames(newResults)[colnames(newResults)=="Row.names"] <- "names"

newResults$v <- newResults$DFS_STATUS

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.txt")))))))))
            , sep = "\t", row.names = T)

pdf(file = paste0("Figures/RNA/Beta.Diversity.Bray_PCoA_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.pdf"))))))))),
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
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_analysis, v)

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
write.csv(res, file = paste0("Results/RNA/edgeR.results_",paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.csv"))))))))))

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR.results_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group,  paste0("_", paste0(alpha, paste0("_age_cutoff_70.pdf"))))))))))), 
    height = 8 , width = 6)
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





############ Young subset#### 

#young tumor comparing recu to no rec 
dds.analysis <- dds_young[,dds_young$definition=="Primary solid Tumor"]
vsd.analysis <- vsd_young[,vsd_young$definition=="Primary solid Tumor"]



#change design 
design(dds.analysis) <- ~ DFS_STATUS

#define variables 
dds.analysis$DFS_STATUS
v= "DFS_STATUS"
tissue <- "tumor"
stage <- "stage_I_TCGA"
age_group <- "young"

#convert to factor 
dds.analysis[[v]] <- factor(dds.analysis[[v]])
vsd.analysis[[v]] <- factor(vsd.analysis[[v]])

#set reference levels 
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "No.Recurrence")
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "No.Recurrence")

#define metadata 
metadata_analysis <- data.frame(colData(dds.analysis))

#define colors 
# rec
col1 <- "blue2"

#no rec
col2 <- "skyblue2"
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
colnames(newResults)[colnames(newResults)=="Row.names"] <- "names"

newResults$v <- newResults$DFS_STATUS

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.txt")))))))))
            , sep = "\t", row.names = T)

pdf(file = paste0("Figures/RNA/Beta.Diversity.Bray_PCoA_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.pdf"))))))))),
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
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_analysis, v)

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
write.csv(res, file = paste0("Results/RNA/edgeR.results_",paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.csv"))))))))))

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR.results_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group,  paste0("_", paste0(alpha, paste0("_age_cutoff_70.pdf"))))))))))), 
    height = 8 , width = 6)
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



################ stage I Lung tissue, comparing recurrence in age groups (according to median)#####


############# Old subset#### 
#old tumor comparing recu to no rec 
dds.analysis <- dds_old[,dds_old$definition=="Solid Tissue Normal"]
vsd.analysis <- vsd_old[,vsd_old$definition=="Solid Tissue Normal"]

#change design 
design(dds.analysis) <- ~ DFS_STATUS

#define variables 
dds.analysis$DFS_STATUS
v= "DFS_STATUS"
tissue <- "Lung"
stage <- "stage_I_TCGA"
age_group <- "old"

#convert to factor 
dds.analysis[[v]] <- factor(dds.analysis[[v]])
vsd.analysis[[v]] <- factor(vsd.analysis[[v]])

#set reference levels 
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "No.Recurrence")
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "No.Recurrence")

#define metadata 
metadata_analysis <- data.frame(colData(dds.analysis))

#define colors 
# rec
col1 <- "darkred"

#no rec
col2 <- "orangered"
######################Plot PCoA


#Create Distance Matrix
vegdist = vegdist(t(assay(vsd.analysis)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =8)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = colData(vsd.analysis), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "names"

newResults$v <- newResults$DFS_STATUS

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.txt")))))))))
            , sep = "\t", row.names = T)

pdf(file = paste0("Figures/RNA/Beta.Diversity.Bray_PCoA_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.pdf"))))))))),
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
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_analysis, v)

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
write.csv(res, file = paste0("Results/RNA/edgeR.results_",paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.csv"))))))))))

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.2
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR.results_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group,  paste0("_", paste0(alpha, paste0("_age_cutoff_70.pdf"))))))))))), 
    height = 8 , width = 6)
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





################ Young subset#### 

#young tumor comparing recu to no rec 
dds.analysis <- dds_young[,dds_young$definition=="Solid Tissue Normal"]
vsd.analysis <- vsd_young[,vsd_young$definition=="Solid Tissue Normal"]



#change design 
design(dds.analysis) <- ~ DFS_STATUS

#define variables 
dds.analysis$DFS_STATUS
v= "DFS_STATUS"
tissue <- "Lung"
stage <- "stage_I_TCGA"
age_group <- "young"

#convert to factor 
dds.analysis[[v]] <- factor(dds.analysis[[v]])
vsd.analysis[[v]] <- factor(vsd.analysis[[v]])

#set reference levels 
dds.analysis[[v]] <- relevel(dds.analysis[[v]], ref = "No.Recurrence")
vsd.analysis[[v]] <- relevel(vsd.analysis[[v]], ref = "No.Recurrence")

#define metadata 
metadata_analysis <- data.frame(colData(dds.analysis))

#define colors 
# rec
col1 <- "darkblue"

#no rec
col2 <- "skyblue"
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
colnames(newResults)[colnames(newResults)=="Row.names"] <- "names"

newResults$v <- newResults$DFS_STATUS

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ v,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="v",suffixes=c("",".centroid"))


#stats 
x <- adonis2(vegdist ~ v, data = newResults)
write.table(x, file = paste0("Results/RNA/Beta.Diversity.Bray_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.txt")))))))))
            , sep = "\t", row.names = T)

pdf(file = paste0("Figures/RNA/Beta.Diversity.Bray_PCoA_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.pdf"))))))))),
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
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("No Recurrence", "Recurrence"), size=10)) +
  labs(title = "Beta Diversity, Bray", subtitle = paste0("Adonis = ", paste0(x$`Pr(>F)`)))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

##############Differential Analysis
mycounts.analysis <- assay(dds.analysis)

#get genes table and add 1
x= as(mycounts.analysis, "matrix")
x= x+1 

# get your group variable (the condition you want to analyze according to it)
group = get_variable(metadata_analysis, v)

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
write.csv(res, file = paste0("Results/RNA/edgeR.results_",paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group, paste0("_age_cutoff_70.csv"))))))))))

###### volcano plot
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
alpha = 0.25
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 0 & res$FDR < alpha ] <- col1
cols[res$logFC < 0 & res$FDR < alpha ] <- col2


#plot with labels 
pdf(file = paste0("Figures/RNA/edgeR.results_", paste0(tissue, paste0("_", paste0(stage, paste0("_", paste0(v, paste0("_", paste0(age_group,  paste0("_", paste0(alpha, paste0("_age_cutoff_70.pdf"))))))))))), 
    height = 8 , width = 6)
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



######## plot IPA heatmap for TCGA stage I , rec vs no rec, age cutoff of 70 



#plot IPA heatmap 
#### replace wtih your IPA results file 
#subset for tumor only 

IPA_res <- read.csv(file = "IPA__results_file")

#get max and min values 
max(IPA_res$TCGA_Tumor_Old, na.rm = TRUE)
min(IPA_res$TCGA_Tumor_Old, na.rm = TRUE)
max(IPA_res$TCGA_Tumor_Young, na.rm = TRUE)
min(IPA_res$TCGA_Tumor_Young, na.rm = TRUE)

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)
col_fun = colorRamp2(c(-2.668, 0, 9.065), c("blue", "white", "orange"))
col_fun(seq(-2.668, 0, 9.065))

#convert data to matrix 
IPA_res_mat <- IPA_res

IPA_res_mat <- IPA_res_mat %>% arrange(desc(TCGA_Tumor_Old))

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


#####lung subset 

IPA_res <- read.csv(file = "IPA__results_file")


#get max and min values 
max(IPA_res$TCGA_Lung_Old, na.rm = TRUE)
min(IPA_res$TCGA_Lung_Old, na.rm = TRUE)
max(IPA_res$TCGA_Lung_Young, na.rm = TRUE)
min(IPA_res$TCGA_Lung_Young, na.rm = TRUE)

#####plot heatmap using complex heatmaps 

# Define colors for each levels of qualitative variables
library(circlize)

# Create a custom color scale using only orange
col_fun <- colorRamp2(c(0, 2.887), c("white", "orange"))
col_fun(seq(0, 2.887))

#convert data to matrix 
IPA_res_mat <- IPA_res

IPA_res_mat <- IPA_res_mat %>% arrange(desc(TCGA_Lung_Old))

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








