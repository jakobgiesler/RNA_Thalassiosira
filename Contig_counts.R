### RNA Contig presence-absence analysis of 9 T. rotula strains ###

# Housekeeping
library(ape)
library(readr)
setwd("/Users/jgiesler/Library/Mobile Documents/com~apple~CloudDocs/PhD/R Directory/Transcriptomics/RNA_Thalassiosira/Contig counts")

# Data Import
file_list <- list.files() #create list of files to import
ldf <- lapply(file_list, read_csv) #check format, if txt use read.delim

# Data Edit
names <- ldf[[1]][["Name"]] #store all contig names in separate df
ldf_sub <- lapply(ldf, function(x) { x[c(1,2,4:10)] <- NULL; x }) #only keep the counts for each sample
df <- as.data.frame(do.call(cbind, ldf_sub)) #bind samples together (do.call because of list object)
colnames(df) <- c('J1','J2','J3','J4','J5','J6','J7','J8','J9') #set sample names (strains) as colnames
rownames(df) <- names #set contig names as row names
df<- ifelse(df[,c(1:9)]>= 10,1,0) #create presences (1) - absences (0)
df_matrix <- t(as.matrix(df)) #transform to matrix object
mat.df <- nj(dist.gene(df_matrix)) # neighbor joining tree construction

# Plotting Data
plot(mat.df, "phylo") # first glimpse but looks ugly
nodelabels() 
plot(mat.df,type = "unrooted", cex = 0.4) # remove root structure
