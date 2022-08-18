#_______________________________________________________________________________
### RNA Analysis of 9 Thalassiosira rotula strains (arctic vs. temperate) ###


#_______________________________________________________________________________
#### RNA Contig presence-absence analysis ####

# Housekeeping
rm(list = ls())
library(ape)
library(readr)
library(UpSetR)
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
#-> make the plot more beautiful with 'ggtree'

#_______________________________________________________________________________
#### SNP Analysis ####

setwd("../SNP files/") #change working directory to import files

#Data Import
file_list <- list.files() #create list of files to import
ldf <- lapply(file_list, read.delim) #check format, if txt use read.delim
colnames(ldf[[1]]) #get the colnames and sort out what you need (1:5,31,35)
str(ldf[[1]]) #check that you don't have stored the characters as factors, if so past options(StringsAsFactor=FALSE) in the beginning of the script.
ldf_sub <- lapply(ldf, function(x) { x[c(6:30,32:34,36,37)] <- NULL; x })
str(ldf_sub[[1]])

#remove SNPs of contigs that are not present in all strains (get from presence/absence analysis above)
head(df)
dfx <- as.data.frame(df)
str(dfx)
dfx$sum <- rowSums(dfx)
dfx <- subset(dfx, dfx$sum==9)
dfx <- rownames(dfx) #store row names
dfx <- gsub(' mapping','',dfx) #remove the "mapping" and blank space from the name to make them match
ldf_sub <- lapply(ldf_sub, function(x) {subset(x,x$Chromosome %in% dfx)}) # subset SNP lists using the present contig names
str(ldf_sub[[1]])

#filter for amino-acid changing SNPs
ldf_sub <- lapply(ldf_sub, function(x) {subset(x,x$Non.synonymous=="Yes")})
str(ldf_sub[[1]])

#create SNP IDs based on the SNP information and contig name
ldf_sub <- lapply(ldf_sub, function(x) {cbind(x,x$id <- paste(x$Chromosome,x$Region,x$Type,x$Reference,x$Allele,x$Non.synonymous,sep="_"))})

#rename column 
ldf_sub <- lapply(ldf_sub, function(x) {colnames(x)[8]<-'id'; x})
str(ldf_sub[[1]])

#remove everything except the IDs (we dont need counts anymore, only presence absence of respective SNPs)
ldf_sub <- lapply(ldf_sub, function(x) { x[1:7] <- NULL; x })

#shorten names in the list
names_ldf <- c('J1','J2','J3','J4','J5','J6','J7','J8','J9')
names(ldf_sub) <- names_ldf 
str(ldf_sub[[1]])

#change the format
listInput <- lapply(ldf_sub,function(x) x[[1]])

#Use UpSetR to give all intersections of SNPs between strains. Every SNP has a 
#unique name containing information (contig, position,which base replaced) and 
#can be identified across the different strains (if present). The package also
#stores the data as binary matrix (can be copied for plots later on)
SNPs <- upset(fromList(listInput), order.by = "freq", nsets = 9)
SNPs <- SNPs$New_data

#since the matrix does not have rownames, those will be attached here
x1 <- unlist(listInput, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
rownames(SNPs) <- x1
head(SNPs)

#kick out false positives (all SNPs that are present in all samples due to mistakes in the reference assembly)
SNPs$Sums <- rowSums(SNPs)
table(SNPs$Sums) #334 present in all strains -> remove
SNP_A <- SNPs[SNPs$Sums<9,]
colnames(SNP_A)
table(SNP_A$Sums)
SNP_A$Sums <- NULL

#create a tree from SNPs
snp_matrix <- t(as.matrix(SNP_A)) #transform to matrix object
mat.snp <- nj(dist.gene(snp_matrix)) # neighbor joining tree construction

# Plotting Data
plot(mat.snp, "phylo") # first glimpse but looks ugly
nodelabels() 
plot(mat.snp,type = "unrooted", cex = 0.4) #remove roots
#-> make the plot more beautiful with 'ggtree'

#_______________________________________________________________________________
#### SNP Analysis of BUSCO genes (only one version of the specific gene present)

#Change working directory
setwd("/Users/jgiesler/Library/Mobile Documents/com~apple~CloudDocs/PhD/R Directory/Transcriptomics/RNA_Thalassiosira")
bs <- read.delim("busco.txt") #import busco genes (already matched and named after our contigs)
head(bs)

#add a column to the SNP_A df with the same names as in the bs file
SNP_A$contig <- substring(rownames(SNP_A), 1, 20)
SNP_busco <- subset(SNP_A,SNP_A$contig %in% unique(bs$contig))
SNP_busco$contig <- NULL

#create a tree from SNPs (only BUSCO genes)
snp_busco_matrix <- t(as.matrix(SNP_busco)) #transform to matrix object
mat.snp.busco <- nj(dist.gene(snp_busco_matrix)) # neighbor joining tree construction

# Plotting Data
plot(mat.snp.busco, "phylo") # first glimpse but looks ugly
nodelabels() 
plot(mat.snp.busco,type = "unrooted", cex = 0.75) #remove roots
#-> make the plot more beautiful with 'ggtree'
#mat.snp.busco$tip.label <- c(jj1,jj2,jj3,jj4,jj5,jj6)#relabel 

