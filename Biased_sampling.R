
#Function to combine two rows in a dataframe
combine_rows <- function(data, row1, row2) {
  data[row2, ] <- data[row1, ] + data[row2, ]
  data[-row1, ]
}

#Read the fish dataset in R
fish <- read.table("final_alignment.phylip")

# install.packages("tidyverse")
library(tidyverse)
library(dplyr)
#install.packages("rfishbase")                   
library(rfishbase)


#Rename columns
colnames(fish) <- c("species_name","seq")

#Extract CO1 nucleotides (This step is to extract sequences with CO1 gene. So that we can remove sequences without CO1 gene)
fishCO1 <- substr(fish$seq,2292,2973)

#Check the class of fishCO1 dataset
class(fishCO1)

#Convert to dataframe
fishCO1 <- data.frame(fishCO1)

#Check the converted dataset
class(fishCO1)

#Bind the column "species name" to the fishCO1 dataframe
fishCO1_with_speciesname <- cbind(fish$species_name,fishCO1)

#Replace "-" with nothing 
df_0 <- gsub('-','',fishCO1_with_speciesname$fishCO1)

#Convert to dataframe
df_0 <- as.data.frame(df_0)

#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(fish$species_name,df_0)

#Remove rows with empty CO1 regions (For the phylogenetic placement we do need CO1 region. Therefore the empty CO1 regions are useless for further analysis)
df_2 <- subset(df_1, df_1$df_0 != "")

#Rename column names
colnames(df_2) <- c("species_name","seq")

#Keep only rows with sequences
fish_multigene <- merge(fish,df_2, by= "species_name")

#New dataframe with first two columns (Get the dataframe which we need for the down stream analysis)
fish_multigene <- fish_multigene[, 1:2]

#Rename columns
colnames(fish_multigene) <- c("species_name","seq")

#Read the fish taxonomy spread sheet in R
fish_taxonomy <- read.csv("/home/thanu/Desktop/FishData/Stratified/PFCtaxonomy .csv")

#Extract family column and genus_species column for downstream analysis
df_fish_family_species <- fish_taxonomy[,c("family","genus.species")]

#rename columns
colnames(df_fish_family_species) <- c("family","species_name")

#Replace the space with "_" and covert the resultant list into a df
df_fish_family_species1 <- as.data.frame(gsub(' ','_',df_fish_family_species$species_name))

#Bind the colum "family" of fish_family_species dataset with new dataframe df_fish_family_species1
newdf_fish_taxon <- cbind(df_fish_family_species1,df_fish_family_species$family)

#rename columns
colnames(newdf_fish_taxon) <- c("species_name","family")

#merge fish taxon and sequences
#df_intersect <- intersect(newdf_fish_taxon,fish_multigene)
df_fish <- merge(newdf_fish_taxon,fish_multigene,by="species_name")

#Checking the families with the number of species in each family and convert to a dataframe
df_tt <- data.frame(table(df_fish$family))

#Rename columns
colnames(df_tt) <- c("family","species_count")

#Check the umber of families in fish sample
length(unique(df_fish$family))

# Fishbase dataset matched against the information obtained from FishTreeofLife dataset 
# Extract all of the species names that are available on FishBase
AllFish <- fishbase

#Paste Genus and Species columns together to get the species name 
AllFish$FullName <- paste(AllFish$Genus, AllFish$Species) 
FishBaseSpecies <- AllFish$FullName # 33104 species names
#Match the species labels from FishTreeofLife dataset with the species names from FishBase

#Make FishBaseSpecies into a dataframe first so it can be merged with the family column of AllFish dataset

dfFishBaseSpecies <- data.frame(FishBaseSpecies)

#Checking the columns of fishBase dataset
names(fishbase)

#combine the species name and family columns together in fishbase data
FB_family_species <- cbind(AllFish$Family,dfFishBaseSpecies)

#Rename columns
colnames(FB_family_species) <- c("family","species_name")

#Checking the families with the number of species in each family
tt_FB <- table(FB_family_species$family)
df_tt_FB <- data.frame(tt_FB <- table(FB_family_species$family))

#Rename columns
colnames(df_tt_FB) <- c("family","species_count")

#Check the number of families in fish sample
length(unique(FB_family_species$family))

#Merge the two dataframes of rfish base and fishTreeofLife dataset by "family" column
df_merged <- merge(df_tt, df_tt_FB, by="family")

#Rename columns
colnames(df_merged) <- c("family","FishTreeofLife_species_count","FishBase_species_count")
#There are 364 familes common in fishTree and FishBase
#By looking at this dataset, I found Psychrolutidae has more species count in fishTree than rFishbase. 
#When checking the literature deeply, I was able to identify that Psychrolutidae is a sub set of Cottilae in old work.
#Likewise Alepisaundridae has more species count in fishTree than fishBase. 
#According to literature, there is a close relationship between Alepisauridae and Paralepididae.
#Usually, fishTree dataset should be a subset of Fishbase.
#Hence, this count shows (either the fishbase is incomplete or a descripency in taxonomy between 2 datasets)
#As the mentioned families are correlated, I have combined the two families together for down stream analysis

#this is to get meaningful values (which implies fishTree is a subset of FishBase)
#Combine Psychrolutidae and cottidae family group (row numbers for Psychrolutidae and cottidae are 249 and 102 respectively)
df_merged_com1 <- combine_rows(df_merged, 294, 102)
#Combine Alepisauridae and Paralepididae familes (row numbers for Psychrolutidae and cottidae are 10 and 257 respectively)
df_merged_com2 <- combine_rows(df_merged_com1, 10, 257)
tail(df_merged_com2)

#Calculation
#Get the species count percentage per family by dividing the FishTreeofLife species count by FishBase species count
df_merged_1<- df_merged_com2 %>% mutate(present_species_count_percentage = (FishTreeofLife_species_count/ FishBase_species_count)*100)
#df_merged_1
#missing data
df_merged_2<- df_merged_1 %>% mutate(missing_percentage=100-(FishTreeofLife_species_count/ FishBase_species_count)*100)

#Let's check a venn diagram,
# Chart
library(grid)
library(futile.logger)
library(VennDiagram)
venn.diagram(
  x = list(df_fish$family,FB_family_species$family),
  fill = c("green", "red"),
  category.names = c("FishTree", "FishBase"),
  filename = 'venn_diagram',
  output=TRUE
)
#According to this 364 families are common in both datasets.
#14 familes are there only in fishTree and not in fishbase(Which cannot be acceptable according to our assumption: fishTree is a subset of FishBase).
#I checked these families deeply and found descripency of taxonomy between two datasets
#Hence, we can exclude these 14 for downstream analysis.
#185 familes can be seen in fishbase but not in fishTree (Which can be acceptable)
#So, I'm going to include these 185 families for downstream analysis.
#Now there are 547 families together for the missing data analysis 
#(present species percentage for these 185 families are 0% therefore, the missing species percentage is 100%)
#Number of families in the genetic dataset is 378 (in fish tree of life dataset).

#Union of the two dataset
unio <- union(df_fish$family,FB_family_species$family)
#get the families in fishbase but not in fishTree using the setdiff function.
setdif_FishBase <- data.frame(setdiff(unio,df_fish$family))
#Add the colum the missing species percentage which is 100%
setdif_FishBase$missing_percentage <- 100

#563-378 = 185 (fish families in fishTree but not in fish base)
#assign column names for the setdif_FishBase dataframe
colnames(setdif_FishBase) <- c("family","missing_percentage")
#make a new dataframe by extracting family and missing percentage columns
df_merged_3 <- df_merged_2[,c("family","missing_percentage")]
#Adding 185 families with 100% missing species percentage to 362 families, 
#the percentages are zeros here for present species
new <- rbind(df_merged_3,setdif_FishBase)
#checking the column names of the new dataframe
names(new)

#get the families in fishbase but not in fishTree using the setdiff function.
setdif_FishBase <- data.frame(setdiff(unio,df_fish$family))
#Add the colum the present species percentage which is 0%
setdif_FishBase$present_species_count_percentage <- 0

#563-378 = 185 (fish families in fishTree but not in fish base)
#assign column names for the setdif_FishBase dataframe
colnames(setdif_FishBase) <- c("family","present_species_count_percentage")
#make a new dataframe by extracting family and missing percentage columns
df_merged_4 <- df_merged_2[,c("family","present_species_count_percentage")]
#Adding 185 families with 100% missing species percentage to 362 families, 
#the percentages are zeros here for present species
new_present <- rbind(df_merged_4,setdif_FishBase)
#check the bottom of the dataset
tail(new_present)

attach(new)
#Histogram
hist(missing_percentage,main = "Histogram of present missing species per family", xlab = "Present missing species per family (%)")
#Statistics
sd(missing_percentage)
mean(missing_percentage)
var(missing_percentage)
length(missing_percentage)

detach(new)
attach(new_present)
#Histogram
hist(present_species_count_percentage ,main = "Histogram of present present species per family", xlab = "Present present species per family (%)")
#Statistics
sd(present_species_count_percentage)
mean(present_species_count_percentage)
var(present_species_count_percentage)
length(present_species_count_percentage)

#I'm using exponential to simulate missingness
#random exponential distribution (rate=(1/mean))
# expo <- rexp(n=547, rate=(1/74.49158))
# hist(expo, main = "Histogram of exponential distribution")
# expo
set.seed(1002)
#get 10 sets of 547 values
sample_10 <- data.frame(replicate(10,rexp(n=547, rate=(1/74.49158))))
head(sample_10)
#now using the 10 sets of random numbers do the distribution
for (col in 2:ncol(sample_10)) {
  hist((sample_10[,col]),main = "Exponential distribution for missing species", xlab = "Exponential random numbers")
}

set.seed(1004)
#get 10 sets of 547 values
sample_10_present <- data.frame(replicate(10,rexp(n=547, rate=(1/25.50842))))
head(sample_10_present)

#change values over 100 to 100
sample_10_present[sample_10_present>100] <- 100

#now using the 10 sets of random numbers do the distribution
for (col in 2:ncol(sample_10_present)) {
  hist((sample_10_present[,col]),main = "Exponential distribution for present species", xlab = "Exponential random numbers")
}

#Bind the dataframe with random numbers (sample_10_present) to the dataframe consists of family name and species count
binded_df <- cbind(new_present,sample_10_present)
tail(binded_df)

#get species count to the dataframe (Here,the families with zero species are left out)
merged_df <- merge(df_tt,binded_df,by="family")
detach(new_present)
attach(merged_df)
#Calculate the sample species to extract from the dataset
#merged_df1<- merged_df %>% mutate(sampled_species = (round(species_count*X1/100)))

merged_df1<- merged_df %>% mutate(sampled_species = (ceiling(species_count*X1/100)))

head(merged_df1)
