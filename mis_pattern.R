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


combine_rows <- function(data, row1, row2) {
  data[row2, ] <- data[row1, ] + data[row2, ]
  data[-row1, ]
}
#combine Psychrolutidae and cottidae family groups and Alepisauridae and Paralepididae familes as they have correlated together
#this is to get meaningful values
df_merged_com1 <- combine_rows(df_merged, 294, 102)
df_merged_com2 <- combine_rows(df_merged_com1, 10, 257)
tail(df_merged_com2)

#Calculation
#Get the species count percentage per family by dividing the FishTreeofLife species count by FishBase species count
df_merged_1<- df_merged_com2 %>% mutate(present_species_count_percentage = (FishTreeofLife_species_count/ FishBase_species_count)*100)
#df_merged_1
#missing data
df_merged_2<- df_merged_1 %>% mutate(missing_percentage=100-(FishTreeofLife_species_count/ FishBase_species_count)*100)
unio <- union(df_fish$family,FB_family_species$family)
#attach dataframe "df_merged_1"
setdif_FishBase <- data.frame(setdiff(unio,df_fish$family))
#563-**= 185 (fish families in fishTree but not in fish base)
colnames(setdif_FishBase) <- family


#the percentages are zeros here
new <- merge(df_merged_2,setdif_FishBase,all=TRUE)
names(new)
new <- new[,c("family","missing_percentage")]
new[is.na(new)] <- 100
tail(new)
attach(new)
hist(missing_percentage)
sd(missing_percentage)
mean(missing_percentage)
var(missing_percentage)
length(missing_percentage)
expo2 <- rexp(n=547, rate=(1/74.49158))
hist(expo, main = "Histogram of exponential distribution")
expo2

#Merge by species name
df_merged_species <- merge(df_fish_family_species, FB_family_species, by="species_name")
#rename columns
colnames(df_merged_species) <- c("species","FishTree_fam_name","FishBase_fam_name")

species_diff_family <- df_merged_species %>% 
  filter(FishTree_fam_name!=FishBase_fam_name)
#Summary statistics

attach(df_merged_1) 
hist(present_species_count_percentage)

max(present_species_count_percentage)
min(present_species_count_percentage)
sd(present_species_count_percentage)
mean(present_species_count_percentage)
var(present_species_count_percentage)
#for missing data
attach(df_merged_2)
hist(missing_percentage)
sd(missing_percentage)
mean(missing_percentage)
var(missing_percentage)

#Filter the families with a present species count percentage more than 100; fishTreeofLife dataset has more species than fishbase species
df_fish_filtered <- df_merged_1 %>% 
  filter(present_species_count_percentage>100)
#According to this,  Alepisauridae and Psychrolutidae families have more speciesin fish Treeof Life dataset than fishbase dataset

df_fish_Alepisauridae <- df_fish %>% 
  filter(family=="Alepisauridae")

df_fish_Alepisauridae1 <- FB_family_species %>% 
  filter(family=="Alepisauridae")

df_fish_Psychrolutidae <- df_fish %>% 
  filter(family=="Psychrolutidae")

df_fish_Psychrolutidae1 <- FB_family_species %>% 
  filter(family=="Psychrolutidae")

#install.packages("VennDiagram")
library(VennDiagram)
library(grid)
library(futile.logger)
# Chart
venn.diagram(
  x = list(df_fish$family,FB_family_species$family),
  fill = c("green", "red"),
  category.names = c("FishTree", "FishBase"),
  filename = 'venn_diagram',
  output=TRUE
)
interse <- intersect(df_fish$family,FB_family_species$family)
#364
length(unique(interse))
length(interse)
unio <- union(df_fish$family,FB_family_species$family)
#563
data.frame(unio)
setdif_FishTree <- setdiff(unio,FB_family_species$family)
FishTree14 <- data.frame(setdif_FishTree)
colnames(FishTree14) <- family
ab <- merge(FishTree14,df_fish,by="family")
#now I'm checking whether the species in these 14 families available in Fishbase dataset
#for that I'm removing "_" in species_name column indf_fish dataframe
head(df_fish$species_name)
df_a <- as.data.frame(gsub('_',' ',df_fish$species_name))
new_df_fish <- cbind(df_fish,df_a)
names(new_df_fish)
colnames(new_df_fish) <- c("speciesname","family","seq","species_name")
df_merg_fishtree <- merge(new_df_fish,FishTree14,by="family")
#let's see wheter thwew are common species in fb_family_species dataframe and df_merg_fishtree dataframe
df_merg_ft_fb <- merge(FB_family_species,df_merg_fishtree,by="species_name")
#we found misclassification in family names in fishtree of life data set and Fishbase data set
#so we remove these 14 families
fam_in_FB_not_FT <- setdiff(unio,setdif_FishTree)
fam_in_FB_not_FT <- data.frame(fam_in_FB_not_FT)
colnames(fam_in_FB_not_FT) <- family

df_fish_14fam_reduced <- data.frame(setdiff(df_tt$family,FishTree14$family))
colnames(df_fish_14fam_reduced) <- family
df_fish_14_fam_red <- merge(df_tt,df_fish_14fam_reduced,by="family" )
names(df_fish_14_fam_red)
final_merge <- union_all(df_fish_14_fam_red,df_tt_FB,by="family")
#563-549= 14 (fish families in fishTree but not in fish base)


#check the species in 14 families of FishTree
setdif_FishTree
df_fish_Cottidae1 <- FB_family_species %>% 
  filter(family=="Cottidae")

df_fish_Cottidae <- df_fish %>% 
  filter(family=="Cottidae")

df_fish_Anotopteridae1 <- FB_family_species %>% 
  filter(family=="Anotopteridae")

df_fish_Anotopteridae <- df_fish %>% 
  filter(family=="Anotopteridae")

#364 families included in both the datasets

val <- rnorm(n=364, mean=39.15796, sd=31.12749)
hist(val)

expo <- rexp(n=364, rate=(1/39.15796))
hist(expo)
