# Metageno-Guerrero

This repository contains the automatized pipeline used in the article [...] for the metagenomic analysis of the 208 samples collected from Guerrero, Mexico

# Content table
[Conda and the environment](#conda-and-the-environment)

[Read Quality](#read-quality)

[Trimming and Filtering](#trimming-and-filtering)

[Metagenome Assembly](#metagenome-assembly)

[Metagenome Binning](#metagenome-binning)

[Taxonomic assignment](#taxonomic-assignment)

[Rstudio](#rstudio)
  - [Exploring Taxonomy](#exploring-taxonomy)
  - [Diveristy](#diversity)
  - [Taxonomic Analysis](#taxonomic-analysis)

---

## Conda and the environment
Conda is an open-source program that provides package, dependency, and environment management. Environments are a fundamental part of bioinformatics and allow us to make reproducible research

## Read Quality
To check on the quality of the sample we use the tool FASTQC, this program helps us detect any quality issue our data might have, this process is done across all selected samples at the same time, rather than one by one, giving us as an end result a graphic chart where we are able to see the quality of each run and we can compare it with each other

## Trimming and Filtering
In this step we will cut and filter low quality information of our samples that didn´t pass the quality check of FASTQC, in a way we can see it as removing the imperfections of our samples. For this whole process we are using a tool named Trimmomatic. During this process bases will be removed if their Phred score is below 20, as well as, those reads that end with less than 25 paired bases after the cleaning is done. It is worth mentioning that the adapters of the lectures will be eliminated as well.

## Metagenome Assembly

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Metagenome Binning
In the binning of the metagenome, the original genomes of the samples are separated. This allows the individual analysis of the species that have enough reads to reconstruct their genome, these reconstructions are known as MAGS. To check MAGS qualities the program CheckM is used, this program checks that the genome has the complete information and that it is not cross-contaminated with info of other genomes.

## Taxonomic Asignment
For the taxonomical assignment, the genomes of the samples are compared with data bases that contain complete genomes of different species. The database used in this case is Kraken2. This taxonomical classification program delivers high quality and high speed classifications. Afterwards, the results of Kraken are put together in a graphic of the Krona program, which allows us to explore the different abundances of the many microorganisms of the sample being analyzed.

## Rstudio
RStudio is an integrated development environment (IDE), which is widely used in statistical computing, data analysis, and data visualization. It provides a user-friendly interface that simplifies the process of writing and executing R code, making it especially helpful for both beginners and experienced data scientists. RStudio can be run locally as a desktop application or accessed remotely via a web browser through RStudio Server, which is particularly useful in collaborative or high-performance computing environments, which is the way in which we worked due to the high computational requirements needed for this research

  ### Exploring Taxonomy
For the taxonomical assignment, the genomes of the samples are compared with data bases that contain complete genomes of different species. The database used in this case is Kraken2. This taxonomical classification program delivers high quality and high speed classifications. Afterwards, the results of Kraken are put together in a graphic of the Krona program, which allows us to explore the different abundances of the many microorganisms of the sample being analyzed.
 
  ### Diversity
In Rstudio, we can determine the levels of taxonomic diversity of the metagenome, in this case, we are going to use two different measures, these being alpha and beta diversities. Alpha diversity determines the different species in an environment. Beta diversity is the difference among two or more environments. Afterwards both results are going to be represented in a graphic
  
  ### Taxonomic Analysis
Finally we are going to create graphics using Rstudio to be able to visualize the taxonomic diversity of the relative and absolute abundance. This can also be used to explore diversities and abundances of specific families





 
