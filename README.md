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

Figures
- [PCoA in the anthropometric-demographic space and its association with BMI percentile, age group, and lifestyle](#pcoa)

- [Alpha and beta diversity of Rural and Urban metagenomes.](#alpha-and-beta-diversity-of-rural-and-urban-metagenomes)
- [Differentially abundant bacterial and eukaryotic taxa between Rural and Urban groups](#differentially-abundant-bacterial-and-eukaryotic-taxa-between-rural-and-urban-groups)
- [Mean relative abundance of the most prevalent bacterial and eukaryotic genera in rural and urban samples](#mean-relative-abundance-of-the-most-prevalent-bacterial-and-eukaryotic-genera-in-rural-and-urban-samples)
- Volcano plots
  - [Setup](#setup)
  - [Rural vs Urban](#rural-vs-urban)
  - [BMI<25 vs BMI≥25](#bmi)
  - [BMI≥25 Rural vs BMI≥25 Urban](#bmiru)
- [BUSCO completeness assessment of MAGs](#busco-completeness-assessment-of-mags)
- [Functional annotation profiles of high-quality MAGs](#fun-anno)
   - [Part A](#part-a)
   - [Part B](#part-b)
   - [Part C](#part-c)
- [Functional annotation of genes from selected taxa](#functional-annotation-of-genes-from-selected-taxa)
- [Differential enrichment of COG categories in Blautia and Clostridia between Rural and Urban groups](#diff-enr)
- [Co-occurrence networks of bacterial and eukaryotic taxa based on Spearman correlations](#occurr)
- [Differential distribution of functional annotations reconstructed from high-quality metagenome-assembled genomes (MAGs)](#diff-dis)

---

## Conda and the environment
Conda is an open-source program that provides package, dependency, and environment management. Environments are a fundamental part of bioinformatics and allow us to make reproducible research

## Read Quality
To check on the quality of the sample we use the tool FASTQC, this program helps us detect any quality issue our data might have, this process is done across all selected samples at the same time, rather than one by one, giving us as an end result a graphic chart where we are able to see the quality of each run and we can compare it with each other

## Trimming and Filtering
In this step we will cut and filter low quality information of our samples that didn´t pass the quality check of FASTQC, in a way we can see it as removing the imperfections of our samples. For this whole process we are using a tool named Trimmomatic. During this process bases will be removed if their Phred score is below 20, as well as, those reads that end with less than 25 paired bases after the cleaning is done. It is worth mentioning that the adapters of the lectures will be eliminated as well.

## Metagenome Assembly

During the assembly of the metagenome the final objective is to ideally obtain the whole sequence of a chromosome. There are several methods and programs in which this can be done, for example the greedy extension, overlap layout consensus, De Bruijns graphs, etc.  

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

# Figures

<a id="pcoa"></a>
## PCoA in the anthropometric-demographic space and its association with BMI percentile, age group, and lifestyle

###  REQUIRED LIBRARIES 
```{r}
library(ggplot2)
library(vegan)
library(dplyr)
library(ggcorrplot)
library(patchwork)
```
### --- 1. DATA LOADING AND CLEANING ---
``` {r}
data <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE,
  sep = ",",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
```
### Remove empty columns and rows (only those completely empty)

``` {r}
data <- data[, colSums(!is.na(data)) > 0]
data <- data[rowSums(is.na(data)) < ncol(data), ]
```
### Numeric variables
```
data$BMI <- as.numeric(data$BMI)
data$Age <- as.numeric(data$Age)
```
### --- 2. DERIVED VARIABLES ---
```
data <- data %>%
  mutate(
    Percentil_group = case_when(
      Percentil_formulas %in% c("Underweight", "Malnutrition") ~ "Underweight",
      Percentil_formulas %in% c("Overweight", "Obesity") ~ "Overweight/Obesity",
      Percentil_formulas == "Normal Weight" ~ "Normal Weight",
      TRUE ~ NA_character_
    ),
    Lifestyle_num = ifelse(Lifestyle == "Urban", 1, 0)
  )

```
### --- 3. FILTERING for analysis ---

```
data_complete <- data %>%
  filter(!is.na(BMI), !is.na(Age), !is.na(Percentil_group), !is.na(Age_group), !is.na(Lifestyle))
```
### --- 4. DISTANCE MATRIX AND PCoA ---
```
dist_matrix <- vegdist(data_complete[, c("BMI", "Age")], method = "bray")
pcoa <- cmdscale(dist_matrix, k = 2, eig = TRUE)
scores_pcoa <- as.data.frame(pcoa$points)
colnames(scores_pcoa) <- c("Dim1", "Dim2")
```
### Add categorical variables to color points

```
scores_pcoa$Percentil_group <- data_complete$Percentil_group
scores_pcoa$Age_group <- data_complete$Age_group
scores_pcoa$Lifestyle <- data_complete$Lifestyle
```
### Labels with counts for Percentil_group
```
percentil_counts <- data_complete %>%
  count(Percentil_group)
percentil_labels <- setNames(paste0(percentil_counts$Percentil_group, " (", percentil_counts$n, ")"), percentil_counts$Percentil_group)

```
### --- 5. PERMANOVA ---
```
set.seed(123)
adonis_Percentil <- adonis2(dist_matrix ~ Percentil_group, data = data_complete)
adonis_Age <- adonis2(dist_matrix ~ Age_group, data = data_complete)
adonis_Lifestyle <- adonis2(dist_matrix ~ Lifestyle, data = data_complete)

cat("PERMANOVA RESULTS:\n")
print(adonis_Percentil)
print(adonis_Age)
print(adonis_Lifestyle)
```
### --- 6. UNIFIED THEME --
```
custom_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.title = element_text(face = "plain", size = 10),
    axis.text = element_text(face = "plain", size = 10),
    legend.title = element_text(face = "plain", size = 11),
    legend.text = element_text(size = 10),
    aspect.ratio = 1,
    plot.margin = margin(8,8,8,8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
```
### --- 7.1 PLOTS ---
#### Labels with counts for Age_group
```
age_counts <- data_complete %>%
  count(Age_group)
age_labels <- setNames(paste0(age_counts$Age_group, " (", age_counts$n, ")"), age_counts$Age_group)

```
#### Labels with counts for Lifestyle
```
lifestyle_counts <- data_complete %>%
  count(Lifestyle)
lifestyle_labels <- setNames(paste0(lifestyle_counts$Lifestyle, " (", lifestyle_counts$n, ")"), lifestyle_counts$Lifestyle)

```
### --- 7. 2 PLOTS ---
```
p_percentil <- ggplot(scores_pcoa, aes(Dim1, Dim2, color = Percentil_group)) +
  geom_point(size = 1.8, alpha = 0.6) +
  labs(title = "PCoA: Percentil Group", color = "Percentil Group") +
  scale_color_manual(values = c(
    "Overweight/Obesity" = "#D55E00",
    "Normal Weight" = "#0072B2",
    "Underweight" = "#E69F00"
  ), labels = percentil_labels) +
  custom_theme
```
### Counts across the whole dataset for Age_group and Lifestyle
```
age_counts_full <- data %>%
  filter(!is.na(Age_group)) %>%
  count(Age_group)

age_labels_full <- setNames(paste0(age_counts_full$Age_group, " (", age_counts_full$n, ")"), age_counts_full$Age_group)

lifestyle_counts_full <- data %>%
  filter(!is.na(Lifestyle)) %>%
  count(Lifestyle)

lifestyle_labels_full <- setNames(paste0(lifestyle_counts_full$Lifestyle, " (", lifestyle_counts_full$n, ")"), lifestyle_counts_full$Lifestyle)

```
### Use these complete labels in the plots
```
p_age <- ggplot(scores_pcoa, aes(Dim1, Dim2, color = Age_group)) +
  geom_point(size = 1.8, alpha = 0.6) +
  labs(title = "PCoA: Age Group", color = "Age Group") +
  scale_color_brewer(palette = "Set1", labels = age_labels_full) +
  custom_theme

p_lifestyle <- ggplot(scores_pcoa, aes(Dim1, Dim2, color = Lifestyle)) +
  geom_point(size = 1.8, alpha = 0.6) +
  labs(title = "PCoA: Lifestyle", color = "Lifestyle") +
  scale_color_manual(values = c(
    "Rural" = "#E69F00",
    "Urban" = "#0072B2"
  ), labels = lifestyle_labels_full) +
  custom_theme

```
### Create data_filtered for correlation with Percentil_num and Lifestyle_num
```
data_filtered <- data_complete %>%
  mutate(
    Lifestyle_num = ifelse(Lifestyle == "Urban", 1, 0),
    Percentil_num = case_when(
      Percentil_group == "Underweight" ~ 1,
      Percentil_group == "Normal Weight" ~ 3,
      Percentil_group == "Overweight/Obesity" ~ 5, # aquí se juntaron Overweight y Obesity
      TRUE ~ NA_real_
    )
  )

```
### Correlation with Percentil_num, Age, and Lifestyle_num
```
cor_data <- data_filtered %>% 
  select(Percentil_num, Age, Lifestyle_num) %>% 
  na.omit()

cor_matrix <- cor(cor_data)

p_corr <- ggcorrplot(
  cor_matrix,
  lab = TRUE,
  colors = c("#0072B2", "white", "#D55E00"),
  title = "Correlation Matrix",
  type = "upper"
) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.title = element_blank(),
    axis.text = element_text(size = 10),
    aspect.ratio = 1,
    plot.margin = margin(8,8,8,8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

```
### --- 9. COMBINE PLOTS ---
```
library(patchwork)
```
### Themes to adjust margins by row
```
custom_theme_top <- custom_theme + theme(plot.margin = margin(t = 8, r = 8, b = 0, l = 8))
custom_theme_bottom <- custom_theme + theme(plot.margin = margin(t = 0, r = 8, b = 8, l = 8))

p_percentil <- p_percentil + custom_theme_top
p_age <- p_age + custom_theme_top
p_lifestyle <- p_lifestyle + custom_theme_bottom
p_corr <- p_corr + custom_theme_bottom

final_plot <- (p_percentil + p_age) / (p_lifestyle + p_corr) +
  plot_layout(guides = "keep") +
  plot_annotation(tag_levels = "A")

print(final_plot)
```
### --- Load required libraries ---
```
library(ggplot2)
library(dplyr)
library(ggpubr)
library(grid)
library(gtable)
library(patchwork)
library(tibble)

```
### --- Load data ---
```
data <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE,
  sep = ",",
  fileEncoding = "latin1", # prueba latin1 para acentos
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

```
### Remove empty columns and rows
```
data <- data[, colSums(!is.na(data)) > 0]
data <- data[rowSums(is.na(data)) < ncol(data), ]
```
### Create derived variables
```
data <- data %>%
  mutate(
    Age_group_recode = case_when(
      Age_group %in% c("Infant", "Child", "Teenager") ~ "Child",
      Age_group %in% c("Adult", "Elderly") ~ "Adult",
      TRUE ~ NA_character_
    ),
    Age_group_recode = factor(Age_group_recode, levels = c("Child", "Adult")),
    Gender_Lifestyle = paste(Gender, Lifestyle, sep = " "),
    Gender_Lifestyle = factor(Gender_Lifestyle,
                              levels = c("Male Rural", "Female Rural", "Male Urban", "Female Urban")),
    BMI_group = case_when(
      BMI >= 30 ~ "Overweight/Obesity",
      BMI >= 18.5 & BMI < 30 ~ "Normal Weight",
      BMI < 18.5 ~ "Underweight",
      TRUE ~ NA_character_
    ),
    BMI_group = factor(BMI_group, levels = c("Underweight", "Normal Weight", "Overweight/Obesity")),
    Percentil_group = case_when(
      Percentil_formulas %in% c("Malnutrition", "Underweight") ~ "Underweight",
      Percentil_formulas == "Normal Weight" ~ "Normal Weight",
      Percentil_formulas %in% c("Overweight", "Obesity") ~ "Overweight/Obesity",
      TRUE ~ NA_character_
    ),
    Percentil_group = factor(Percentil_group,
                             levels = c("Underweight", "Normal Weight", "Overweight/Obesity"))
  )

```
### Filter for analysis (only rows without NA in key variables)
```
data_filtered <- data %>%
  filter(!is.na(BMI), !is.na(Age), !is.na(Percentil_formulas), !is.na(Age_group), !is.na(Lifestyle))

```
### --- Define comparisons for the statistical test ---
```
comparisons_list <- list(
  c("Male Rural", "Male Urban"),
  c("Female Rural", "Female Urban"),
  c("Male Rural", "Female Rural"),
  c("Male Urban", "Female Urban")
)

comparisons_child <- list(
  c("Male Rural", "Male Urban"),
  c("Female Rural", "Female Urban"),
  c("Male Rural", "Female Rural"),
  c("Male Urban", "Female Urban")
)

comparisons_adult <- comparisons_child  # mismo set para ambos si quieres

```
### Create labels with total counts from the original dataset for Lifestyle and Percentil_group
```
labels_lifestyle <- data %>%
  filter(!is.na(Lifestyle)) %>%
  count(Lifestyle) %>%
  mutate(label = paste0(Lifestyle, " (n = ", n, ")")) %>%
  select(Lifestyle, label) %>%
  deframe()

labels_percentil_group <- data %>%
  filter(!is.na(Percentil_group)) %>%
  count(Percentil_group) %>%
  mutate(label = paste0(Percentil_group, " (n = ", n, ")")) %>%
  select(Percentil_group, label) %>%
  deframe()

```
### --- Plot 1: BMI density by Lifestyle (filtered data, full legend) ---
```
plot1 <- ggplot(data_filtered, aes(x = BMI, fill = Lifestyle)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(
    values = c("Rural" = "#FF7F00", "Urban" = "#1E90FF"),
    labels = labels_lifestyle
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0, face = "plain", size = 12),
    axis.title = element_text(face = "plain", size = 10),
    axis.text = element_text(size = 9, face = "plain"),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  labs(
    title = "BMI Density by Lifestyle",
    x = "BMI",
    y = "Density",
    fill = "Lifestyle"
  )

```
### --- Plot 2: Age density by Percentil_group ---
```
plot2 <- ggplot(data_filtered, aes(x = Age, fill = Percentil_group)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(
    values = c(
      "Underweight" = "#F4A261",
      "Normal Weight" = "#457B9D",
      "Overweight/Obesity" = "#E63946"
    ),
    labels = labels_percentil_group
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0, face = "plain", size = 12),
    axis.title = element_text(face = "plain", size = 10),
    axis.text = element_text(size = 9, face = "plain"),
    plot.margin = margin(10, 15, 10, 15)
  ) +
  labs(
    title = "Age Distribution by Percentil Group",
    x = "Age",
    y = "Density",
    fill = "Percentil Group"
  )
```
### Calculate n by group and age (already done)
```
group_counts <- data_filtered %>%
  count(Age_group_recode, Gender_Lifestyle)

```
### Position to place the texts (n) below the X axis
```
y_pos_n <- min(data_filtered$BMI, na.rm = TRUE) - 2  # ajusta -2 según tu rango de BMI

```
### Create text label for n
```
group_counts <- group_counts %>%
  mutate(n_label = paste0("n = ", n))

```
### Plot with geom_text for n
```
plot3 <- ggplot(data_filtered, aes(x = Gender_Lifestyle, y = BMI, fill = Gender_Lifestyle)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(position = position_jitter(width = 0.15), size = 0.6, alpha = 0.3) +
  facet_wrap(~Age_group_recode, nrow = 1, strip.position = "top") +
  scale_fill_manual(values = c(
    "Male Rural" = "#D62828",
    "Female Rural" = "coral",
    "Male Urban" = "skyblue",
    "Female Urban" = "#F2A104"
  )) +
```
### Add n below each bar
```
  geom_text(
    data = group_counts,
    aes(x = Gender_Lifestyle, y = y_pos_n, label = n_label),
    inherit.aes = FALSE,
    size = 3,
    angle = 0,
    vjust = -3,
    hjust = 0.5
  ) +
  stat_compare_means(
    data = subset(data_filtered, Age_group_recode == "Child"),
    aes(x = Gender_Lifestyle, y = BMI, fill = Gender_Lifestyle),
    comparisons = comparisons_child,
    method = "wilcox.test",
    label = "p.signif",
    exact = FALSE,
    inherit.aes = FALSE,
    label.y = c(34, 28, 22, 24)
  ) +
  stat_compare_means(
    data = subset(data_filtered, Age_group_recode == "Adult"),
    aes(x = Gender_Lifestyle, y = BMI, fill = Gender_Lifestyle),
    comparisons = comparisons_adult,
    method = "wilcox.test",
    label = "p.signif",
    exact = FALSE,
    inherit.aes = FALSE,
    label.y = c(40, 34, 30, 30)
  ) +
  labs(
    title = "BMI by Gender & Lifestyle in Children and Adults",
    y = "BMI",
    x = NULL,
    fill = "Gender & Lifestyle"
  ) +
  guides(fill = guide_legend(title = "Gender & Lifestyle")) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain", size = 10),
    axis.text.x = element_blank(),  # Ocultar etiquetas X originales
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(face = "plain", size = 9),
    axis.title.y = element_text(face = "plain", size = 9),
    plot.title = element_text(face = "plain", size = 10, hjust = 0),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 15)
  )


library(patchwork)


combined_plot <- (plot1 | plot2)  / plot3 + 
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.02, 0.95)  # X, Y posición de la etiqueta, Y más bajo (de 1 a 0)
    )
  )


print(combined_plot)

```
## FINAL FIG. 1 COMPOSED
```
plot1 <- plot1 + theme(plot.title = element_blank())
plot2 <- plot2 + theme(plot.title = element_blank())
plot33 <- plot3 + theme(plot.title = element_blank())
p_lifestyle <- p_lifestyle + theme(plot.title = element_blank())


combined_plot <- (plot1 | plot2)  / (plot3 |p_lifestyle) + 
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.02, 0.95)  # X, Y posición de la etiqueta, Y más bajo (de 1 a 0)
    )
  )

print(combined_plot)
```
## Alpha and beta diversity of Rural and Urban metagenomes.

### --- Libraries ---
```{r}
library(dplyr)
library(tidyr)
 library(ggplot2)

```
### --- 0. Standardize column and file names ---
```{r}
cat("\n--- Renombrando columnas de Kraken2 ---\n")
kraken2_aligned <- kraken2_data_full %>%
   rename(
     Domain      = domain,
   Supergroup  = supergroup,
     Kingdom     = subgroup,
     Phylum      = phylum,
     Class       = class,
     Order       = order,
     Family      = family,
     Genus       = genus,
     Species     = species
   )
```
### Assign key numeric columns
```{r}
kraken2_aligned$reads <- as.numeric(kraken2_data_full$reads_clade)
 kraken2_aligned$percent <- as.numeric(kraken2_data_full$percent)

 kaiju_data$percent <- as.numeric(kaiju_data$percent)
 kaiju_data$reads   <- as.numeric(kaiju_data$reads)
 
```
### Set common file names
```{r}
 cat("\n--- Estandarizando nombres de archivo ---\n")
 kaiju_data$file <- gsub("_kaiju\\.out$", "", kaiju_data$file)
 kraken2_aligned$file <- gsub("\\.report$", "", kraken2_aligned$file)
 
```
### Confirm number of unique files
```{r}
 cat("\nArchivos únicos antes de unir:\n")
 cat("  Kaiju:   ", length(unique(kaiju_data$file)), "\n")
 cat("  Kraken2: ", length(unique(kraken2_aligned$file)), "\n")
 
```
### --- 1. Define taxonomic hierarchy ---
```{r}
 tax_levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Supergroup", "Domain")
 
```
### --- 2. Function to detect the most specific taxonomic level ---
```{r}
get_deepest_tax <- function(row) {
  for (level in tax_levels) {
     value <- row[[level]]
     if (!is.na(value) && value != "") {
       return(paste0(level, ":", value))
     }
   }
   return(NA_character_)
 }
 
```
### --- 3. Apply function to both datasets ---
```{r}
 cat("\nAplicando nivel taxonómico más específico...\n")
 kaiju_data$deepest_tax <- apply(kaiju_data[, tax_levels], 1, get_deepest_tax)
 kraken2_aligned$deepest_tax <- apply(kraken2_aligned[, tax_levels], 1, get_deepest_tax)
 
```
### Check example
```{r}
 cat("  Ejemplo Kaiju:     ", head(kaiju_data$deepest_tax, 1), "\n")
 cat("  Ejemplo Kraken2:   ", head(kraken2_aligned$deepest_tax, 1), "\n")
 
```
### --- 4. Create unique keys: file + deepest_tax ---
```{r}
 kaiju_keys <- paste(kaiju_data$file, kaiju_data$deepest_tax)
 kraken2_keys <- paste(kraken2_aligned$file, kraken2_aligned$deepest_tax)
 
 cat("\nLongitudes de claves únicas:\n")
 cat("  Kaiju:     ", length(unique(kaiju_keys)), "\n")
 cat("  Kraken2:   ", length(unique(kraken2_keys)), "\n")

```
### --- 5. Filter Kraken2: only what is NOT in Kaiju ---
```{r}
 kraken2_filtered <- kraken2_aligned[!(kraken2_keys %in% kaiju_keys), ]
 
 cat("\nLecturas después de filtrar Kraken2:\n")
 cat("  Kraken2 original: ", nrow(kraken2_aligned), "\n")
 cat("  Kraken2 filtrado: ", nrow(kraken2_filtered), "\n")
 
```
### --- 6. Merge data without logical duplicates ---
```{r}
 cat("\n--- Uniendo datasets (Kaiju + Kraken2 filtrado) ---\n")
 kaiju_merged <- bind_rows(kaiju_data, kraken2_filtered)

```
### --- 7. Final validations ---
```{r}
 cat("\nValidaciones finales:\n")
cat("  Kaiju        - Total reads: ", sum(kaiju_data$reads, na.rm = TRUE), "\n")
cat("  Kraken2 FILT - Total reads: ", sum(kraken2_filtered$reads, na.rm = TRUE), "\n")
 cat("  Fusionado    - Total reads: ", sum(kaiju_merged$reads, na.rm = TRUE), "\n")
 
cat("  Kaiju        - Nrows: ", nrow(kaiju_data), "\n")
 cat("  Kraken2 FILT - Nrows: ", nrow(kraken2_filtered), "\n")
 cat("  Fusionado    - Nrows: ", nrow(kaiju_merged), "\n")
 
cat("  Archivos fusionados: ", length(unique(kaiju_merged$file)), "\n\n")
  
```
### --- Rename Kraken2 columns ---
```{r}
 kraken2_aligned <- kraken2_data_full %>%
   rename(
     Domain      = domain,
     Supergroup  = supergroup,
     Kingdom     = subgroup,
     Phylum      = phylum,
     Class       = class,
     Order       = order,
     Family      = family,
     Genus       = genus,
    Species     = species
   )
 
```
### --- Use file_base from Kraken2 as reference ---
```{r}
 kraken2_aligned$file <- kraken2_data_full$file_base
 
```
### --- Extract file_base from Kaiju ---
```{r}
 kaiju_data$file <- gsub("_kaiju\\.out$", "", basename(kaiju_data$file))
kaiju_data$file <- gsub("\\.out$", "", kaiju_data$file)  # por si hay archivos erróneos
 
```
### Validate if codes match
```{r}
 cat("\n--- Comparación de archivos únicos ---\n")
 cat("Kaiju archivos únicos:   ", length(unique(kaiju_data$file)), "\n")
 cat("Kraken2 archivos únicos: ", length(unique(kraken2_aligned$file)), "\n")
 cat("Coincidencias exactas:   ", length(intersect(unique(kaiju_data$file), unique(kraken2_aligned$file))), "\n")
 cat("Diferencias:             ", length(setdiff(unique(kraken2_aligned$file), unique(kaiju_data$file))), "\n")
 
```
### --- Establish taxonomic levels ---
```{r}
 tax_levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Supergroup", "Domain")
 
```
### --- Detect most specific level ---
```{r}
 get_deepest_tax <- function(row) {
   for (level in tax_levels) {
     value <- row[[level]]
     if (!is.na(value) && value != "") return(paste0(level, ":", value))
   }
   return(NA_character_)
 }
 
 kaiju_data$deepest_tax <- apply(kaiju_data[, tax_levels], 1, get_deepest_tax)
 kraken2_aligned$deepest_tax <- apply(kraken2_aligned[, tax_levels], 1, get_deepest_tax)
 
```
### --- Unique keys ---
```{r}
 kaiju_keys   <- paste(kaiju_data$file, kaiju_data$deepest_tax)
 kraken2_keys <- paste(kraken2_aligned$file, kraken2_aligned$deepest_tax)
 
```
### --- Conversion of numeric values ---
```{r}
kaiju_data$percent   <- as.numeric(kaiju_data$percent)
 kaiju_data$reads     <- as.numeric(kaiju_data$reads)
 kraken2_aligned$reads   <- as.numeric(kraken2_aligned$reads_clade)
kraken2_aligned$percent <- as.numeric(kraken2_aligned$percent)
 
```
### --- Filter duplicates ---
```{r}
kraken2_filtered <- kraken2_aligned[!(kraken2_keys %in% kaiju_keys), ]
 
```
### --- Merge ---
```{r}
 kaiju_merged <- bind_rows(kaiju_data, kraken2_filtered)
 
```
### --- Assign file_base and Group to merged object ---
```{r}
 kaiju_merged <- kaiju_merged %>%
  mutate(
     file_base = gsub("_kaiju\\.out$", "", basename(file)),  # extraer nombre base
     Group = case_when(
       file_base %in% rural_files_kaiju ~ "Rural",
       file_base %in% urban_files_kaiju ~ "Urban",
       TRUE ~ NA_character_
    )
   )
 
```
### Ensure reads are numeric
```{r}
 kaiju_merged$reads <- as.numeric(kaiju_merged$reads)
 
```
### Replace NA or empty in Order and Genus with "Unclassified"
```{r}
kaiju_merged$Order[is.na(kaiju_merged$Order) | kaiju_merged$Order == ""] <- "Unclassified"
 kaiju_merged$Genus[is.na(kaiju_merged$Genus) | kaiju_merged$Genus == ""] <- "Unclassified"
 
 
```
### --- Urban and Rural lists adapted with suffix "_kaiju.out" ---
```{r}
 urban_files_kaiju <- paste0(c("37082_2#1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                               "37082_2#14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                              "37082_1#26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                               "37082_1#13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                               "37035_2#23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                               "37082_1#10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                               "37035_2#3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                               "36703_3#2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                               "37082_2#13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                               "37082_1#16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                               "37082_2#11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                               "37035_2#1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                               "37035_2#11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                               "37082_1#21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                               "37035_2#20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                             "36703_3#10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                               "37035_2#26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22"), 
                             "_kaiju.out")
  rural_files_kaiju <- paste0(c("37082_3#17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                               "36703_3#26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                              "37082_3#7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                               "37035_1#18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                               "37035_7#7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                               "37082_2#16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                               "37082_2#19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                               "37035_7#5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                              "37035_1#10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                               "37082_3#19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                               "37035_1#15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                               "37082_3#24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                               "37082_3#26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                               "37082_2#20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                               "37082_3#3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                               "37035_7#22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                               "37035_1#3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                               "37082_2#21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                               "37082_2#26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                               "36703_3#28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                               "37082_3#28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                               "37082_3#23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                               "37082_3#14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                               "37035_7#17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                               "37035_7#18", "37035_7#23", "37035_7#25"),
                             "_kaiju.out")
 
```
### --- Correctly assign file_base and Group ---
```{r}
 kaiju_merged <- kaiju_merged %>%
   mutate(
     file_base = paste0(gsub("\\.out$", "", basename(file)), "_kaiju.out"),
     Group = case_when(
       file_base %in% urban_files_kaiju ~ "Urban",
       file_base %in% rural_files_kaiju ~ "Rural",
       TRUE ~ NA_character_
     )
   ) %>%
   filter(!is.na(Group))
 
 
```
### Save as CSV
```{r}
write.csv(kaiju_merged, "kaiju_merged_final.csv", row.names = FALSE)
 
```
### Save as TXT (tab-delimited)
```{r}
write.table(kaiju_merged, "kaiju_merged_final.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)
 
```

### GENERAL ALPHA DIVERSITY

#### Ensure independence
```{r}

rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv); invisible(gc())

```
#### --- 1. DATA LOADING AND CLEANING ---
```{r}
kaiju_merged <- read.csv(
  file = "/home/alumno21/axel/files/kaiju_merged_final.csv",
  header = TRUE,
  sep = ",",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

```
#### --- Urban and Rural lists adapted with suffix "_kaiju.out" ---
```{r}
urban_files_kaiju <- paste0(c("37082_2#1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                              "37082_2#14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                              "37082_1#26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                              "37082_1#13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                              "37035_2#23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                              "37082_1#10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                              "37035_2#3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                              "36703_3#2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                              "37082_2#13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                              "37082_1#16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                              "37082_2#11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                              "37035_2#1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                              "37035_2#11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                              "37082_1#21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                              "37035_2#20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                              "36703_3#10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                              "37035_2#26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22"), 
                            "_kaiju.out")

rural_files_kaiju <- paste0(c("37082_3#17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                              "36703_3#26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                              "37082_3#7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                              "37035_1#18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                              "37035_7#7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                              "37082_2#16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                              "37082_2#19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                              "37035_7#5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                              "37035_1#10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                              "37082_3#19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                              "37035_1#15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                              "37082_3#24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                              "37082_3#26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                              "37082_2#20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                              "37082_3#3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                              "37035_7#22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                              "37035_1#3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                              "37082_2#21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                              "37082_2#26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                              "36703_3#28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                              "37082_3#28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                              "37082_3#23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                              "37082_3#14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                              "37035_7#17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                              "37035_7#18", "37035_7#23", "37035_7#25"),
                            "_kaiju.out")

```
#### --- Correctly assign file_base and Group ---
```{r}
kaiju_merged <- kaiju_merged %>%
  mutate(
    file_base = paste0(gsub("\\.out$", "", basename(file)), "_kaiju.out"),
    Group = case_when(
      file_base %in% urban_files_kaiju ~ "Urban",
      file_base %in% rural_files_kaiju ~ "Rural",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

```
#### ---------- 0) Finest available taxonomic key ----------
#### (kept for compatibility, but not used because we work at Genus level)
```{r}
tax_key <- if ("taxon_id" %in% names(kaiju_merged)) "taxon_id" else "Organism"

```
#### ---------- 1) Abundances per sample at GENUS level ----------
##### - We do NOT filter Domain/Kingdom.
##### - Consolidated by sample + Genus.
```{r}
alltaxa_by_sample <- kaiju_merged %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)
  ) %>%
  group_by(file_base, Group, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  group_by(file_base, Group) %>%
  mutate(rel_abund = reads / sum(reads, na.rm = TRUE)) %>%
  ungroup()

```
#### (Optional check)
```{r}
message("Rango de reads totales por muestra (nivel Genus):")
print(
  alltaxa_by_sample %>%
    group_by(file_base) %>%
    summarise(TotalReads = sum(reads), .groups = "drop") %>%
    summarise(min = min(TotalReads), q1 = quantile(TotalReads, .25),
              median = median(TotalReads), q3 = quantile(TotalReads, .75),
              max = max(TotalReads))
)

```
#### ---------- 2) Alpha diversity metrics per sample ----------
```{r}
alpha_per_sample <- alltaxa_by_sample %>%
  group_by(file_base, Group) %>%
  summarise(
    observed_taxa = sum(reads > 0),                                   # riqueza total (nº de géneros)
    H = { p <- rel_abund[rel_abund > 0]; -sum(p * log(p)) },          # Shannon (ln)
    simpson = 1 - sum(rel_abund^2),                                   # Gini–Simpson
    pielou = ifelse(observed_taxa > 1, H / log(observed_taxa), NA),   # Pielou
    berger_parker = max(rel_abund),                                   # Dominancia
    .groups = "drop"
  )

```
#### ---------- 3) Long format for faceting ----------
```{r}
alpha_long <- alpha_per_sample %>%
  pivot_longer(
    cols = c(observed_taxa, H, simpson, pielou, berger_parker),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("observed_taxa","H","simpson","pielou","berger_parker"),
      labels = c("Observed Taxa","Shannon Index","Simpson Index",
                 "Pielou's Evenness","Berger–Parker")
    )
  )

```
#### ---------- 4) n labels per group ----------
```{r}
n_labels <- alpha_long %>%
  group_by(metric, Group) %>%
  summarise(n = dplyr::n(), ymax = max(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = paste0("n = ", n))

```
#### ---------- 5) Wilcoxon and stars ----------
```{r}
p_stars <- function(p){
  if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**"
  else if (p < 0.05) "*" else ""
}

wilcox_df <- alpha_long %>%
  group_by(metric) %>%
  summarise(
    p = tryCatch(wilcox.test(value ~ Group)$p.value, error = function(e) NA_real_),
    ymax = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(star = vapply(p, p_stars, character(1)),
         x = 1.5, y = ymax * 1.02)

```
#### ---------- 6) Style and plot ----------
```{r}
pal <- c(Rural = "#E9B44C", Urban = "#4F86C6")

library(grid)  # por unit()

```
#### --- Clean base theme ---
```{r}
base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
```
#### Only bottom and left borders (axis lines)
```{r}
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.4, colour = "black"),
    axis.line.y = element_line(linewidth = 0.4, colour = "black"),
    
```
#### Facets
```{r}
    strip.placement = "outside",
    strip.clip = "off",
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text.x = element_text(
      size = 9,
      hjust = 0, vjust = 0.5,
      margin = margin(t = 6, r = 6, b = 6, l = 6)
    ),
    
    panel.spacing.x = unit(16, "pt"),
    panel.spacing.y = unit(8, "pt"),
    plot.title = element_text(hjust = 0.5, face = "plain"),
    plot.title.position = "plot",
    plot.margin = margin(12, 16, 12, 12),
    legend.position = "right"
  )

p_alpha <- ggplot(alpha_long, aes(Group, value, fill = Group)) +
  geom_boxplot(width = 0.65, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.6) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1, strip.position = "top") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = NULL,
       title = "Alpha diversity analysis across lifestyle groups (all taxa)") +
  base_theme +
  coord_cartesian(clip = "off") +
  geom_text(
    data = wilcox_df,
    aes(x = x, y = y, label = star),
    inherit.aes = FALSE, size = 8
  )

p_alpha


```
### ALPHA DIVERSITY OF BACTERIA
#### -------- PALETTE AND THEME --------
```{r}
pal <- c(Rural = "#E9B44C", Urban = "#4F86C6")

base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(linewidth = 0.4, colour = "black"),
    axis.line.y = element_line(linewidth = 0.4, colour = "black"),
    strip.placement = "outside", strip.clip = "off",
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text.x = element_text(size = 9, hjust = 0, vjust = 0.5,
                                margin = margin(t = 6, r = 6, b = 6, l = 6)),
    panel.spacing.x = unit(16, "pt"),
    panel.spacing.y = unit(8, "pt"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(12, 16, 12, 12)
  )

p_stars <- function(p){
  if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**"
  else if (p < 0.05) "*" else "ns"
}

```
#### -------- 1) FILTER BACTERIA AND WORK AT GENUS LEVEL --------
```{r}
bact_genus <- kaiju_merged %>%
  filter(Domain == "Bacteria") %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>%
  group_by(file_base, Group, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  group_by(file_base, Group) %>%
  mutate(rel_abund = reads / sum(reads, na.rm = TRUE)) %>%
  ungroup()

```
#### --- Optional checks ---
```{r}
message("Lecturas TOT por muestra (Bacteria, nivel Genus):")
print(bact_genus %>%
        group_by(file_base) %>%
        summarise(TotalReads = sum(reads), .groups = "drop") %>%
        summarise(min = min(TotalReads), q1 = quantile(TotalReads, .25),
                  median = median(TotalReads), q3 = quantile(TotalReads, .75),
                  max = max(TotalReads)))

message("Top 10 géneros por lecturas (global, Bacteria):")
print(bact_genus %>% group_by(Genus) %>%
        summarise(reads = sum(reads), .groups = "drop") %>%
        arrange(desc(reads)) %>% head(10))

```
#### -------- 2) ALPHA METRICS PER SAMPLE --------
```{r}
alpha_bact_genus <- bact_genus %>%
  group_by(file_base, Group) %>%
  summarise(
    observed_taxa = sum(reads > 0),                          # riqueza (nº de géneros)
    H = { p <- rel_abund[rel_abund > 0]; -sum(p * log(p)) }, # Shannon (ln)
    simpson = 1 - sum(rel_abund^2),                          # Gini–Simpson
    pielou = ifelse(observed_taxa > 1, H / log(observed_taxa), NA_real_),
    berger_parker = max(rel_abund),                          # Dominancia
    .groups = "drop"
  )

alpha_long <- alpha_bact_genus %>%
  pivot_longer(
    cols = c(observed_taxa, H, simpson, pielou, berger_parker),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric,
                    levels = c("observed_taxa","H","simpson","pielou","berger_parker"),
                    labels = c("Observed","Shannon","Simpson",
                               "Pielou","Berger–Parker"))
  )

```
#### -------- 3) WILCOXON AND STARS --------
```{r}
wilcox_df <- alpha_long %>%
  group_by(metric) %>%
  summarise(
    p = tryCatch(wilcox.test(value ~ Group)$p.value, error = function(e) NA_real_),
    ymax = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(star = vapply(p, p_stars, character(1)),
         x = 1.5, y = ymax * 1.02)

```
#### -------- 4) PLOT --------
```{r}
p_alpha_bact_genus <- ggplot(alpha_long, aes(Group, value, fill = Group)) +
  geom_boxplot(width = 0.65, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.6) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1, strip.position = "top") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = NULL,
       title = "Alpha diversity across lifestyle groups (Bacteria at Genus level)") +
  base_theme +
  coord_cartesian(clip = "off") +
  geom_text(data = wilcox_df, aes(x = x, y = y, label = star),
            inherit.aes = FALSE, size = 8)

p_alpha_bact_genus


```
### ALPHA DIVERSITY OF EUKARYOTA
#### -------- PALETTE AND THEME --------
```{r}
pal <- c(Rural = "#E9B44C", Urban = "#4F86C6")

base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(linewidth = 0.4, colour = "black"),
    axis.line.y = element_line(linewidth = 0.4, colour = "black"),
    strip.placement = "outside", strip.clip = "off",
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text.x = element_text(size = 9, hjust = 0, vjust = 0.5,
                                margin = margin(t = 6, r = 6, b = 6, l = 6)),
    panel.spacing.x = unit(16, "pt"),
    panel.spacing.y = unit(8, "pt"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(12, 16, 12, 12)
  )

p_stars <- function(p){
  if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**"
  else if (p < 0.05) "*" else ""
}

```
#### -------- 1) FILTER EUKARYOTA AND WORK AT GENUS LEVEL --------
```{r}
euk_genus <- kaiju_merged %>%
  filter(Domain == "Eukaryota") %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>%
  group_by(file_base, Group, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  group_by(file_base, Group) %>%
  mutate(rel_abund = reads / sum(reads, na.rm = TRUE)) %>%
  ungroup()

```
#### --- Optional checks ---
```{r}
message("Lecturas TOT por muestra (Eukaryota, nivel Genus):")
print(euk_genus %>%
        group_by(file_base) %>%
        summarise(TotalReads = sum(reads), .groups = "drop") %>%
        summarise(min = min(TotalReads), q1 = quantile(TotalReads, .25),
                  median = median(TotalReads), q3 = quantile(TotalReads, .75),
                  max = max(TotalReads)))

message("Top 10 géneros por lecturas (global, Eukaryota):")
print(euk_genus %>% group_by(Genus) %>%
        summarise(reads = sum(reads), .groups = "drop") %>%
        arrange(desc(reads)) %>% head(10))

```
#### -------- 2) ALPHA METRICS PER SAMPLE --------
```{r}
alpha_euk_genus <- euk_genus %>%
  group_by(file_base, Group) %>%
  summarise(
    observed_taxa = sum(reads > 0),                          # riqueza (nº de géneros)
    H = { p <- rel_abund[rel_abund > 0]; -sum(p * log(p)) }, # Shannon (ln)
    simpson = 1 - sum(rel_abund^2),                          # Gini–Simpson
    pielou = ifelse(observed_taxa > 1, H / log(observed_taxa), NA_real_),
    berger_parker = max(rel_abund),                          # Dominancia
    .groups = "drop"
  )

alpha_long <- alpha_euk_genus %>%
  pivot_longer(
    cols = c(observed_taxa, H, simpson, pielou, berger_parker),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric,
                    levels = c("observed_taxa","H","simpson","pielou","berger_parker"),
                    labels = c("Observed","Shannon","Simpson ",
                               "Pielou","Berger–Parker"))
  )

```
#### -------- 3) WILCOXON AND STARS --------
```{r}
wilcox_df <- alpha_long %>%
  group_by(metric) %>%
  summarise(
    p = tryCatch(wilcox.test(value ~ Group)$p.value, error = function(e) NA_real_),
    ymax = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(star = vapply(p, p_stars, character(1)),
         x = 1.5, y = ymax * 1.02)

```
#### -------- 4) PLOT --------
```{r}
p_alpha_euk_genus <- ggplot(alpha_long, aes(Group, value, fill = Group)) +
  geom_boxplot(width = 0.65, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.6) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1, strip.position = "top") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = NULL,
       title = "Alpha diversity across lifestyle groups (Eukaryota at Genus level)") +
  base_theme +
  coord_cartesian(clip = "off") +
  geom_text(data = wilcox_df, aes(x = x, y = y, label = star),
            inherit.aes = FALSE, size = 8)

p_alpha_euk_genus
```
### BETA DIVERSITY
### --- Libraries ---
```{r}
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(patchwork)
library(uwot)
library(tibble)
library(grid)
library(gtable)
library(gridExtra)

set.seed(42)  # Para reproducibilidad

kaiju_merged$reads <- as.numeric(kaiju_merged$reads)
kaiju_merged$Order[is.na(kaiju_merged$Order) | kaiju_merged$Order == ""] <- "Unclassified"

```
#### --- Create abundance matrix at Order level ---
```{r}
abundance_matrix <- kaiju_merged %>%
  group_by(file_base, Order) %>%
  summarise(reads = sum(reads, na.rm=TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = reads, values_fill = 0) %>%
  column_to_rownames("file_base")

```
#### --- Normalize to relative abundance ---
```{r}
abundance_matrix_rel <- abundance_matrix / rowSums(abundance_matrix)

```
#### --- Metadata for PERMANOVA and plots ---
```{r}
metadata_df <- data.frame(Sample = rownames(abundance_matrix_rel)) %>%
  mutate(Group = case_when(
    Sample %in% urban_files_kaiju ~ "Urban",   # ← usa las listas CON sufijo
    Sample %in% rural_files_kaiju ~ "Rural",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("Urban","Rural")))

abundance_matrix_rel <- abundance_matrix_rel[metadata_df$Sample, , drop = FALSE]


```
#### --- Calculate Bray-Curtis distances ---
```{r}
dist_bray <- vegdist(abundance_matrix_rel, method = "bray")

```
#### --- PERMANOVA ---
```{r}
permanova_res <- adonis2(dist_bray ~ Group, data = metadata_df)
print(permanova_res)

r2_value <- round(permanova_res$R2[1], 3)
p_value <- permanova_res$`Pr(>F)`[1]

```
#### --- NMDS ---
```{r}
nmds_res <- metaMDS(dist_bray, k=2, trymax=100)
nmds_df <- as.data.frame(nmds_res$points) %>%
  mutate(Sample = rownames(nmds_res$points)) %>%
  left_join(metadata_df, by = "Sample")

```
#### --- PCoA ---
```{r}
pcoa_res <- cmdscale(dist_bray, k=2, eig=TRUE)
pcoa_df <- as.data.frame(pcoa_res$points) %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  mutate(Sample = rownames(abundance_matrix_rel)) %>%
  left_join(metadata_df, by = "Sample")

```
#### --- UMAP ---
```{r}
umap_res <- umap(abundance_matrix_rel, n_neighbors=15, metric="euclidean")
umap_df <- as.data.frame(umap_res) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(Sample = rownames(abundance_matrix_rel)) %>%
  left_join(metadata_df, by = "Sample")

```
#### --- PCA (on relative abundance matrix) ---
```{r}
pca_res <- prcomp(abundance_matrix_rel, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x[,1:2]) %>%
  mutate(Sample = rownames(abundance_matrix_rel)) %>%
  left_join(metadata_df, by = "Sample")

```
#### --- Base function for ggplot ---
```{r}
base_plot <- function(df, x, y, title) {
  ggplot(df, aes_string(x=x, y=y, color="Group")) +
    geom_point(size=1.5, alpha=0.5) +
    scale_color_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
    theme_minimal() +
    labs(title = title, color = "Group") +
    theme(plot.title = element_text(hjust = 0.15, face = "plain", size=10),
          axis.title = element_text(face = "plain", size=9),   # Títulos de ejes sin negrita
          axis.text = element_text(face = "plain"))
}

```
#### --- Create plots without legend ---
```{r}
p_nmds <- base_plot(nmds_df, "MDS1", "MDS2", "NMDS (Bray-Curtis)") + theme(legend.position = "none")
p_pcoa <- base_plot(pcoa_df, "PCoA1", "PCoA2", "PCoA (Bray-Curtis)") + theme(legend.position = "none")
p_umap <- base_plot(umap_df, "UMAP1", "UMAP2", "UMAP") + theme(legend.position = "none")
p_pca  <- base_plot(pca_df, "PC1", "PC2", "PCA") + theme(legend.position = "none")

```
#### --- Extract legend from p_pcoa to place on the right ---
```{r}
get_legend <- function(a_gplot) {
  tmp <- ggplotGrob(a_gplot)
  legend_index <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[legend_index]]
  return(legend)
}

legend_plot <- base_plot(pcoa_df, "PCoA1", "PCoA2", "PCoA (Bray-Curtis)") + theme(legend.position = "right")

legend_grob <- get_legend(legend_plot)

```
#### --- Create statistical text ---
```{r}
annot_nms <- paste(
  sprintf("PERMANOVA R² = %.3f", r2_value),
  sprintf("p = %.3g", p_value),
  sprintf("Stress NMDS = %.2f", nmds_res$stress),
  sprintf("n (Rural) = %d", sum(metadata_df$Group == "Rural")),
  sprintf("n (Urban) = %d", sum(metadata_df$Group == "Urban")),
  sep = "\n"
)

caption_title <- textGrob(
  "Statistics",
  gp = gpar(fontsize = 8, fontface = "bold"),
  x = unit(0, "npc"), hjust = 0
)
space <- nullGrob()
space_height <- unit(4, "mm")
caption_body <- textGrob(
  annot_nms,
  gp = gpar(fontsize = 8, fontface = "plain"),
  x = unit(0, "npc"), hjust = 0
)
caption_grob <- arrangeGrob(
  grobs = list(caption_title, space, caption_body),
  ncol = 1,
  heights = unit.c(grobHeight(caption_title), space_height, grobHeight(caption_body))
)

space_between_legend_and_text <- unit(15, "mm")
legend_table <- gtable(
  widths = unit(0.3, "npc"),
  heights = unit.c(grobHeight(legend_grob), space_between_legend_and_text, grobHeight(caption_grob))
)
legend_table <- gtable_add_grob(legend_table, legend_grob, t = 1, l = 1)
legend_table <- gtable_add_grob(legend_table, caption_grob, t = 3, l = 1)

```
#### 2x2 plots grid without legend
```{r}
plots_grid <- (p_nmds + p_pcoa) / (p_umap + p_pca) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    plot.margin = margin(2, 2, 0, 2)
  )

```
#### Combine plots and legend + statistics closer, with adjusted widths
```{r}
combined_plot <- wrap_elements(plots_grid) + wrap_elements(legend_table) +
  plot_layout(widths = c(8, 5))

print(combined_plot)

library(ggplot2)
library(patchwork)
library(gridExtra)

```
#### 1) Comfortable margins per subplot (neither tight nor loose)
```{r}
tema_solo_dos_lineas_comfy <- theme(
  panel.background   = element_rect(fill = "white", colour = NA),
  plot.background    = element_rect(fill = "white", colour = NA),
  panel.grid.major   = element_blank(),
  panel.grid.minor   = element_blank(),
  axis.text          = element_text(colour = "black"),
  axis.title         = element_text(colour = "black"),
  axis.ticks         = element_line(colour = "black"),
  axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
  axis.line.x.top    = element_blank(),
  axis.line.y.left   = element_line(colour = "black", linewidth = 0.5),
  axis.line.y.right  = element_blank(),
  plot.margin        = margin(6, 10, 6, 10)  # ↑ dale un poco de aire lateral
)

p_nmds <- p_nmds + tema_solo_dos_lineas_comfy
p_pcoa <- p_pcoa + tema_solo_dos_lineas_comfy
p_umap <- p_umap + tema_solo_dos_lineas_comfy
p_pca  <- p_pca  + tema_solo_dos_lineas_comfy

```
#### 2) Spacing between grid panels
```{r}
plots_grid <- (p_nmds + p_pcoa) / (p_umap + p_pca) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    panel.spacing   = unit(6, "pt"),     # ← pequeño pero evita choques
    plot.margin     = margin(0, 0, 0, 0)
  )

```
#### 3) Tags without borders, but with minimal spacing
```{r}
plots_grid_tagged <- plots_grid +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag          = element_text(size = 12),
    plot.tag.position = c(0.02, 0.98),
    plot.margin       = margin(2, 2, 2, 2)
  )

```
#### 4) Right column narrower with moderate spacing
```{r}
space_between_legend_and_text <- unit(15, "mm")  # ↓ antes 15 mm
legend_table <- gtable(
  widths  = unit(0.24, "npc"),
  heights = unit.c(grobHeight(legend_grob),
                   space_between_legend_and_text,
                   grobHeight(caption_grob))
)
legend_table <- gtable_add_grob(legend_table, legend_grob, t = 1, l = 1)
legend_table <- gtable_add_grob(legend_table, caption_grob, t = 3, l = 1)

```
#### 5) EXTERNAL margin of the set (title not stuck or clashing)
```{r}
combined_plot <-
  wrap_elements(plots_grid_tagged) + wrap_elements(legend_table) +
  plot_layout(widths = c(9.5, 3.5)) +
  plot_annotation(
    title = "Rural and urban beta diversity",
    theme = theme(
      plot.title   = element_text(hjust = 0.5, size = 14,
                                  margin = margin(b = 8)),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin  = margin(8, 12, 8, 12)  # ← colchón externo uniforme
    )
  )

```
### FINAL WELL-FORMATTED BETA DIVERSITY GRAPH
```{r}
print(combined_plot)

```
### PLOT ALPHA AND BETA DIVERSITY TOGETHER
```{r}

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(patchwork)

```
### ========= PALETTE AND THEME =========
```{r}
pal <- c(Rural = "#E9B44C", Urban = "#4F86C6")

base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line.x       = element_line(linewidth = 0.4, colour = "black"),
    axis.line.y       = element_line(linewidth = 0.4, colour = "black"),
    strip.placement   = "outside", strip.clip = "off",
    strip.background  = element_rect(fill = "white", colour = NA),
    strip.text.x      = element_text(size = 10, hjust = 0.5),
    panel.spacing.x   = unit(12, "pt"),
    panel.spacing.y   = unit(8,  "pt"),
    legend.position   = "right",
    plot.margin       = margin(10, 16, 10, 12)
  )

```
### ========= ALPHA PANEL FUNCTION (reusable) =========
```{r}
alpha_panel <- function(alpha_long_df, titulo = NULL) {
```
### shorter labels and line breaks where helpful
```{r}
  alpha_long_df <- alpha_long_df %>%
    mutate(metric = factor(metric,
                           levels = c("Observed","Shannon","Simpson","Pielou","Berger–Parker"),
                           labels = c("Observed","Shannon","Simpson","Pielou","Berger–\nParker")
    ))
  
```
### p-value and stars
```{r}
  wilcox_df <- alpha_long_df %>%
    group_by(metric) %>%
    summarise(
      p    = tryCatch(wilcox.test(value ~ Group)$p.value, error = function(e) NA_real_),
      ymax = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      star = case_when(
        is.na(p)      ~ "",
        p < 0.001     ~ "***",
        p < 0.01      ~ "**",
        p < 0.05      ~ "*",
        TRUE          ~ ""
      ),
      x = 1.5, y = ymax * 1.04
    )
  
  ggplot(alpha_long_df, aes(Group, value, fill = Group, color = Group)) +
    geom_boxplot(width = 0.6, alpha = 0.9, outlier.shape = NA, color = NA) +
    geom_jitter(width = 0.12, size = 0.7, alpha = 0.55, show.legend = FALSE) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_color_manual(values = pal, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
    labs(x = NULL, y = NULL, title = titulo) +
    base_theme +
    coord_cartesian(clip = "off") +
    geom_text(data = wilcox_df, aes(x = x, y = y, label = star),
              inherit.aes = FALSE, size = 4) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
}

```
### ========= YOUR TWO ALPHA PANELS =========
```{r}
p1 <- alpha_panel(alpha_euk_genus %>%
                    pivot_longer(c(observed_taxa, H, simpson, pielou, berger_parker),
                                 names_to = "metric", values_to = "value") %>%
                    mutate(metric = recode(metric,
                                           observed_taxa = "Observed", H = "Shannon",
                                           simpson = "Simpson", pielou = "Pielou",
                                           berger_parker = "Berger–Parker")),
                  "") +
  theme(legend.position = "none")

p2 <- alpha_panel(alpha_bact_genus %>%
                    pivot_longer(c(observed_taxa, H, simpson, pielou, berger_parker),
                                 names_to = "metric", values_to = "value") %>%
                    mutate(metric = recode(metric,
                                           observed_taxa = "Observed", H = "Shannon",
                                           simpson = "Simpson", pielou = "Pielou",
                                           berger_parker = "Berger–Parker")),
                  "")

```
### ========= PCA (without grid, more space) =========
```{r}
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(x = "PC1", y = "PC2", title = "") +
  scale_color_manual(values = c(Rural = "#E69F00", Urban = "#0072B2")) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "plain"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(hjust = 0, face = "bold"),
```
### This adds the full border that closes the box:
```{r}
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

p3 <- p_pca +
  scale_color_manual(values = pal, name = "Group") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left   = element_line(),
    legend.position = "none",
    plot.margin = margin(6, 10, 6, 8)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)))

p3 <- p_pca +
  scale_color_manual(values = pal, name = "Group") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left   = element_line(),
    legend.position = "none",
    plot.margin = margin(8, 14, 8, 14)  # más aire en todos los lados
  ) +
```
### more space at the edges
```{r}
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08)))
```
### ========= FINAL COMPOSITION (vertical, breathing) =========
```{r}
final_plot <- (p1 / p2 / p3) +
  plot_layout(heights = c(1, 1, 1.3), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "plain"),
    plot.tag        = element_text(size = 14, face = "plain"),
    panel.spacing   = unit(10, "pt")
  )

final_plot
```

### Differentially abundant bacterial and eukaryotic taxa between Rural and Urban groups

#### Requiered libraries: dplyr, tidyr, stringr, ggplot2, purrr, broom
```{r}

rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv); invisible(gc())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})


```
### --- 1. DATA LOADING AND CLEANING ---
```{r}
kaiju_merged <- read.csv(
  file = "/home/alumno21/axel/files/kaiju_merged_final.csv",
  header = TRUE,
  sep = ",",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

```
### --- Urban and Rural lists adapted with suffix "_kaiju.out" ---
```{r}
urban_files_kaiju <- paste0(c("37082_2 # 1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                              "37082_2 # 14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                              "37082_1 # 26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                              "37082_1 # 13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                              "37035_2 # 23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                              "37082_1 # 10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                              "37035_2 # 3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                              "36703_3 # 2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                              "37082_2 # 13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                              "37082_1 # 16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                              "37082_2 # 11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                              "37035_2 # 1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                              "37035_2 # 11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                              "37082_1 # 21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                              "37035_2 # 20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                              "36703_3 # 10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                              "37035_2 # 26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22"),
                            "_kaiju.out")

rural_files_kaiju <- paste0(c("37082_3 # 17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                              "36703_3 # 26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                              "37082_3 # 7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                              "37035_1 # 18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                              "37035_7 # 7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                              "37082_2 # 16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                              "37082_2 # 19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                              "37035_7 # 5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                              "37035_1 # 10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                              "37082_3 # 19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                              "37035_1 # 15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                              "37082_3 # 24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                              "37082_3 # 26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                              "37082_2 # 20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                              "37082_3 # 3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                              "37035_7 # 22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                              "37035_1 # 3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                              "37082_2 # 21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                              "37082_2 # 26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                              "36703_3 # 28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                              "37082_3 # 28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                              "37082_3 # 23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                              "37082_3 # 14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                              "37035_7 # 17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                              "37035_7 # 18", "37035_7#23", "37035_7#25"),
                            "_kaiju.out")

```
### --- Correctly assign file_base and Group ---
```{r}
kaiju_merged <- kaiju_merged %>%
  mutate(
    file_base = paste0(gsub("\\.out$", "", basename(file)), "_kaiju.out"),
    Group = case_when(
      file_base %in% urban_files_kaiju ~ "Urban",
      file_base %in% rural_files_kaiju ~ "Rural",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(purrr); library(broom)
})

```
### ---------- Parameter: genera of interest ----------
```{r}
target_genera <- c("Clostridium", "Blautia", "Ruminococcus", "Prevotella", "Streptococcus")

```
### ---------- Utilities ----------
```{r}
p_to_stars <- function(p){
  dplyr::case_when(
    is.na(p)   ~ "",
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ ""
  )
}

```
### n labels per group
```{r}
n_by_group <- kaiju_merged %>%
  distinct(file_base, Group) %>%
  count(Group) %>%
  mutate(label = paste0(Group, " (n=", n, ")"))
lab_urban <- n_by_group %>% filter(Group=="Urban") %>% pull(label)
lab_rural <- n_by_group %>% filter(Group=="Rural") %>% pull(label)

```
### ---------- Per-sample summary: reads per Genus (only Bacteria, only selected) ----------
```{r}
per_sample_genus <- kaiju_merged %>%
  filter(Domain == "Bacteria",
         !is.na(Genus), Genus != "", Genus != "Unclassified",
         Genus %in% target_genera) %>%
  group_by(file_base, Group, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
```
### Guarantee the presence of all selected genus (0 if absent)
```{r}
  complete(file_base, Group, Genus = target_genera, fill = list(reads = 0))

```
### Maintain th order of the genus like in target_genera 
```{r}
per_sample_genus$Genus <- factor(per_sample_genus$Genus, levels = target_genera)

```
### ---------- Prevalence and average reads per Group (only selected) ----------
### Prevalence: % of samples of the group with reads > 0 for that genus
```{r}
preval_mean <- per_sample_genus %>%
  group_by(Group) %>%
  mutate(n_group = n_distinct(file_base)) %>%
  group_by(Group, Genus, n_group) %>%
  summarise(
    prevalence = 100 * mean(reads > 0, na.rm = TRUE),
    mean_reads = mean(reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Group, match(Genus, target_genera))

```
### (Optional) print summary table
```{r}
cat("\n=== Géneros seleccionados — prevalencia y media de lecturas ===\n")
print(preval_mean %>% mutate(prevalence = round(prevalence,1),
mean_reads = round(mean_reads,1)))

```
### --- Data for the plot ---
```{r}
plot_df_all <- per_sample_genus %>%
  mutate(log_reads = log10(reads + 1),
         Group_lab = ifelse(Group=="Urban", lab_urban, lab_rural))

```
### For boxplot: exclude 0 (presence only)
```{r}
plot_df_nz <- plot_df_all %>% filter(reads > 0)

```
### --- Wilcoxon per genus (you can keep all points, including zeros) ---
```{r}
wilcox_tbl <- plot_df_all %>%
  group_by(Genus) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(log_reads ~ Group, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(stars = p_to_stars(p_value)) %>%
```
### Position fo the text: use the maximum of the NO-CERO per genus
```{r}
  left_join(plot_df_nz %>% group_by(Genus) %>%
              summarise(y_pos = max(log_reads, na.rm = TRUE) + 0.15, .groups="drop"),
            by = "Genus")

```
### --- Theme (same) ---
```{r}
base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left  = element_line(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(hjust = 0.5, size = 8),
    plot.title = element_text(hjust = 0.5, size = 14, face = "plain"),
    plot.margin = margin(10, 10, 10, 10),
    panel.spacing = unit(6, "pt") # un pelín menos de espacio entre facet
  )

```
### --- Plot: boxplot with NON-ZEROS + jitter with ALL ---
```{r}
plot_boxbactpred <- ggplot() +
  geom_boxplot(
    data = plot_df_nz,
    aes(x = Group_lab, y = log_reads, fill = Group),
    width = 0.65, alpha = 0.9, outlier.shape = NA
  ) +
  geom_jitter(
    data = plot_df_nz, # <- sin ceros
    aes(x = Group_lab, y = log_reads),
    width = 0.15, size = 0.3, alpha = 0.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Genus, nrow = 1, scales = "free_y") +
  labs(
    title = "",
    x = "Group", y = "Reads (log10)"
  ) +
  scale_fill_manual(values = c("Urban" = " # 5DA5DA", "Rural" = "#F4A460")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.07))) +
  base_theme +
  geom_text(
    data = wilcox_tbl,
    aes(x = 1.5, y = y_pos, label = stars),
    inherit.aes = FALSE, size = 8
  )

print(plot_boxbactpred)




```
### EUKARYOTA

#### ---------- Parameters ----------
```{r}
target_taxa <- c("Agaricales", "Chorda", "Eimeriidae", "Halteriidae", "Saccharomycetales")

p_to_stars <- function(p){
  dplyr::case_when(
    is.na(p)   ~ "",
    p < 1e-4   ~ "****",
    p < 1e-3   ~ "***",
    p < 1e-2   ~ "**",
    p < 5e-2   ~ "*",
    TRUE       ~ ""
  )
}

```
### ---------- n labels per group ----------
```{r}
n_by_group <- kaiju_merged %>% distinct(file_base, Group) %>%
  count(Group) %>%
  mutate(label = paste0(Group, " (n=", n, ")"))
lab_urban <- n_by_group %>% filter(Group=="Urban") %>% pull(label)
lab_rural <- n_by_group %>% filter(Group=="Rural") %>% pull(label)

```
### ---------- Per-sample summary (EXCLUDING CHORDATA) ONLY target_taxa ----------
```{r}
per_sample_genus <- kaiju_merged %>%
  filter(
    Domain == "Eukaryota",
    is.na(Phylum) | !grepl("^Chordata$", Phylum, ignore.case = TRUE),
    !is.na(Genus), Genus != "", Genus != "Unclassified",
    Genus %in% target_taxa
  ) %>%
  group_by(file_base, Group, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
```
### Make sure of the presence of al the taxons of all the samples
```{r}
  complete(file_base, Group, Genus = target_taxa, fill = list(reads = 0))

```
### Order of facets according to target_taxa
```{r}
per_sample_genus$Genus <- factor(per_sample_genus$Genus, levels = target_taxa)

```
### ---------- Data for plot ----------
```{r}
plot_df_all <- per_sample_genus %>%
  mutate(
    log_reads = log10(reads + 1),
    Group_lab = ifelse(Group=="Urban", lab_urban, lab_rural)
  )

```
### For the jitter (without 0 to avoid the formation of "bands")
```{r}
plot_df_nz <- plot_df_all %>% filter(reads > 0)

```
### ---------- Wilcoxon per genus ----------
```{r}
wilcox_tbl <- plot_df_all %>%
  group_by(Genus) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(log_reads ~ Group, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = max(plot_df_nz$log_reads[plot_df_nz$Genus==first(Genus)], na.rm = TRUE) + 0.15,
    .groups = "drop"
  ) %>%
  mutate(stars = p_to_stars(p_value))

```
### ---------- Theme ----------
```{r}
base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left  = element_line(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(hjust = 0.5, size = 9),
    plot.title = element_text(hjust = 0.5, size = 14, face = "plain"),
    plot.margin = margin(10, 10, 10, 10)
  )

```
### ---------- Plot ----------
```{r}
p_euk_sel <- ggplot() +
  geom_boxplot(
    data = plot_df_nz,
    aes(x = Group_lab, y = log_reads, fill = Group),
    width = 0.65, alpha = 0.9, outlier.shape = NA
  ) +
  geom_jitter(
    data = plot_df_nz, # sin ceros
    aes(x = Group_lab, y = log_reads),
    width = 0.15, size = 0.3, alpha = 0.5, inherit.aes = FALSE
  ) +
  facet_wrap(~ Genus, nrow = 1, scales = "free_y") +
  labs(
    title = "",
    x = "Group", y = "Reads (log10)"
  ) +
  scale_fill_manual(values = c("Urban" = " # 5DA5DA", "Rural" = "#F4A460")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.07))) +
  base_theme +
  geom_text(
    data = wilcox_tbl,
    aes(x = 1.5, y = y_pos, label = stars),
    inherit.aes = FALSE, size = 8
  )

print(p_euk_sel)


library(patchwork)

final_fig5 <- plot_boxbactpred / (p_euk_sel+ theme(legend.position = "none")) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 14, face = "plain") # estilo de las letras
  )

final_fig5
```
## Mean relative abundance of the most prevalent bacterial and eukaryotic genera in rural and urban samples.

### Graph boxplots of EUKARYOTA
```{r}
 
 
 rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv); invisible(gc())
 
 suppressPackageStartupMessages({
   library(dplyr); library(tidyr); library(ggplot2); library(scales)
 })
 
 
```
### --- 1. Upload and clean the data---
```{r}
 kaiju_merged <- read.csv(
   file = "/home/alumno21/axel/files/kaiju_merged_final.csv",
   header = TRUE,
   sep = ",",
   fileEncoding = "latin1",
   stringsAsFactors = FALSE,
   na.strings = c("", "NA")
 )
 
```
### --- Urban and Rural list adapted with the suffix "_kaiju.out" ---
```{r}
 urban_files_kaiju <- paste0(c("37082_2 # 1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                               "37082_2 # 14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                               "37082_1 # 26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                              "37082_1 # 13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                               "37035_2 # 23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                               "37082_1 # 10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                               "37035_2 # 3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                               "36703_3 # 2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                               "37082_2 # 13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                               "37082_1 # 16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                               "37082_2 # 11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                               "37035_2 # 1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                               "37035_2 # 11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                               "37082_1 # 21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                               "37035_2 # 20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                               "36703_3 # 10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                               "37035_2 # 26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22"),
                             "_kaiju.out")
 
 rural_files_kaiju <- paste0(c("37082_3 # 17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                               "36703_3 # 26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                               "37082_3 # 7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                               "37035_1 # 18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                               "37035_7 # 7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                               "37082_2 # 16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                               "37082_2 # 19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                               "37035_7 # 5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                               "37035_1 # 10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                               "37082_3 # 19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                               "37035_1 # 15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                               "37082_3 # 24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                               "37082_3 # 26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                               "37082_2 # 20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                               "37082_3 # 3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                               "37035_7 # 22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                               "37035_1 # 3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                               "37082_2 # 21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                               "37082_2 # 26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                               "36703_3 # 28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                               "37082_3 # 28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                               "37082_3 # 23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                               "37082_3 # 14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                               "37035_7 # 17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                               "37035_7 # 18", "37035_7#23", "37035_7#25"),
                             "_kaiju.out")
 
```
### --- Assign correctly flie_base and Group ---
```{r}
 kaiju_merged <- kaiju_merged %>%
   mutate(
     file_base = paste0(gsub("\\.out$", "", basename(file)), "_kaiju.out"),
     Group = case_when(
       file_base %in% urban_files_kaiju ~ "Urban",
       file_base %in% rural_files_kaiju ~ "Rural",
       TRUE ~ NA_character_
     )
   ) %>%
   filter(!is.na(Group))
 
 
```
### --- Selected taxa (level Genus) ---
```{r}
 target_taxa <- c("Agaricales", "Chorda", "Eimeriidae", "Halteriidae", "Saccharomycetales")
 
```
### --- Totals per sample (denominator): total of Eukaryota EXCLUDING Chordata ---
```{r}
 totals_sample <- kaiju_merged %>%
   filter(Domain == "Eukaryota",
          is.na(Phylum) | !grepl("^Chordata$", Phylum, ignore.case = TRUE)) %>%
   group_by(file_base, Group) %>%
   summarise(total_reads = sum(reads, na.rm = TRUE), .groups = "drop")
 
```
### --- Reads of the selected taxa (Genus), also EXCLUDING Chordata ---
```{r}
 sel_reads <- kaiju_merged %>%
   filter(Domain == "Eukaryota",
          is.na(Phylum) | !grepl("^Chordata$", Phylum, ignore.case = TRUE),
          !is.na(Genus), Genus != "", Genus != "Unclassified",
          Genus %in% target_taxa) %>%
   group_by(file_base, Group, Genus) %>%
   summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
   complete(file_base, Group, Genus = target_taxa, fill = list(reads = 0)) %>%
   left_join(totals_sample, by = c("file_base","Group")) %>%
   mutate(pct = if_else(total_reads > 0, 100 * reads / total_reads, NA_real_))
 
```
### Maintain the order of the taxons
```{r}
 sel_reads$Genus <- factor(sel_reads$Genus, levels = target_taxa)
 
```
### --- Average of % per group (for the bars)---
```{r}
 group_avg <- sel_reads %>%
   group_by(Group, Genus) %>%
   summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop")
 
```
### ====================== Wilcoxon by Genus ======================
```{r}
 get_pval_per_genus <- function(df) {
   df <- df %>% filter(!is.na(pct))
   if (n_distinct(df$Group) < 2) return(NA_real_)
   tryCatch(wilcox.test(pct ~ Group, data = df)$p.value, error = function(e) NA_real_)
 }
 
 wilcox_res <- sel_reads %>%
   group_by(Genus) %>%
   summarise(p = get_pval_per_genus(cur_data()), .groups = "drop") %>%
   mutate(stars = case_when(
     is.na(p)  ~ "",
     p <= 1e-4 ~ "****",
     p <= 1e-3 ~ "***",
     p <= 0.01 ~ "**",
     p <= 0.05 ~ "*",
     TRUE      ~ ""
   ))
 
```
### ----  Height of stars/brackets (with padding) ----
```{r}
 overall_max <- max(group_avg$mean_pct, na.rm = TRUE)
 y_pad <- overall_max * 0.06
 
 anno_y <- group_avg %>%
   group_by(Genus) %>%
   summarise(ypos = max(mean_pct, na.rm = TRUE) + y_pad, .groups = "drop")
 
 anno_df <- wilcox_res %>%
   left_join(anno_y, by = "Genus") %>%
   filter(stars != "")

```
### --- Palette for Group ---
```{r}
 pal_group <- c(Rural = " # E9B44C", Urban = "#4F86C6")
 
```
### ---------- Theme ----------
```{r}
 base_theme <- theme_minimal(base_size = 11) +
   theme(
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     axis.line.x.bottom = element_line(),
     axis.line.y.left  = element_line(),
     strip.text.x = element_text(hjust = 0.5, size = 9),
     plot.title = element_text(hjust = 0, size = 14, face = "plain",
                               margin = margin(t = 8, b = 6)),
     plot.margin = margin(15, 12, 12, 12)
   )
 
```
### ---------- Data for BRACKETS drawn by hand ----------
### Numerical position of the genus in the x axis
```{r}
 lvl <- levels(group_avg$Genus)
 brackets_manual <- anno_df %>%
   mutate(
     gx = match(Genus, lvl),
```
### Horizontal offset: center above the bars (adjust if width/dodge is changed)
```{r}
     off = 0.23, # ~ Half of the space between the two bars
     x_rural = gx - off,
     x_urban = gx + off,
     y = ypos
   )
 
```
### ---  Final graphic with Brackets + stars---
```{r}
 p_sep <- ggplot(group_avg, aes(x = Genus, y = mean_pct, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8),
            width = 0.7, color = "black", linewidth = 0.2) +
```
### % labels
```{r}
   geom_text(aes(label = sprintf("%.3f%%", mean_pct)),
             position = position_dodge(width = 0.8),
             vjust = -0.3, size = 3, color = "black") +
```
### Horizontal BRACKET (30% higher)
```{r}
   geom_segment(data = brackets_manual,
                aes(x = x_rural, xend = x_urban,
                    y = y * 1.7, yend = y * 1.7),
               inherit.aes = FALSE, linewidth = 0.6) +
```
### Bracket legs (2% below the new level)
```{r}
   geom_segment(data = brackets_manual,
                aes(x = x_rural, xend = x_rural,
                    y = y * 1.7, yend = y * 1.7 - (y * 1.7 * 0.02)),
               inherit.aes = FALSE, linewidth = 0.6) +
   geom_segment(data = brackets_manual,
                aes(x = x_urban, xend = x_urban,
                    y = y * 1.7, yend = y * 1.7 - (y * 1.7 * 0.02)),
                inherit.aes = FALSE, linewidth = 0.6) +
```
### Stars (also in the new level)
```{r}
   geom_text(data = brackets_manual,
             aes(x = gx, y = y * 1.7, label = stars),
             inherit.aes = FALSE, vjust = -0.2, size = 5) +
   scale_fill_manual(values = pal_group) +
   scale_y_continuous(
     labels = label_percent(scale = 1, accuracy = 0.1),
     expand = expansion(mult = c(0, 0.30)) # Superior air for brackets/stars
   ) +
   labs(
     title = "Mean relative abundance per sample of the most prevalent eukaryotic taxa",
     x = "Genus", y = "Mean per sample"
   ) +
   base_theme

```
### Corrected graph
```{r}
 print(p_sep)
 
 
 
 suppressPackageStartupMessages({
   library(dplyr); library(tidyr); library(ggplot2); library(scales)
 })
 
```
### --- Selected Ggenus (Domain = Bacteria) ---
```{r}
 target_genera <- c("Clostridium", "Faecalibacterium", "Ruminococcus", "Prevotella", "Streptococcus")
 
```
### --- Totals per sample: total of Bacteria per sample ---
```{r}
 totals_sample <- kaiju_merged %>%
   filter(Domain == "Bacteria") %>%
   group_by(file_base, Group) %>%
   summarise(total_reads = sum(reads, na.rm = TRUE), .groups = "drop")
 
```
### --- Reads of the selected genus ---
```{r}
 sel_reads <- kaiju_merged %>%
   filter(Domain == "Bacteria",
          !is.na(Genus), Genus != "", Genus != "Unclassified",
          Genus %in% target_genera) %>%
   group_by(file_base, Group, Genus) %>%
   summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
   complete(file_base, Group, Genus = target_genera, fill = list(reads = 0)) %>%
   left_join(totals_sample, by = c("file_base","Group")) %>%
   mutate(pct = if_else(total_reads > 0, 100 * reads / total_reads, NA_real_))
 
```
### Order of genus on the x axis
```{r}
 sel_reads$Genus <- factor(sel_reads$Genus, levels = target_genera)
 
```
### --- Average per group (for the bars)---
```{r}
 group_avg <- sel_reads %>%
   group_by(Group, Genus) %>%
   summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop")
 
```
### ====================== Wilcoxon by Genus ======================
```{r}
 get_pval_per_genus <- function(df) {
   df <- df %>% filter(!is.na(pct))
   if (n_distinct(df$Group) < 2) return(NA_real_)
   tryCatch(wilcox.test(pct ~ Group, data = df)$p.value, error = function(e) NA_real_)
 }
 
 wilcox_res <- sel_reads %>%
   group_by(Genus) %>%
   summarise(p = get_pval_per_genus(cur_data()), .groups = "drop") %>%
   mutate(stars = case_when(
     is.na(p)  ~ "",
     p <= 1e-4 ~ "****",
     p <= 1e-3 ~ "***",
     p <= 0.01 ~ "**",
     p <= 0.05 ~ "*",
     TRUE      ~ ""
   ))

```
### ---- Height for stars/brackets (with padding)----
```{r}
 overall_max <- max(group_avg$mean_pct, na.rm = TRUE)
 y_pad <- overall_max * 0.06
 
 anno_y <- group_avg %>%
   group_by(Genus) %>%
   summarise(ypos = max(mean_pct, na.rm = TRUE) + y_pad, .groups = "drop")
 
 anno_df <- wilcox_res %>%
   left_join(anno_y, by = "Genus") %>%
   filter(stars != "")
 
```
### --- Colores per Group (consistent) ---
```{r}
 pal_group <- c(Rural = " # E9B44C", Urban = "#4F86C6")
 
```
### ---------- Theme ----------
```{r}
 base_theme <- theme_minimal(base_size = 11) +
   theme(
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     axis.line.x.bottom = element_line(),
     axis.line.y.left  = element_line(),
     strip.text.x = element_text(hjust = 0.5, size = 9),
     plot.title = element_text(hjust = 0, size = 14, face = "plain",
                               margin = margin(t = 8, b = 6)),
     plot.margin = margin(15, 12, 12, 12)
   )
 
```
### ---------- Data for BRACKETS drawn by hand ----------
```{r}
 lvl <- levels(group_avg$Genus)
 brackets_manual <- anno_df %>%
   mutate(
     gx  = match(Genus, lvl), # Numerical position of the genus on x
     off = 0.23, # Horizontal separation to each bar (adjust if width/dodge changes)
     x_rural = gx - off,
     x_urban = gx + off,
     y = ypos
   )
```
### --- Gráfico final (barras en paralelo + brackets + estrellas) Final graphic (barrs in parallel + brackets + stars ) ---
```{r}
 p_sep2 <- ggplot(group_avg, aes(x = Genus, y = mean_pct, fill = Group)) +
   geom_col(position = position_dodge(width = 0.8),
            width = 0.7, color = "black", linewidth = 0.2) +
   geom_text(aes(label = sprintf("%.3f%%", mean_pct)),
             position = position_dodge(width = 0.8),
             vjust = -0.3, size = 3, color = "black") +
```
### Horizontal BRACKET (30% higher)
```{r}
   geom_segment(data = brackets_manual,
                aes(x = x_rural, xend = x_urban,
                    y = y * 1.6, yend = y * 1.6),
                inherit.aes = FALSE, linewidth = 0.6) +
```
### Bracket legs (2% below the new level)
```{r}
   geom_segment(data = brackets_manual,
                aes(x = x_rural, xend = x_rural,
                    y = y * 1.6, yend = y * 1.6 - (y * 1.6 * 0.02)),
                inherit.aes = FALSE, linewidth = 0.6) +
   geom_segment(data = brackets_manual,
                aes(x = x_urban, xend = x_urban,
                    y = y * 1.6, yend = y * 1.6 - (y * 1.6 * 0.02)),
                inherit.aes = FALSE, linewidth = 0.6) +
```
### Stars (also in the new level)
```{r}
   geom_text(data = brackets_manual,
             aes(x = gx, y = y * 1.6, label = stars),
             inherit.aes = FALSE, vjust = -0.2, size = 5) +
   scale_fill_manual(values = pal_group, name = "Group") +
   scale_y_continuous(
     labels = label_percent(scale = 1, accuracy = 0.1),
     expand = expansion(mult = c(0, 0.30)) # Air for brackets/stars
   ) +
   labs(
     title = "Mean relative abundance per group of the most prevalent bacterial taxa",
     x = "Genus", y = "Mean per sample"
   ) +
   base_theme
 print(p_sep2)
 
 
 
 p_sep2 <- p_sep2 + theme(legend.position = "none") # Hide the legend from the top part
 
 
```
### Graph both boxplots for EUK and BACT
```{r}
 
 final_fig <- (p_sep2 + labs(title = NULL)) /
   (p_sep  + labs(title = NULL)) +
   plot_layout(guides = "collect") +
   plot_annotation(tag_levels = "A") & # ← Add letters A, B...
   theme(
     legend.position = "right",
     plot.tag = element_text(size = 14, face = "plain")
   )
 
 final_fig
 
```

# Volcano graphs

## Setup


```{r}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(readr); library(purrr); library(tibble)
})

set.seed(1234)

```
### ================== PARAMETERS ==================
```{r}
PSEUDO      <- 1e-8 # to avoid log2(0) and division by 0
ALPHA_Q     <- 0.05 # FDR threshold
MIN_PREV    <- 0.10 # minimum prevalence (≥10% of samples with >0) to keep a genus
LABEL_TOP_N <- 15 # how many labels to put on the volcano (per signal)
OUT_PREFIX  <- "volcano_Rural_vs_Urban"

```
### ============== CHECKS & PREPARATION ==============
```{r}
stopifnot(exists("kaiju_merged"))

dat0 <- kaiju_merged %>% as_tibble()

```
### Ensure minimum types/columns
```{r}
need_cols <- c("file_base","Group","Domain","Genus","reads")
if (!all(need_cols %in% names(dat0))) {
```
### Try to split taxon_name if Domain/Genus/Phylum missing
```{r}
  if ("taxon_name" %in% names(dat0)) {
    dat0 <- dat0 %>%
      tidyr::separate(
        taxon_name,
        into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
                 "Subclass","Order","Family","Genus","Species"),
        sep = ";", fill = "right", extra = "drop", remove = FALSE
      )
  }
}
stopifnot(all(need_cols %in% names(dat0)))

```
### Normalize basic content
```{r}
dat0 <- dat0 %>%
  mutate(
    Group  = case_when(Group %in% c("Rural","Urban") ~ Group, TRUE ~ NA_character_),
    Domain = str_trim(Domain),
    Genus  = str_trim(Genus),
    reads  = suppressWarnings(as.numeric(reads))
  ) %>%
  filter(!is.na(Group), !is.na(Domain), !is.na(Genus), !is.na(reads))

```
### Exclude Homo and empty entries
```{r}
dat0 <- dat0 %>%
  filter(!(Domain == "Eukaryota" & Genus %in% c("Homo","Homo sapiens"))) %>%
  filter(Genus != "")

```
### If Phylum exists, exclude Chordata in Eukaryota
```{r}
if ("Phylum" %in% names(dat0)) {
  dat0 <- dat0 %>%
    mutate(Phylum = str_trim(Phylum)) %>%
    filter(!(Domain == "Eukaryota" & Phylum == "Chordata"))
} else if ("taxon_name" %in% names(dat0)) {
```
### Extract Phylum from taxon_name if not already separated above (for safety)
```{r}
  if (!"Phylum" %in% names(dat0)) {
    dat0 <- dat0 %>%
      tidyr::separate(
        taxon_name,
        into = c("Organism","Domain2","Supergroup","Kingdom","Phylum","Class",
                 "Subclass","Order","Family","Genus2","Species"),
        sep = ";", fill = "right", extra = "drop", remove = FALSE
      ) %>%
      mutate(Phylum = ifelse(is.na(Phylum), "", Phylum)) %>%
      select(-Domain2, -Genus2)
  }
  dat0 <- dat0 %>%
    filter(!(Domain == "Eukaryota" & Phylum == "Chordata"))
}

```
### Keep only Bacteria and Eukaryota
```{r}
dat0 <- dat0 %>% filter(Domain %in% c("Bacteria","Eukaryota"))

```
### =========== Relative abundance within DOMAIN per sample ===========
### Total reads per sample and domain
```{r}
totals_domain <- dat0 %>%
  group_by(file_base, Domain) %>%
  summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")

genus_reads <- dat0 %>%
  group_by(file_base, Group, Domain, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  left_join(totals_domain, by = c("file_base","Domain")) %>%
  mutate(rel_abund = if_else(total_reads_domain > 0, reads / total_reads_domain, 0)) %>%
  ungroup()

```
### =========== Prevalence filter for robustness ===========
### Prevalence = fraction of samples with rel_abund > 0 (in any of the groups)
```{r}
prev_tbl <- genus_reads %>%
  group_by(Domain, Genus) %>%
  summarise(prev = mean(rel_abund > 0), .groups = "drop")

keepers <- prev_tbl %>%
  filter(prev >= MIN_PREV) %>%
  select(Domain, Genus)

genus_reads_f <- genus_reads %>%
  inner_join(keepers, by = c("Domain","Genus"))

```
### =========== Statistical function by domain ===========
```{r}
test_by_domain <- function(df_domain, domain_label) {
```
### Long table per sample
### Wilcoxon by genus: Rural vs Urban on relative abundances (within domain)
```{r}
  test_tbl <- df_domain %>%
    group_by(Genus) %>%
    summarise(
      n_Rural = sum(Group == "Rural"),
      n_Urban = sum(Group == "Urban"),
      mean_Rural = mean(rel_abund[Group == "Rural"], na.rm = TRUE),
      mean_Urban = mean(rel_abund[Group == "Urban"], na.rm = TRUE),
      log2FC = log2((mean_Rural + PSEUDO) / (mean_Urban + PSEUDO)),
      p = tryCatch(
        {
          x <- rel_abund[Group == "Rural"]
          y <- rel_abund[Group == "Urban"]
          if (length(x) >= 2 && length(y) >= 2 && (sd(x) > 0 || sd(y) > 0)) {
            wilcox.test(x, y, exact = FALSE)$p.value
          } else { NA_real_ }
        },
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      q = p.adjust(p, method = "BH"),
      negLog10Q = -log10(pmax(q, .Machine$double.xmin)),
      Domain = domain_label
    ) %>%
    arrange(q, desc(abs(log2FC)))
  
  test_tbl
}

```
### Split by domain and test
```{r}
res_bact <- genus_reads_f %>% filter(Domain == "Bacteria")  %>% test_by_domain("Bacteria")
res_euk  <- genus_reads_f %>% filter(Domain == "Eukaryota") %>% test_by_domain("Eukaryota")

```
### =========== (OPTIONAL) ZicoSeq if available ===========
```{r}
run_zicoseq_safe <- function(df_domain, domain_label) {
  if (!requireNamespace("ZicoSeq", quietly = TRUE)) return(NULL)
```
### Build count matrix per sample × genus and metadata with Group
```{r}
  mat <- df_domain %>%
    select(file_base, Genus, reads, Group) %>%
    group_by(file_base, Genus, Group) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    pivot_wider(names_from = Genus, values_from = reads, values_fill = 0)
  if (nrow(mat) < 4) return(NULL)
  meta <- mat %>% select(file_base, Group)
  counts <- mat %>% select(-file_base, -Group) %>% as.data.frame()
  rownames(counts) <- mat$file_base
```
### Run ZicoSeq (simple Group model)
```{r}
  zres <- tryCatch(
    ZicoSeq::ZicoSeq(
      feature.dat = t(counts),
      meta.dat    = meta,
      grp        = "Group",
      adj.method = "BH",
      n.perm.max = 1000,
      msg        = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(zres)) return(NULL)
  out <- tibble(
    Genus    = rownames(zres$raw.pval),
    p_zico   = as.numeric(zres$raw.pval[,1]),
    q_zico   = as.numeric(zres$adj.pval[,1]),
    stat_zico= as.numeric(zres$stat[,1]),
    Domain   = domain_label
  )
  out
}

z_bact <- run_zicoseq_safe(dat0 %>% filter(Domain=="Bacteria"),  "Bacteria")
z_euk  <- run_zicoseq_safe(dat0 %>% filter(Domain=="Eukaryota"), "Eukaryota")

```
### Merge (if present) for reference; does not change main volcano (Wilcoxon)
```{r}
if (!is.null(z_bact)) {
  res_bact <- res_bact %>%
    left_join(z_bact %>% select(Genus, p_zico, q_zico, stat_zico), by = "Genus")
}
if (!is.null(z_euk)) {
  res_euk <- res_euk %>%
    left_join(z_euk %>% select(Genus, p_zico, q_zico, stat_zico), by = "Genus")
}

```
### =========== Save tables ===========
```{r}
readr::write_tsv(res_bact, paste0(OUT_PREFIX, "_Bacteria_stats.tsv"))
readr::write_tsv(res_euk,  paste0(OUT_PREFIX, "_Eukaryota_stats.tsv"))

```
### =========== Volcano plot helper ===========
```{r}
volcano_plot <- function(df_stats, title_txt, out_png) {
  df_stats <- df_stats %>%
    mutate(sig = ifelse(q <= ALPHA_Q, "FDR ≤ 0.05", "NS"))
  
```
### Choose labels: the most significant + highest |log2FC|
```{r}
  lab_df <- df_stats %>%
    arrange(q, desc(abs(log2FC))) %>%
    slice_head(n = LABEL_TOP_N)
  
  p <- ggplot(df_stats, aes(x = log2FC, y = negLog10Q)) +
    geom_point(aes(shape = sig), alpha = 0.8, size = 2.2) +
    geom_hline(yintercept = -log10(ALPHA_Q), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = Genus),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4, min.segment.length = 0
    ) +
    labs(
      title = title_txt,
      x = "log2 Fold-Change (Rural / Urban)",
      y = expression(-log[10]("FDR (BH)")),
      shape = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
  
  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
  message("Volcano guardado: ", out_png)
  invisible(p)
}

```
### =========== Draw volcano plots ===========
```{r}
volcano_plot(res_bact, "Bacteria — Rural vs Urban", paste0(OUT_PREFIX, "_Bacteria.png"))
volcano_plot(res_euk,  "Eukaryota (sin Chordata) — Rural vs Urban", paste0(OUT_PREFIX, "_Eukaryota.png"))

message("Listo. Archivos generados:\n - ", OUT_PREFIX, "_Bacteria_stats.tsv",
        "\n - ", OUT_PREFIX, "_Eukaryota_stats.tsv",
        "\n - ", OUT_PREFIX, "_Bacteria.png",
        "\n - ", OUT_PREFIX, "_Eukaryota.png")




volcano_plot <- function(df_stats, title_txt, out_png) {
  df_stats <- df_stats %>%
    mutate(sig = ifelse(q <= ALPHA_Q, "FDR ≤ 0.05", "NS"))
  
```
### Top labels
```{r}
  lab_df <- df_stats %>%
    arrange(q, desc(abs(log2FC))) %>%
    slice_head(n = LABEL_TOP_N)
  
  p <- ggplot(df_stats, aes(x = log2FC, y = negLog10Q)) +
    geom_point(aes(shape = sig), alpha = 0.8, size = 2.2, na.rm = TRUE) +
    geom_hline(yintercept = -log10(ALPHA_Q), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = Genus),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4, min.segment.length = 0
    ) +
    labs(
      title = title_txt,
      x = "log2 Fold-Change (Rural / Urban)",
      y = expression(-log[10]("FDR (BH)")),
      shape = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")
  
```
### ATTEMPT 1: ragg (if installed)
```{r}
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(out_png, width = 8, height = 6, units = "in", res = 300)
    print(p)
    dev.off()
  } else {
```
### ATTEMPT 2: base png
```{r}
    png(out_png, width = 8, height = 6, units = "in", res = 300)
    print(p)
    dev.off()
  }
  
  message("Volcano guardado: ", out_png)
  invisible(p)
}


```
## Rural vs Urban
```{r}

rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv); invisible(gc())

```
### --- 1. DATA LOADING AND CLEANING ---
```{r}
kaiju_merged <- read.csv(
  file = "/home/alumno21/axel/files/kaiju_merged_final.csv",
  header = TRUE,
  sep = ",",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

```
### --- Urban and Rural lists adapted with suffix "_kaiju.out" ---
```{r}
urban_files_kaiju <- paste0(c("37082_2 # 1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                              "37082_2 # 14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                              "37082_1 # 26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                              "37082_1 # 13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                              "37035_2 # 23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                              "37082_1 # 10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                              "37035_2 # 3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                              "36703_3 # 2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                              "37082_2 # 13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                              "37082_1 # 16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                              "37082_2 # 11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                              "37035_2 # 1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                              "37035_2 # 11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                              "37082_1 # 21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                              "37035_2 # 20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                              "36703_3 # 10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                              "37035_2 # 26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22"),
                            "_kaiju.out")

rural_files_kaiju <- paste0(c("37082_3 # 17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                              "36703_3 # 26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                              "37082_3 # 7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                              "37035_1 # 18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                              "37035_7 # 7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                              "37082_2 # 16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                              "37082_2 # 19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                              "37035_7 # 5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                              "37035_1 # 10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                              "37082_3 # 19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                              "37035_1 # 15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                              "37082_3 # 24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                              "37082_3 # 26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                              "37082_2 # 20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                              "37082_3 # 3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                              "37035_7 # 22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                              "37035_1 # 3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                              "37082_2 # 21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                              "37082_2 # 26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                              "36703_3 # 28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                              "37082_3 # 28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                              "37082_3 # 23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                              "37082_3 # 14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                              "37035_7 # 17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                              "37035_7 # 18", "37035_7#23", "37035_7#25"),
                            "_kaiju.out")

```
### --- Correctly assign file_base and Group ---
```{r}
kaiju_merged <- kaiju_merged %>%
  mutate(
    file_base = paste0(gsub("\\.out$", "", basename(file)), "_kaiju.out"),
    Group = case_when(
      file_base %in% urban_files_kaiju ~ "Urban",
      file_base %in% rural_files_kaiju ~ "Rural",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

```
### Default device depending on environment
```{r}
if (Sys.getenv("RSTUDIO") == "1") {
  options(device = "RStudioGD")
} else {
  options(bitmapType = "cairo")
}
graphics.off()

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(readr); library(purrr); library(tibble)
  library(ggrepel)
})

set.seed(1234)

```
### ================== PARAMETERS ==================
```{r}
PSEUDO      <- 1e-8 # para evitar log2(0) y divisiones por 0
ALPHA_Q     <- 0.05 # umbral FDR
MIN_PREV    <- 0.20 # prevalencia mínima (≥10% de muestras con >0)
LABEL_TOP_N <- 15 # cuántas etiquetas en el volcán
OUT_PREFIX  <- "volcano_Rural_vs_Urban"

```
### ============== CHECKS & PREPARATION ==============
```{r}
stopifnot(exists("kaiju_merged"))

dat0 <- kaiju_merged %>% as_tibble()
colnames(dat0)
```
### Ensure types/columns; if Domain/Genus/Phylum is missing, try from taxon_name
```{r}
need_cols <- c("file_base","Group","Domain","Genus","reads")
if (!all(need_cols %in% names(dat0))) {
  if ("taxon_name" %in% names(dat0)) {
    dat0 <- dat0 %>%
      tidyr::separate(
        taxon_name,
        into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
                 "Subclass","Order","Family","Genus","Species"),
        sep = ";", fill = "right", extra = "drop", remove = FALSE
      )
  }
}
stopifnot(all(need_cols %in% names(dat0)))

```
### Minimal normalization
```{r}
dat0 <- dat0 %>%
  mutate(
    Group  = case_when(Group %in% c("Rural","Urban") ~ Group, TRUE ~ NA_character_),
    Domain = str_trim(Domain),
    Genus  = str_trim(Genus),
    reads  = suppressWarnings(as.numeric(reads))
  ) %>%
  filter(!is.na(Group), !is.na(Domain), !is.na(Genus), !is.na(reads))

```
### Exclude Homo and empty entries
```{r}
dat0 <- dat0 %>%
  filter(!(Domain == "Eukaryota" & Genus %in% c("Homo","Homo sapiens"))) %>%
  filter(Genus != "")

```
### Exclude Chordata in Eukaryota
```{r}
if ("Phylum" %in% names(dat0)) {
  dat0 <- dat0 %>% mutate(Phylum = str_trim(Phylum)) %>%
    filter(!(Domain == "Eukaryota" & Phylum == "Chordata"))
} else if ("taxon_name" %in% names(dat0)) {
  if (!"Phylum" %in% names(dat0)) {
    dat0 <- dat0 %>%
      tidyr::separate(
        taxon_name,
        into = c("Organism","Domain2","Supergroup","Kingdom","Phylum","Class",
                 "Subclass","Order","Family","Genus2","Species"),
        sep = ";", fill = "right", extra = "drop", remove = FALSE
      ) %>%
      mutate(Phylum = ifelse(is.na(Phylum), "", Phylum)) %>%
      select(-Domain2, -Genus2)
  }
  dat0 <- dat0 %>% filter(!(Domain == "Eukaryota" & Phylum == "Chordata"))
}

```
### Keep only Bacteria and Eukaryota
```{r}
dat0 <- dat0 %>% filter(Domain %in% c("Bacteria","Eukaryota"))

```
### ===== Relative abundance within DOMAIN per sample =====
```{r}
totals_domain <- dat0 %>%
  group_by(file_base, Domain) %>%
  summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")

genus_reads <- dat0 %>%
  group_by(file_base, Group, Domain, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  left_join(totals_domain, by = c("file_base","Domain")) %>%
  mutate(rel_abund = if_else(total_reads_domain > 0, reads / total_reads_domain, 0)) %>%
  ungroup()

```
### ===== Global prevalence filter =====
```{r}
prev_tbl <- genus_reads %>%
  group_by(Domain, Genus) %>%
  summarise(prev = mean(rel_abund > 0), .groups = "drop")

keepers <- prev_tbl %>% filter(prev >= MIN_PREV) %>% select(Domain, Genus)

genus_reads_f <- genus_reads %>% inner_join(keepers, by = c("Domain","Genus"))

```
### ===== Test by domain (Wilcoxon R vs U on rel_abund) =====
```{r}
test_by_domain <- function(df_domain, domain_label) {
  df_domain %>%
    group_by(Genus) %>%
    summarise(
      n_Rural    = sum(Group == "Rural"),
      n_Urban    = sum(Group == "Urban"),
      mean_Rural = mean(rel_abund[Group == "Rural"], na.rm = TRUE),
      mean_Urban = mean(rel_abund[Group == "Urban"], na.rm = TRUE),
      log2FC     = log2((mean_Rural + PSEUDO) / (mean_Urban + PSEUDO)),
      p = tryCatch({
        x <- rel_abund[Group == "Rural"]; y <- rel_abund[Group == "Urban"]
        if (length(x) >= 2 && length(y) >= 2 && (sd(x) > 0 || sd(y) > 0)) {
          wilcox.test(x, y, exact = FALSE)$p.value
        } else NA_real_
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      q         = p.adjust(p, method = "BH"),
      negLog10Q = -log10(pmax(q, .Machine$double.xmin)),
      Domain    = domain_label
    ) %>%
    arrange(q, desc(abs(log2FC)))
}

res_bact <- genus_reads_f %>% filter(Domain == "Bacteria")  %>% test_by_domain("Bacteria")
res_euk  <- genus_reads_f %>% filter(Domain == "Eukaryota") %>% test_by_domain("Eukaryota")

```
### ===== (Optional) ZicoSeq if available: add Zico p/q columns =====
```{r}
run_zicoseq_safe <- function(df_all, domain_label) {
  if (!requireNamespace("ZicoSeq", quietly = TRUE)) return(NULL)
  mat <- df_all %>%
    select(file_base, Genus, reads, Group) %>%
    group_by(file_base, Genus, Group) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    pivot_wider(names_from = Genus, values_from = reads, values_fill = 0)
  if (nrow(mat) < 4) return(NULL)
  meta   <- mat %>% select(file_base, Group)
  counts <- mat %>% select(-file_base, -Group) %>% as.data.frame()
  rownames(counts) <- mat$file_base
  
  zres <- tryCatch(ZicoSeq::ZicoSeq(
    feature.dat = t(counts), meta.dat = meta, grp = "Group",
    adj.method = "BH", n.perm.max = 1000, msg = FALSE
  ), error = function(e) NULL)
  if (is.null(zres)) return(NULL)
  
  tibble(
    Genus     = rownames(zres$raw.pval),
    p_zico    = as.numeric(zres$raw.pval[,1]),
    q_zico    = as.numeric(zres$adj.pval[,1]),
    stat_zico = as.numeric(zres$stat[,1]),
    Domain    = domain_label
  )
}

z_bact <- run_zicoseq_safe(dat0 %>% filter(Domain=="Bacteria"),  "Bacteria")
z_euk  <- run_zicoseq_safe(dat0 %>% filter(Domain=="Eukaryota"), "Eukaryota")

if (!is.null(z_bact)) res_bact <- res_bact %>% left_join(z_bact %>% select(Genus, p_zico, q_zico, stat_zico), by = "Genus")
if (!is.null(z_euk))  res_euk  <- res_euk  %>% left_join(z_euk  %>% select(Genus, p_zico, q_zico, stat_zico), by = "Genus")

```
### ===== Save tables =====
```{r}
readr::write_tsv(res_bact, paste0(OUT_PREFIX, "_Bacteria_stats.tsv"))
readr::write_tsv(res_euk,  paste0(OUT_PREFIX, "_Eukaryota_stats.tsv"))

volcano_plot <- function(df_stats, title_txt, out_png) {
  df_stats <- df_stats %>%
    mutate(
      sig = case_when(
        is.na(q) ~ "NA",
        q <= ALPHA_Q ~ "FDR ≤ 0.05",
        TRUE ~ "NS"
      ),
      group_color = ifelse(log2FC >= 0, "Rural", "Urban")
    )
  
  lab_df <- df_stats %>%
    arrange(q, desc(abs(log2FC))) %>%
    slice_head(n = LABEL_TOP_N)
  
  p <- ggplot(df_stats, aes(x = log2FC, y = negLog10Q)) +
    geom_point(aes(shape = sig, color = group_color), alpha = 0.6, size = 2.5, na.rm = TRUE) +
    scale_color_manual(values = c("Rural" = " # E69F00", "Urban" = "#0072B2")) +
    geom_hline(yintercept = -log10(ALPHA_Q), linetype = "dashed", size = 0.1) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.1) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = Genus, color = group_color),
      size = 3.5, max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0, show.legend = FALSE
    ) +
    labs(
      title = title_txt,
      x = "log2 Fold-Change (Rural / Urban)",
      y = expression(-log[10]("FDR (BH)")),
      shape = NULL,
      color = "Higher in"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")
  
```
### Show in the Plots panel
```{r}
  print(p)
  
```
### Save robustly with cairo
```{r}
  tp <- if (.Platform$OS.type == "windows") "cairo-png" else "cairo"
  png(out_png, width = 8, height = 6, units = "in", res = 300, type = tp)
  print(p); dev.off()
  
  message("Volcano guardado: ", out_png)
  invisible(p)
}

```
### Calls
```{r}
volcano_plot(res_bact, "", paste0(OUT_PREFIX, "_Bacteria.png"))
volcano_plot(res_euk,  "", paste0(OUT_PREFIX, "_Eukaryota.png"))

---------------------------------
```
<a id="bmi"></a>
## BMI<25 vs BMI≥25
### Default device depending on environment
```{r}
if (Sys.getenv("RSTUDIO") == "1") {
  options(device = "RStudioGD")
} else {
  options(bitmapType = "cairo")
}
graphics.off()

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(readr); library(purrr); library(tibble)
  library(ggrepel)
})

set.seed(1234)

```
### ================== PARAMETERS ==================
```{r}
PSEUDO      <- 1e-8 # to avoid log2(0) and divisions by 0
ALPHA_Q     <- 0.05 # FDR threshold
MIN_PREV    <- 0.20 # minimum prevalence (≥10% of samples with >0)
LABEL_TOP_N <- 15 # how many labels on the volcano plot
OUT_PREFIX  <- "volcano_BMI_groups"

```
### ============== CHECKS & PREPARATION ==============
```{r}
stopifnot(exists("kaiju_merged"))

dat0 <- kaiju_merged %>% as_tibble()

```
### Ensure types/columns; if Domain/Genus/Phylum is missing, try from taxon_name
```{r}
need_cols <- c("file_base","BMI_group","Domain","Genus","reads")
if (!all(need_cols %in% names(dat0))) {
  if ("taxon_name" %in% names(dat0)) {
    dat0 <- dat0 %>%
      tidyr::separate(
        taxon_name,
        into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
                 "Subclass","Order","Family","Genus","Species"),
        sep = ";", fill = "right", extra = "drop", remove = FALSE
      )
  }
}
stopifnot(all(need_cols %in% names(dat0)))

```
### Minimal normalization
```{r}
dat0 <- dat0 %>%
  mutate(
    BMI_group = case_when(
      BMI_group %in% c("BMI<25","BMI≥25") ~ BMI_group,
      TRUE ~ NA_character_
    ),
    Domain = str_trim(Domain),
    Genus  = str_trim(Genus),
    reads  = suppressWarnings(as.numeric(reads))
  ) %>%
  filter(!is.na(BMI_group), !is.na(Domain), !is.na(Genus), !is.na(reads))

```
### Exclude Homo and empty entries
```{r}
dat0 <- dat0 %>%
  filter(!(Domain == "Eukaryota" & Genus %in% c("Homo","Homo sapiens"))) %>%
  filter(Genus != "")

```
### Exclude Chordata in Eukaryota
```{r}
if ("Phylum" %in% names(dat0)) {
  dat0 <- dat0 %>% mutate(Phylum = str_trim(Phylum)) %>%
    filter(!(Domain == "Eukaryota" & Phylum == "Chordata"))
} else if ("taxon_name" %in% names(dat0)) {
  if (!"Phylum" %in% names(dat0)) {
    dat0 <- dat0 %>%
      tidyr::separate(
        taxon_name,
        into = c("Organism","Domain2","Supergroup","Kingdom","Phylum","Class",
                 "Subclass","Order","Family","Genus2","Species"),
        sep = ";", fill = "right", extra = "drop", remove = FALSE
      ) %>%
      mutate(Phylum = ifelse(is.na(Phylum), "", Phylum)) %>%
      select(-Domain2, -Genus2)
  }
  dat0 <- dat0 %>% filter(!(Domain == "Eukaryota" & Phylum == "Chordata"))
}

```
### Keep only Bacteria and Eukaryota
```{r}
dat0 <- dat0 %>% filter(Domain %in% c("Bacteria","Eukaryota"))

```
### ===== Relative abundance within DOMAIN per sample =====
```{r}
totals_domain <- dat0 %>%
  group_by(file_base, Domain) %>%
  summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")

genus_reads <- dat0 %>%
  group_by(file_base, BMI_group, Domain, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  left_join(totals_domain, by = c("file_base","Domain")) %>%
  mutate(rel_abund = if_else(total_reads_domain > 0, reads / total_reads_domain, 0)) %>%
  ungroup()

```
### ===== Global prevalence filter =====
```{r}
prev_tbl <- genus_reads %>%
  group_by(Domain, Genus) %>%
  summarise(prev = mean(rel_abund > 0), .groups = "drop")

keepers <- prev_tbl %>% filter(prev >= MIN_PREV) %>% select(Domain, Genus)

genus_reads_f <- genus_reads %>% inner_join(keepers, by = c("Domain","Genus"))

```
### ===== Test by domain (Wilcoxon BMI<25 vs BMI≥25 on rel_abund) =====
```{r}
test_by_domain <- function(df_domain, domain_label) {
  df_domain %>%
    group_by(Genus) %>%
    summarise(
      n_low     = sum(BMI_group == "BMI<25"),
      n_high    = sum(BMI_group == "BMI≥25"),
      mean_low  = mean(rel_abund[BMI_group == "BMI<25"],  na.rm = TRUE),
      mean_high = mean(rel_abund[BMI_group == "BMI≥25"], na.rm = TRUE),
      log2FC    = log2((mean_low + PSEUDO) / (mean_high + PSEUDO)), # BMI<25 / BMI≥25
      p = tryCatch({
        x <- rel_abund[BMI_group == "BMI<25"]
        y <- rel_abund[BMI_group == "BMI≥25"]
        if (length(x) >= 2 && length(y) >= 2 && (sd(x) > 0 || sd(y) > 0)) {
          wilcox.test(x, y, exact = FALSE)$p.value
        } else NA_real_
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      q         = p.adjust(p, method = "BH"),
      negLog10Q = -log10(pmax(q, .Machine$double.xmin)),
      Domain    = domain_label
    ) %>%
    arrange(q, desc(abs(log2FC)))
}

res_bact <- genus_reads_f %>% filter(Domain == "Bacteria")  %>% test_by_domain("Bacteria")
res_euk  <- genus_reads_f %>% filter(Domain == "Eukaryota") %>% test_by_domain("Eukaryota")

```
### ===== (Optional) ZicoSeq if available: add Zico p/q columns =====
```{r}
run_zicoseq_safe <- function(df_all, domain_label) {
  if (!requireNamespace("ZicoSeq", quietly = TRUE)) return(NULL)
  mat <- df_all %>%
    select(file_base, Genus, reads, BMI_group) %>%
    group_by(file_base, Genus, BMI_group) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    pivot_wider(names_from = Genus, values_from = reads, values_fill = 0)
  if (nrow(mat) < 4) return(NULL)
  meta   <- mat %>% select(file_base, BMI_group)
  counts <- mat %>% select(-file_base, -BMI_group) %>% as.data.frame()
  rownames(counts) <- mat$file_base
  
  zres <- tryCatch(ZicoSeq::ZicoSeq(
    feature.dat = t(counts), meta.dat = meta, grp = "BMI_group",
    adj.method = "BH", n.perm.max = 1000, msg = FALSE
  ), error = function(e) NULL)
  if (is.null(zres)) return(NULL)
  
  tibble(
    Genus     = rownames(zres$raw.pval),
    p_zico    = as.numeric(zres$raw.pval[,1]),
    q_zico    = as.numeric(zres$adj.pval[,1]),
    stat_zico = as.numeric(zres$stat[,1]),
    Domain    = domain_label
  )
}

z_bact <- run_zicoseq_safe(dat0 %>% filter(Domain=="Bacteria"),  "Bacteria")
z_euk  <- run_zicoseq_safe(dat0 %>% filter(Domain=="Eukaryota"), "Eukaryota")

if (!is.null(z_bact)) res_bact <- res_bact %>% left_join(z_bact %>% select(Genus, p_zico, q_zico, stat_zico), by = "Genus")
if (!is.null(z_euk))  res_euk  <- res_euk  %>% left_join(z_euk  %>% select(Genus, p_zico, q_zico, stat_zico), by = "Genus")

```
### ===== Save tables =====
```{r}
readr::write_tsv(res_bact, paste0(OUT_PREFIX, "_Bacteria_stats.tsv"))
readr::write_tsv(res_euk,  paste0(OUT_PREFIX, "_Eukaryota_stats.tsv"))

```
### ===== Volcano plot (color by group with higher abundance) =====
```{r}
volcano_plot <- function(df_stats, title_txt, out_png) {
  df_stats <- df_stats %>%
    mutate(
      sig = case_when(
        is.na(q) ~ "NA",
        q <= ALPHA_Q ~ "FDR ≤ 0.05",
        TRUE ~ "NS"
      ),
      group_color = ifelse(log2FC >= 0, "BMI<25", "BMI≥25") # log2FC = low/high
    )
  
  lab_df <- df_stats %>% arrange(q, desc(abs(log2FC))) %>% slice_head(n = LABEL_TOP_N)
  
  p <- ggplot(df_stats, aes(x = log2FC, y = negLog10Q)) +
    geom_point(aes(shape = sig, color = group_color), alpha = 0.6, size = 2.5, na.rm = TRUE) +
    scale_color_manual(values = c("BMI<25" = "lightblue", "BMI≥25" = "coral")) + # blue <25, orange ≥25
    geom_hline(yintercept = -log10(ALPHA_Q), linetype = "dashed", size = 0.1) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.1) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = Genus, color = group_color),
      size = 3.5, max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0, show.legend = FALSE
    ) +
    labs(
      title = title_txt,
      x = "log2 Fold-Change (BMI<25 / BMI≥25)",
      y = expression(-log[10]("FDR (BH)")),
      shape = NULL,
      color = "Higher in"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")
  
```
### Show in the Plots panel
```{r}
  print(p)
  
```
### Save robustly with cairo
```{r}
  tp <- if (.Platform$OS.type == "windows") "cairo-png" else "cairo"
  png(out_png, width = 8, height = 6, units = "in", res = 300, type = tp)
  print(p); dev.off()
  
  message("Volcano guardado: ", out_png)
  invisible(p)
}

```
### ===== Draw volcano plots =====
```{r}
volcano_plot(res_bact, "", paste0(OUT_PREFIX, "_Bacteria.png"))
volcano_plot(res_euk,  "", paste0(OUT_PREFIX, "_Eukaryota.png"))

```
### ===== Summary =====
```{r}
message("\nResumen Bacteria: n=", nrow(res_bact),
        " | q<=0.05: ", sum(res_bact$q <= ALPHA_Q, na.rm = TRUE))
message("Resumen Eukaryota: n=", nrow(res_euk),
        " | q<=0.05: ", sum(res_euk$q <= ALPHA_Q, na.rm = TRUE))
message("\nListo. Archivos generados:\n - ", OUT_PREFIX, "_Bacteria_stats.tsv",
        "\n - ", OUT_PREFIX, "_Eukaryota_stats.tsv",
        "\n - ", OUT_PREFIX, "_Bacteria.png",
        "\n - ", OUT_PREFIX, "_Eukaryota.png")
```

<a id="bmiru"></a>
##  BMI≥25 Rural vs BMI≥25 Urban 
### Default device depending on environment
```{r}
if (Sys.getenv("RSTUDIO") == "1") {
  options(device = "RStudioGD")
} else {
  options(bitmapType = "cairo")
}
graphics.off()

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(readr); library(purrr); library(tibble)
  library(ggrepel)
})

set.seed(1234)

```
### ================== PARAMETERS ==================
```{r}
PSEUDO      <- 1e-8 # to avoid log2(0) and divisions by 0
ALPHA_Q     <- 0.05 # FDR threshold
MIN_PREV    <- 0.2 # minimum prevalence (≥10% of samples with >0)
LABEL_TOP_N <- 8 # how many labels on the volcano plot
OUT_PREFIX  <- "volcano_BMI26_Rural_vs_Urban"

```
### ============== CHECKS & PREPARATION ==============
```{r}
stopifnot(exists("kaiju_merged"))

dat0 <- kaiju_merged %>% as_tibble()

```
### Ensure types/columns; if Domain/Genus/Phylum is missing, try from taxon_name
```{r}
need_cols <- c("file_base","context_group","Domain","Genus","reads")
if (!all(need_cols %in% names(dat0))) {
  if ("taxon_name" %in% names(dat0)) {
    dat0 <- dat0 %>%
      tidyr::separate(
        taxon_name,
        into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
                 "Subclass","Order","Family","Genus","Species"),
        sep = ";", fill = "right", extra = "drop", remove = FALSE
      )
  }
}
stopifnot(all(need_cols %in% names(dat0)))

```
### Minimal normalization
```{r}
dat0 <- dat0 %>%
  mutate(
    context_group = case_when(
      context_group %in% c("BMI≥25 Rural","BMI≥25 Urban") ~ context_group,
      TRUE ~ NA_character_
    ),
    Domain = str_trim(Domain),
    Genus  = str_trim(Genus),
    reads  = suppressWarnings(as.numeric(reads))
  ) %>%
  filter(!is.na(context_group), !is.na(Domain), !is.na(Genus), !is.na(reads))

```
### Exclude Homo and empty entries
```{r}
dat0 <- dat0 %>%
  filter(!(Domain == "Eukaryota" & Genus %in% c("Homo","Homo sapiens"))) %>%
  filter(Genus != "")

```
### Exclude Chordata in Eukaryota
```{r}
if ("Phylum" %in% names(dat0)) {
  dat0 <- dat0 %>% mutate(Phylum = str_trim(Phylum)) %>%
    filter(!(Domain == "Eukaryota" & Phylum == "Chordata"))
} else if ("taxon_name" %in% names(dat0)) {
  if (!"Phylum" %in% names(dat0)) {
    dat0 <- dat0 %>%
      tidyr::separate(
        taxon_name,
        into = c("Organism","Domain2","Supergroup","Kingdom","Phylum","Class",
                 "Subclass","Order","Family","Genus2","Species"),
        sep = ";", fill = "right", extra = "drop", remove = FALSE
      ) %>%
      mutate(Phylum = ifelse(is.na(Phylum), "", Phylum)) %>%
      select(-Domain2, -Genus2)
  }
  dat0 <- dat0 %>% filter(!(Domain == "Eukaryota" & Phylum == "Chordata"))
}

```
### Keep only Bacteria and Eukaryota
```{r}
dat0 <- dat0 %>% filter(Domain %in% c("Bacteria","Eukaryota"))

```
### ===== Relative abundance within DOMAIN per sample =====
```{r}
totals_domain <- dat0 %>%
  group_by(file_base, Domain) %>%
  summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")

genus_reads <- dat0 %>%
  group_by(file_base, context_group, Domain, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  left_join(totals_domain, by = c("file_base","Domain")) %>%
  mutate(rel_abund = if_else(total_reads_domain > 0, reads / total_reads_domain, 0)) %>%
  ungroup()

```
### ===== Global prevalence filter =====
```{r}
prev_tbl <- genus_reads %>%
  group_by(Domain, Genus) %>%
  summarise(prev = mean(rel_abund > 0), .groups = "drop")

keepers <- prev_tbl %>% filter(prev >= MIN_PREV) %>% select(Domain, Genus)

genus_reads_f <- genus_reads %>% inner_join(keepers, by = c("Domain","Genus"))

```
### ===== Test by domain (Wilcoxon BMI≥25 Rural vs Urban) =====
```{r}
test_by_domain <- function(df_domain, domain_label) {
  df_domain %>%
    group_by(Genus) %>%
    summarise(
      n_rural = sum(context_group == "BMI≥25 Rural"),
      n_urban = sum(context_group == "BMI≥25 Urban"),
      mean_rural = mean(rel_abund[context_group == "BMI≥25 Rural"], na.rm = TRUE),
      mean_urban = mean(rel_abund[context_group == "BMI≥25 Urban"], na.rm = TRUE),
      log2FC     = log2((mean_rural + PSEUDO) / (mean_urban + PSEUDO)), # Rural / Urban
      p = tryCatch({
        x <- rel_abund[context_group == "BMI≥25 Rural"]
        y <- rel_abund[context_group == "BMI≥25 Urban"]
        if (length(x) >= 2 && length(y) >= 2 && (sd(x) > 0 || sd(y) > 0)) {
          wilcox.test(x, y, exact = FALSE)$p.value
        } else NA_real_
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      q         = p.adjust(p, method = "BH"),
      negLog10Q = -log10(pmax(q, .Machine$double.xmin)),
      Domain    = domain_label
    ) %>%
    arrange(q, desc(abs(log2FC)))
}

res_bact <- genus_reads_f %>% filter(Domain == "Bacteria")  %>% test_by_domain("Bacteria")
res_euk  <- genus_reads_f %>% filter(Domain == "Eukaryota") %>% test_by_domain("Eukaryota")

```
### ===== (Optional) ZicoSeq if available: add Zico p/q columns =====
```{r}
run_zicoseq_safe <- function(df_all, domain_label) {
  if (!requireNamespace("ZicoSeq", quietly = TRUE)) return(NULL)
  mat <- df_all %>%
    select(file_base, Genus, reads, context_group) %>%
    group_by(file_base, Genus, context_group) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    pivot_wider(names_from = Genus, values_from = reads, values_fill = 0)
  if (nrow(mat) < 4) return(NULL)
  meta   <- mat %>% select(file_base, context_group)
  counts <- mat %>% select(-file_base, -context_group) %>% as.data.frame()
  rownames(counts) <- mat$file_base
  
  zres <- tryCatch(ZicoSeq::ZicoSeq(
    feature.dat = t(counts), meta.dat = meta, grp = "context_group",
    adj.method = "BH", n.perm.max = 1000, msg = FALSE
  ), error = function(e) NULL)
  if (is.null(zres)) return(NULL)
  
  tibble(
    Genus     = rownames(zres$raw.pval),
    p_zico    = as.numeric(zres$raw.pval[,1]),
    q_zico    = as.numeric(zres$adj.pval[,1]),
    stat_zico = as.numeric(zres$stat[,1]),
    Domain    = domain_label
  )
}

z_bact <- run_zicoseq_safe(dat0 %>% filter(Domain=="Bacteria"),  "Bacteria")
z_euk  <- run_zicoseq_safe(dat0 %>% filter(Domain=="Eukaryota"), "Eukaryota")

if (!is.null(z_bact)) res_bact <- res_bact %>% left_join(z_bact %>% select(Genus, p_zico, q_zico, stat_zico), by = "Genus")
if (!is.null(z_euk))  res_euk  <- res_euk  %>% left_join(z_euk  %>% select(Genus, p_zico, q_zico, stat_zico), by = "Genus")

```
### ===== Save tables =====
```{r}
readr::write_tsv(res_bact, paste0(OUT_PREFIX, "_Bacteria_stats.tsv"))
readr::write_tsv(res_euk,  paste0(OUT_PREFIX, "_Eukaryota_stats.tsv"))

```
### ===== Volcano plot (color by group with higher abundance) =====
```{r}
volcano_plot <- function(df_stats, title_txt, out_png) {
  df_stats <- df_stats %>%
    mutate(
      sig = case_when(
        is.na(q) ~ "NA",
        q <= ALPHA_Q ~ "FDR ≤ 0.05",
        TRUE ~ "NS"
      ),
      group_color = ifelse(log2FC >= 0, "BMI≥25 Rural", "BMI≥25 Urban") # Rural / Urban
    )
  
  lab_df <- df_stats %>% arrange(q, desc(abs(log2FC))) %>% slice_head(n = LABEL_TOP_N)
  
  p <- ggplot(df_stats, aes(x = log2FC, y = negLog10Q)) +
    geom_point(aes(shape = sig, color = group_color), alpha = 0.6, size = 2.5, na.rm = TRUE) +
    scale_color_manual(values = c("BMI≥25 Rural" = "darkgreen", "BMI≥25 Urban" = "darkred")) +  
    geom_hline(yintercept = -log10(ALPHA_Q), linetype = "dashed", size = 0.1) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.1) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = Genus, color = group_color),
      size = 3.5, max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0, show.legend = FALSE
    ) +
    labs(
      title = title_txt,
      x = "log2 Fold-Change (BMI≥25 Rural / BMI≥25 Urban)",
      y = expression(-log[10]("FDR (BH)")),
      shape = NULL,
      color = "Higher in"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")
  
```
### Show in the Plots panel
```{r}
  print(p)
  
```
### Save robustly with cairo
```{r}
  tp <- if (.Platform$OS.type == "windows") "cairo-png" else "cairo"
  png(out_png, width = 8, height = 6, units = "in", res = 300, type = tp)
  print(p); dev.off()
  
  message("Volcano guardado: ", out_png)
  invisible(p)
}

```
### ===== Draw volcano plots =====
```{r}
volcano_plot(res_bact, "", paste0(OUT_PREFIX, "_Bacteria.png"))
volcano_plot(res_euk,  "", paste0(OUT_PREFIX, "_Eukaryota.png"))

```
### ===== Summary =====
```{r}
message("\nResumen Bacteria: n=", nrow(res_bact),
        " | q<=0.05: ", sum(res_bact$q <= ALPHA_Q, na.rm = TRUE))
message("Resumen Eukaryota: n=", nrow(res_euk),
        " | q<=0.05: ", sum(res_euk$q <= ALPHA_Q, na.rm = TRUE))
message("\nListo. Archivos generados:\n - ", OUT_PREFIX, "_Bacteria_stats.tsv",
        "\n - ", OUT_PREFIX, "_Eukaryota_stats.tsv",
        "\n - ", OUT_PREFIX, "_Bacteria.png",
        "\n - ", OUT_PREFIX, "_Eukaryota.png")

```

## BUSCO completeness assessment of MAGs

#### For the comprehension of the specific demographic signatures, we recovered 713 MAGs in total by focusing on the 10 highest-read individuals for each of the ten differential lineages identified previously, five bacterial and five eukaryotic taxa (Fig. 10). Rural and Urban cohorts contributed roughly equal shares of these genomes (Fig. 10A–B). Across MAGs, single-copy BUSCOs dominated the profiles, whereas duplicated, fragmented, and missing categories formed minor fractions, yielding near-identical distributions between lifestyles (Fig. 10A vs. 10B).

### Summary figure of MAGs
### We need the BUSCO file "tabla_busco.csv"
```{r}
buscos <- read.csv("/home/alumno21/axel/files/tabla_busco.csv", sep = ",", stringsAsFactors = FALSE)
colnames(buscos)

library(dplyr)
library(tidyr)
library(ggplot2)

urban_files <- c("37082_2 # 1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                 "37082_2 # 14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                 "37082_1 # 26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                 "37082_1 # 13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                 "37035_2 # 23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                 "37082_1 # 10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                 "37035_2 # 3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                 "36703_3 # 2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                 "37082_2 # 13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                 "37082_1 # 16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                 "37082_2 # 11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                 "37035_2 # 1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                 "37035_2 # 11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                 "37082_1 # 21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                 "37035_2 # 20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                 "36703_3 # 10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                 "37035_2 # 26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22")

rural_files <- c("37082_3 # 17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                 "36703_3 # 26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                 "37082_3 # 7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                 "37035_1 # 18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                 "37035_7 # 7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                 "37082_2 # 16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                 "37082_2 # 19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                 "37035_7 # 5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                 "37035_1 # 10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                 "37082_3 # 19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                 "37035_1 # 15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                 "37082_3 # 24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                 "37082_3 # 26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                 "37082_2 # 20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                 "37082_3 # 3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                 "37035_7 # 22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                 "37035_1 # 3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                 "37082_2 # 21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                 "37082_2 # 26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                 "36703_3 # 28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                 "37082_3 # 28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                 "37082_3 # 23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                 "37082_3 # 14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                 "37035_7 # 17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                 "37035_7 # 18", "37035_7#23", "37035_7#25")

```
### Assign a unique name to the MAG
```{r}
buscos$id <- paste0(buscos$File, "_", buscos$MAG)

```
### Ensure that S, D, F, M are numeric
```{r}
buscos <- buscos %>%
  mutate(across(c(S, D, F, M), as.numeric))

```
### Assign Rural or Urban group according to external vectors
```{r}
buscos$Group <- ifelse(buscos$File %in% rural_files, "Rural",
                       ifelse(buscos$File %in% urban_files, "Urban", NA))

```
### Transform to long format
```{r}
buscos_long <- buscos %>%
  select(id, Group, S, D, F, M) %>%
  pivot_longer(cols = c(S, D, F, M),
               names_to = "Status",
               values_to = "Percent")

```
### More readable labels
```{r}
buscos_long$Status <- recode(buscos_long$Status,
                             "S" = "Single-copy",
                             "D" = "Duplicated",
                             "F" = "Fragmented",
                             "M" = "Missing")

```
### Create a consecutive label ID per group
```{r}
buscos_long <- buscos_long %>%
  group_by(Group) %>%
  mutate(label_id = dense_rank(id)) %>%
  ungroup()

```
### Improved color palette
```{r}
busco_colors <- c(
  "Single-copy" = " # 1f77b4",
  "Duplicated" = " # ff7f0e",
  "Fragmented" = " # 2ca02c",
  "Missing" = " # d62728"
)

```
### --- Individual plots ---
```{r}

plot_rural <- ggplot(buscos_long %>% filter(Group == "Rural"),
                     aes(x = label_id, y = Percent, fill = Status)) +
  geom_bar(stat = "identity", width = 0.8, alpha=0.85) +
  scale_fill_manual(values = busco_colors) +
  labs(
    title = "",
    x = "MAGs (consecutive IDs)",
    y = "Percentage",
    fill = "BUSCO Category"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 12, face = "plain", hjust = 0.5)
  )

plot_urban <- ggplot(buscos_long %>% filter(Group == "Urban"),
                     aes(x = label_id, y = Percent, fill = Status)) +
  geom_bar(stat = "identity", width = 0.8, alpha=0.85) +
  scale_fill_manual(values = busco_colors) +
  labs(
    title = "",
    x = "MAGs (consecutive IDs)",
    y = "Percentage",
    fill = "BUSCO Category"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "plain", hjust = 0.5)
  )



```
### --- Combine plots with patchwork ---
```{r}


combined_busco_plot <- plot_rural + plot_urban +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "",
    tag_levels = "A", # Letras automáticas: A, B...
    tag_prefix = "", # Sin prefijos, solo letras
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20), # márgenes en todos lados
      plot.tag = element_text(size = 16, face = "bold"),
      plot.tag.position = c(0.1, 0.85)
    )
  )
```
### Final BUSCO table wit colors
```{r}
combined_busco_plot
```

<a id="fun-anno"></a>
## Functional annotation profiles of high-quality MAGs
#### (>95% completeness) from Rural and Urban groups.

## Part A
### We need the functional annotations file "all_annotations_trimmed.csv"

### --- Libraries ---
```{r}
library(dplyr)
library(ggplot2)
library(scales) # Para percent_format()

```
### --- Load data ---
```{r}
anot <- read.csv("/home/alumno21/axel/files/all_annotations_trimmed.csv",
                 sep = ",", stringsAsFactors = FALSE)

```
### --- Define files by group ---
```{r}
urban_files <- c("37082_2 # 1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                 "37082_2 # 14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                 "37082_1 # 26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                 "37082_1 # 13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                 "37035_2 # 23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                 "37082_1 # 10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                 "37035_2 # 3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                 "36703_3 # 2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                 "37082_2 # 13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                 "37082_1 # 16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                 "37082_2 # 11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                 "37035_2 # 1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                 "37035_2 # 11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                 "37082_1 # 21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                 "37035_2 # 20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                 "36703_3 # 10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                 "37035_2 # 26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22")
rural_files <- c("37082_3 # 17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                 "36703_3 # 26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                 "37082_3 # 7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                 "37035_1 # 18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                 "37035_7 # 7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                 "37082_2 # 16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                 "37082_2 # 19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                 "37035_7 # 5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                 "37035_1 # 10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                 "37082_3 # 19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                 "37035_1 # 15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                 "37082_3 # 24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                 "37082_3 # 26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                 "37082_2 # 20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                 "37082_3 # 3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                 "37035_7 # 22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                 "37035_1 # 3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                 "37082_2 # 21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                 "37082_2 # 26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                 "36703_3 # 28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                 "37082_3 # 28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                 "37082_3 # 23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                 "37082_3 # 14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                 "37035_7 # 17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                 "37035_7 # 18", "37035_7#23", "37035_7#25")

```
### --- Function to classify simplified domain ---
```{r}
assign_domain <- function(term) {
  term_lower <- tolower(term)
  if (grepl("virus|myoviridae|caudovirales|podoviridae|siphoviridae|ssdna|dsdna", term_lower)) {
    return("Viral")
  } else if (grepl("bacteria", term_lower)) {
    return("Bacteria")
  } else if (grepl("archaea", term_lower)) {
    return("Archaea")
  } else if (grepl("eukaryota", term_lower)) {
    return("Eukaryota")
  } else if (term_lower == "root") {
    return("Other")
  } else {
    return("Other")
  }
}

```
### --- Step 1: Prepare data ---
```{r}
anot <- anot %>%
  filter(grepl("\\|", max_annot_lvl)) %>%
  mutate(
    max_annot_lvl_clean = sub(".*\\|", "", max_annot_lvl),
    Domain_simplified = sapply(max_annot_lvl_clean, assign_domain),
    Grupo = case_when(
      File %in% rural_files ~ "Rural",
      File %in% urban_files ~ "Urban",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Grupo))

```
### --- Step 2: Relative abundance per MAG ---
```{r}
score_total_por_file <- anot %>%
  group_by(File) %>%
  summarise(score_total = sum(score, na.rm = TRUE), .groups = "drop")

score_por_file_nivel <- anot %>%
  group_by(File, max_annot_lvl_clean, Domain_simplified) %>%
  summarise(score_nivel = sum(score, na.rm = TRUE), .groups = "drop")

abundancia_relativa <- score_por_file_nivel %>%
  left_join(score_total_por_file, by = "File") %>%
  mutate(abund_rel = score_nivel / score_total) %>%
  left_join(anot %>% select(File, Grupo) %>% distinct(), by = "File") %>%
  filter(!is.na(abund_rel))

```
### --- Step 3: Average by group and domain ---
```{r}
abundancia_promedio <- abundancia_relativa %>%
  group_by(Grupo, Domain_simplified) %>%
  summarise(promedio_abund = mean(abund_rel), .groups = "drop")

```
### --- Step 4: Count unique MAGs per group ---
```{r}
n_mags_por_grupo <- anot %>%
  select(Grupo, MAG) %>%
  distinct() %>%
  group_by(Grupo) %>%
  summarise(n_MAGs = n(), .groups = "drop")
```
### --- Extra Step 4: Count annotated genes/ORFs per group ---
```{r}
n_genes_por_grupo <- anot %>%
  group_by(Grupo) %>%
  summarise(n_genes = n(), .groups = "drop")


```
### --- Step 5: Count unique microorganisms per group ---
```{r}
n_microbes_por_grupo <- anot %>%
  select(Grupo, max_annot_lvl_clean) %>%
  distinct() %>%
  group_by(Grupo) %>%
  summarise(n_microbes = n(), .groups = "drop")

```
### --- Required libraries ---
```{r}
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)
library(RColorBrewer)


```
### --- Step 6: Text of annotated genes/ORFs per group ---
```{r}
texto_genes <- paste0(
  "Annotated genes (Rural: ", n_genes_por_grupo$n_genes[n_genes_por_grupo$Grupo == "Rural"],
  " | Urban: ", n_genes_por_grupo$n_genes[n_genes_por_grupo$Grupo == "Urban"],")"
)

```
### --- Step 7: Percentages in legend ---
```{r}
porcentaje_por_dominio <- abundancia_promedio %>%
  group_by(Domain_simplified) %>%
  summarise(total_abund = sum(promedio_abund), .groups = "drop") %>%
  mutate(label = paste0(Domain_simplified, " (", percent(total_abund / sum(total_abund)), ")"))

labels_legend <- porcentaje_por_dominio$label
names(labels_legend) <- porcentaje_por_dominio$Domain_simplified

```
### --- Step 8: Bar plot ---
```{r}
plot_abundancia <- ggplot(abundancia_promedio, aes(x = Grupo, y = promedio_abund, fill = Domain_simplified)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(
    values = brewer.pal(n = length(labels_legend), name = "Set1"),
    labels = labels_legend
  ) +
  labs(
    title = "",
    x = "Group",
    y = "Relative abundance (%)",
    fill = "Domain"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 12, face = "plain", hjust = 0.5),
    plot.margin = margin(t = 20, r = 20, b = 0, l = 20)
  )

```
### --- Step 9: Text box below the plot ---
```{r}
plot_texto <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = texto_genes, size = 3, hjust = 0.5) +
  theme_void()

```
### --- Step 10: Combine both with patchwork ---
### plot_funanondoms <- plot_abundancia / plot_texto + plot_layout(heights = c(10, 1))
```{r}
plot_funanondoms <- plot_abundancia 

```
### Show results
```{r}
plot_funanondoms


```
## Part B

### --- Libraries ---
```{r}
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(RColorBrewer)
library(taxize)

```
### --- Configure Entrez API Key for taxize ---
```{r}
Sys.setenv(ENTREZ_KEY = "99c57399ef66d8906cbc244b43bbf0239008")

```
### --- Load data ---
```{r}
anot <- read.csv("/home/alumno21/axel/files/all_annotations_trimmed.csv",
                 sep = ",", stringsAsFactors = FALSE)

```
### --- Define files by group ---
```{r}
urban_files <- c("37082_2 # 1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                 "37082_2 # 14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                 "37082_1 # 26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                 "37082_1 # 13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                 "37035_2 # 23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                 "37082_1 # 10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                 "37035_2 # 3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                 "36703_3 # 2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                 "37082_2 # 13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                 "37082_1 # 16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                 "37082_2 # 11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                 "37035_2 # 1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                 "37035_2 # 11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                 "37082_1 # 21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                 "37035_2 # 20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                 "36703_3 # 10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                 "37035_2 # 26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22")
rural_files <- c("37082_3 # 17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                 "36703_3 # 26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                 "37082_3 # 7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                 "37035_1 # 18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                 "37035_7 # 7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                 "37082_2 # 16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                 "37082_2 # 19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                 "37035_7 # 5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                 "37035_1 # 10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                 "37082_3 # 19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                 "37035_1 # 15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                 "37082_3 # 24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                 "37082_3 # 26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                 "37082_2 # 20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                 "37082_3 # 3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                 "37035_7 # 22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                 "37035_1 # 3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                 "37082_2 # 21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                 "37082_2 # 26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                 "36703_3 # 28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                 "37082_3 # 28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                 "37082_3 # 23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                 "37082_3 # 14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                 "37035_7 # 17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                 "37035_7 # 18", "37035_7#23", "37035_7#25")

```
### --- Function to classify simplified domain ---
```{r}
assign_domain <- function(term) {
  term_lower <- tolower(term)
  if (grepl("virus|myoviridae|caudovirales|podoviridae|siphoviridae|ssdna|dsdna", term_lower)) {
    return("Viral")
  } else if (grepl("bacteria", term_lower)) {
    return("Bacteria")
  } else if (grepl("archaea", term_lower)) {
    return("Archaea")
  } else if (grepl("eukaryota", term_lower)) {
    return("Eukaryota")
  } else if (term_lower == "root") {
    return("Other")
  } else {
    return("Other")
  }
}

```
### --- Step 1: Prepare data ---
```{r}
anot <- anot %>%
  filter(grepl("\\|", max_annot_lvl)) %>%
  mutate(
    max_annot_lvl_clean = sub(".*\\|", "", max_annot_lvl),
    Domain_simplified = sapply(max_annot_lvl_clean, assign_domain),
    Grupo = case_when(
      File %in% rural_files ~ "Rural",
      File %in% urban_files ~ "Urban",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Grupo))

```
### --- Step 2: Filter only bacteria and clean Preferred_name ---
```{r}
anot_bacteria <- anot %>%
  filter(Domain_simplified == "Bacteria") %>%
  mutate(
    Preferred_name_clean = ifelse(
      grepl("\\|", Preferred_name),
      sub(".*\\|", "", Preferred_name),
      Preferred_name
    )
  ) %>%
  filter(Preferred_name_clean != "")

```
### --- Step 3: Calculate relative abundance per file and Preferred_name_clean ---
```{r}
score_total_por_file_bact <- anot_bacteria %>%
  group_by(File) %>%
  summarise(score_total = sum(score, na.rm = TRUE), .groups = "drop")

score_por_file_name <- anot_bacteria %>%
  group_by(File, Preferred_name_clean) %>%
  summarise(score_name = sum(score, na.rm = TRUE), .groups = "drop")

abundancia_relativa_name <- score_por_file_name %>%
  left_join(score_total_por_file_bact, by = "File") %>%
  mutate(abund_rel = score_name / score_total) %>%
  left_join(anot_bacteria %>% select(File, Grupo) %>% distinct(), by = "File") %>%
  filter(!is.na(abund_rel), !is.na(Preferred_name_clean))

```
### --- Step 4: Average abundance by group and Preferred_name_clean ---
```{r}
abundancia_promedio_name <- abundancia_relativa_name %>%
  group_by(Grupo, Preferred_name_clean) %>%
  summarise(promedio_abund = mean(abund_rel), .groups = "drop")

```
### --- Step 5: Calculate total abundance for each taxon ---
```{r}
porcentaje_por_name <- abundancia_promedio_name %>%
  group_by(Preferred_name_clean) %>%
  summarise(total_abund = sum(promedio_abund), .groups = "drop") %>%
  arrange(desc(total_abund))

```
### --- Step 6: Query NCBI taxonomic hierarchies for top 50 names ---
```{r}
top_raw_names <- porcentaje_por_name %>%
  slice_head(n = 50) %>%
  pull(Preferred_name_clean)

```
### Clean for taxize (eliminate "unclassified" to avoid any errors)
```{r}
cleaned_top_names <- gsub("unclassified ", "", top_raw_names)
```
### Filter empty entries, dashes or NA
```{r}
cleaned_top_names <- cleaned_top_names[cleaned_top_names != "" & cleaned_top_names != "-" & !is.na(cleaned_top_names)]

```
### Replace "Bacteroidetes" for "Bacteroidota" (NCBI approved name)
```{r}
cleaned_top_names <- gsub("^Bacteroidetes$", "Bacteroidota", cleaned_top_names)

```
###  Now consult with "ask" = TRUE for those ambiguous , or directly FLASE if already cleaned
```{r}
tax_data <- classification(cleaned_top_names, db = "ncbi", ask = FALSE)

```
### --- Step 7 (modified): Extract the most specific taxon for each name ---
```{r}
tax_levels <- c("species", "genus", "family", "order", "class", "phylum", "superkingdom")

```
### We start tax_info with the size andoriginal names to avoid issues
```{r}
tax_info <- data.frame(Name = top_raw_names,
                       Taxon_especifico = NA_character_,
                       Rank = NA_character_,
                       stringsAsFactors = FALSE)

for (i in seq_along(tax_data)) {
  taxon <- tax_data[[i]]
  
```
### Skip if taxon is NULL or from logic class (as FALSE)
```{r}
  if (is.null(taxon) || is.logical(taxon)) next
  
```
### Filter levels of interest
```{r}
  taxon_filtered <- taxon %>%
    filter(rank %in% tax_levels)
  
  if (nrow(taxon_filtered) > 0) {
    taxon_filtered$level_num <- match(taxon_filtered$rank, tax_levels)
    most_specific <- taxon_filtered %>%
      arrange(level_num) %>%
      slice(1)
    tax_info$Taxon_especifico[i] <- most_specific$name
    tax_info$Rank[i] <- most_specific$rank
  }
}

```
### Add the column Original_Name (the complete name, with possible "unclassified")
```{r}
tax_info$Original_Name <- top_raw_names

```
### Create column Taxon_final with prefix "unclassified" if the original name had it
```{r}
tax_info$Taxon_final <- ifelse(grepl("^unclassified ", tax_info$Original_Name),
                               paste0("unclassified ", tax_info$Taxon_especifico),
                               tax_info$Taxon_especifico)

```
### --- Step 9: Merge abundance with tax_info to use Taxon_final ---
```{r}
abundancia_promedio_name <- abundancia_promedio_name %>%
  left_join(tax_info %>% select(Original_Name, Taxon_final),
            by = c("Preferred_name_clean" = "Original_Name")) %>%
  mutate(Taxon_final = ifelse(is.na(Taxon_final), Preferred_name_clean, Taxon_final))

```
### --- Step 10: Select top 15 most specific taxa to plot ---
```{r}
top_taxa_final <- abundancia_promedio_name %>%
  group_by(Taxon_final) %>%
  summarise(total_abund = sum(promedio_abund), .groups = "drop") %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = 15) %>%
  pull(Taxon_final)

abundancia_promedio_taxon_top <- abundancia_promedio_name %>%
  filter(Taxon_final %in% top_taxa_final)

```
### --- Step 11: Labels and colors for plot ---
```{r}
labels_taxon_top <- abundancia_promedio_taxon_top %>%
  group_by(Taxon_final) %>%
  summarise(total_abund = sum(promedio_abund)) %>%
  ungroup() %>%
  mutate(label = paste0(Taxon_final, " (", percent(total_abund / sum(total_abund)), ")")) %>%
  arrange(match(Taxon_final, top_taxa_final)) %>%
  pull(label)
names(labels_taxon_top) <- top_taxa_final

colors_top <- colorRampPalette(brewer.pal(12, "Set3"))(15)

```
### --- Step 12: Text of annotated genes for bacteria ---
```{r}
texto_genes_bact <- paste0(
  "Bacterial annotated genes (Rural: ", sum(anot_bacteria$Grupo == "Rural"),
  " | Urban: ", sum(anot_bacteria$Grupo == "Urban"), ")"
)

```
### --- Step 13: Bar plot ---
```{r}
plot_bacteria <- ggplot(abundancia_promedio_taxon_top,
                        aes(x = Grupo, y = promedio_abund, fill = Taxon_final)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = colors_top, labels = labels_taxon_top) +
  labs(
    title = "",
    x = "Group",
    y = "Relative abundance (%)",
    fill = "Taxon"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 12, face = "plain", hjust = 0),
    plot.margin = margin(t = 20, r = 20, b = 0, l = 20)
  )

```
### --- Step 14: Text below the plot ---
```{r}
plot_texto_bact <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = texto_genes_bact, size = 3, hjust = 0.5) +
  theme_void()

```
### --- Step 15: Combine plot and text ---
### plot_bacteriafunano <- plot_bacteria / plot_texto_bact + patchwork::plot_layout(heights = c(10, 1))
```{r}
plot_bacteriafunano <- plot_bacteria 

```
### --- Show result ---
```{r}
plot_bacteriafunano




```
### Part C
```{r}

anot_eukaryota <- anot %>%
  filter(Domain_simplified == "Eukaryota") %>%
  filter(grepl("\\|", Preferred_name)) %>% # ❗️ only those with '|'
  mutate(
    Preferred_name_clean = sub(".*\\|", "", Preferred_name)
  ) %>%
  filter(Preferred_name_clean != "")


score_total_por_file_euk <- anot_eukaryota %>%
  group_by(File) %>%
  summarise(score_total = sum(score, na.rm = TRUE), .groups = "drop")

score_por_file_name_euk <- anot_eukaryota %>%
  group_by(File, Preferred_name_clean) %>%
  summarise(score_name = sum(score, na.rm = TRUE), .groups = "drop")

abundancia_relativa_name_euk <- score_por_file_name_euk %>%
  left_join(score_total_por_file_euk, by = "File") %>%
  mutate(abund_rel = score_name / score_total) %>%
  left_join(anot_eukaryota %>% select(File, Grupo) %>% distinct(), by = "File") %>%
  filter(!is.na(abund_rel), !is.na(Preferred_name_clean))


abundancia_promedio_name_euk <- abundancia_relativa_name_euk %>%
  group_by(Grupo, Preferred_name_clean) %>%
  summarise(promedio_abund = mean(abund_rel), .groups = "drop")

top_raw_names_euk <- abundancia_promedio_name_euk %>%
  group_by(Preferred_name_clean) %>%
  summarise(total_abund = sum(promedio_abund)) %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = 50) %>%
  pull(Preferred_name_clean)

cleaned_top_names_euk <- gsub("unclassified ", "", top_raw_names_euk)
cleaned_top_names_euk <- cleaned_top_names_euk[cleaned_top_names_euk != "" & cleaned_top_names_euk != "-" & !is.na(cleaned_top_names_euk)]

```
### Query taxonomy
```{r}
tax_data_euk <- classification(cleaned_top_names_euk, db = "ncbi", ask = FALSE)

tax_info_euk <- data.frame(Name = top_raw_names_euk,
                           Taxon_especifico = NA_character_,
                           Rank = NA_character_,
                           stringsAsFactors = FALSE)

for (i in seq_along(tax_data_euk)) {
  taxon <- tax_data_euk[[i]]
  if (is.null(taxon) || is.logical(taxon)) next
  taxon_filtered <- taxon %>% filter(rank %in% tax_levels)
  if (nrow(taxon_filtered) > 0) {
    taxon_filtered$level_num <- match(taxon_filtered$rank, tax_levels)
    most_specific <- taxon_filtered %>% arrange(level_num) %>% slice(1)
    tax_info_euk$Taxon_especifico[i] <- most_specific$name
    tax_info_euk$Rank[i] <- most_specific$rank
  }
}

tax_info_euk$Original_Name <- top_raw_names_euk
tax_info_euk$Taxon_final <- ifelse(grepl("^unclassified ", tax_info_euk$Original_Name),
                                   paste0("unclassified ", tax_info_euk$Taxon_especifico),
                                   tax_info_euk$Taxon_especifico)
abundancia_promedio_name_euk <- abundancia_promedio_name_euk %>%
  left_join(tax_info_euk %>% select(Original_Name, Taxon_final),
            by = c("Preferred_name_clean" = "Original_Name")) %>%
  mutate(Taxon_final = ifelse(is.na(Taxon_final), Preferred_name_clean, Taxon_final))

top_taxa_final_euk <- abundancia_promedio_name_euk %>%
  group_by(Taxon_final) %>%
  summarise(total_abund = sum(promedio_abund), .groups = "drop") %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = 15) %>%
  pull(Taxon_final)

abundancia_promedio_taxon_top_euk <- abundancia_promedio_name_euk %>%
  filter(Taxon_final %in% top_taxa_final_euk)

```
### Labels
```{r}
labels_taxon_top_euk <- abundancia_promedio_taxon_top_euk %>%
  group_by(Taxon_final) %>%
  summarise(total_abund = sum(promedio_abund)) %>%
  ungroup() %>%
  mutate(label = paste0(Taxon_final, " (", percent(total_abund / sum(total_abund)), ")")) %>%
  arrange(match(Taxon_final, top_taxa_final_euk)) %>%
  pull(label)
names(labels_taxon_top_euk) <- top_taxa_final_euk

colors_top_euk <- colorRampPalette(brewer.pal(12, "Paired"))(15)

texto_genes_euk <- paste0(
  "Eukaryotic annotated genes (Rural: ", sum(anot_eukaryota$Grupo == "Rural"),
  " | Urban: ", sum(anot_eukaryota$Grupo == "Urban"), ")"
)

```
### Plot
```{r}
plot_euk <- ggplot(abundancia_promedio_taxon_top_euk,
                   aes(x = Grupo, y = promedio_abund, fill = Taxon_final)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = colors_top_euk, labels = labels_taxon_top_euk) +
  labs(
    title = "",
    x = "Group",
    y = "Relative abundance (%)",
    fill = "Taxon"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, face = "plain", hjust = 0),
    plot.margin = margin(t = 20, r = 20, b = 0, l = 20)
  )

plot_texto_euk <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = texto_genes_euk, size = 3, hjust = 0.5) +
  theme_void()
plot_euk_final <- plot_euk / plot_texto_euk + patchwork::plot_layout(heights = c(10, 1))
plot_eukfunanon <- plot_euk 

```
### This
```{r}
plot_eukfunanon
library(patchwork)
```
## ABC UNITED
```{r}
final_plotfunaon <- ( plot_funanondoms|plot_bacteriafunano|plot_eukfunanon
) +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "") &
  theme(
    plot.tag = element_text(size = 14, face = "plain")
  )

final_plotfunaon

```


## Functional annotation of genes from selected taxa
#### In high-quality MAGs (>95% completeness). Stacked bar plots display the relative abundance of COG functional categories assigned to
annotated genes from representative bacterial groups (Bacilli,Bacteroidota, Blautia, Clostridia, and Clostridiaceae) in Rural and Urban communities.

### =========================
### Libraries
### =========================
```{r}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(forcats); library(taxize); library(readr)
  library(scales)
})

```
### =========================
### User config (provide these)
### =========================
### Sys.setenv(ENTREZ_KEY = "YOUR_NCBI_KEY")   # already set in your session
```{r}
stopifnot(exists("rural_files"), exists("urban_files"))

```
### If anot is not yet loaded:
```{r}
 anot <- read.csv("/home/alumno21/axel/files/all_annotations_trimmed.csv",
 sep = ",", stringsAsFactors = FALSE)
```
### =========================
### Helpers
### =========================
```{r}
assign_domain <- function(term) {
  term_lower <- tolower(term)
  if (grepl("virus|myoviridae|caudovirales|podoviridae|siphoviridae|ssdna|dsdna", term_lower)) "Viral"
  else if (grepl("bacteria", term_lower)) "Bacteria"
  else if (grepl("archaea", term_lower)) "Archaea"
  else if (grepl("eukaryota", term_lower)) "Eukaryota"
  else if (term_lower == "root") "Other" else "Other"
}

clean_pref_name <- function(x) {
  ifelse(grepl("\\|", x), sub(".*\\|", "", x), x)
}

```
### Robust, UID-based NCBI classification with cleaning and synonym fixes
```{r}
ncbi_classification_cached <- function(names_vec,
                                       ranks_keep = c("species","genus","family","order","class","phylum","superkingdom"),
                                       cache_path = "taxize_ncbi_cache_uid.rds") {
  canonize <- function(x) {
    x <- stringr::str_trim(x)
    x <- sub("^unclassified\\s+", "", x, ignore.case = TRUE)
    syn <- c(
      "Bacteroidetes"   = "Bacteroidota",
      "Planctomycetes"  = "Planctomycetota",
      "Spirochaetes"    = "Spirochaetota",
      "Fusobacteria"    = "Fusobacteriota",
      "Aquificae"       = "Aquificota",
      "Tenericutes"     = "Mycoplasmatota"
    )
    if (!is.na(x) && nzchar(x) && x %in% names(syn)) x <- syn[[x]]
    x
  }
  
  nm_clean <- unique(vapply(names_vec, canonize, character(1)))
  nm_clean <- nm_clean[!is.na(nm_clean) & nm_clean != "" & nm_clean != "-"]
  nm_clean <- nm_clean[grepl("[A-Za-z]", nm_clean)]
  if (length(nm_clean) == 0) {
    return(tibble::tibble(query_name = character(), lineage = character()))
  }
  
  cache <- if (file.exists(cache_path)) readRDS(cache_path) else list()
  to_query <- setdiff(nm_clean, names(cache))
  
  if (length(to_query) > 0) {
    for (nm in to_query) {
      uid <- suppressWarnings(taxize::get_uid(nm, ask = FALSE, messages = FALSE, rows = 1))
      if (is.null(uid) || is.na(uid)) { cache[[nm]] <- NA_character_; next }
      cl <- try(suppressWarnings(taxize::classification(uid, db = "ncbi")[[1]]), silent = TRUE)
      if (inherits(cl, "try-error") || is.null(cl)) { cache[[nm]] <- NA_character_; next }
      cl2 <- cl[cl$rank %in% ranks_keep, , drop = FALSE]
      lineage <- if (nrow(cl2) > 0) paste(cl2$name, collapse = ";") else NA_character_
      cache[[nm]] <- lineage
```
### Sys.sleep(0.34) # be nice to NCBI (optional)
```{r}
    }
    saveRDS(cache, cache_path)
  }
  
  tibble::tibble(
    query_name = names(cache),
    lineage    = unname(unlist(cache))
  )
}

```
### Infer COG letters from multiple columns and add COG_inferred + COG_primary
```{r}
infer_cog_from_columns <- function(df,
                                   cols_to_scan = c("COG_category","EC","KEGG_ko","KEGG_Pathway",
                                                    "Description","GOs","Preferred_name",
                                                    "eggNOG_OGs","BRITE","PFAMs","CAZy",
                                                    "BiGG_Reaction","KEGG_Module",
                                                    "KEGG_Reaction","KEGG_rclass","KEGG_TC",
                                                    "max_annot_lvl","max_annot_lvl_clean",
                                                    "Preferred_name_clean"),
                                   id_cols = c("File","MAG","query")) {
  cols_to_scan <- intersect(cols_to_scan, names(df))
  stopifnot(all(id_cols %in% names(df)))
  if (length(cols_to_scan) == 0) stop("No candidate columns to scan were found in the data frame.")
  
  valid_cogs <- c("J","A","K","L","B","D","Y","V","T","M","N","Z","W","U","O",
                  "C","G","E","F","H","I","P","Q","R","S")
  
  delim   <- "[\\s,;|/()\\[\\]{}:_-]"
  pattern <- paste0("(?:(?<=^)|(?<=", delim, "))(", paste(valid_cogs, collapse="|"),
                    ")(?=(?:$|", delim, "))")
  
  long_scan <- df %>%
    mutate(across(all_of(cols_to_scan), ~as.character(.))) %>%
    select(all_of(id_cols), everything()) %>%
    pivot_longer(cols = all_of(cols_to_scan), names_to = "field", values_to = "value") %>%
    filter(!is.na(value), value != "-") %>%
    mutate(value_up = toupper(value))
  
  cogs_raw <- long_scan %>%
    mutate(matches = stringr::str_extract_all(value_up, pattern)) %>%
    select(all_of(id_cols), matches) %>%
    tidyr::unnest(matches, keep_empty = FALSE) %>%
    rename(cog = matches) %>%
    filter(!is.na(cog))
  
  cogs_per_gene <- cogs_raw %>% distinct(across(all_of(id_cols)), cog)
  
  cog_mode <- cogs_raw %>%
    group_by(across(all_of(id_cols)), cog) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    arrange(desc(n), cog) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    rename(COG_primary = cog)
  
  cog_set <- cogs_per_gene %>%
    group_by(across(all_of(id_cols))) %>%
    summarise(COG_inferred = paste(sort(unique(cog)), collapse = ","), .groups = "drop")
  
  df %>% left_join(cog_set, by = id_cols) %>% left_join(cog_mode, by = id_cols)
}

```
### Match lineage with focus taxa (returns vector of matches per lineage)
```{r}
match_focus_taxa <- function(lineage_str, focus_taxa) {
  if (is.na(lineage_str) || lineage_str == "") return(character(0))
  hits <- focus_taxa[vapply(focus_taxa, function(tx) {
    grepl(paste0("(^|;)", tx, "(;|$)"), lineage_str)
  }, logical(1))]
  hits
}

```
### English labels for COG letters
```{r}
cog_labels <- c(
  J="Translation, ribosome, biogenesis", A="RNA/nuclear processes (rare)",
  K="Transcription", L="Replication, recombination & repair",
  B="Chromatin structure (rare)", D="Cell cycle & division",
  Y="Other uncommon", V="Intracellular trafficking & secretion",
  T="Chaperones & folding", M="Cell wall/membrane/envelope",
  N="Cell motility", Z="Cytoskeleton", W="Organelles",
  U="Cofactors & vitamins", O="Post-translational modification",
  C="Energy production & conversion", G="Carbohydrate metabolism",
  E="Amino acid metabolism", F="Nucleotide metabolism",
  H="Coenzyme metabolism", I="Lipid metabolism",
  P="Inorganic ion transport & metabolism",
  Q="Secondary metabolites biosynthesis",
  R="General function prediction only", S="Function unknown"
)

```
### =========================
### 1) Prepare base table
### =========================
```{r}
anot_prep <- anot %>%
  filter(grepl("\\|", max_annot_lvl)) %>%
  mutate(
    max_annot_lvl_clean = sub(".*\\|", "", max_annot_lvl),
    Domain_simplified   = vapply(max_annot_lvl_clean, assign_domain, character(1)),
    Grupo = case_when(
      File %in% rural_files ~ "Rural",
      File %in% urban_files ~ "Urban",
      TRUE ~ NA_character_
    ),
    Preferred_name_clean = clean_pref_name(Preferred_name)
  ) %>%
  filter(!is.na(Grupo), Preferred_name_clean != "")

```
### Add inferred COGs
```{r}
anot_prep <- infer_cog_from_columns(anot_prep)

```
### =========================
### 2) NCBI lineages for all Preferred_name_clean (cached)
### =========================
```{r}
all_names <- unique(anot_prep$Preferred_name_clean)
lin_tbl <- ncbi_classification_cached(all_names)
anot_prep <- anot_prep %>%
  left_join(lin_tbl, by = c("Preferred_name_clean" = "query_name"))

```
### =========================
### 3) Focus taxa (from your Figure 1 legends)
### =========================
```{r}
focus_bacteria <- c("Clostridia","Clostridiaceae","Blautia","Bacteroidota","Bacilli")
focus_euk      <- c("Ascomycota","Bilateria","Magnoliopsida","Opisthokonta","Streptophyta")

```
### Keep only focus taxa that actually appear in the lineages
```{r}
present_in_lineages <- function(focus, lineages) {
  hits <- vapply(focus, function(tx) any(grepl(paste0("(^|;)", tx, "(;|$)"), lineages, useBytes = TRUE)), logical(1))
  focus[hits]
}
present_bact <- present_in_lineages(focus_bacteria, na.omit(anot_prep$lineage))
present_euk  <- present_in_lineages(focus_euk,      na.omit(anot_prep$lineage))

```
### =========================
### 4) Expand to long by focus taxon
### =========================
```{r}
expand_focus <- function(df, domain_filter, focus_taxa) {
  df %>%
    filter(Domain_simplified == domain_filter, !is.na(lineage)) %>%
    mutate(focus_list = purrr::map(lineage, match_focus_taxa, focus_taxa = focus_taxa)) %>%
    filter(lengths(focus_list) > 0) %>%
    tidyr::unnest(focus_list) %>%
    rename(Focus_taxon = focus_list)
}
bact_focus_long <- expand_focus(anot_prep, "Bacteria", present_bact)
euk_focus_long  <- expand_focus(anot_prep, "Eukaryota", present_euk)

```
### =========================
### 5) Relative abundance of COG by file (counts -> proportions -> group mean)
### =========================
```{r}
summarize_cog_rel <- function(df, case_col = "Grupo") {
  df %>%
    filter(!is.na(COG_primary), COG_primary != "") %>%
    group_by(File, .data[[case_col]], Focus_taxon, COG_primary) %>%
    summarise(n_genes = n(), .groups = "drop_last") %>%
    group_by(File, .data[[case_col]], Focus_taxon) %>%
    mutate(prop = n_genes / sum(n_genes)) %>%
    ungroup() %>%
    group_by(.data[[case_col]], Focus_taxon, COG_primary) %>%
    summarise(prop_mean = mean(prop, na.rm = TRUE), .groups = "drop")
}
bact_cog_rel <- summarize_cog_rel(bact_focus_long)
euk_cog_rel  <- summarize_cog_rel(euk_focus_long)

bact_cog_rel <- bact_cog_rel %>%
  mutate(COG_lab = factor(paste0(COG_primary, " — ", cog_labels[COG_primary]),
                          levels = paste0(names(cog_labels), " — ", cog_labels)))
euk_cog_rel  <- euk_cog_rel %>%
  mutate(COG_lab = factor(paste0(COG_primary, " — ", cog_labels[COG_primary]),
                          levels = paste0(names(cog_labels), " — ", cog_labels)))

```
### =========================
### 6) Plot functions
### =========================
```{r}
plot_cog_stacks <- function(df_rel, title_txt = "") {
  ggplot(df_rel, aes(x = .data[[ "Grupo" ]], y = prop_mean, fill = COG_lab)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "Group", y = "Relative abundance (%)", fill = "COG category", title = title_txt) +
    facet_wrap(~ Focus_taxon, ncol = 2, scales = "free_x") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold")
    )
}

plot_bact_cog_by_taxa <- plot_cog_stacks(bact_cog_rel, "Bacteria: COG distribution by focus taxa")
plot_euk_cog_by_taxa  <- plot_cog_stacks(euk_cog_rel,  "Eukaryota: COG distribution by focus taxa")

```
### =========================
### 7) Show plots
### =========================
```{r}
plot_bact_cog_by_taxa
plot_euk_cog_by_taxa
```
### =========================
### 8) Improving format
### =========================

### ---- Facet order (keep your chosen order) ----
```{r}
bact_cog_rel <- bact_cog_rel %>%
  mutate(Focus_taxon = factor(Focus_taxon, levels = unique(Focus_taxon)))
euk_cog_rel  <- euk_cog_rel %>%
  mutate(Focus_taxon = factor(Focus_taxon, levels = unique(Focus_taxon)))

```
### ---- Colors (same palettes you used) ----
```{r}
present_levels_bact <- levels(droplevels(bact_cog_rel$COG_lab))
present_levels_euk  <- levels(droplevels(euk_cog_rel$COG_lab))

cols_bact <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(present_levels_bact))
names(cols_bact) <- present_levels_bact
cols_euk  <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(present_levels_euk))
names(cols_euk) <- present_levels_euk

```
### ---- Styled plotting function (thin bars, single-row facets, clean strips, no % in legend) ----
```{r}
plot_cog_single_row <- function(df_rel, present_levels, cols, title_txt = "") {
  ggplot(df_rel, aes(x = Grupo, y = prop_mean, fill = COG_lab)) +
    geom_bar(stat = "identity", position = "fill", width = 0.55) +   # thinner bars
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(
      breaks = present_levels,
      values = cols[present_levels],
      labels = present_levels,   # <-- no percentages in legend
      drop   = TRUE
    ) +
    labs(
      title = title_txt,
      x = "Group",
      y = "Relative abundance (%)",
      fill = "COG category"
    ) +
    facet_wrap(~ Focus_taxon, ncol = 3, scales = "free_x", strip.position = "top") +  # single row
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 11),
      legend.text     = element_text(size = 10),
      axis.text.x     = element_text(size = 12),
      axis.text.y     = element_text(size = 12),
      plot.title      = element_text(size = 12, face = "plain", hjust = 0),
      plot.margin     = margin(t = 12, r = 16, b = 6, l = 12),
      panel.spacing.x = unit(1.0, "lines"),
      strip.background = element_blank(),           # remove the rectangle
      strip.text       = element_text(face = "plain", size = 12, margin = margin(b = 6))
    )
}

```
### ---- Final plots in your format ----
```{r}
plot_bacteria_cog <- plot_cog_single_row(
  df_rel = bact_cog_rel,
  present_levels = present_levels_bact,
  cols = cols_bact,
  title_txt = ""  # keep empty like your 2nd figure
)

plot_eukaryota_cog <- plot_cog_single_row(
  df_rel = euk_cog_rel,
  present_levels = present_levels_euk,
  cols = cols_euk,
  title_txt = ""
)

```
### ---- Show ----
```{r}
plot_bacteria_cog
plot_eukaryota_cog
```
### =========================
### 9) Corrections
### =========================

### ---- Complete missing groups with zeros (so both Rural/Urban always appear) ----
```{r}
complete_groups_for_plot <- function(df_rel) {
  all_cogs <- sort(unique(df_rel$COG_primary))
  out <- df_rel %>%
    dplyr::group_by(Focus_taxon) %>%
    tidyr::complete(Grupo = c("Rural","Urban"),
                    COG_primary = all_cogs,
                    fill = list(prop_mean = 0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      COG_lab = factor(paste0(COG_primary, " — ", cog_labels[COG_primary]),
                       levels = paste0(names(cog_labels), " — ", cog_labels))
    )
  out
}

bact_cog_rel_full <- complete_groups_for_plot(bact_cog_rel)
euk_cog_rel_full  <- complete_groups_for_plot(euk_cog_rel)

```
### ---- Identify panels/groups with total 0 to annotate "no data" ----
```{r}
make_nodata_df <- function(df_rel_full) {
  df_rel_full %>%
    dplyr::group_by(Focus_taxon, Grupo) %>%
    dplyr::summarise(total = sum(prop_mean), .groups = "drop") %>%
    dplyr::filter(total == 0) %>%
    dplyr::mutate(y = 0.03, label = "no data")  # text near baseline
}
bact_nodata <- make_nodata_df(bact_cog_rel_full)
euk_nodata  <- make_nodata_df(euk_cog_rel_full)

```
### ---- Plot function: thin bars, single row, clean strips; add optional "no data" ----
```{r}
plot_cog_single_row <- function(df_rel, present_levels, cols, title_txt = "", nodata_df = NULL) {
  p <- ggplot(df_rel, aes(x = Grupo, y = prop_mean, fill = COG_lab)) +
    geom_bar(stat = "identity", position = "fill", width = 0.55) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(breaks = present_levels, values = cols[present_levels],
                      labels = present_levels, drop = TRUE) +
    labs(title = title_txt, x = "Group", y = "Relative abundance (%)", fill = "COG category") +
    facet_wrap(~ Focus_taxon, nrow = 1, scales = "free_x", strip.position = "top") +
    theme_classic(base_size = 12) +
    theme(
      legend.position  = "right",
      legend.title     = element_text(size = 11),
      legend.text      = element_text(size = 10),
      axis.text.x      = element_text(size = 12),
      axis.text.y      = element_text(size = 12),
      plot.title       = element_text(size = 12, face = "plain", hjust = 0),
      plot.margin      = margin(t = 12, r = 16, b = 6, l = 12),
      panel.spacing.x  = unit(1.0, "lines"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "plain", size = 12, margin = margin(b = 6))
    )
  if (!is.null(nodata_df) && nrow(nodata_df) > 0) {
    p <- p + geom_text(data = nodata_df,
                       aes(x = Grupo, y = y, label = label),
                       inherit.aes = FALSE, size = 3, color = "grey30")
  }
  p
}

```
### ---- Colors (reuse your palettes, now with completed data) ----
```{r}
present_levels_bact <- levels(droplevels(bact_cog_rel_full$COG_lab))
present_levels_euk  <- levels(droplevels(euk_cog_rel_full$COG_lab))
cols_bact <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(present_levels_bact)); names(cols_bact) <- present_levels_bact
cols_euk  <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(present_levels_euk)); names(cols_euk) <- present_levels_euk

```
### ---- Plots (now always show both bars; "no data" where appropriate) ----
```{r}
plot_bacteria_cog <- plot_cog_single_row(bact_cog_rel_full, present_levels_bact, cols_bact, "", bact_nodata)
plot_eukaryota_cog <- plot_cog_single_row(euk_cog_rel_full, present_levels_euk, cols_euk, "", euk_nodata)

```
### FINAL PLOT PRESENTED AT PAPER
```{r}
plot_bacteria_cog
plot_eukaryota_cog *not presented for missing values in one of the groups
```

<a id="diff-enr"></a>
## Differential enrichment of COG categories in Blautia and Clostridia between Rural and Urban groups
```{r}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales); library(patchwork)
})

```
### Optional: palette for Group (same as your example)
```{r}
pal_group <- c(Rural = " # E9B44C", Urban = "#4F86C6")

```
### --- Load data ---
```{r}
anot <- read.csv("/home/alumno21/axel/files/all_annotations_trimmed.csv",
                 sep = ",", stringsAsFactors = FALSE)

```
### --- Define files by group ---
```{r}
urban_files <- c("37082_2 # 1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                 "37082_2 # 14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                 "37082_1 # 26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                 "37082_1 # 13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                 "37035_2 # 23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                 "37082_1 # 10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                 "37035_2 # 3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                 "36703_3 # 2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                 "37082_2 # 13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                 "37082_1 # 16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                 "37082_2 # 11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                 "37035_2 # 1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                 "37035_2 # 11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                 "37082_1 # 21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                 "37035_2 # 20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                 "36703_3 # 10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                 "37035_2 # 26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22")
rural_files <- c("37082_3 # 17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                 "36703_3 # 26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                 "37082_3 # 7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                 "37035_1 # 18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                 "37035_7 # 7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                 "37082_2 # 16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                 "37082_2 # 19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                 "37035_7 # 5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                 "37035_1 # 10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                 "37082_3 # 19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                 "37035_1 # 15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                 "37082_3 # 24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                 "37082_3 # 26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                 "37082_2 # 20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                 "37082_3 # 3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                 "37035_7 # 22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                 "37035_1 # 3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                 "37082_2 # 21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                 "37082_2 # 26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                 "36703_3 # 28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                 "37082_3 # 28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                 "37082_3 # 23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                 "37082_3 # 14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                 "37035_7 # 17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                 "37035_7 # 18", "37035_7#23", "37035_7#25")

assign_domain <- function(term) {
  term_lower <- tolower(term)
  if (grepl("virus|myoviridae|caudovirales|podoviridae|siphoviridae|ssdna|dsdna", term_lower)) "Viral"
  else if (grepl("bacteria", term_lower)) "Bacteria"
  else if (grepl("archaea", term_lower)) "Archaea"
  else if (grepl("eukaryota", term_lower)) "Eukaryota"
  else if (term_lower == "root") "Other" else "Other"
}

clean_pref_name <- function(x) {
  ifelse(grepl("\\|", x), sub(".*\\|", "", x), x)
}
```
### Robust, UID-based NCBI classification with cleaning and synonym fixes
```{r}
ncbi_classification_cached <- function(names_vec,
                                       ranks_keep = c("species","genus","family","order","class","phylum","superkingdom"),
                                       cache_path = "taxize_ncbi_cache_uid.rds") {
  canonize <- function(x) {
    x <- stringr::str_trim(x)
    x <- sub("^unclassified\\s+", "", x, ignore.case = TRUE)
    syn <- c(
      "Bacteroidetes"   = "Bacteroidota",
      "Planctomycetes"  = "Planctomycetota",
      "Spirochaetes"    = "Spirochaetota",
      "Fusobacteria"    = "Fusobacteriota",
      "Aquificae"       = "Aquificota",
      "Tenericutes"     = "Mycoplasmatota"
    )
    if (!is.na(x) && nzchar(x) && x %in% names(syn)) x <- syn[[x]]
    x
  }
  
  nm_clean <- unique(vapply(names_vec, canonize, character(1)))
  nm_clean <- nm_clean[!is.na(nm_clean) & nm_clean != "" & nm_clean != "-"]
  nm_clean <- nm_clean[grepl("[A-Za-z]", nm_clean)]
  if (length(nm_clean) == 0) {
    return(tibble::tibble(query_name = character(), lineage = character()))
  }
  
  cache <- if (file.exists(cache_path)) readRDS(cache_path) else list()
  to_query <- setdiff(nm_clean, names(cache))
  
  if (length(to_query) > 0) {
    for (nm in to_query) {
      uid <- suppressWarnings(taxize::get_uid(nm, ask = FALSE, messages = FALSE, rows = 1))
      if (is.null(uid) || is.na(uid)) { cache[[nm]] <- NA_character_; next }
      cl <- try(suppressWarnings(taxize::classification(uid, db = "ncbi")[[1]]), silent = TRUE)
      if (inherits(cl, "try-error") || is.null(cl)) { cache[[nm]] <- NA_character_; next }
      cl2 <- cl[cl$rank %in% ranks_keep, , drop = FALSE]
      lineage <- if (nrow(cl2) > 0) paste(cl2$name, collapse = ";") else NA_character_
      cache[[nm]] <- lineage
```
### Sys.sleep(0.34) # be nice to NCBI (optional)
```{r}
    }
    saveRDS(cache, cache_path)
  }
  
  tibble::tibble(
    query_name = names(cache),
    lineage    = unname(unlist(cache))
  )
}
infer_cog_from_columns <- function(df,
                                   cols_to_scan = c("COG_category","EC","KEGG_ko","KEGG_Pathway",
                                                    "Description","GOs","Preferred_name",
                                                    "eggNOG_OGs","BRITE","PFAMs","CAZy",
                                                    "BiGG_Reaction","KEGG_Module",
                                                    "KEGG_Reaction","KEGG_rclass","KEGG_TC",
                                                    "max_annot_lvl","max_annot_lvl_clean",
                                                    "Preferred_name_clean"),
                                   id_cols = c("File","MAG","query")) {
  cols_to_scan <- intersect(cols_to_scan, names(df))
  stopifnot(all(id_cols %in% names(df)))
  if (length(cols_to_scan) == 0) stop("No candidate columns to scan were found in the data frame.")
  
  valid_cogs <- c("J","A","K","L","B","D","Y","V","T","M","N","Z","W","U","O",
                  "C","G","E","F","H","I","P","Q","R","S")
  
  delim   <- "[\\s,;|/()\\[\\]{}:_-]"
  pattern <- paste0("(?:(?<=^)|(?<=", delim, "))(", paste(valid_cogs, collapse="|"),
                    ")(?=(?:$|", delim, "))")
  
  long_scan <- df %>%
    mutate(across(all_of(cols_to_scan), ~as.character(.))) %>%
    select(all_of(id_cols), everything()) %>%
    pivot_longer(cols = all_of(cols_to_scan), names_to = "field", values_to = "value") %>%
    filter(!is.na(value), value != "-") %>%
    mutate(value_up = toupper(value))
  
  cogs_raw <- long_scan %>%
    mutate(matches = stringr::str_extract_all(value_up, pattern)) %>%
    select(all_of(id_cols), matches) %>%
    tidyr::unnest(matches, keep_empty = FALSE) %>%
    rename(cog = matches) %>%
    filter(!is.na(cog))
  
  cogs_per_gene <- cogs_raw %>% distinct(across(all_of(id_cols)), cog)
  
  cog_mode <- cogs_raw %>%
    group_by(across(all_of(id_cols)), cog) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    arrange(desc(n), cog) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    rename(COG_primary = cog)
  
  cog_set <- cogs_per_gene %>%
    group_by(across(all_of(id_cols))) %>%
    summarise(COG_inferred = paste(sort(unique(cog)), collapse = ","), .groups = "drop")
  
  df %>% left_join(cog_set, by = id_cols) %>% left_join(cog_mode, by = id_cols)
}

```
### Match lineage with focus taxa (returns vector of matches per lineage)
```{r}
match_focus_taxa <- function(lineage_str, focus_taxa) {
  if (is.na(lineage_str) || lineage_str == "") return(character(0))
  hits <- focus_taxa[vapply(focus_taxa, function(tx) {
    grepl(paste0("(^|;)", tx, "(;|$)"), lineage_str)
  }, logical(1))]
  hits
}
```
### --- Step 1: Prepare data ---
```{r}
anot <- anot %>%
  filter(grepl("\\|", max_annot_lvl)) %>%
  mutate(
    max_annot_lvl_clean = sub(".*\\|", "", max_annot_lvl),
    Domain_simplified = sapply(max_annot_lvl_clean, assign_domain),
    Grupo = case_when(
      File %in% rural_files ~ "Rural",
      File %in% urban_files ~ "Urban",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Grupo))

anot_prep <- anot %>%
  filter(grepl("\\|", max_annot_lvl)) %>%
  mutate(
    max_annot_lvl_clean = sub(".*\\|", "", max_annot_lvl),
    Domain_simplified   = vapply(max_annot_lvl_clean, assign_domain, character(1)),
    Grupo = case_when(
      File %in% rural_files ~ "Rural",
      File %in% urban_files ~ "Urban",
      TRUE ~ NA_character_
    ),
    Preferred_name_clean = clean_pref_name(Preferred_name)
  ) %>%
  filter(!is.na(Grupo), Preferred_name_clean != "")

```
### Add inferred COGs
```{r}
anot_prep <- infer_cog_from_columns(anot_prep)

```
### =========================
### 2) NCBI lineages for all Preferred_name_clean (cached)
### =========================
```{r}
all_names <- unique(anot_prep$Preferred_name_clean)
lin_tbl <- ncbi_classification_cached(all_names)
anot_prep <- anot_prep %>%
  left_join(lin_tbl, by = c("Preferred_name_clean" = "query_name"))

```
### =========================
### 3) Focus taxa (from your Figure 1 legends)
### =========================
```{r}
focus_bacteria <- c("Clostridia","Clostridiaceae","Blautia","Bacteroidota","Bacilli")
focus_euk      <- c("Ascomycota","Bilateria","Magnoliopsida","Opisthokonta","Streptophyta")

```
### Keep only focus taxa that actually appear in the lineages
```{r}
present_in_lineages <- function(focus, lineages) {
  hits <- vapply(focus, function(tx) any(grepl(paste0("(^|;)", tx, "(;|$)"), lineages, useBytes = TRUE)), logical(1))
  focus[hits]
}
present_bact <- present_in_lineages(focus_bacteria, na.omit(anot_prep$lineage))
present_euk  <- present_in_lineages(focus_euk,      na.omit(anot_prep$lineage))

```
### =========================
### 4) Expand to long by focus taxon
### =========================
```{r}
expand_focus <- function(df, domain_filter, focus_taxa) {
  df %>%
    filter(Domain_simplified == domain_filter, !is.na(lineage)) %>%
    mutate(focus_list = purrr::map(lineage, match_focus_taxa, focus_taxa = focus_taxa)) %>%
    filter(lengths(focus_list) > 0) %>%
    tidyr::unnest(focus_list) %>%
    rename(Focus_taxon = focus_list)
}
bact_focus_long <- expand_focus(anot_prep, "Bacteria", present_bact)
euk_focus_long  <- expand_focus(anot_prep, "Eukaryota", present_euk)

```
### If you don't have cog_labels in your env, define a minimal mapping
```{r}
if (!exists("cog_labels")) {
  cog_labels <- c(
    J="Translation, ribosome, biogenesis", A="RNA/nuclear processes (rare)",
    K="Transcription", L="Replication, recombination & repair",
    B="Chromatin structure (rare)", D="Cell cycle & division",
    Y="Other uncommon", V="Intracellular trafficking & secretion",
    T="Chaperones & folding", M="Cell wall/membrane/envelope",
    N="Cell motility", Z="Cytoskeleton", W="Organelles",
    U="Cofactors & vitamins", O="Post-translational modification",
    C="Energy production & conversion", G="Carbohydrate metabolism",
    E="Amino acid metabolism", F="Nucleotide metabolism",
    H="Coenzyme metabolism", I="Lipid metabolism",
    P="Inorganic ion transport & metabolism",
    Q="Secondary metabolites biosynthesis",
    R="General function prediction only", S="Function unknown"
  )
}

```
### ===============================
### Main plot: grouped bars + brackets + stars (no SE, no q)
### ===============================

### ---- Build per-file compositions (reusable) ----
```{r}
per_file_cog <- function(focus_long_df) {
  focus_long_df %>%
    dplyr::filter(!is.na(COG_primary), COG_primary != "") %>%
    dplyr::group_by(File, Grupo, Focus_taxon, COG_primary) %>%
    dplyr::summarise(n_genes = dplyr::n(), .groups = "drop_last") %>%
    dplyr::group_by(File, Grupo, Focus_taxon) %>%
    dplyr::mutate(prop = n_genes / sum(n_genes)) %>%
    dplyr::ungroup()
}

bact_perfile <- per_file_cog(bact_focus_long)
euk_perfile  <- per_file_cog(euk_focus_long)

plot_cog_groupbars_style <- function(perfile_df, focus_taxon,
                                     top_n = 6, alpha = 0.05, min_diff = 0.03, min_n = 3,
                                     sig_only = TRUE, show_values = TRUE) {
  df <- perfile_df %>% dplyr::filter(Focus_taxon == !!focus_taxon)
  
```
### Choose top COGs by overall mean
```{r}
  top_cogs <- df %>%
    dplyr::group_by(COG_primary) %>%
    dplyr::summarise(mean_all = mean(prop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_all)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(COG_primary)
  df <- df %>% dplyr::filter(COG_primary %in% top_cogs)
  
```
### Group means (%)
```{r}
  group_avg <- df %>%
    dplyr::group_by(Grupo, COG_primary) %>%
    dplyr::summarise(mean_pct = 100*mean(prop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      COG_lab = factor(paste0(COG_primary, " — ", cog_labels[COG_primary]),
                       levels = paste0(top_cogs, " — ", cog_labels[top_cogs]))
    )
  
```
### Wilcoxon per COG (BH within panel)
```{r}
  tests <- df %>%
    dplyr::group_by(COG_primary) %>%
    dplyr::summarise(
      n_R = sum(Grupo == "Rural"),
      n_U = sum(Grupo == "Urban"),
      p = {
        x <- prop[Grupo == "Rural"]; y <- prop[Grupo == "Urban"]
        if (length(x) >= min_n && length(y) >= min_n)
          suppressWarnings(stats::wilcox.test(x, y)$p.value)
        else NA_real_
      },
      mean_R = 100*mean(prop[Grupo == "Rural"], na.rm = TRUE),
      mean_U = 100*mean(prop[Grupo == "Urban"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      q = p.adjust(p, method = "BH"),
      delta = mean_R - mean_U,
      stars = dplyr::case_when(
        is.na(q) ~ "",
        q < 1e-4 ~ "****",
        q < 1e-3 ~ "***",
        q < 1e-2 ~ "**",
        q < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  tests_sig <- tests %>% dplyr::filter(stars != "", abs(delta) >= 100*min_diff)
  
```
### Keep only significant COGs if requested
```{r}
  if (sig_only) {
    if (nrow(tests_sig) == 0) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "No significant COGs",
                   size = 5, fontface = "bold") +
          labs(title = focus_taxon) +
          theme_void(base_size = 12)
      )
    }
    keep_cogs <- tests_sig$COG_primary
    group_avg <- group_avg %>% dplyr::filter(COG_primary %in% keep_cogs)
    tests     <- tests_sig
  } else {
    tests <- tests_sig
  }
  
  
  
```
### --- DROP unused levels after the sig_only filter ---
```{r}
  group_avg <- group_avg %>%
    dplyr::mutate(COG_lab = forcats::fct_drop(COG_lab))
  
```
### Y-position of the bracket: a little bit higher than the highest barr of each COG
```{r}
  anno_y <- group_avg %>%
    dplyr::group_by(COG_lab) %>%
    dplyr::summarise(y = max(mean_pct, na.rm = TRUE) * 1.10, .groups = "drop")
  
```
###  Map of the real positions in x for what has to be drawn
```{r}
  x_map <- group_avg %>%
    dplyr::distinct(COG_lab) %>%
    dplyr::arrange(COG_lab) %>%
    dplyr::mutate(gx = dplyr::row_number(),
                  off = 0.23, # mismo offset que usas en las barras
                  x_rural = gx - off,
                  x_urban = gx + off)
  
```
### Brackets + stars using the map of real positions
```{r}
  brackets_manual <- tests %>%
    dplyr::mutate(
      COG_lab = factor(paste0(COG_primary, " — ", cog_labels[COG_primary]),
                       levels = levels(group_avg$COG_lab))
    ) %>%
    dplyr::inner_join(x_map,  by = "COG_lab") %>%
    dplyr::left_join(anno_y, by = "COG_lab")
  
  
  
  p <- ggplot(group_avg, aes(x = COG_lab, y = mean_pct, fill = Grupo)) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    {if (show_values)
      geom_text(aes(label = sprintf("%.1f%%", mean_pct)),
                position = position_dodge(width = 0.8),
                vjust = 0, size = 3)} +
```
### bracket & stars (single panel, grouped bars)
```{r}
    geom_segment(data = brackets_manual,
                 aes(x = x_rural, xend = x_urban, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(data = brackets_manual,
                 aes(x = x_rural, xend = x_rural, y = y, yend = y - 0.02*max(group_avg$mean_pct, na.rm=TRUE)),
                 inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(data = brackets_manual,
                 aes(x = x_urban, xend = x_urban, y = y, yend = y - 0.02*max(group_avg$mean_pct, na.rm=TRUE)),
                 inherit.aes = FALSE, linewidth = 0.6) +
    geom_text(data = brackets_manual,
              aes(x = gx, y = y + 0.01*max(group_avg$mean_pct, na.rm=TRUE), label = stars),
              inherit.aes = FALSE, vjust = -0.2, size = 5) +
    scale_fill_manual(values = pal_group, name = "Group") +
    scale_x_discrete(labels = function(x) substr(x, 1, 1)) +
    scale_y_continuous(labels = label_percent(scale = 1, accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.18))) +
    labs(title = focus_taxon, x = "COG", y = "Mean per sample") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x.bottom = element_line(),
      axis.line.y.left   = element_line(),
      axis.text.x = element_text(angle = 0, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 11, face = "plain",
                                margin = margin(t = 8, b = 6)),
      plot.margin = margin(15, 12, 12, 12)
    )
  p
}

```
### ===============================
### 
### ===============================
```{r}
make_taxon_plots_groupstyle <- function(perfile_df, taxa_vec,
                                        top_n = 6, alpha = 0.05, min_diff = 0.03, min_n = 3,
                                        sig_only = TRUE, show_values = TRUE) {
  pls <- lapply(taxa_vec, function(tx)
    plot_cog_groupbars_style(
      perfile_df, focus_taxon = tx,
      top_n = top_n, alpha = alpha, min_diff = min_diff, min_n = min_n,
      sig_only = sig_only, show_values = show_values
    )
  )
  names(pls) <- taxa_vec
  pls
}

```
### ===============================
### PLOTS
### ===============================
### Choose taxa present in your data:
```{r}
bact_taxa <- intersect(c("Blautia","Clostridia"),
                       unique(bact_perfile$Focus_taxon))

euk_taxa  <- intersect(c("Ascomycota","Bilateria","Magnoliopsida","Opisthokonta","Streptophyta"),
                       unique(euk_perfile$Focus_taxon))

```
### Build lists of plots
```{r}
bact_plots <- make_taxon_plots_groupstyle(bact_perfile, bact_taxa,
                                          top_n = 6, min_diff = 0.03, sig_only = TRUE)
euk_plots  <- make_taxon_plots_groupstyle(euk_perfile,  euk_taxa,
                                          top_n = 6, min_diff = 0.03, sig_only = TRUE)

```
### Combine with patchwork (one row; shared legend)
```{r}
bact_combined <- wrap_plots(bact_plots, nrow = 1, guides = "collect") &
  theme(legend.position = "right")
euk_combined  <- wrap_plots(euk_plots,  nrow = 1, guides = "collect") &
  theme(legend.position = "right")

```
### Show
```{r}
bact_combined
```
##### euk_combined no significant data


<a id="occurr"></a>
## Co-occurrence networks of bacterial and eukaryotic taxa based on Spearman correlations.

### --- Required libraries ---
```{r}
library(dplyr)
library(tidyr)
library(tibble)
library(boot)
library(stringr)

kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

```
### If separate taxonomic columns do not exist but 'taxon_name' does
```{r}
if (!all(c("Domain","Genus") %in% names(kaiju_merged)) && "taxon_name" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    separate(
      taxon_name,
      into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
               "Subclass","Order","Family","Genus","Species"),
      sep = ";", fill = "right", extra = "drop"
    )
}

```
### Fill minimal gaps
```{r}
kaiju_merged <- kaiju_merged %>%
  mutate(across(c(Genus, Family, Order), ~ ifelse(is.na(.) | . == "", "Unclassified", .)))

```
### file_base if needed
```{r}
if (!"file_base" %in% names(kaiju_merged)) {
  if ("file" %in% names(kaiju_merged)) {
    kaiju_merged <- kaiju_merged %>% mutate(file_base = gsub("^.*/", "", file))
  } else {
    stop("No 'file_base' or 'file' column found in kaiju_merged.")
  }
}

```
### --- Assign groups (adjust your 'urban_files' and 'rural_files' lists) ---
```{r}
kaiju_merged <- kaiju_merged %>%
  mutate(Group = case_when(
    file_base %in% paste0(urban_files, "_kaiju.out") ~ "Urban",
    file_base %in% paste0(rural_files, "_kaiju.out") ~ "Rural",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))
cat("Número de muestras asignadas a grupos:\n"); print(table(kaiju_merged$Group))

```
### --- Filter ONLY Bacteria/Eukaryota and keep GENUS ---
```{r}
full_data <- kaiju_merged %>%
  filter(Domain %in% c("Bacteria","Eukaryota")) %>%
  filter(!is.na(Genus) & Genus != "" & Genus != "Unclassified") %>%
  mutate(Genus = str_trim(Genus))

```
### ------------------------------------------------------------------
### 2) PRE-FILTERING OF ABUNDANT GENERA (accelerates bootstrap)
### (same thresholds for both domains)
### ------------------------------------------------------------------
```{r}
threshold_genus <- 8 # lecturas promedio mínimas

```
### Bacteria
```{r}
genus_bact_urban <- full_data %>% filter(Group == "Urban",  Domain == "Bacteria") %>%
  group_by(Genus) %>% summarise(mean_reads = mean(reads), .groups = "drop") %>%
  filter(mean_reads > threshold_genus) %>% pull(Genus)

genus_bact_rural <- full_data %>% filter(Group == "Rural",  Domain == "Bacteria") %>%
  group_by(Genus) %>% summarise(mean_reads = mean(reads), .groups = "drop") %>%
  filter(mean_reads > threshold_genus) %>% pull(Genus)

common_genus_bact <- intersect(genus_bact_urban, genus_bact_rural)

```
### Eukaryota
```{r}
genus_euk_urban <- full_data %>% filter(Group == "Urban",  Domain == "Eukaryota") %>%
  group_by(Genus) %>% summarise(mean_reads = mean(reads), .groups = "drop") %>%
  filter(mean_reads > threshold_genus) %>% pull(Genus)

genus_euk_rural <- full_data %>% filter(Group == "Rural",  Domain == "Eukaryota") %>%
  group_by(Genus) %>% summarise(mean_reads = mean(reads), .groups = "drop") %>%
  filter(mean_reads > threshold_genus) %>% pull(Genus)

common_genus_euk <- intersect(genus_euk_urban, genus_euk_rural)

cat("Common genera after filter — Bacteria:", length(common_genus_bact),
    " | Eukaryota:", length(common_genus_euk), "\n")

```
### Keep only those abundant genera per domain
```{r}
full_data_filtered <- full_data %>%
  filter((Domain == "Bacteria"  & Genus %in% common_genus_bact) |
           (Domain == "Eukaryota" & Genus %in% common_genus_euk))

```
### ------------------------------------------------------------------
### 3) Function to create matrix (GENUS) log2-normalized per sample
### ------------------------------------------------------------------
```{r}
crear_matriz <- function(df, taxon_col = "Genus") {
  df %>%
    group_by(file_base, !!sym(taxon_col)) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(taxon_col), values_from = reads, values_fill = 0) %>%
    column_to_rownames("file_base") %>%
    {
      mat <- sweep(., 1, rowSums(.), "/") # proporción por muestra
      log2(mat + 1e-6)
    }
}

```
### Matrices per group (Bacteria + Eukaryota together, columns = GENUS)
```{r}
mat_urban <- full_data_filtered %>% filter(Group == "Urban") %>% crear_matriz("Genus")
mat_rural <- full_data_filtered %>% filter(Group == "Rural") %>% crear_matriz("Genus")

cat("Dimensiones matriz Urban (genus):", dim(mat_urban), "\n")
cat("Dimensiones matriz Rural (genus):", dim(mat_rural), "\n")

```
### GENUS -> DOMAIN map
```{r}
genus_domain_map <- full_data_filtered %>%
  distinct(Genus, Domain)

```
### ------------------------------------------------------------------
### 4) Bootstrap Spearman
### ------------------------------------------------------------------
```{r}
bootstrap_correlation <- function(vec1, vec2, n = 200) {
  stat <- function(data, indices) {
    tryCatch(cor(data[indices,1], data[indices,2], method="spearman"),
             warning=function(w) NA, error=function(e) NA)
  }
  boot_obj <- boot(data = cbind(vec1, vec2), statistic = stat, R = n)
  if (all(is.na(boot_obj$t))) return(c(mean = NA, lower = NA, upper = NA))
  ci <- tryCatch(boot.ci(boot_obj, type = "perc"), error = function(e) NULL)
  c(mean = mean(boot_obj$t, na.rm = TRUE),
    lower = ifelse(is.null(ci), NA, ci$percent[4]),
    upper = ifelse(is.null(ci), NA, ci$percent[5]))
}

```
### ------------------------------------------------------------------
### 5) Correlations (all at GENUS level) and relationship classification
### ------------------------------------------------------------------
```{r}
correlacion_spearman_boot <- function(matrix, genus_domain_map) {
  out <- tibble()
  pairs <- combn(colnames(matrix), 2, simplify = FALSE)
  for (pp in pairs) {
    g1 <- pp[1]; g2 <- pp[2]
    x <- matrix[[g1]]; y <- matrix[[g2]]
    if (any(is.na(x)) || any(is.na(y))) next
    if (sd(x) == 0 || sd(y) == 0) next
    
    d1 <- genus_domain_map %>% filter(Genus == g1) %>% pull(Domain)
    d2 <- genus_domain_map %>% filter(Genus == g2) %>% pull(Domain)
    if (length(d1) != 1 || length(d2) != 1) next
    
    relation <- if (d1 == "Bacteria"   & d2 == "Bacteria")   "Bacteria-Bacteria" else
      if (d1 == "Eukaryota" & d2 == "Eukaryota") "Eukaryota-Eukaryota" else
        "Bacteria-Eukaryota"
    boot_res <- bootstrap_correlation(x, y, n = 200)
    out <- bind_rows(out, tibble(
      Genus_1 = g1, Genus_2 = g2, Relation = relation,
      Spearman_rho_mean = boot_res["mean"],
      CI_lower = boot_res["lower"], CI_upper = boot_res["upper"]
    ))
  }
  out
}

```
### Run by group
```{r}
cat("Calculando correlaciones (Urban)...\n")
edges_urban <- correlacion_spearman_boot(mat_urban, genus_domain_map)
cat("Calculando correlaciones (Rural)...\n")
edges_rural <- correlacion_spearman_boot(mat_rural, genus_domain_map)

```
### ------------------------------------------------------------------
### 6) Filters by relationship and export
### ------------------------------------------------------------------
```{r}
filtrar_por_relacion <- function(edges_df, relation_type) {
  rho_thresh <- dplyr::case_when(
    relation_type == "Bacteria-Bacteria"     ~ 0.7,
    relation_type == "Bacteria-Eukaryota"    ~ 0.3,
    relation_type == "Eukaryota-Eukaryota"   ~ 0.3,
    TRUE ~ 0.3
  )
  edges_df %>%
    filter(
      Relation == relation_type,
      abs(Spearman_rho_mean) > rho_thresh,
      (CI_lower > 0 & CI_upper > 0) | (CI_lower < 0 & CI_upper < 0)
    )
}

crear_nodos_de_edges <- function(edges_df, nodes_df) {
  genes_in_edges <- unique(c(edges_df$Genus_1, edges_df$Genus_2))
  nodes_df %>%
    filter(Genus %in% genes_in_edges) %>%
    transmute(ID = gsub(" ", "_", Genus), Domain)
}

relaciones <- c("Bacteria-Bacteria", "Bacteria-Eukaryota", "Eukaryota-Eukaryota")
grupos <- list(Urban = edges_urban, Rural = edges_rural)

for (grupo_nombre in names(grupos)) {
  edges_total <- grupos[[grupo_nombre]]
  for (rel in relaciones) {
    edges_filtradas <- filtrar_por_relacion(edges_total, rel)
    if (nrow(edges_filtradas) == 0) {
      cat("Warning: No correlations in", grupo_nombre, "-", rel, "\n")
      next
    }
    nodes_filtrados <- crear_nodos_de_edges(edges_filtradas, full_data_filtered)
    
```
### Output files
```{r}
    edges_sif <- edges_filtradas %>%
      mutate(Source = gsub(" ", "_", Genus_1),
             Target = gsub(" ", "_", Genus_2),
             interaction = "correlation") %>%
      select(Source, interaction, Target)
    
    edges_attr <- edges_filtradas %>%
      mutate(Source = gsub(" ", "_", Genus_1),
             Target = gsub(" ", "_", Genus_2)) %>%
      select(Source, Target, Relation, Spearman_rho_mean)
    
    nodes_attr <- nodes_filtrados %>%
      transmute(name = ID, Domain)
    
    write.table(edges_sif,
                paste0("network_", tolower(grupo_nombre), "_", rel, ".sif"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(edges_attr,
                paste0("edges_", tolower(grupo_nombre), "_", rel, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    write.table(nodes_attr,
                paste0("nodes_", tolower(grupo_nombre), "_", rel, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
}

cat("Done! Networks (GENUS) exported by relationship and group.\n")



```
### ===================== Adjustments (much stricter) =====================
```{r}




suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(boot)
})
set.seed(1234)
```
### Filters by group (prevalence and mean %, NORMALIZED BY DOMAIN)
```{r}
MIN_PREV_BACT <- 0.20
MIN_PREV_EUK  <- 0.10
MIN_MEAN_BACT <- 0.06 # % dentro de Bacteria
MIN_MEAN_EUK  <- 0.015 # % dentro de Eukaryota
TOP_N_BACT    <- 250
TOP_N_EUK     <- 120

```
### Correlación / Bootstrap
```{r}
BOOT_R               <- 300
PSEUDO               <- 1e-6
PRESCREEN_ABS_RHO    <- 0.35 # global (antes 0.20)
PRESCREEN_ABS_RHO_WL <- 0.20 # si toca whitelist
MIN_SAMPLES_PER_G    <- 8

```
### Minimum co-occurrence (samples where BOTH > 0)
```{r}
MIN_COOC_FRAC <- 0.35 # ≥35% de las muestras del grupo
MIN_COOC_ABS  <- 20 # y al menos 20 muestras

```
### Final edge filter (only BE)
```{r}
THRESH_BE <- 0.45 # |rho_mean| mínimo
ALPHA_Q   <- 0.05 # BH FDR

```
### ===================== Whitelist of key genera =====================
```{r}
force_include_euk  <- c("Saccharomyces")
force_include_bact <- character(0)

```
### ======================== Datos esperados en 'kaiju_merged' ========================
```{r}
kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

if (!all(c("Domain","Genus") %in% names(kaiju_merged)) && "taxon_name" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    tidyr::separate(
      taxon_name,
      into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
               "Subclass","Order","Family","Genus","Species"),
      sep = ";", fill = "right", extra = "drop"
    )
}

if (!"file_base" %in% names(kaiju_merged)) {
  if ("file" %in% names(kaiju_merged)) kaiju_merged <- kaiju_merged %>% mutate(file_base = gsub("^.*/", "", file))
  else stop("kaiju_merged needs 'file_base' or 'file'.")
}

if (!"Group" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    mutate(Group = case_when(
      exists("urban_files") & file_base %in% paste0(urban_files, "_kaiju.out") ~ "Urban",
      exists("rural_files") & file_base %in% paste0(rural_files, "_kaiju.out") ~ "Rural",
      TRUE ~ NA_character_
    ))
}
kaiju_merged <- kaiju_merged %>% filter(!is.na(Group))

```
### Only Bacteria/Eukaryota at GENUS level; exclude Homo
```{r}
dat <- kaiju_merged %>%
  filter(Domain %in% c("Bacteria","Eukaryota")) %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus=="" | Genus=="Unclassified", NA, str_trim(Genus))) %>%
  filter(!is.na(Genus)) %>%
  filter(!(Domain=="Eukaryota" & Genus=="Homo"))

```
### ------------------------ % within DOMAIN per sample ------------------------
```{r}
totals_domain <- dat %>%
  group_by(file_base, Group, Domain) %>%
  summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")

genus_reads <- dat %>%
  group_by(file_base, Group, Domain, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  left_join(totals_domain, by = c("file_base","Group","Domain")) %>%
  mutate(pct_domain = if_else(total_reads_domain > 0, 100 * reads / total_reads_domain, 0))

```
### ------------------------ prevalence and means per group -------------------
```{r}
prev_mean <- genus_reads %>%
  group_by(Group, Domain, Genus) %>%
  summarise(prevalence = mean(reads > 0),
            mean_pct_dom = mean(pct_domain),
            .groups = "drop")

```
### Keepers per domain + enforce whitelist
```{r}
keepers_bact <- prev_mean %>%
  filter(Domain=="Bacteria") %>%
  group_by(Genus) %>%
  summarise(
    ok = all(prevalence[Group=="Rural"] >= MIN_PREV_BACT,
             prevalence[Group=="Urban"] >= MIN_PREV_BACT,
             mean_pct_dom[Group=="Rural"] >= MIN_MEAN_BACT,
             mean_pct_dom[Group=="Urban"] >= MIN_MEAN_BACT),
    overall = mean(mean_pct_dom),
    .groups = "drop"
  ) %>% filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_BACT) %>% pull(Genus)

keepers_euk <- prev_mean %>%
  filter(Domain=="Eukaryota") %>%
  group_by(Genus) %>%
  summarise(
    ok = all(prevalence[Group=="Rural"] >= MIN_PREV_EUK,
             prevalence[Group=="Urban"] >= MIN_PREV_EUK,
             mean_pct_dom[Group=="Rural"] >= MIN_MEAN_EUK,
             mean_pct_dom[Group=="Urban"] >= MIN_MEAN_EUK),
    overall = mean(mean_pct_dom),
    .groups = "drop"
  ) %>% filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_EUK) %>% pull(Genus)

keepers_bact <- union(keepers_bact,
                      intersect(force_include_bact, unique(dat$Genus[dat$Domain=="Bacteria"])))
keepers_euk  <- union(keepers_euk,
                      intersect(force_include_euk,  unique(dat$Genus[dat$Domain=="Eukaryota"])))

message("Kept genera (domain-normalized) — Bacteria=", length(keepers_bact),
        " | Eukaryota=", length(keepers_euk))

dat_filt <- genus_reads %>%
  filter((Domain=="Bacteria"  & Genus %in% keepers_bact) |
           (Domain=="Eukaryota" & Genus %in% keepers_euk))

```
### ======================= matrices per group (counts + CLR) ======================
```{r}
counts_from_long <- function(df_long, group) {
  df_long %>%
    filter(Group == group) %>%
    select(file_base, Genus, reads) %>%
    group_by(file_base, Genus) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    pivot_wider(names_from = Genus, values_from = reads, values_fill = 0) %>%
    column_to_rownames("file_base") %>% as.matrix()
}

clr_from_counts <- function(counts_mat) {
  if (nrow(counts_mat) < MIN_SAMPLES_PER_G) return(NULL)
  prop <- sweep(counts_mat, 1, rowSums(counts_mat), "/"); prop[is.na(prop)] <- 0
  logx <- log(prop + PSEUDO); sweep(logx, 1, rowMeans(logx), "-")
}

cnt_rural <- counts_from_long(dat_filt, "Rural")
cnt_urban <- counts_from_long(dat_filt, "Urban")
stopifnot(nrow(cnt_rural) >= MIN_SAMPLES_PER_G, nrow(cnt_urban) >= MIN_SAMPLES_PER_G)

mat_clr_rural <- clr_from_counts(cnt_rural)
mat_clr_urban <- clr_from_counts(cnt_urban)

```
### Alinear columnas
```{r}
common_cols <- Reduce(intersect, list(colnames(mat_clr_rural), colnames(mat_clr_urban)))
mat_clr_rural <- mat_clr_rural[, common_cols, drop = FALSE]
mat_clr_urban <- mat_clr_urban[, common_cols, drop = FALSE]
cnt_rural     <- cnt_rural[,     common_cols, drop = FALSE]
cnt_urban     <- cnt_urban[,     common_cols, drop = FALSE]

```
### Genus -> Domain map (only common)
```{r}
genus_domain_map <- dat_filt %>% distinct(Genus, Domain) %>% filter(Genus %in% common_cols)

```
### =================== Prescreen + co-occurrence + bootstrap =================
```{r}
fast_pairwise_corr_BE <- function(mat_clr, counts_mat, genus_domain_map,
                                  prescreen_global = PRESCREEN_ABS_RHO,
                                  prescreen_wl     = PRESCREEN_ABS_RHO_WL,
                                  min_cooc_frac = MIN_COOC_FRAC,
                                  min_cooc_abs  = MIN_COOC_ABS,
                                  boot_R = BOOT_R,
                                  whitelist = character(0)) {
  
  cols <- colnames(mat_clr)
  dom_vec <- setNames(as.character(genus_domain_map$Domain),
                      as.character(genus_domain_map$Genus))
  bact_cols <- cols[dom_vec[cols] == "Bacteria"]
  euk_cols  <- cols[dom_vec[cols] == "Eukaryota"]
  if (length(bact_cols) == 0 || length(euk_cols) == 0) return(tibble())
  
```
### Co-occurrence
```{r}
  pres_abs <- max(min_cooc_abs, ceiling(min_cooc_frac * nrow(counts_mat)))
  pres_mat <- counts_mat > 0
  
```
### Spearman masivo = Pearson sobre rangos
```{r}
  mat_rank <- apply(mat_clr, 2, rank, ties.method = "average")
  R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                  mat_rank[, euk_cols,  drop = FALSE],
                  method = "pearson", use = "pairwise.complete.obs")
  
```
### Candidates: global prescreen
```{r}
  sel_global <- which(abs(R) >= prescreen_global, arr.ind = TRUE)
  pairs_df <- tibble(
    Genus_1 = bact_cols[sel_global[,"row"]],
    Genus_2 = euk_cols[ sel_global[,"col"]],
    from_whitelist = FALSE
  )
  
```
### Extra candidates: if whitelist applies, looser prescreen
```{r}
  if (length(whitelist) > 0) {
    wl_b <- intersect(bact_cols, whitelist)
    wl_e <- intersect(euk_cols,  whitelist)
    if (length(wl_b) > 0) {
      R_wlb <- R[match(wl_b, bact_cols), , drop = FALSE]
      sel_wlb <- which(abs(R_wlb) >= prescreen_wl, arr.ind = TRUE)
      if (length(sel_wlb))
        pairs_df <- bind_rows(pairs_df,
                              tibble(Genus_1 = wl_b[sel_wlb[,"row"]],
                                     Genus_2 = euk_cols[sel_wlb[,"col"]],
                                     from_whitelist = TRUE))
    }
    if (length(wl_e) > 0) {
      R_wle <- R[, match(wl_e, euk_cols), drop = FALSE]
      sel_wle <- which(abs(R_wle) >= prescreen_wl, arr.ind = TRUE)
      if (length(sel_wle))
        pairs_df <- bind_rows(pairs_df,
                              tibble(Genus_1 = bact_cols[sel_wle[,"row"]],
                                     Genus_2 = wl_e[sel_wle[,"col"]],
                                     from_whitelist = TRUE))
    }
  }
  
```
### Unique + minimum co-occurrence
```{r}
  pairs_df <- pairs_df %>%
    distinct(Genus_1, Genus_2, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(cooc = sum(pres_mat[, Genus_1] & pres_mat[, Genus_2])) %>%
    ungroup() %>%
    filter(cooc >= pres_abs)
  
  message("Prescreen pairs (post co-occur): ", nrow(pairs_df),
          " (", sum(pairs_df$from_whitelist), " via whitelist)")
  
  if (nrow(pairs_df) == 0) return(tibble())
  
  boot_one <- function(g1, g2) {
    x <- mat_clr[, g1]; y <- mat_clr[, g2]
    if (sd(x) == 0 || sd(y) == 0) return(NULL)
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    rho0 <- unname(ct$estimate); p0 <- ct$p.value
    stat <- function(dd, idx) suppressWarnings(cor(dd[idx,1], dd[idx,2], method = "spearman"))
    bt <- boot(cbind(x, y), statistic = stat, R = boot_R)
    if (all(is.na(bt$t))) return(NULL)
    ci <- tryCatch(boot.ci(bt, type = "perc"), error = function(e) NULL)
```
### Sign consistency in bootstrap
```{r}
    scons <- mean(sign(bt$t) == sign(rho0), na.rm = TRUE)
    tibble(
      Genus_1 = g1, Genus_2 = g2, Relation = "Bacteria-Eukaryota",
      rho = rho0, p = p0,
      rho_mean = mean(bt$t, na.rm = TRUE),
      CI_lower = ifelse(is.null(ci), NA, ci$percent[4]),
      CI_upper = ifelse(is.null(ci), NA, ci$percent[5]),
      sign_consistency = scons
    )
  }
  
  out_list <- lapply(seq_len(nrow(pairs_df)), function(i) boot_one(pairs_df$Genus_1[i], pairs_df$Genus_2[i]))
  bind_rows(out_list)
}

message("Computing BE correlations (Rural)…")
edges_rural <- fast_pairwise_corr_BE(
  mat_clr_rural, cnt_rural, genus_domain_map,
  prescreen_global = PRESCREEN_ABS_RHO,
  prescreen_wl     = PRESCREEN_ABS_RHO_WL,
  whitelist = union(force_include_euk, force_include_bact),
  boot_R = BOOT_R
)

message("Computing BE correlations (Urban)…")
edges_urban <- fast_pairwise_corr_BE(
  mat_clr_urban, cnt_urban, genus_domain_map,
  prescreen_global = PRESCREEN_ABS_RHO,
  prescreen_wl     = PRESCREEN_ABS_RHO_WL,
  whitelist = union(force_include_euk, force_include_bact),
  boot_R = BOOT_R
)

```
### BH (todo es BE)
```{r}
add_q <- function(df) if (nrow(df)) df %>% mutate(q = p.adjust(p, method = "BH")) else df
edges_rural <- add_q(edges_rural)
edges_urban <- add_q(edges_urban)

```
### --------- Very strict final filter ----------
```{r}
filter_BE <- function(df) {
  if (nrow(df)==0) return(df[0,])
  df %>%
    mutate(abs_rho = abs(rho_mean)) %>%
```
### CI must exceed threshold (not just same sign)
```{r}
    mutate(ci_strong = (rho_mean > 0 & CI_lower >= THRESH_BE) |
             (rho_mean < 0 & CI_upper <= -THRESH_BE)) %>%
    filter(
      abs_rho >= THRESH_BE,
      abs(rho) >= THRESH_BE, # efecto puntual también grande
      ci_strong,
      sign_consistency >= 0.90, # ≥90% de réplicas con el mismo signo
      q <= ALPHA_Q
    ) %>%
    select(-abs_rho, -ci_strong)
}

edges_rural_f <- filter_BE(edges_rural)
edges_urban_f <- filter_BE(edges_urban)

message("BE edges kept | Rural: ", nrow(edges_rural_f), " | Urban: ", nrow(edges_urban_f))

```
### ============================= Nodes & export =============================
```{r}
make_nodes <- function(edges_df) {
  if (nrow(edges_df)==0) return(tibble(name=character(), Domain=character()))
  nodes <- unique(c(edges_df$Genus_1, edges_df$Genus_2))
  tibble(name = nodes) %>% left_join(genus_domain_map, by = c("name"="Genus"))
}

write_net <- function(edges, nodes, tag) {
  if (nrow(edges) == 0) { message("No edges for ", tag); return(invisible(NULL)) }
  write.table(edges %>% transmute(Source = Genus_1, interaction = "corr", Target = Genus_2),
              paste0("network_", tag, ".sif"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(edges %>% transmute(Source = Genus_1, Target = Genus_2, Relation,
                                  rho_mean, CI_lower, CI_upper, p, q, sign_consistency),
              paste0("edges_",   tag, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(nodes %>% mutate(name = gsub(" ", "_", name)) %>% select(name, Domain),
              paste0("nodes_",   tag, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

nodes_rural <- make_nodes(edges_rural_f)
nodes_urban <- make_nodes(edges_urban_f)

write_net(edges_rural_f, nodes_rural, "rural_BE_genus_tight")
write_net(edges_urban_f, nodes_urban, "urban_BE_genus_tight")

message("Ready. BE networks VERY tight (domain-normalized, strong CI, co-occurrence, sign-consistent).")



```





















































<a id="diff-dis"></a>
### Differential distribution of functional annotations reconstructed from high-quality metagenome-assembled genomes (MAGs)

### =========================================
### 1. Load base annotations
### =========================================
```{r}
anot <- read.csv("/home/alumno21/axel/files/all_annotations_trimmed.csv",
                 sep = ",", stringsAsFactors = FALSE)

```
### --- Define files by group ---
```{r}
urban_files <- c("37082_2 # 1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
                 "37082_2 # 14", "37035_2#12", "37035_2#6", "37082_1#17", "37035_2#14",
                 "37082_1 # 26", "37035_1#30", "37035_1#32", "37082_1#15", "37082_2#15",
                 "37082_1 # 13", "37035_2#10", "37082_1#31", "37035_2#17", "37035_2#8",
                 "37035_2 # 23", "37035_2#31", "37035_2#24", "37082_2#5", "36703_3#5",
                 "37082_1 # 10", "36703_3#7", "37082_2#9", "37082_2#3", "37035_2#2",
                 "37035_2 # 3", "37035_2#19", "37035_2#21", "36703_3#1", "37082_1#24",
                 "36703_3 # 2", "37035_2#4", "37035_2#15", "37035_2#18", "37035_2#28",
                 "37082_2 # 13", "37082_1#22", "37082_1#29", "37082_1#19", "37035_2#30",
                 "37082_1 # 16", "37035_1#31", "37035_2#7", "37082_1#30", "37035_2#16",
                 "37082_2 # 11", "37082_1#14", "37035_2#5", "37082_2#4", "37082_1#18",
                 "37035_2 # 1", "37082_1#23", "37082_2#12", "37082_1#11", "37082_1#12",
                 "37035_2 # 11", "37035_2#25", "37082_1#32", "37082_1#9", "37035_2#29",
                 "37082_1 # 21", "37082_2#2", "37035_2#27", "36703_3#3", "37082_2#6",
                 "37035_2 # 20", "37082_2#7", "37082_2#8", "37082_2#10", "37082_1#28",
                 "36703_3 # 10", "37035_2#9", "37082_1#25", "36703_3#8", "36703_3#9",
                 "37035_2 # 26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22")
rural_files <- c("37082_3 # 17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
                 "36703_3 # 26", "37035_7#10", "36703_3#21", "37082_2#22", "37035_7#2",
                 "37082_3 # 7", "37035_7#6", "37035_1#7", "37035_7#9", "37082_2#30",
                 "37035_1 # 18", "37035_7#4", "37082_3#13", "37082_3#32", "37035_1#8",
                 "37035_7 # 7", "37035_1#19", "37082_3#29", "37035_7#13", "37035_7#12",
                 "37082_2 # 16", "36703_3#25", "37082_3#27", "37082_3#5", "37082_3#21",
                 "37082_2 # 19", "37082_3#16", "37035_1#5", "37082_3#1", "37035_7#11",
                 "37035_7 # 5", "36703_3#13", "37035_7#14", "37035_1#1", "37082_3#11",
                 "37035_1 # 10", "37035_1#12", "37082_3#4", "36703_3#17", "36703_3#27",
                 "37082_3 # 19", "37082_2#18", "36703_3#29", "36703_3#12", "36703_3#32",
                 "37035_1 # 15", "37035_1#27", "37035_1#13", "37035_7#8", "37035_1#6",
                 "37082_3 # 24", "36703_3#30", "37035_7#1", "37035_1#16", "37035_7#15",
                 "37082_3 # 26", "37035_1#23", "37035_1#2", "37082_2#27", "37035_7#3",
                 "37082_2 # 20", "36703_3#16", "37082_3#8", "37035_1#25", "36703_3#14",
                 "37082_3 # 3", "37035_1#4", "37082_2#29", "37082_3#30", "37082_2#31",
                 "37035_7 # 22", "37035_7#16", "37082_2#17", "36703_3#18", "37035_1#11",
                 "37035_1 # 3", "37035_1#14", "37082_3#9", "36703_3#23", "37082_2#28",
                 "37082_2 # 21", "37082_3#31", "36703_3#20", "37082_2#25", "36703_3#19",
                 "37082_2 # 26", "37082_3#6", "37035_1#17", "37082_2#23", "36703_3#15",
                 "36703_3 # 28", "37082_3#12", "37082_2#32", "37082_3#10", "36703_3#22",
                 "37082_3 # 28", "36703_3#24", "37082_3#18", "37082_3#20", "37035_1#24",
                 "37082_3 # 23", "37082_3#2", "37035_1#20", "37082_3#22", "37082_3#25",
                 "37082_3 # 14", "37035_1#9", "36703_3#11", "37035_1#21", "37035_7#20",
                 "37035_7 # 17", "37035_7#21", "37035_7#19", "37035_1#26", "37035_7#24",
                 "37035_7 # 18", "37035_7#23", "37035_7#25")

```
### --- Function to classify simplified domain ---
```{r}
assign_domain <- function(term) {
  term_lower <- tolower(term)
  if (grepl("virus|myoviridae|caudovirales|podoviridae|siphoviridae|ssdna|dsdna", term_lower)) {
    return("Viral")
  } else if (grepl("bacteria", term_lower)) {
    return("Bacteria")
  } else if (grepl("archaea", term_lower)) {
    return("Archaea")
  } else if (grepl("eukaryota", term_lower)) {
    return("Eukaryota")
  } else if (term_lower == "root") {
    return("Other")
  } else {
    return("Other")
  }
}

```
### --- Step 1: Prepare data ---
```{r}
anot <- anot %>%
  filter(grepl("\\|", max_annot_lvl)) %>%
  mutate(
    max_annot_lvl_clean = sub(".*\\|", "", max_annot_lvl),
    Domain_simplified = sapply(max_annot_lvl_clean, assign_domain),
    Grupo = case_when(
      File %in% rural_files ~ "Rural",
      File %in% urban_files ~ "Urban",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Grupo))

```
### =========================================
### 2. Crear anot_with_cog
### =========================================
### Add COG_primary column from IDs/descriptions
```{r}
anot_with_cog <- infer_cog_from_columns(anot)
```
### --- 1) Ensure File column ---
```{r}
if (!"File" %in% names(anot_with_cog)) {
  cand_file <- intersect(c("sample","Sample","archivo","Archivo"), names(anot_with_cog))
  if (length(cand_file) == 1) {
    anot_with_cog <- anot_with_cog %>% rename(File = !!sym(cand_file))
  } else {
    stop("Column 'File' not found. Rename the sample column to 'File'.")
  }
}


```
### Normalize normal values
```{r}
anot_with_cog <- anot_with_cog %>%
  mutate(Grupo = ifelse(Grupo %in% c("Rural","Urban"), Grupo, NA_character_))

```
### --- 3) Ensure COG_primary (if your 'infer_cog_from_columns' has not set it yet) ---
```{r}
if (!"COG_primary" %in% names(anot_with_cog)) {
  anot_with_cog <- infer_cog_from_columns(anot_with_cog)
}
if (!"COG_primary" %in% names(anot_with_cog)) {
  stop("'COG_primary' still missing. Check 'infer_cog_from_columns'.")
}

```
### --- 4) Ensure pathway column (metabolic_pathways or kegg_pathway_ids) ---
### If you don't have 'metabolic_pathways' but you do have standalone KEGG IDs, build 'kegg_pathway_ids'
```{r}
if (!("metabolic_pathways" %in% names(anot_with_cog) || "kegg_pathway_ids" %in% names(anot_with_cog))) {
```
### Try extracting 'mapXXXXX' from multiple candidate columns
```{r}
  cand_cols <- intersect(
    c("KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_ko","BRITE",
      "MetaCyc_Pathway","BiGG_Reaction","GOs","eggNOG_OGs","Description"),
    names(anot_with_cog)
  )
  if (length(cand_cols) == 0) {
    stop("No pathway columns found. Provide 'metabolic_pathways' or some KEGG information.")
  }
  extract_ids <- function(x, pattern, normalize = function(v) v) {
    m <- stringr::str_extract_all(x, pattern)
    vapply(m, function(v) {
      if (length(v) == 0) return(NA_character_)
      v <- normalize(v); v <- unique(v[nzchar(v)])
      ifelse(length(v) > 0, paste(v, collapse = ";"), NA_character_)
    }, character(1))
  }
  anot_with_cog <- anot_with_cog %>%
    tidyr::unite(".all_text", dplyr::all_of(cand_cols), sep = " | ", remove = FALSE, na.rm = TRUE) %>%
    mutate(
      kegg_pathway_ids = extract_ids(.all_text, "(?i)\\b(?:map|ko)\\d{5}\\b",
                                     function(v){ v <- tolower(v); sub("^ko","map",v) })
    ) %>%
    select(-.all_text)
```
### if nothing is available, at least create an empty column so the code does not fail
```{r}
  if (!"kegg_pathway_ids" %in% names(anot_with_cog)) {
    anot_with_cog$kegg_pathway_ids <- NA_character_
  }
}

```
### --- 5) lineage_ncbi (robust and cached) ---
```{r}
if (!"lineage_ncbi" %in% names(anot_with_cog)) {
```
### column with the name to query in NCBI
```{r}
  if (!"Preferred_name_ncbi" %in% names(anot_with_cog)) {
    anot_with_cog <- anot_with_cog %>%
      mutate(Preferred_name_ncbi = ifelse(grepl("\\|", Preferred_name %||% ""),
                                          sub(".*\\|", "", Preferred_name),
                                          Preferred_name))
  }
  
```
### name cleaner to avoid empty queries
```{r}
  .clean_ncbi_name <- function(x) {
    x <- trimws(x)
    x <- sub("^unclassified\\s+", "", x, ignore.case = TRUE)
    x <- sub("^uncultured\\s+",   "", x, ignore.case = TRUE)
    x <- sub("^metagenome\\s+",   "", x, ignore.case = TRUE)
    x[!nzchar(x)] <- NA_character_
    x
  }
  
```
### cached and error-tolerant function
```{r}
  ncbi_classification_cached <- function(names_vec,
                                         ranks_keep = c("species","genus","family","order","class","phylum","superkingdom"),
                                         cache_path = "taxize_ncbi_cache_uid.rds") {
```
### Try loading Taxize; if it's not there, leave with NA without breaking
```{r}
    if (!requireNamespace("taxize", quietly = TRUE)) {
      warning("Paquete 'taxize' no disponible; 'lineage_ncbi' se dejará como NA.")
      return(tibble(query_name = unique(names_vec), lineage_ncbi = NA_character_))
    }
    library(taxize)
    
    nm <- unique(.clean_ncbi_name(names_vec))
```
### Filter empty spaces, dashes and entries without letter
```{r}
    nm <- nm[!is.na(nm) & nm != "-" & grepl("[A-Za-z]", nm)]
    
    cache <- if (file.exists(cache_path)) readRDS(cache_path) else list()
    to_query <- setdiff(nm, names(cache))
    
    if (length(to_query) > 0) {
      for (n in to_query) {
```
### Protection: do not send empty queries
```{r}
        if (is.na(n) || !nzchar(n) || !grepl("[A-Za-z]", n)) { cache[[n]] <- NA_character_; next }
        
        uid <- trySuppressWarnings(get_uid(n, ask = FALSE, messages = FALSE, rows = 1))
        if (inherits(uid, "try-error") || is.null(uid) || length(uid) == 0 || is.na(uid)) {
          cache[[n]] <- NA_character_; next
        }
        
        cl <- trySuppressWarnings(classification(uid, db = "ncbi")[[1]])
        if (inherits(cl, "try-error") || is.null(cl)) {
          cache[[n]] <- NA_character_; next
        }
        
        cl2 <- cl[cl$rank %in% ranks_keep, , drop = FALSE]
        cache[[n]] <- if (nrow(cl2) > 0) paste(cl2$name, collapse = ";") else NA_character_
        
      }
      saveRDS(cache, cache_path)
    }
    
    tibble(query_name = names(cache), lineage_ncbi = unname(unlist(cache)))
  }
  
```
### silent helpers
```{r}
  trySuppressWarnings <- function(expr) suppressWarnings(try(expr, silent = TRUE))
  
```
### run the function (with global failure handling)
```{r}
  lin_tbl <- trySuppressWarnings(
    ncbi_classification_cached(anot_with_cog$Preferred_name_ncbi)
  )
  
  if (inherits(lin_tbl, "try-error") || is.null(lin_tbl)) {
    warning("Could not retrieve NCBI taxonomy; 'lineage_ncbi' will remain as NA.")
    anot_with_cog$lineage_ncbi <- NA_character_
  } else {
    anot_with_cog <- anot_with_cog %>%
      left_join(lin_tbl, by = c("Preferred_name_ncbi" = "query_name"))
  }
}

```
### --- 6) Domain_simplified (if missing) ---
```{r}
if (!"Domain_simplified" %in% names(anot_with_cog)) {
  anot_with_cog <- anot_with_cog %>%
    mutate(Domain_simplified = case_when(
      str_detect(lineage_ncbi %||% "", "(?i)(^|;)Bacteria(;|$)")   ~ "Bacteria",
      str_detect(lineage_ncbi %||% "", "(?i)(^|;)Eukaryota(;|$)")  ~ "Eukaryota",
      TRUE ~ NA_character_
    ))
}

```
### --- 7) Minimal final cleaning ---
### Ensure that Group only has Rural/Urban and remove rows without these groups
```{r}
anot_with_cog <- anot_with_cog %>% filter(Grupo %in% c("Rural","Urban"))

```
### =========================================
### Librerías
### =========================================
```{r}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(ggrepel); library(forcats)
})

```
### =========================================
### Helpers básicos + núcleo de análisis
### =========================================
```{r}

.col_or_stop <- function(df, candidates, must = TRUE) {
  hit <- intersect(candidates, names(df))
  if (!length(hit)) { if (must) stop("No encuentro: ", paste(candidates, collapse=", ")); return(NA_character_) }
  hit[[1]]
}

.detect_taxa_any <- function(lineage_vec, taxa) {
  if (!length(lineage_vec)) return(logical(0))
  safe <- ifelse(is.na(lineage_vec), "", lineage_vec)
  pat  <- paste0("(^|;)", taxa, "(;|$)")
  stringr::str_detect(safe, regex(pat, ignore_case = TRUE))
}

.separate_ids <- function(df, colname) {
  df %>%
    filter(!is.na(.data[[colname]]), .data[[colname]] != "", .data[[colname]] != "-") %>%
    mutate(.raw = .data[[colname]]) %>%
    tidyr::separate_rows(.raw, sep = "\\s*;\\s*") %>%
    mutate(path_id = str_trim(.raw)) %>%
    filter(path_id != "") %>%
    select(-.raw)
}

.fisher_p <- function(a, A, b, B) {
  if (any(is.na(c(a,A,b,B)))) return(NA_real_)
  mat <- matrix(c(a, A-a, b, B-b), nrow=2, byrow=TRUE)
  suppressWarnings(stats::fisher.test(mat)$p.value)
}

```
### ---- Pathway table per microorganism (presence/absence per file) ----
```{r}
analyze_paths_for_micro <- function(anot,
                                    taxa,
                                    cog_keep = NULL, # p.ej. "V" o "LA"
                                    domain_hint = c("Bacteria","Eukaryota"),
                                    microorganism_label = taxa,
                                    min_files_per_group = 1,
                                    drop_overview = c("map01100","map00000")) {
  
  col_file    <- .col_or_stop(anot, c("File"))
  col_group   <- .col_or_stop(anot, c("Grupo","Group"))
  col_domain  <- .col_or_stop(anot, c("Domain_simplified","domain"))
  col_lineage <- .col_or_stop(anot, c("lineage_ncbi","lineage"))
  col_cog     <- .col_or_stop(anot, c("COG_primary","primary_cog"), must = FALSE)
  col_path    <- .col_or_stop(anot, c("metabolic_pathways","kegg_pathway_ids","path_id"))
  
  df0 <- anot %>%
    filter(.data[[col_domain]] %in% domain_hint) %>%
    filter(.detect_taxa_any(.data[[col_lineage]], taxa))
  
  if (!is.null(cog_keep)) {
    if (is.na(col_cog)) stop("No COG column (COG_primary) found to apply the filter.")
    keep_letters <- strsplit(cog_keep, "")[[1]]
    df0 <- df0 %>% filter(.data[[col_cog]] %in% keep_letters)
  }
  
  sub_paths <- df0 %>%
    .separate_ids(col_path) %>%
    distinct(File = .data[[col_file]], Grupo = .data[[col_group]], path_id)
  
  if (length(drop_overview)) {
    sub_paths <- sub_paths %>% filter(!path_id %in% drop_overview)
  }
  
  if (nrow(sub_paths) == 0) {
    warning("No pathways for: ", microorganism_label)
    return(tibble())
  }
  
  files_per_group <- anot %>%
    distinct(File = .data[[col_file]], Grupo = .data[[col_group]]) %>%
    count(Grupo, name = "n_files_group")
  
  counts <- sub_paths %>%
    count(Grupo, path_id, name = "files_with") %>%
    tidyr::complete(Grupo = unique(files_per_group$Grupo), path_id, fill = list(files_with = 0L)) %>%
    left_join(files_per_group, by = "Grupo")
  
  wide <- counts %>%
    pivot_wider(names_from = Grupo, values_from = c(files_with, n_files_group), values_fill = 0)
  
  for (nm in c("files_with_Rural","files_with_Urban","n_files_group_Rural","n_files_group_Urban"))
    if (!nm %in% names(wide)) wide[[nm]] <- 0L
  
  res <- wide %>%
    transmute(
      path_id,
      files_with_Rural = files_with_Rural,
      files_with_Urban = files_with_Urban,
      n_files_Rural    = n_files_group_Rural,
      n_files_Urban    = n_files_group_Urban,
      prev_Rural       = ifelse(n_files_Rural>0, files_with_Rural/n_files_Rural, NA_real_),
      prev_Urban       = ifelse(n_files_Urban>0, files_with_Urban/n_files_Urban, NA_real_),
      prev_Rural_pct   = 100*prev_Rural,
      prev_Urban_pct   = 100*prev_Urban,
      delta_pp         = prev_Urban_pct - prev_Rural_pct,
      p = mapply(.fisher_p,
                 a = files_with_Urban, A = n_files_Urban,
                 b = files_with_Rural, B = n_files_Rural),
      Microorganism = microorganism_label,
      domain_type   = domain_hint[[1]]
    ) %>%
    filter(files_with_Rural >= min_files_per_group | files_with_Urban >= min_files_per_group) %>%
    mutate(
      q = p.adjust(p, method = "BH"),
      enriched = case_when(
        is.na(delta_pp) ~ "Tie",
        delta_pp > 0    ~ "Urban",
        delta_pp < 0    ~ "Rural",
        TRUE            ~ "Tie"
      )
    )
  
  res
}

```
### =========================================
### Build subsets (Blautia, Clostridia, Ascomycota) and merge into volcano_df
### =========================================
```{r}

res_blautia <- analyze_paths_for_micro(
  anot_with_cog, taxa = "Blautia",   cog_keep = "V",
  domain_hint = "Bacteria", microorganism_label = "Blautia · COG V",
  min_files_per_group = 1
)

res_clostridia <- analyze_paths_for_micro(
  anot_with_cog, taxa = "Clostridia", cog_keep = "LA",
  domain_hint = "Bacteria", microorganism_label = "Clostridia · COG L/A",
  min_files_per_group = 1
)

res_ascomycota <- analyze_paths_for_micro(
  anot_with_cog, taxa = "Ascomycota", cog_keep = NULL,
  domain_hint = "Eukaryota", microorganism_label = "Ascomycota",
  min_files_per_group = 1
)

volcano_df <- bind_rows(res_blautia, res_clostridia, res_ascomycota) %>%
  mutate(q_global = p.adjust(p, method = "BH"))

```
### =========================================
### Plotting function and final plot
### =========================================
```{r}

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel); library(forcats)
})

plot_volcano_all <- function(df_all,
                             label_top_frac = 0.15, # etiqueta solo el top 15%
                             p_for_y = c("q_global","q_micro","p_raw"),
                             point_size = 3, label_size = 1.6) {
  
  p_for_y <- match.arg(p_for_y)
  
  dfp <- df_all %>%
    mutate(
      p_raw    = p,
      q_micro  = q,
      q_global = ifelse(is.na(q_global), p.adjust(p_raw, method = "BH"), q_global),
      yval = dplyr::case_when(
        p_for_y == "q_global" ~ -log10(q_global),
        p_for_y == "q_micro"  ~ -log10(q_micro),
        TRUE                  ~ -log10(p_raw)
      ),
```
### Maintain original labels for the tables, but create a simple name for the legend
```{r}
      Microorganism = factor(Microorganism,
                             levels = c("Blautia · COG V","Clostridia · COG L/A","Ascomycota")),
      Micro_simple  = fct_recode(Microorganism,
                                 "Blautia"    = "Blautia · COG V",
                                 "Clostridia" = "Clostridia · COG L/A",
                                 "Ascomycota" = "Ascomycota"),
      enriched = factor(enriched, levels = c("Rural","Urban","Tie"))
    )
  
```
### Top 15% por evidencia (yval alto = q pequeño) Top 15% by evidence (yval high = q low)
```{r}
  n_lab  <- max(1, floor(nrow(dfp) * label_top_frac))
  lab_df <- dfp %>% arrange(desc(yval)) %>% slice_head(n = n_lab)
  
```
### Palettes (shapes with fillers; "Tie" in gray but outside the legend)
```{r}
  fill_cols <- c(Rural = " # E9B44C", Urban = "#4F86C6", Tie = "grey70")
  shp_vals  <- c("Blautia" = 21, "Clostridia" = 22, "Ascomycota" = 24)
  
  ggplot(dfp, aes(x = delta_pp, y = yval)) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
    geom_vline(xintercept = 0,         linetype = 2, linewidth = 0.3) +
    geom_point(aes(shape = Micro_simple, fill = enriched),
               size = point_size, stroke = 0.3, color = "black", alpha = 0.8) +
    scale_shape_manual(values = shp_vals, name = "Taxa") +
    scale_fill_manual(
      values = fill_cols,
      breaks = c("Rural","Urban"), # oculta "Tie" en la leyenda
      name   = "Group",
      guide  = guide_legend(override.aes = list(shape = 21, color = "black"))
    )+
    labs(
      x = "Data distribuition",
      y = "−log10(q) (Benjamini–Hochberg)",
      title = ""
    ) +
    theme_classic(base_size = 9) +
    theme(legend.position = "right",
          plot.margin = margin(8,12,8,8)) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = path_id),
      size = label_size,
      min.segment.length = 0,
      segment.size = 0.2,
      box.padding = 0.12,
      point.padding = 0.15,
      force = 3, force_pull = 2,
      max.overlaps = Inf,
      seed = 1,
      direction = "both"
    )
}

```
### === Plot ===
```{r}
p <- plot_volcano_all(
  volcano_df,
  label_top_frac = 0.15, p_for_y = "q_global",
  point_size = 2.6, label_size = 1.4
)
print(p)

```


