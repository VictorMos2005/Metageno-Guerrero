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
```{r}

```
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
### TABLA FINAL BUSCO COLORES
```{r}
combined_busco_plot
```
