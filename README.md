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

## Fig. 1

### FIGURE 1. PCoA in the anthropometric–demographic space and its association with BMI percentile, age group, and lifestyle
```{r}


```
```
### --- REQUIRED LIBRARIES ---
```{r}
library(ggplot2)
library(vegan)
library(dplyr)
library(ggcorrplot)
library(patchwork)

```
```
### --- 1. 1. DATA LOADING AND CLEANING ---
```{r}
data <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE,
  sep = ",",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

```
```
### Remove empty columns and rows (only those completely empty)
```{r}
data <- data[, colSums(!is.na(data)) > 0]
data <- data[rowSums(is.na(data)) < ncol(data), ]

```
```
### Numeric variables
```{r}
data$BMI <- as.numeric(data$BMI)
data$Age <- as.numeric(data$Age)

```
```
### --- 2. 2. DERIVED VARIABLES ---
```{r}
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
```
### --- 3. 3. FILTERING for analysis ---
```{r}
data_complete <- data %>%
  filter(!is.na(BMI), !is.na(Age), !is.na(Percentil_group), !is.na(Age_group), !is.na(Lifestyle))

```
```
### --- 4. 4. DISTANCE MATRIX AND PCoA ---
```{r}
dist_matrix <- vegdist(data_complete[, c("BMI", "Age")], method = "bray")
pcoa <- cmdscale(dist_matrix, k = 2, eig = TRUE)
scores_pcoa <- as.data.frame(pcoa$points)
colnames(scores_pcoa) <- c("Dim1", "Dim2")

```
```
### Add categorical variables to color points
```{r}
scores_pcoa$Percentil_group <- data_complete$Percentil_group
scores_pcoa$Age_group <- data_complete$Age_group
scores_pcoa$Lifestyle <- data_complete$Lifestyle

```
```
### Labels with counts for Percentil_group
```{r}
percentil_counts <- data_complete %>%
  count(Percentil_group)
percentil_labels <- setNames(paste0(percentil_counts$Percentil_group, " (", percentil_counts$n, ")"), percentil_counts$Percentil_group)

```
```
### --- 5. 5. PERMANOVA ---
```{r}
set.seed(123)
adonis_Percentil <- adonis2(dist_matrix ~ Percentil_group, data = data_complete)
adonis_Age <- adonis2(dist_matrix ~ Age_group, data = data_complete)
adonis_Lifestyle <- adonis2(dist_matrix ~ Lifestyle, data = data_complete)

cat("PERMANOVA RESULTS:\n")
print(adonis_Percentil)
print(adonis_Age)
print(adonis_Lifestyle)

```
```
### --- 6. 6. UNIFIED THEME ---
```{r}
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
```
### --- 7. 7. PLOTS ---
```
### Labels with counts for Age_group
```{r}
age_counts <- data_complete %>%
  count(Age_group)
age_labels <- setNames(paste0(age_counts$Age_group, " (", age_counts$n, ")"), age_counts$Age_group)

```
```
### Labels with counts for Lifestyle
```{r}
lifestyle_counts <- data_complete %>%
  count(Lifestyle)
lifestyle_labels <- setNames(paste0(lifestyle_counts$Lifestyle, " (", lifestyle_counts$n, ")"), lifestyle_counts$Lifestyle)

```
```
### --- 7. 7. PLOTS ---
```{r}
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
```
### Counts across the whole dataset for Age_group and Lifestyle
```{r}
age_counts_full <- data %>%
  filter(!is.na(Age_group)) %>%
  count(Age_group)

age_labels_full <- setNames(paste0(age_counts_full$Age_group, " (", age_counts_full$n, ")"), age_counts_full$Age_group)

lifestyle_counts_full <- data %>%
  filter(!is.na(Lifestyle)) %>%
  count(Lifestyle)

lifestyle_labels_full <- setNames(paste0(lifestyle_counts_full$Lifestyle, " (", lifestyle_counts_full$n, ")"), lifestyle_counts_full$Lifestyle)

```
```
### Use these complete labels in the plots
```{r}
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
```
### Create data_filtered for correlation with Percentil_num and Lifestyle_num
```{r}
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
```
### Correlation with Percentil_num, Age, and Lifestyle_num
```{r}
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
```
### --- 9. 9. COMBINE PLOTS ---
```{r}
library(patchwork)

```
```
### Themes to adjust margins by row
```{r}
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
```
### --- Load required libraries ---
```{r}
library(ggplot2)
library(dplyr)
library(ggpubr)
library(grid)
library(gtable)
library(patchwork)
library(tibble)

```
```
### --- Load data ---
```{r}
data <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE,
  sep = ",",
  fileEncoding = "latin1", # prueba latin1 para acentos
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

```
```
### Remove empty columns and rows
```{r}
data <- data[, colSums(!is.na(data)) > 0]
data <- data[rowSums(is.na(data)) < ncol(data), ]

```
```
### Create derived variables
```{r}
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
```
### Filter for analysis (only rows without NA in key variables)
```{r}
data_filtered <- data %>%
  filter(!is.na(BMI), !is.na(Age), !is.na(Percentil_formulas), !is.na(Age_group), !is.na(Lifestyle))

```
```
### --- Define comparisons for the statistical test ---
```{r}
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
```
### Create labels with total counts from the original dataset for Lifestyle and Percentil_group
```{r}
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
```
### --- Plot 1: BMI density by Lifestyle (filtered data, full legend) ---
```{r}
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
```
### --- Plot 2: Age density by Percentil_group ---
```{r}
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
```
### Calculate n by group and age (already done)
```{r}
group_counts <- data_filtered %>%
  count(Age_group_recode, Gender_Lifestyle)

```
```
### Position to place the texts (n) below the X axis
```{r}
y_pos_n <- min(data_filtered$BMI, na.rm = TRUE) - 2  # ajusta -2 según tu rango de BMI

```
```
### Create text label for n
```{r}
group_counts <- group_counts %>%
  mutate(n_label = paste0("n = ", n))

```
```
### Plot with geom_text for n
```{r}
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
```
### Add n below each bar
```{r}
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
```
### FINAL FIG. 1 COMPOSED
```{r}

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
```



 
