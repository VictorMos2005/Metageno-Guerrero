##FIG4. Average percent reads and unique order counts across classification methods


### GENERAR LOS GRAFICOS DE KRAKEN2

# Librerías necesarias
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1. Leer archivo Kraken2 ---
input_path <- "/home/alumno21/axel/files/kraken207tags.txt"  # Modifica si cambia
kraken2_data <- read.delim(input_path, header = FALSE, stringsAsFactors = FALSE)
colnames(kraken2_data) <- c("percent", "reads_clade", "reads_direct", "rank", "taxid", "name", "file")
kraken2_data$name <- trimws(kraken2_data$name)
kraken2_data$file_base <- basename(kraken2_data$file)

# --- 2. Mapear niveles taxonómicos ---
rank_map <- c(
  "U" = "unclassified", "R" = "root", "R1" = "cellular_organisms",
  "D" = "domain", "D1" = "supergroup", "D2" = "subgroup",
  "P" = "phylum", "C" = "class", "O" = "order",
  "F" = "family", "G" = "genus", "S" = "species", "S1" = "subspecies"
)


kraken2_data$tax_level <- rank_map[kraken2_data$rank]

# --- 5. Definir archivos urbanos y rurales ---
urban_files <- c("37082_2#1", "37082_1#20", "37082_1#27", "37035_2#13", "37035_1#29",
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
                 "37035_2#26", "36703_3#6", "37035_2#32", "36703_3#4", "37035_2#22")
rural_files <- c("37082_3#17", "37082_3#15", "37035_1#22", "36703_3#31", "37082_2#24",
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
                 "37035_7#18", "37035_7#23", "37035_7#25")
kraken2_data <- kraken2_data %>%
  mutate(Group = case_when(
    file_base %in% urban_files ~ "Urban",
    file_base %in% rural_files ~ "Rural",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))


# Instalar taxonomizr SOLO UNA VEZ
# install.packages("taxonomizr")  # Ya no es necesario si ya está instalado

library(taxonomizr)
library(dplyr)

# ----------------------------------
# 1. Crear base de datos SQLite con taxonomía
#    ⚠️ Solo necesitas hacerlo una vez
# ----------------------------------

sqlFile <- "/home/alumno21/axel/ncbi_taxonomy.sqlite"
nodesFile <- "/home/alumno21/axel/nodes.dmp"
namesFile <- "/home/alumno21/axel/names.dmp"

# ⚠️ Estos dos comandos SÓLO se usan si no has creado la base o quieres sobreescribirla:
# read.names.sql(namesFile, sqlFile)
# read.nodes.sql(nodesFile, sqlFile)

# Si ya hiciste estos pasos y el archivo .sqlite ya tiene las tablas, NO necesitas correrlos de nuevo.

# ----------------------------------
# 2. Obtener jerarquía taxonómica a partir de taxid
#    ✅ Esto sí debes correr cada vez que cambie tu archivo `kraken2_data`
# ----------------------------------

# Extraer taxids únicos del dataframe original
unique_taxids <- unique(kraken2_data$taxid)

# Obtener jerarquía para cada taxid
taxa_hierarchy <- getTaxonomy(
  unique_taxids,
  sqlFile = sqlFile,
  desiredTaxa = unname(rank_map)
)


# Convertir la jerarquía a data.frame y añadir taxid como columna
taxa_hierarchy_df <- as.data.frame(taxa_hierarchy)
taxa_hierarchy_df$taxid <- as.integer(rownames(taxa_hierarchy_df))

# ----------------------------------
# 3. Unir jerarquía con tu tabla original de Kraken2
# ----------------------------------

kraken2_data_full <- kraken2_data %>%
  left_join(taxa_hierarchy_df, by = "taxid")


# Calcular porcentaje de lecturas clasificadas por dominio y grupo
domain_percent <- kraken2_data_full %>%
  filter(!is.na(domain), domain != "", !is.na(Group)) %>%
  group_by(Group, domain) %>%
  summarise(total_reads = sum(reads_clade, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(percent = 100 * total_reads / sum(total_reads)) %>%
  ungroup() %>%
  mutate(label = paste0(round(percent, 1), "%"))


# (Opcional) Ver el resultado
#View(kraken2_data_full)

library(dplyr)

tax_levels <- c("domain", "phylum", "class", "subgroup", "order", "family", "genus", "species")

# Para cada nivel, calcular promedio de lecturas totales por grupo
avg_reads_level <- lapply(tax_levels, function(level) {
  kraken2_data_full %>%
    filter(!is.na(.data[[level]]), .data[[level]] != "", !is.na(Group)) %>%
    group_by(file_base, Group, taxon = .data[[level]]) %>%
    summarise(total_reads = sum(reads_clade, na.rm = TRUE), .groups = "drop") %>%
    group_by(file_base, Group) %>%
    summarise(total_reads_level = sum(total_reads), .groups = "drop") %>%
    group_by(Group) %>%
    summarise(avg_reads = mean(total_reads_level), .groups = "drop") %>%
    mutate(tax_level = level)
}) %>% bind_rows()

kraken_g <- ggplot(domain_percent, aes(x = percent, y = reorder(domain, percent), fill = Group)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7), color = "white") +
  
  # Labels inside if percent > 5, outside otherwise
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            hjust = ifelse(domain_percent$percent > 5, 1.1, -0.1),
            color = "black", size = 3) +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Domain - Classified Reads Percentage per Group (Kraken2)",
    x = "Classified Reads Percentage (%)",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 9),
    axis.text = element_text(face = "plain", size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    aspect.ratio = 1
  ) +
  
  xlim(0, max(domain_percent$percent) + 5)

kraken_g

library(dplyr)
library(ggplot2)

# --- Datos % de lecturas clasificadas ---
# domain_percent: columnas Group, domain, percent, label

# --- Calcular microorganismos únicos a nivel orden ---
unique_orders <- kraken2_data_full %>%
  filter(!is.na(domain), domain != "", !is.na(order), order != "") %>%
  group_by(Group, domain) %>%
  summarise(unique_orders = n_distinct(order), .groups = "drop")

# Crear etiquetas para el gráfico
unique_orders <- unique_orders %>%
  mutate(label = as.character(unique_orders))

# --- Gráfico de % de lecturas clasificadas (kraken_g) ---
kraken_g <- ggplot(domain_percent, aes(x = percent, y = reorder(domain, percent), fill = Group)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7), color = "white") +
  
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            hjust = ifelse(domain_percent$percent > 5, 1.1, -0.1),
            color = "black", size = 3) +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Domain - Classified Reads Percentage per Group (Kraken2)",
    x = "Classified Reads Percentage (%)",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 9),
    axis.text = element_text(face = "plain", size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    aspect.ratio = 1
  ) +
  
  xlim(0, max(domain_percent$percent) + 5)

# --- Gráfico de microorganismos únicos a nivel orden ---
unique_orders_g <- ggplot(unique_orders, aes(x = unique_orders, y = reorder(domain, unique_orders), fill = Group)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7), color = "white") +
  
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            hjust = ifelse(unique_orders$unique_orders > 10, 1.1, -0.1),  # Ajusta umbral para etiqueta dentro/fuera
            color = "black", size = 3) +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Domain - Unique Orders Count per Group (Kraken2)",
    x = "Unique Orders Count",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 9),
    axis.text = element_text(face = "plain", size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    aspect.ratio = 1
  ) +
  
  xlim(0, max(unique_orders$unique_orders) + 5)

# Mostrar ambos gráficos
kraken_g
unique_orders_g

library(dplyr)
library(ggplot2)
library(tidyr)


# 1. Preparar dataset para porcentaje de lecturas
domain_percent2 <- domain_percent %>%
  rename(domain = domain, value = percent) %>%  # Ajusta 'percent' si es necesario
  mutate(Metric = "Average Percent Reads")

# 2. Preparar dataset para número de órdenes únicos
unique_orders2 <- unique_orders %>%
  rename(domain = domain, value = unique_orders) %>%
  mutate(Metric = "Unique Order Count")

# 3. Unir datasets
combined_data <- bind_rows(
  domain_percent2 %>% select(Group, domain, value, Metric),
  unique_orders2 %>% select(Group, domain, value, Metric)
)

# 4. Controlar orden de factores
combined_data$domain <- factor(combined_data$domain, levels = unique(domain_percent$domain))
combined_data$Metric <- factor(combined_data$Metric, levels = c("Average Percent Reads", "Unique Order Count"))

# 5. Gráfico corregido para evitar que esté "pegado"
combined_g <- ggplot(combined_data, aes(x = value, y = domain, fill = Group)) +
  geom_col(position = position_dodge2(width = 1.2, preserve = "single", padding = 0.25), color = "white", width = 0.5) +
  
  facet_wrap(~ Metric, scales = "free_x") +
  
  geom_text(
    aes(label = round(value, 1), group = Group),   # ← agrega group=Group
    position = position_dodge2(width = 1.2, preserve = "single", padding = 0.25),  # ← dodge estable
    hjust = -0.2,
    color = "black", size = 3
  )+
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  coord_cartesian(clip = "off") +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Percentage of domain reads and unique genera count (Kraken2)",
    x = "Value",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    aspect.ratio = NULL
  )
#
#
#
#
#FINAL ACTUALIZADO PARA BARRAS EN KRAKEN##
##########################################
combined_g
##########################################
























# --- Librerías ---
library(dplyr)
library(tidyr)
library(ggplot2)

# --- Cargar datos Kaiju ---
input_path <- "/home/alumno21/axel/files/207_kaijusG.txt"
kaiju_data <- read.table(input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "M")

# Asegurar que reads sea numérico
kaiju_data$reads <- as.numeric(kaiju_data$reads)

# --- Separar taxonomía ---
kaiju_data <- kaiju_data %>%
  separate(taxon_name, into = c("Organism", "Domain", "Supergroup", "Kingdom", 
                                "Phylum", "Class", "Subclass", "Order", 
                                "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "drop")

# Reemplazar NA o vacíos en Order y Genus por "Unclassified"
kaiju_data$Order[is.na(kaiju_data$Order) | kaiju_data$Order == ""] <- "Unclassified"
kaiju_data$Genus[is.na(kaiju_data$Genus) | kaiju_data$Genus == ""] <- "Unclassified"

# Extraer nombre base del archivo
kaiju_data <- kaiju_data %>% mutate(file_base = gsub("^.*/", "", file))

# --- Listas Urban y Rural adaptadas con sufijo "_kaiju.out" ---
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

# --- Asignar grupo ---
kaiju_data <- kaiju_data %>%
  mutate(Group = case_when(
    file_base %in% urban_files_kaiju ~ "Urban",
    file_base %in% rural_files_kaiju ~ "Rural",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))

# --- Filtrar dominios válidos (remover NA o vacíos) ---
kaiju_data_filtered <- kaiju_data %>%
  filter(!is.na(Domain), Domain != "")

# --- Calcular abundancia relativa por archivo y dominio ---
# Sumar reads por archivo y dominio
reads_by_file_domain <- kaiju_data_filtered %>%
  group_by(file_base, Group, Domain) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE), .groups = "drop")

# Calcular total reads por archivo (para normalizar)
total_reads_by_file <- reads_by_file_domain %>%
  group_by(file_base) %>%
  summarise(file_reads = sum(total_reads), .groups = "drop")

# Añadir total_reads al dataframe por archivo
reads_norm <- reads_by_file_domain %>%
  left_join(total_reads_by_file, by = "file_base") %>%
  mutate(percent = 100 * total_reads / file_reads)

# --- Promediar abundancia relativa por grupo y dominio ---
domain_avg <- reads_norm %>%
  group_by(Group, Domain) %>%
  summarise(avg_percent = mean(percent), .groups = "drop") %>%
  mutate(label = paste0(round(avg_percent, 1), "%"))

# --- Graficar ---
kaiju_g <-ggplot(domain_avg, aes(x = avg_percent, y = reorder(Domain, avg_percent), fill = Group)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7), color = "white") +
  
  # Etiquetas: dentro si hay suficiente espacio, fuera si es pequeño
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            hjust = ifelse(domain_avg$avg_percent > 5, 1.1, -0.1),
            color = "black", size = 3) +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Average Domain Composition per Group (Kaiju)",
    x = "Average percentage",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 9),
    axis.text = element_text(face = "plain", size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    aspect.ratio = 1
  ) +
  
  xlim(0, max(domain_avg$avg_percent) + 5)

unique_microbes_kaiju <- kaiju_data_filtered %>%
  filter(!is.na(Domain), Domain != "", !is.na(Order), Order != "") %>%
  group_by(Group, Domain) %>%
  summarise(unique_genera = n_distinct(Order), .groups = "drop") %>%
  mutate(label = as.character(unique_genera))
unique_kaiju_g <- ggplot(unique_microbes_kaiju, aes(x = unique_genera, y = reorder(Domain, unique_genera), fill = Group)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7), color = "white") +
  
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            hjust = ifelse(unique_microbes_kaiju$unique_genera > 20, 1.1, -0.1),  # Ajusta umbral para etiqueta dentro/fuera
            color = "black", size = 3) +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Unique Order Count per Domain and Group (Kaiju)",
    x = "Unique Order Count",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 9),
    axis.text = element_text(face = "plain", size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    aspect.ratio = 1
  ) +
  
  xlim(0, max(unique_microbes_kaiju$unique_genera) + 5)
kaiju_g  # Gráfico original de abundancia relativa promedio

unique_kaiju_g  # Nuevo gráfico de microorganismos únicos



##Para unir ambas graficas

# --- Datos de abundancia relativa promedio (ya calculados) ---
domain_avg2 <- domain_avg %>%
  rename(Domain = Domain) %>%  # por si acaso
  mutate(value = avg_percent,
         Metric = "Average Percent Reads") %>%
  select(Group, Domain, value, Metric)

# --- Calcular géneros únicos por dominio y grupo ---
unique_microbes_kaiju <- kaiju_data_filtered %>%
  filter(!is.na(Domain), Domain != "", !is.na(Order), Order != "") %>%
  group_by(Group, Domain) %>%
  summarise(unique_genera = n_distinct(Order), .groups = "drop") %>%
  mutate(label = as.character(unique_genera))

unique_microbes_kaiju2 <- unique_microbes_kaiju %>%
  mutate(value = unique_genera,
         Metric = "Unique Order Count") %>%
  select(Group, Domain, value, Metric)

# --- Unir ambos datasets ---
combined_kaiju <- bind_rows(domain_avg2, unique_microbes_kaiju2)

# Para controlar el orden de los factores en y
combined_kaiju$Domain <- factor(combined_kaiju$Domain, levels = unique(domain_avg$Domain))
combined_kaiju$Metric <- factor(combined_kaiju$Metric, levels = c("Average Percent Reads", "Unique Order Count"))

# --- Gráfico combinado ---
combined_kaiju_g <- ggplot(combined_kaiju, aes(x = value, y = Domain, fill = Group)) +
  geom_col(position = position_dodge2(width = 1.2, preserve = "single", padding = 0.25), color = "white", width = 0.6) +
  
  facet_wrap(~ Metric, scales = "free_x") +
  
  geom_text(
    aes(label = round(value, 1), group = Group),   # ← agrega group=Group
    position = position_dodge2(width = 1.2, preserve = "single", padding = 0.25),  # ← dodge estable
    hjust = -0.2,
    color = "black", size = 3
  )+
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  coord_cartesian(clip = "off") +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Percentage of domain reads and unique genera count (Kaiju)",
    x = "Value",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.title = element_text(face = "plain", size = 10),
    axis.text = element_text(face = "plain", size = 9),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    aspect.ratio = NULL
  )

########################################
########################################
#GRAFICA FINAL DE STACKED BARS DE KAIJU 
combined_kaiju_g
########################################







### Kaiju-merged


################################################################################
################################################################################
#KAIJU COMPLEMENTADO CON KRAKEN2
################################################################################

# --- Librerías ---
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 0. Estandarizar nombres de columnas y archivos ---
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

# Asignar columnas numéricas clave
kraken2_aligned$reads <- as.numeric(kraken2_data_full$reads_clade)
kraken2_aligned$percent <- as.numeric(kraken2_data_full$percent)

kaiju_data$percent <- as.numeric(kaiju_data$percent)
kaiju_data$reads   <- as.numeric(kaiju_data$reads)

# Establecer nombres de archivo comunes
cat("\n--- Estandarizando nombres de archivo ---\n")
kaiju_data$file <- gsub("_kaiju\\.out$", "", kaiju_data$file)
kraken2_aligned$file <- gsub("\\.report$", "", kraken2_aligned$file)

# Confirmar cantidad de archivos únicos
cat("\nArchivos únicos antes de unir:\n")
cat("  Kaiju:   ", length(unique(kaiju_data$file)), "\n")
cat("  Kraken2: ", length(unique(kraken2_aligned$file)), "\n")

# --- 1. Definir jerarquía taxonómica ---
tax_levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Supergroup", "Domain")

# --- 2. Función para detectar el nivel taxonómico más específico ---
get_deepest_tax <- function(row) {
  for (level in tax_levels) {
    value <- row[[level]]
    if (!is.na(value) && value != "") {
      return(paste0(level, ":", value))
    }
  }
  return(NA_character_)
}

# --- 3. Aplicar función a ambos datasets ---
cat("\nAplicando nivel taxonómico más específico...\n")
kaiju_data$deepest_tax <- apply(kaiju_data[, tax_levels], 1, get_deepest_tax)
kraken2_aligned$deepest_tax <- apply(kraken2_aligned[, tax_levels], 1, get_deepest_tax)

# Verificar ejemplo
cat("  Ejemplo Kaiju:     ", head(kaiju_data$deepest_tax, 1), "\n")
cat("  Ejemplo Kraken2:   ", head(kraken2_aligned$deepest_tax, 1), "\n")

# --- 4. Crear claves únicas: file + deepest_tax ---
kaiju_keys <- paste(kaiju_data$file, kaiju_data$deepest_tax)
kraken2_keys <- paste(kraken2_aligned$file, kraken2_aligned$deepest_tax)

cat("\nLongitudes de claves únicas:\n")
cat("  Kaiju:     ", length(unique(kaiju_keys)), "\n")
cat("  Kraken2:   ", length(unique(kraken2_keys)), "\n")

# --- 5. Filtrar Kraken2: solo lo que NO está en Kaiju ---
kraken2_filtered <- kraken2_aligned[!(kraken2_keys %in% kaiju_keys), ]

cat("\nLecturas después de filtrar Kraken2:\n")
cat("  Kraken2 original: ", nrow(kraken2_aligned), "\n")
cat("  Kraken2 filtrado: ", nrow(kraken2_filtered), "\n")

# --- 6. Unir datos sin duplicados lógicos ---
cat("\n--- Uniendo datasets (Kaiju + Kraken2 filtrado) ---\n")
kaiju_merged <- bind_rows(kaiju_data, kraken2_filtered)

# --- 7. Validaciones finales ---
cat("\nValidaciones finales:\n")
cat("  Kaiju        - Total reads: ", sum(kaiju_data$reads, na.rm = TRUE), "\n")
cat("  Kraken2 FILT - Total reads: ", sum(kraken2_filtered$reads, na.rm = TRUE), "\n")
cat("  Fusionado    - Total reads: ", sum(kaiju_merged$reads, na.rm = TRUE), "\n")

cat("  Kaiju        - Nrows: ", nrow(kaiju_data), "\n")
cat("  Kraken2 FILT - Nrows: ", nrow(kraken2_filtered), "\n")
cat("  Fusionado    - Nrows: ", nrow(kaiju_merged), "\n")

cat("  Archivos fusionados: ", length(unique(kaiju_merged$file)), "\n\n")


library(dplyr)

# --- Renombrar columnas Kraken2 ---
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

# --- Usamos file_base de Kraken2 como referencia ---
kraken2_aligned$file <- kraken2_data_full$file_base

# --- Extraer file_base desde Kaiju ---
kaiju_data$file <- gsub("_kaiju\\.out$", "", basename(kaiju_data$file))
kaiju_data$file <- gsub("\\.out$", "", kaiju_data$file)  # por si hay archivos erróneos

# Validar si coinciden los códigos
cat("\n--- Comparación de archivos únicos ---\n")
cat("Kaiju archivos únicos:   ", length(unique(kaiju_data$file)), "\n")
cat("Kraken2 archivos únicos: ", length(unique(kraken2_aligned$file)), "\n")
cat("Coincidencias exactas:   ", length(intersect(unique(kaiju_data$file), unique(kraken2_aligned$file))), "\n")
cat("Diferencias:             ", length(setdiff(unique(kraken2_aligned$file), unique(kaiju_data$file))), "\n")

# --- Establecer niveles taxonómicos ---
tax_levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Supergroup", "Domain")

# --- Detectar nivel más específico ---
get_deepest_tax <- function(row) {
  for (level in tax_levels) {
    value <- row[[level]]
    if (!is.na(value) && value != "") return(paste0(level, ":", value))
  }
  return(NA_character_)
}

kaiju_data$deepest_tax <- apply(kaiju_data[, tax_levels], 1, get_deepest_tax)
kraken2_aligned$deepest_tax <- apply(kraken2_aligned[, tax_levels], 1, get_deepest_tax)

# --- Claves únicas ---
kaiju_keys   <- paste(kaiju_data$file, kaiju_data$deepest_tax)
kraken2_keys <- paste(kraken2_aligned$file, kraken2_aligned$deepest_tax)

# --- Conversión de valores numéricos ---
kaiju_data$percent   <- as.numeric(kaiju_data$percent)
kaiju_data$reads     <- as.numeric(kaiju_data$reads)
kraken2_aligned$reads   <- as.numeric(kraken2_aligned$reads_clade)
kraken2_aligned$percent <- as.numeric(kraken2_aligned$percent)

# --- Filtrar duplicados ---
kraken2_filtered <- kraken2_aligned[!(kraken2_keys %in% kaiju_keys), ]

# --- Unir ---
kaiju_merged.gz <- bind_rows(kaiju_data, kraken2_filtered)

# --- Asignar file_base y Group al objeto fusionado ---
kaiju_merged <- kaiju_merged %>%
  mutate(
    file_base = gsub("_kaiju\\.out$", "", basename(file)),  # extraer nombre base
    Group = case_when(
      file_base %in% rural_files_kaiju ~ "Rural",
      file_base %in% urban_files_kaiju ~ "Urban",
      TRUE ~ NA_character_
    )
  )

# Asegurar que reads sea numérico
kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

# Reemplazar NA o vacíos en Order y Genus por "Unclassified"
kaiju_merged$Order[is.na(kaiju_merged$Order) | kaiju_merged$Order == ""] <- "Unclassified"
kaiju_merged$Genus[is.na(kaiju_merged$Genus) | kaiju_merged$Genus == ""] <- "Unclassified"


# --- Listas Urban y Rural adaptadas con sufijo "_kaiju.out" ---
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

# --- Asignar correctamente file_base y Group ---
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


# Guardar como CSV
write.csv(kaiju_merged, "kaiju_merged_final.csv", row.names = FALSE)

# Guardar como TXT (tabulado)
write.table(kaiju_merged, "kaiju_merged_final.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)




# --- Filtrar dominios válidos (remover NA o vacíos) ---
kaiju_merged_filtered <- kaiju_merged %>%
  filter(!is.na(Domain), Domain != "")

# --- Calcular abundancia relativa por archivo y dominio ---
# Sumar reads por archivo y dominio
reads_by_file_domain <- kaiju_merged_filtered %>%
  group_by(file_base, Group, Domain) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE), .groups = "drop")

# Calcular total reads por archivo (para normalizar)
total_reads_by_file <- reads_by_file_domain %>%
  group_by(file_base) %>%
  summarise(file_reads = sum(total_reads), .groups = "drop")

# Añadir total_reads al dataframe por archivo
reads_norm <- reads_by_file_domain %>%
  left_join(total_reads_by_file, by = "file_base") %>%
  mutate(percent = 100 * total_reads / file_reads)

# --- Promediar abundancia relativa por grupo y dominio ---
domain_avg <- reads_norm %>%
  group_by(Group, Domain) %>%
  summarise(avg_percent = mean(percent), .groups = "drop") %>%
  mutate(label = paste0(round(avg_percent, 1), "%"))

# --- Graficar ---
kaiju_merged_g <-ggplot(domain_avg, aes(x = avg_percent, y = reorder(Domain, avg_percent), fill = Group)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7), color = "white") +
  
  # Etiquetas: dentro si hay suficiente espacio, fuera si es pequeño
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            hjust = ifelse(domain_avg$avg_percent > 5, 1.1, -0.1),
            color = "black", size = 3) +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Average Domain Composition per Group (Kaiju merged)",
    x = "Average percentage",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 9),
    axis.text = element_text(face = "plain", size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    aspect.ratio = 1
  ) +
  
  xlim(0, max(domain_avg$avg_percent) + 5)
unique_microbes_kaiju_merged <- kaiju_merged_filtered %>%
  filter(!is.na(Domain), Domain != "", !is.na(Order), Order != "") %>%
  group_by(Group, Domain) %>%
  summarise(unique_genera = n_distinct(Order), .groups = "drop") %>%
  mutate(label = as.character(unique_genera))
unique_kaiju_merged_g <- ggplot(unique_microbes_kaiju_merged, aes(x = unique_genera, y = reorder(Domain, unique_genera), fill = Group)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7), color = "white") +
  
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            hjust = ifelse(unique_microbes_kaiju_merged$unique_genera > 20, 1.1, -0.1),  # Ajusta umbral para etiqueta dentro/fuera
            color = "black", size = 3) +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Unique Order Count per Domain and Group (Kaiju merged)",
    x = "Unique Order Count",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 9),
    axis.text = element_text(face = "plain", size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    aspect.ratio = 1
  ) +
  
  xlim(0, max(unique_microbes_kaiju_merged$unique_genera) + 5)
kaiju_merged_g  # Gráfico original de abundancia relativa promedio DE BARRAS

unique_kaiju_merged_g  # Nuevo gráfico de microorganismos únicos DE BARRAS

library(dplyr)
library(ggplot2)

# --- Datos de abundancia relativa promedio (ya calculados) ---
domain_avg2 <- domain_avg %>%
  rename(Domain = Domain) %>%  # por si acaso
  mutate(value = avg_percent,
         Metric = "Percent Reads") %>%
  select(Group, Domain, value, Metric)

# --- Calcular géneros únicos por dominio y grupo ---
unique_microbes_kaiju_merged <- kaiju_merged_filtered %>%
  filter(!is.na(Domain), Domain != "", !is.na(Order), Order != "") %>%
  group_by(Group, Domain) %>%
  summarise(unique_genera = n_distinct(Order), .groups = "drop") %>%
  mutate(label = as.character(unique_genera))

unique_microbes_kaiju_merged2 <- unique_microbes_kaiju_merged %>%
  mutate(value = unique_genera,
         Metric = "Unique Order") %>%
  select(Group, Domain, value, Metric)

# --- Unir ambos datasets ---
combined_kaiju_merged <- bind_rows(domain_avg2, unique_microbes_kaiju_merged2)

# Para controlar el orden de los factores en y
# Ordenar Domain según el promedio de value (mayor a menor)
combined_kaiju_merged <- combined_kaiju_merged %>%
  group_by(Domain) %>%
  mutate(mean_value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Domain = reorder(Domain, mean_value))
combined_kaiju_merged$Metric <- factor(combined_kaiju_merged$Metric, levels = c("Percent Reads", "Unique Order"))

# --- Gráfico combinado ---
combined_kaiju_merged_g <- ggplot(combined_kaiju_merged, aes(x = value, y = Domain, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), color = "white", width = 0.6) +
  
  facet_wrap(~ Metric, scales = "free_x") +
  
  geom_text(aes(label = round(value, 1)),
            position = position_dodge(width = 0.8),
            hjust = -0.1,
            color = "black", size = 3,
            check_overlap = TRUE )+
  
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) + 
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Percentage of domain reads and unique genera count (Kaiju merged)",
    x = "Value",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.title = element_text(face = "plain", size = 10),
    axis.text = element_text(face = "plain", size = 9),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    aspect.ratio = 1
  )

#combined_kaiju_merged_g


# --- Datos de abundancia relativa promedio (ya calculados) ---
domain_avg2 <- domain_avg %>%
  rename(Domain = Domain) %>%  # por si acaso
  mutate(value = avg_percent,
         Metric = "Average Percent Reads") %>%
  select(Group, Domain, value, Metric)

# --- Calcular géneros únicos por dominio y grupo ---
unique_microbes_kaiju_merged <- kaiju_merged_filtered %>%
  filter(!is.na(Domain), Domain != "", !is.na(Order), Order != "") %>%
  group_by(Group, Domain) %>%
  summarise(unique_genera = n_distinct(Order), .groups = "drop") %>%
  mutate(label = as.character(unique_genera))

unique_microbes_kaiju_merged2 <- unique_microbes_kaiju_merged %>%
  mutate(value = unique_genera,
         Metric = "Unique Order Count") %>%
  select(Group, Domain, value, Metric)

# --- Unir ambos datasets ---
combined_kaiju_merged <- bind_rows(domain_avg2, unique_microbes_kaiju_merged2)

# Para controlar el orden de los factores en y
combined_kaiju_merged$Domain <- factor(combined_kaiju_merged$Domain, levels = unique(domain_avg$Domain))
combined_kaiju_merged$Metric <- factor(combined_kaiju_merged$Metric, levels = c("Average Percent Reads", "Unique Order Count"))

# --- Gráfico combinado ---
combined_kaiju_merged_g <- ggplot(combined_kaiju_merged, aes(x = value, y = Domain, fill = Group)) +
  geom_col(position = position_dodge2(width = 1.2, preserve = "single", padding = 0.25), color = "white", width = 0.6) +
  
  facet_wrap(~ Metric, scales = "free_x") +
  
  geom_text(
    aes(label = round(value, 1), group = Group),   # ← agrega group=Group
    position = position_dodge2(width = 1.2, preserve = "single", padding = 0.25),  # ← dodge estable
    hjust = -0.2,
    color = "black", size = 3
  )+
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  coord_cartesian(clip = "off") +
  
  scale_fill_manual(values = c("Urban" = "#0072B2", "Rural" = "#E69F00")) +
  
  labs(
    title = "Percentage of domain reads and unique genera count (Kaiju merged)",
    x = "Value",
    y = "Domain",
    fill = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.title = element_text(face = "plain", size = 10),
    axis.text = element_text(face = "plain", size = 9),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    aspect.ratio = NULL
  )

########################################
########################################
#GRAFICA FINAL DE STACKED BARS DE KAIJU MERGED
combined_kaiju_merged_g
########################################


###############################################################################
###############################################################################
###############################################################################
# Grafica final de comparacion de bases de datos
###############################################################################


# 1) Unificar datasets y estandarizar nombres/orden
kaiju_df <- combined_kaiju %>%
  mutate(Database = "Kaiju") %>%
  rename(domain = Domain)

kraken_df <- combined_data %>%
  mutate(Database = "Kraken2")  # <- este ya tiene columnas Group, domain, value, Metric

ftp_df <- combined_kaiju_merged %>%
  mutate(Database = "FTP") %>%               # <- renombramos aquí (nunca "Kaiju merged")
  rename(domain = Domain)

all_df <- bind_rows(kaiju_df, kraken_df, ftp_df) %>%
  mutate(
    # asegurar que jamás salga "Kaiju merged" en ninguna columna de texto
    across(where(is.character), ~ gsub("Kaiju merged", "FTP", .)),
    domain    = factor(domain, levels = c("Eukaryota","Bacteria","Archaea")),
    Group     = factor(Group,  levels = c("Rural","Urban")),
    Metric    = factor(Metric, levels = c("Average Percent Reads","Unique Order Count")),
    Database  = factor(Database, levels = c("Kaiju","Kraken2","FTP"))
  )

# 2) Función de etiqueta según métrica
lab_value <- function(val, metric) {
  ifelse(metric == "Average Percent Reads",
         sprintf("%.1f", val),                                 # 1 decimal
         label_number(accuracy = 1, big.mark = ",")(val))      # enteros con separador de miles
}

# 3) Plot tipo “lollipop” + puntos con posición en paralelo
p_better <- ggplot(all_df,
                   aes(x = value, y = domain, color = Group)) +
  # línea fina desde 0 hasta el valor (ancla visual)
  geom_segment(aes(x = 0, xend = value, yend = domain),
               linewidth = 0.25, alpha = 0.25, lineend = "round") +
  # puntos (Rural/Urban en paralelo) con shape por Database
  geom_point(aes(shape = Database),
             position = position_dodge2(width = 0.35, preserve = "single"),
             size = 2.6, stroke = 0.2) +
  # etiquetas a la derecha del punto
  geom_text(
    aes(
      # mueve las etiquetas a la derecha según la métrica
      x = value + ifelse(Metric == "Average Percent Reads", 2, 50),
      label = lab_value(value, Metric)
    ),
    position = position_dodge2(width = 0.35, preserve = "single"),
    size = 3, color = "black", hjust = 0
  ) +
  facet_grid(Database ~ Metric, scales = "free_x") +
  # colores exactos solicitados
  scale_color_manual(values = c(Rural = "#E69F00", Urban = "#0072B2"), name = "Group") +
  scale_shape_manual(values = c(Kaiju = 16, Kraken2 = 17, FTP = 8), name = "Database") +
  # más espacio a los lados (evita que se peguen a los márgenes)
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(x = "Value", y = "Domain") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left   = element_line(),
    strip.text = element_text(size = 10, face = "plain"),
    # ocultar los rótulos laterales de Database (se queda solo la leyenda)
    strip.text.y = element_blank(),
    strip.background.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "plain"),
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 9),
    panel.spacing.x = unit(16, "pt"),
    panel.spacing.y = unit(10, "pt"),
    plot.margin = margin(10, 28, 10, 10)  # margen derecho mayor para etiquetas
  )
#####################TABLA FINAL DE COMPARACIONES 2 
p_better




