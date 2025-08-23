install.packages("dplyr")
install.packages("tidyr")
install.packages("tibble")
install.packages("boot")
install.packages("stringr")

# --- Librerías necesarias ---
library(dplyr)
library(tidyr)
library(tibble)
library(boot)
library(stringr)



kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

# Si no existen columnas taxonómicas separadas pero sí 'taxon_name'
if (!all(c("Domain","Genus") %in% names(kaiju_merged)) && "taxon_name" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    separate(
      taxon_name,
      into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
               "Subclass","Order","Family","Genus","Species"),
      sep = ";", fill = "right", extra = "drop"
    )
}

# Rellena vacíos mínimos
kaiju_merged <- kaiju_merged %>%
  mutate(across(c(Genus, Family, Order), ~ ifelse(is.na(.) | . == "", "Unclassified", .)))

# file_base si hace falta
if (!"file_base" %in% names(kaiju_merged)) {
  if ("file" %in% names(kaiju_merged)) {
    kaiju_merged <- kaiju_merged %>% mutate(file_base = gsub("^.*/", "", file))
  } else {
    stop("No 'file_base' or 'file' column found in kaiju_merged.")
  }
}


data <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE,
  sep = ",",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

# Limpieza otra vez
data <- data[, colSums(!is.na(data)) > 0]
data <- data[rowSums(is.na(data)) < ncol(data), ]

# Variables numéricas
data$BMI <- as.numeric(data$BMI)
data$Age <- as.numeric(data$Age)

# Crear la columna BMI_group (como ya lo tenías)
data <- data %>%
  mutate(
    BMI_group = case_when(
      Percentil_formulas %in% c("Overweight", "Obesity") ~ "Overweight/Obesity",
      Percentil_formulas == "Normal Weight" ~ "Normal Weight",
      Percentil_formulas == "Malnutrition" ~ "Malnutrition",
      TRUE ~ NA_character_
    ),
    Lifestyle_num = ifelse(Lifestyle == "Urban", 1, 0)
  )

# Ahora sí, crear la nueva variable agrupada
data <- data %>%
  mutate(
    BMI_index_group = case_when(
      BMI_group == "Overweight/Obesity" ~ ">25 BMI INDEX",
      BMI_group %in% c("Normal Weight", "Malnutrition") ~ "<25 BMI INDEX",
      TRUE ~ NA_character_
    )
  )

# --- Definir categorías ---
data <- data %>%
  mutate(
    BMI_index_group = case_when(
      BMI_group %in% c("Overweight/Obesity") ~ ">25 BMI INDEX",
      BMI_group %in% c("Normal Weight", "Malnutrition") ~ "<25 BMI INDEX",
      TRUE ~ NA_character_
    )
  )

# --- Extraer listas de códigos 'file' ---
files_over25 <- data %>%
  filter(BMI_index_group == ">25 BMI INDEX") %>%
  pull(Lane) %>%
  unique()

files_under25 <- data %>%
  filter(BMI_index_group == "<25 BMI INDEX") %>%
  pull(Lane) %>%
  unique()


# --- Asignar grupos de BMI ---
kaiju_merged <- kaiju_merged %>%
  mutate(BMI_group = case_when(
    file_base %in% paste0(files_over25, "_kaiju.out")  ~ "BMI≥25",
    file_base %in% paste0(files_under25, "_kaiju.out") ~ "BMI<25",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(BMI_group))


# --- Asignar grupos (ajusta tus listas de archivos 'urban_files' y 'rural_files') ---
kaiju_merged <- kaiju_merged %>%
  mutate(Group = case_when(
    file_base %in% paste0(urban_files, "_kaiju.out") ~ "Urban",
    file_base %in% paste0(rural_files, "_kaiju.out") ~ "Rural",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))
cat("Número de muestras asignadas a grupos:\n"); print(table(kaiju_merged$Group))

# --- Filtrar SOLO Bacteria/Eukaryota y quedarnos con GENUS ---
full_data <- kaiju_merged %>%
  filter(Domain %in% c("Bacteria","Eukaryota")) %>%
  filter(!is.na(Genus) & Genus != "" & Genus != "Unclassified") %>%
  mutate(Genus = str_trim(Genus))

# ------------------------------------------------------------------
# 2) FILTRADO PREVIO DE GÉNEROS ABUNDANTES (acelera bootstrap)
#    (mismos umbrales para ambos dominios)
# ------------------------------------------------------------------
threshold_genus <- 8  # lecturas promedio mínimas

# Bacteria
genus_bact_urban <- full_data %>% filter(Group == "Urban",  Domain == "Bacteria") %>%
  group_by(Genus) %>% summarise(mean_reads = mean(reads), .groups = "drop") %>%
  filter(mean_reads > threshold_genus) %>% pull(Genus)

genus_bact_rural <- full_data %>% filter(Group == "Rural",  Domain == "Bacteria") %>%
  group_by(Genus) %>% summarise(mean_reads = mean(reads), .groups = "drop") %>%
  filter(mean_reads > threshold_genus) %>% pull(Genus)

common_genus_bact <- intersect(genus_bact_urban, genus_bact_rural)

# Eukaryota
genus_euk_urban <- full_data %>% filter(Group == "Urban",  Domain == "Eukaryota") %>%
  group_by(Genus) %>% summarise(mean_reads = mean(reads), .groups = "drop") %>%
  filter(mean_reads > threshold_genus) %>% pull(Genus)

genus_euk_rural <- full_data %>% filter(Group == "Rural",  Domain == "Eukaryota") %>%
  group_by(Genus) %>% summarise(mean_reads = mean(reads), .groups = "drop") %>%
  filter(mean_reads > threshold_genus) %>% pull(Genus)

common_genus_euk <- intersect(genus_euk_urban, genus_euk_rural)

cat("Géneros comunes tras filtro — Bacteria:", length(common_genus_bact),
    " | Eukaryota:", length(common_genus_euk), "\n")

# Mantener solo esos géneros abundantes por dominio
full_data_filtered <- full_data %>%
  filter((Domain == "Bacteria"  & Genus %in% common_genus_bact) |
           (Domain == "Eukaryota" & Genus %in% common_genus_euk))

# ------------------------------------------------------------------
# 3) Función para crear matriz (GENUS) log2-normalizada por muestra
# ------------------------------------------------------------------
crear_matriz <- function(df, taxon_col = "Genus") {
  df %>%
    group_by(file_base, !!sym(taxon_col)) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(taxon_col), values_from = reads, values_fill = 0) %>%
    column_to_rownames("file_base") %>%
    {
      mat <- sweep(., 1, rowSums(.), "/")   # proporción por muestra
      log2(mat + 1e-6)
    }
}

# Matrices por grupo (Bacteria + Eukaryota juntos, columnas = GENUS)
mat_urban <- full_data_filtered %>% filter(Group == "Urban") %>% crear_matriz("Genus")
mat_rural <- full_data_filtered %>% filter(Group == "Rural") %>% crear_matriz("Genus")

cat("Dimensiones matriz Urban (genus):", dim(mat_urban), "\n")
cat("Dimensiones matriz Rural (genus):", dim(mat_rural), "\n")

# Mapa GENUS -> DOMAIN
genus_domain_map <- full_data_filtered %>%
  distinct(Genus, Domain)

# ------------------------------------------------------------------
# 4) Bootstrap Spearman
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# 5) Correlaciones (todas a nivel GENUS) y clasificación de relación
# ------------------------------------------------------------------
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

# Ejecutar por grupo
cat("Calculando correlaciones (Urban)...\n")
edges_urban <- correlacion_spearman_boot(mat_urban, genus_domain_map)
cat("Calculando correlaciones (Rural)...\n")
edges_rural <- correlacion_spearman_boot(mat_rural, genus_domain_map)

# ------------------------------------------------------------------
# 6) Filtros por relación y exportación
# ------------------------------------------------------------------
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
      cat("Aviso: No hay correlaciones en", grupo_nombre, "-", rel, "\n")
      next
    }
    nodes_filtrados <- crear_nodos_de_edges(edges_filtradas, full_data_filtered)
    
    # Archivos de salida
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

cat("¡Listo! Redes (GENUS) exportadas por relación y grupo.\n")



# ===================== Ajustes (mucho más estrictos) =====================




suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(boot)
})
set.seed(1234)
# Filtros por grupo (prevalencia y % medio, NORMALIZADO POR DOMINIO)
MIN_PREV_BACT <- 0.20
MIN_PREV_EUK  <- 0.10
MIN_MEAN_BACT <- 0.06     # % dentro de Bacteria
MIN_MEAN_EUK  <- 0.015    # % dentro de Eukaryota
TOP_N_BACT    <- 250
TOP_N_EUK     <- 120

# Correlación / Bootstrap
BOOT_R               <- 300
PSEUDO               <- 1e-6
PRESCREEN_ABS_RHO    <- 0.35   # global (antes 0.20)
PRESCREEN_ABS_RHO_WL <- 0.20   # si toca whitelist
MIN_SAMPLES_PER_G    <- 8

# Co-ocurrencia mínima (muestras donde AMBOS > 0)
MIN_COOC_FRAC <- 0.35          # ≥35% de las muestras del grupo
MIN_COOC_ABS  <- 20            # y al menos 20 muestras

# Filtro final de aristas (solo BE)
THRESH_BE <- 0.45              # |rho_mean| mínimo
ALPHA_Q   <- 0.05              # BH FDR

# ===================== Whitelist de géneros clave =====================
force_include_euk  <- c("Saccharomyces")
force_include_bact <- character(0)

# ======================== Datos esperados en 'kaiju_merged' ========================
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

# Solo Bacteria/Eukaryota a nivel GENUS; excluye Homo
dat <- kaiju_merged %>%
  filter(Domain %in% c("Bacteria","Eukaryota")) %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus=="" | Genus=="Unclassified", NA, str_trim(Genus))) %>%
  filter(!is.na(Genus)) %>%
  filter(!(Domain=="Eukaryota" & Genus=="Homo"))

# ------------------------ % dentro del DOMINIO por muestra ------------------------
totals_domain <- dat %>%
  group_by(file_base, Group, Domain) %>%
  summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")

genus_reads <- dat %>%
  group_by(file_base, Group, Domain, Genus) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  left_join(totals_domain, by = c("file_base","Group","Domain")) %>%
  mutate(pct_domain = if_else(total_reads_domain > 0, 100 * reads / total_reads_domain, 0))

# ------------------------ prevalencia y medias por grupo -------------------
prev_mean <- genus_reads %>%
  group_by(Group, Domain, Genus) %>%
  summarise(prevalence = mean(reads > 0),
            mean_pct_dom = mean(pct_domain),
            .groups = "drop")

# Keepers por dominio + forzar whitelist
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

# ======================= matrices por grupo (counts + CLR) ======================
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

# Alinear columnas
common_cols <- Reduce(intersect, list(colnames(mat_clr_rural), colnames(mat_clr_urban)))
mat_clr_rural <- mat_clr_rural[, common_cols, drop = FALSE]
mat_clr_urban <- mat_clr_urban[, common_cols, drop = FALSE]
cnt_rural     <- cnt_rural[,     common_cols, drop = FALSE]
cnt_urban     <- cnt_urban[,     common_cols, drop = FALSE]

# Genus -> Domain map (solo comunes)
genus_domain_map <- dat_filt %>% distinct(Genus, Domain) %>% filter(Genus %in% common_cols)

# =================== Prescreen + co-ocurrencia + bootstrap =================
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
  
  # Co-ocurrencia
  pres_abs <- max(min_cooc_abs, ceiling(min_cooc_frac * nrow(counts_mat)))
  pres_mat <- counts_mat > 0
  
  # Spearman masivo = Pearson sobre rangos
  mat_rank <- apply(mat_clr, 2, rank, ties.method = "average")
  R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                  mat_rank[, euk_cols,  drop = FALSE],
                  method = "pearson", use = "pairwise.complete.obs")
  
  # Candidatos: prescreen global
  sel_global <- which(abs(R) >= prescreen_global, arr.ind = TRUE)
  pairs_df <- tibble(
    Genus_1 = bact_cols[sel_global[,"row"]],
    Genus_2 = euk_cols[ sel_global[,"col"]],
    from_whitelist = FALSE
  )
  
  # Candidatos extra: si toca whitelist, prescreen más laxo
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
  
  # Únicos + co-ocurrencia mínima
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
    # Consistencia de signo en bootstrap
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

# BH (todo es BE)
add_q <- function(df) if (nrow(df)) df %>% mutate(q = p.adjust(p, method = "BH")) else df
edges_rural <- add_q(edges_rural)
edges_urban <- add_q(edges_urban)

# --------- Filtro final muy estricto ----------
filter_BE <- function(df) {
  if (nrow(df)==0) return(df[0,])
  df %>%
    mutate(abs_rho = abs(rho_mean)) %>%
    # CI debe rebasar el umbral (no sólo mismo signo)
    mutate(ci_strong = (rho_mean > 0 & CI_lower >= THRESH_BE) |
             (rho_mean < 0 & CI_upper <= -THRESH_BE)) %>%
    filter(
      abs_rho >= THRESH_BE,
      abs(rho) >= THRESH_BE,             # efecto puntual también grande
      ci_strong,
      sign_consistency >= 0.90,          # ≥90% de réplicas con el mismo signo
      q <= ALPHA_Q
    ) %>%
    select(-abs_rho, -ci_strong)
}

edges_rural_f <- filter_BE(edges_rural)
edges_urban_f <- filter_BE(edges_urban)

message("BE edges kept | Rural: ", nrow(edges_rural_f), " | Urban: ", nrow(edges_urban_f))

# ============================= Nodos & export =============================
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

























## =========================================================
## Redes GENUS (Bacteria/Eukaryota) por:
##   1) Group (Urban vs Rural)
##   2) BMI_group (BMI≥25 vs BMI<25)
##    → En nombres de archivo: BMI<25 -> BMI24 ;  BMI≥25 -> BMI26
## =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(boot); library(readr)
})

set.seed(1234)

## ======================= PARÁMETROS GLOBALES =======================
# Filtros por grupo (prevalencia y % medio, NORMALIZADO POR DOMINIO)
MIN_PREV_BACT <- 0.20
MIN_PREV_EUK  <- 0.10
MIN_MEAN_BACT <- 0.06     # % dentro de Bacteria
MIN_MEAN_EUK  <- 0.015    # % dentro de Eukaryota
TOP_N_BACT    <- 250
TOP_N_EUK     <- 120

# Correlación / Bootstrap
BOOT_R               <- 300
PSEUDO               <- 1e-6
PRESCREEN_ABS_RHO    <- 0.35   # global
PRESCREEN_ABS_RHO_WL <- 0.20   # si toca whitelist
MIN_SAMPLES_PER_G    <- 8

# Co-ocurrencia mínima (muestras donde AMBOS > 0)
MIN_COOC_FRAC <- 0.35          # ≥35% de las muestras del grupo
MIN_COOC_ABS  <- 20            # y al menos 20 muestras

# Filtro final de aristas (solo BE)
THRESH_BE <- 0.45              # |rho_mean| mínimo
ALPHA_Q   <- 0.05              # BH FDR

# Whitelist de géneros clave
force_include_euk  <- c("Saccharomyces")
force_include_bact <- character(0)

## ======================= ENTRADAS =======================
# Se asume que ya cargaste 'kaiju_merged' en el ambiente.
# Si no, descomenta y pon tu ruta:
# kaiju_merged <- readr::read_tsv("tu_kaiju_merged.tsv", show_col_types = FALSE)

# Asegura tipos:
kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

# Si no existen columnas taxonómicas separadas pero sí 'taxon_name'
if (!all(c("Domain","Genus") %in% names(kaiju_merged)) && "taxon_name" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    tidyr::separate(
      taxon_name,
      into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
               "Subclass","Order","Family","Genus","Species"),
      sep = ";", fill = "right", extra = "drop"
    )
}

# file_base si hace falta
if (!"file_base" %in% names(kaiju_merged)) {
  if ("file" %in% names(kaiju_merged)) {
    kaiju_merged <- kaiju_merged %>% mutate(file_base = gsub("^.*/", "", file))
  } else {
    stop("No 'file_base' or 'file' column found in kaiju_merged.")
  }
}

## ============= METADATOS: construir BMI_group para BMI≥25 / BMI<25 =============
# Lee tu archivo de metadatos (ajusta la ruta si cambia)
meta <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE, sep = ",", fileEncoding = "latin1",
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

# Limpieza mínima
meta <- meta[, colSums(!is.na(meta)) > 0]
meta <- meta[rowSums(is.na(meta)) < ncol(meta), ]
meta$BMI <- suppressWarnings(as.numeric(meta$BMI))
meta$Age <- suppressWarnings(as.numeric(meta$Age))

# Agrupaciones de IMC (como lo tenías + colapsado a 2 niveles)
meta <- meta %>%
  mutate(
    BMI_group = case_when(
      Percentil_formulas %in% c("Overweight", "Obesity") ~ "Overweight/Obesity",
      Percentil_formulas == "Normal Weight" ~ "Normal Weight",
      Percentil_formulas == "Malnutrition" ~ "Malnutrition",
      TRUE ~ NA_character_
    ),
    BMI_index_group = case_when(
      BMI_group == "Overweight/Obesity" ~ ">25 BMI INDEX",
      BMI_group %in% c("Normal Weight", "Malnutrition") ~ "<25 BMI INDEX",
      TRUE ~ NA_character_
    )
  )

# En tus metadatos, el identificador de secuenciación parece estar en 'Lane'
# y los archivos de Kaiju acaban en "_kaiju.out". Ajusta si tu sufijo es distinto.
files_over25 <- meta %>% filter(BMI_index_group == ">25 BMI INDEX") %>% pull(Lane) %>% unique()
files_under25 <- meta %>% filter(BMI_index_group == "<25 BMI INDEX") %>% pull(Lane) %>% unique()

# Asignar BMI_group a kaiju_merged (etiquetas finales BMI≥25 / BMI<25)
kaiju_merged <- kaiju_merged %>%
  mutate(BMI_group = case_when(
    file_base %in% paste0(files_over25, "_kaiju.out")  ~ "BMI≥25",
    file_base %in% paste0(files_under25, "_kaiju.out") ~ "BMI<25",
    TRUE ~ NA_character_
  ))

## ============= (Opcional) asignar Group (Urban/Rural) si no existe =============
if (!"Group" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    mutate(Group = case_when(
      exists("urban_files") & file_base %in% paste0(urban_files, "_kaiju.out") ~ "Urban",
      exists("rural_files") & file_base %in% paste0(rural_files, "_kaiju.out") ~ "Rural",
      TRUE ~ NA_character_
    ))
}

# ================== 1) Lista de Lanes "Adult" ==================
# (meta debe estar cargado; ajusta el nombre de la columna si difiere)
adult_lanes <- meta %>%
  dplyr::filter(!is.na(Age_group), Age_group == "Adult") %>%
  dplyr::pull(Lane) %>%
  unique()

# si tus archivos Kaiju siempre terminan con este sufijo, déjalo así:
KAIJU_SUFFIX <- "_kaiju.out"
adult_files  <- paste0(adult_lanes, KAIJU_SUFFIX)

# ================== 2) Columna context_group en kaiju_merged ==================
# Requisitos previos:
# - kaiju_merged$file_base ya creado (o crea uno desde 'file')
# - kaiju_merged$BMI_group con valores "BMI≥25" / "BMI<25"
# - kaiju_merged$Group con valores "Urban" / "Rural"

kaiju_merged <- kaiju_merged %>%
  dplyr::mutate(
    # marcar si la fila corresponde a archivos de muestras Adult
    is_adult_file = file_base %in% adult_files,
    
    # NUEVA columna 'context_group' SOLO para Adult + BMI≥25
    context_group = dplyr::case_when(
      is_adult_file & BMI_group == "BMI≥25" & Group == "Urban" ~ "BMI≥25 Urban",
      is_adult_file & BMI_group == "BMI≥25" & Group == "Rural" ~ "BMI≥25 Rural",
      TRUE ~ NA_character_
    )
  )


## ======================= Helpers genéricos =======================
# Limpieza base (solo Bacteria/Eukaryota a nivel GENUS; excluye Homo)
prep_dat <- function(df_all) {
  df_all %>%
    filter(Domain %in% c("Bacteria","Eukaryota")) %>%
    mutate(Genus = ifelse(is.na(Genus) | Genus=="" | Genus=="Unclassified", NA, str_trim(Genus))) %>%
    filter(!is.na(Genus)) %>%
    filter(!(Domain=="Eukaryota" & Genus=="Homo"))
}

# Métricas por dominio dentro de cada muestra (% dentro del DOMINIO)
domain_norm_tables <- function(dat, group_var = "Group") {
  stopifnot(all(c("file_base", group_var, "Domain", "Genus", "reads") %in% names(dat)))
  totals_domain <- dat %>%
    group_by(.data[[group_var]], file_base, Domain) %>%
    summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")
  genus_reads <- dat %>%
    group_by(.data[[group_var]], file_base, Domain, Genus) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(totals_domain, by = c(group_var, "file_base", "Domain")) %>%
    mutate(pct_domain = if_else(total_reads_domain > 0, 100 * reads / total_reads_domain, 0))
  list(genus_reads = genus_reads, totals_domain = totals_domain)
}

# Elegir géneros “keepers” por dominio (prevalencia y % medio) usando el group_var
pick_keepers <- function(genus_reads, group_var = "Group") {
  prev_mean <- genus_reads %>%
    group_by(.data[[group_var]], Domain, Genus) %>%
    summarise(prevalence = mean(reads > 0),
              mean_pct_dom = mean(pct_domain),
              .groups = "drop")
  
  # NOTA: asumimos 2 niveles en group_var; se ordenan por aparición
  lvls <- sort(unique(prev_mean[[group_var]]))
  
  keepers_bact <- prev_mean %>%
    filter(Domain=="Bacteria") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_BACT,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_BACT),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    { if (nrow(.)==0) { message("⚠ Sin métricas para Bacteria en pick_keepers()"); . } else . } %>%
    {
      # Depuración: ¿cuántos fallan por prevalencia vs media?
      tmp <- .
      if (nrow(tmp)) {
        if (all(!tmp$ok)) message("⚠ Ningún género bacteriano pasó los filtros: MIN_PREV_BACT=", MIN_PREV_BACT, 
                                  " y/o MIN_MEAN_BACT=", MIN_MEAN_BACT)
      }
      tmp
    } %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_BACT) %>% pull(Genus)
  
  keepers_euk <- prev_mean %>%
    filter(Domain=="Eukaryota") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_EUK,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_EUK),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    { if (nrow(.)==0) { message("⚠ Sin métricas para Eukaryota en pick_keepers()"); . } else . } %>%
    {
      tmp <- .
      if (nrow(tmp)) {
        if (all(!tmp$ok)) message("⚠ Ningún género eucariota pasó los filtros: MIN_PREV_EUK=", MIN_PREV_EUK, 
                                  " y/o MIN_MEAN_EUK=", MIN_MEAN_EUK)
      }
      tmp
    } %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_EUK) %>% pull(Genus)
  
  keepers_bact <- union(keepers_bact,
                        intersect(force_include_bact, unique(genus_reads$Genus[genus_reads$Domain=="Bacteria"])))
  keepers_euk  <- union(keepers_euk,
                        intersect(force_include_euk,  unique(genus_reads$Genus[genus_reads$Domain=="Eukaryota"])))
  
  if (length(keepers_bact)==0) message("⚠ keepers_bact vacío tras aplicar TOP_N_BACT=", TOP_N_BACT)
  if (length(keepers_euk)==0)  message("⚠ keepers_euk vacío tras aplicar TOP_N_EUK=", TOP_N_EUK)
  
  list(keepers_bact = keepers_bact, keepers_euk = keepers_euk)
}

# Matrices de counts por grupo (según group_var) y CLR
counts_from_long <- function(df_long, group_var, group_value) {
  df_long %>%
    filter(.data[[group_var]] == group_value) %>%
    select(file_base, Genus, reads) %>%
    group_by(file_base, Genus) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Genus, values_from = reads, values_fill = 0) %>%
    column_to_rownames("file_base") %>% as.matrix()
}

clr_from_counts <- function(counts_mat) {
  if (is.null(counts_mat) || nrow(counts_mat) < MIN_SAMPLES_PER_G) return(NULL)
  prop <- sweep(counts_mat, 1, rowSums(counts_mat), "/"); prop[is.na(prop)] <- 0
  logx <- log(prop + PSEUDO); sweep(logx, 1, rowMeans(logx), "-")
}

# Correlaciones BE rápidas con bootstrap
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
  
  # Co-ocurrencia
  pres_abs <- max(min_cooc_abs, ceiling(min_cooc_frac * nrow(counts_mat)))
  pres_mat <- counts_mat > 0
  
  # Spearman masivo = Pearson sobre rangos
  mat_rank <- apply(mat_clr, 2, rank, ties.method = "average")
  R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                  mat_rank[, euk_cols,  drop = FALSE],
                  method = "pearson", use = "pairwise.complete.obs")
  
  # Candidatos: prescreen global
  sel_global <- which(abs(R) >= prescreen_global, arr.ind = TRUE)
  pairs_df <- tibble(
    Genus_1 = bact_cols[sel_global[,"row"]],
    Genus_2 = euk_cols[ sel_global[,"col"]],
    from_whitelist = FALSE
  )
  
  # Candidatos extra: whitelist (más laxo)
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
  
  # Únicos + co-ocurrencia mínima
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

# BH helper
add_q <- function(df) if (nrow(df)) df %>% mutate(q = p.adjust(p, method = "BH")) else df

# Filtro final “muy estricto”
filter_BE <- function(df) {
  if (nrow(df)==0) return(df[0,])
  df %>%
    mutate(abs_rho = abs(rho_mean)) %>%
    mutate(ci_strong = (rho_mean > 0 & CI_lower >= THRESH_BE) |
             (rho_mean < 0 & CI_upper <= -THRESH_BE)) %>%
    filter(
      abs_rho >= THRESH_BE,
      abs(rho) >= THRESH_BE,
      ci_strong,
      sign_consistency >= 0.90,
      q <= ALPHA_Q
    ) %>%
    select(-abs_rho, -ci_strong)
}

make_nodes <- function(edges_df, genus_domain_map) {
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

# --- Mapea etiquetas a tags seguros solo para los nombres de archivos ---
label_to_tag <- function(lbl, group_var = NULL) {
  # Caso especial: BMI_group -> forzar distinción clara en nombres
  if (!is.null(group_var) && group_var == "BMI_group") {
    if (grepl("^\\s*BMI\\s*<\\s*25\\s*$", lbl))  return("BMI24")
    if (grepl("^\\s*BMI\\s*[≥>=]\\s*25\\s*$", lbl)) return("BMI26")
  }
  # Fallback general: limpiar caracteres no alfanuméricos
  out <- gsub("[^A-Za-z0-9]+", "", lbl)
  if (nzchar(out)) out else "grp"
}

## ======================= FUNCIÓN “PIPELINE” GENERAL =======================
run_pipeline_for <- function(df, group_var = "Group", labels = NULL, tag_prefix = NULL) {
  # Asegurar grupos válidos
  df <- df %>% filter(!is.na(.data[[group_var]]))
  stopifnot(length(unique(df[[group_var]])) == 2)
  
  if (is.null(labels)) {
    labels <- sort(unique(df[[group_var]]))
  }
  gA <- labels[1]; gB <- labels[2]
  
  message(">>> Ejecutando para ", group_var, " = {", gA, ", ", gB, "} …")
  
  dat <- prep_dat(df)
  
  # % dentro de dominio
  tabs <- domain_norm_tables(dat, group_var = group_var)
  genus_reads <- tabs$genus_reads
  
  # keepers
  k <- pick_keepers(genus_reads, group_var = group_var)
  keepers_bact <- k$keepers_bact; keepers_euk <- k$keepers_euk
  message("Kept genera — Bacteria=", length(keepers_bact), " | Eukaryota=", length(keepers_euk))
  
  # Filtrado a keepers
  dat_filt <- genus_reads %>%
    filter((Domain=="Bacteria"  & Genus %in% keepers_bact) |
             (Domain=="Eukaryota" & Genus %in% keepers_euk))
  
  # Matrices por grupo
  cnt_A <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gA)
  cnt_B <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gB)
  stopifnot(!is.null(cnt_A), !is.null(cnt_B), nrow(cnt_A) >= MIN_SAMPLES_PER_G, nrow(cnt_B) >= MIN_SAMPLES_PER_G)
  
  mat_A <- clr_from_counts(cnt_A)
  mat_B <- clr_from_counts(cnt_B)
  
  # Alinear columnas
  common_cols <- Reduce(intersect, list(colnames(mat_A), colnames(mat_B)))
  mat_A <- mat_A[, common_cols, drop = FALSE]
  mat_B <- mat_B[, common_cols, drop = FALSE]
  cnt_A <- cnt_A[, common_cols, drop = FALSE]
  cnt_B <- cnt_B[, common_cols, drop = FALSE]
  
  # Genus -> Domain map
  genus_domain_map <- dat_filt %>%
    distinct(Genus, Domain) %>%
    filter(Genus %in% common_cols)
  
  # Correlaciones BE con bootstrap
  message("Computing BE correlations (", gA, ")…")
  edges_A <- fast_pairwise_corr_BE(
    mat_A, cnt_A, genus_domain_map,
    prescreen_global = PRESCREEN_ABS_RHO,
    prescreen_wl     = PRESCREEN_ABS_RHO_WL,
    whitelist = union(force_include_euk, force_include_bact),
    boot_R = BOOT_R
  )
  
  message("Computing BE correlations (", gB, ")…")
  edges_B <- fast_pairwise_corr_BE(
    mat_B, cnt_B, genus_domain_map,
    prescreen_global = PRESCREEN_ABS_RHO,
    prescreen_wl     = PRESCREEN_ABS_RHO_WL,
    whitelist = union(force_include_euk, force_include_bact),
    boot_R = BOOT_R
  )
  
  edges_A <- add_q(edges_A); edges_B <- add_q(edges_B)
  
  # Filtro final muy estricto
  edges_A_f <- filter_BE(edges_A)
  edges_B_f <- filter_BE(edges_B)
  message("BE edges kept | ", gA, ": ", nrow(edges_A_f), " | ", gB, ": ", nrow(edges_B_f))
  
  # Nodos y export
  nodes_A <- make_nodes(edges_A_f, genus_domain_map)
  nodes_B <- make_nodes(edges_B_f, genus_domain_map)
  
  # --- Construcción de tags para los archivos de salida (solo afecta a nombres) ---
  prefix <- if (is.null(tag_prefix)) tolower(group_var) else tag_prefix
  tagA <- paste0(prefix, "_", label_to_tag(gA, group_var), "_BE_genus_tight")
  tagB <- paste0(prefix, "_", label_to_tag(gB, group_var), "_BE_genus_tight")
  
  write_net(edges_A_f, nodes_A, tagA)
  write_net(edges_B_f, nodes_B, tagB)
  
  invisible(list(
    edges_A = edges_A_f, edges_B = edges_B_f,
    nodes_A = nodes_A, nodes_B = nodes_B,
    genus_domain_map = genus_domain_map
  ))
}

## ======================= CORRER PIPELINE PARA AMBOS CASOS =======================
# 1) Urban vs Rural (si existe)
if ("Group" %in% names(kaiju_merged) && n_distinct(na.omit(kaiju_merged$Group)) == 2) {
  res_Group <- run_pipeline_for(
    kaiju_merged %>% filter(!is.na(Group)),
    group_var = "Group",
    labels = c("Rural","Urban"),         # fuerza el orden (opcional)
    tag_prefix = "group"
  )
} else {
  message("Saltando Group: no hay dos niveles en 'Group'.")
}

# 2) BMI≥25 vs BMI<25 (nuevo)
if ("BMI_group" %in% names(kaiju_merged) || any(!is.na(kaiju_merged$BMI_group))) {
  res_BMI <- run_pipeline_for(
    kaiju_merged %>% filter(!is.na(BMI_group)),
    group_var = "BMI_group",
    labels = c("BMI<25","BMI≥25"),       # orden consistente
    tag_prefix = "bmi"
  )
} else {
  message("No hay 'BMI_group' asignado a kaiju_merged (revisa mapeo Lane->file_base).")
}

message("Listo. Redes exportadas para los dos esquemas de agrupación (cuando disponibles).")
























## =========================================================
## Redes GENUS (Bacteria/Eukaryota) por:
##   1) Group (Urban vs Rural)
##   2) BMI_group (BMI≥25 vs BMI<25)  -> nombres: BMI<25 -> BMI24 ; BMI≥25 -> BMI26
##   3) context_group (SOLO Adult + BMI≥25) en:
##        - BMI≥25 Urban  -> BMI26Urban (solo en nombres de archivo)
##        - BMI≥25 Rural  -> BMI26Rural (solo en nombres de archivo)
## =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(boot); library(readr)
})

set.seed(1234)

## ======================= PARÁMETROS GLOBALES =======================
# Filtros por grupo (prevalencia y % medio, NORMALIZADO POR DOMINIO)
MIN_PREV_BACT <- 0.20
MIN_PREV_EUK  <- 0.10
MIN_MEAN_BACT <- 0.06     # % dentro de Bacteria
MIN_MEAN_EUK  <- 0.015    # % dentro de Eukaryota
TOP_N_BACT    <- 250
TOP_N_EUK     <- 120

# Correlación / Bootstrap
BOOT_R               <- 300
PSEUDO               <- 1e-6
PRESCREEN_ABS_RHO    <- 0.35   # global
PRESCREEN_ABS_RHO_WL <- 0.20   # si toca whitelist
MIN_SAMPLES_PER_G    <- 8

# Co-ocurrencia mínima (muestras donde AMBOS > 0)
MIN_COOC_FRAC <- 0.35          # ≥35% de las muestras del grupo
MIN_COOC_ABS  <- 20            # y al menos 20 muestras

# Filtro final de aristas (solo BE)
THRESH_BE <- 0.45              # |rho_mean| mínimo
ALPHA_Q   <- 0.05              # BH FDR

# Whitelist de géneros clave
force_include_euk  <- c("Saccharomyces")
force_include_bact <- character(0)

## ======================= ENTRADAS =======================
# Se asume que ya cargaste 'kaiju_merged' en el ambiente.
# Si no, descomenta y pon tu ruta:
# kaiju_merged <- readr::read_tsv("tu_kaiju_merged.tsv", show_col_types = FALSE)

# Asegura tipos:
kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

# Si no existen columnas taxonómicas separadas pero sí 'taxon_name'
if (!all(c("Domain","Genus") %in% names(kaiju_merged)) && "taxon_name" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    tidyr::separate(
      taxon_name,
      into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
               "Subclass","Order","Family","Genus","Species"),
      sep = ";", fill = "right", extra = "drop"
    )
}

# file_base si hace falta
if (!"file_base" %in% names(kaiju_merged)) {
  if ("file" %in% names(kaiju_merged)) {
    kaiju_merged <- kaiju_merged %>% mutate(file_base = gsub("^.*/", "", file))
  } else {
    stop("No 'file_base' or 'file' column found in kaiju_merged.")
  }
}

## ============= METADATOS: construir BMI_group para BMI≥25 / BMI<25 =============
# Lee tu archivo de metadatos (ajusta la ruta si cambia)
meta <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE, sep = ",", fileEncoding = "latin1",
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

# Limpieza mínima
meta <- meta[, colSums(!is.na(meta)) > 0]
meta <- meta[rowSums(is.na(meta)) < ncol(meta), ]
meta$BMI <- suppressWarnings(as.numeric(meta$BMI))
meta$Age <- suppressWarnings(as.numeric(meta$Age))

# Agrupaciones de IMC (como lo tenías + colapsado a 2 niveles)
meta <- meta %>%
  mutate(
    BMI_group = case_when(
      Percentil_formulas %in% c("Overweight", "Obesity") ~ "Overweight/Obesity",
      Percentil_formulas == "Normal Weight" ~ "Normal Weight",
      Percentil_formulas == "Malnutrition" ~ "Malnutrition",
      TRUE ~ NA_character_
    ),
    BMI_index_group = case_when(
      BMI_group == "Overweight/Obesity" ~ ">25 BMI INDEX",
      BMI_group %in% c("Normal Weight", "Malnutrition") ~ "<25 BMI INDEX",
      TRUE ~ NA_character_
    )
  )

# En tus metadatos, el identificador de secuenciación parece estar en 'Lane'
# y los archivos de Kaiju acaban en "_kaiju.out". Ajusta si tu sufijo es distinto.
files_over25 <- meta %>% filter(BMI_index_group == ">25 BMI INDEX") %>% pull(Lane) %>% unique()
files_under25 <- meta %>% filter(BMI_index_group == "<25 BMI INDEX") %>% pull(Lane) %>% unique()

# Asignar BMI_group a kaiju_merged (etiquetas finales BMI≥25 / BMI<25)
kaiju_merged <- kaiju_merged %>%
  mutate(BMI_group = case_when(
    file_base %in% paste0(files_over25, "_kaiju.out")  ~ "BMI≥25",
    file_base %in% paste0(files_under25, "_kaiju.out") ~ "BMI<25",
    TRUE ~ NA_character_
  ))

## ============= (Opcional) asignar Group (Urban/Rural) si no existe =============
if (!"Group" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    mutate(Group = case_when(
      exists("urban_files") & file_base %in% paste0(urban_files, "_kaiju.out") ~ "Urban",
      exists("rural_files") & file_base %in% paste0(rural_files, "_kaiju.out") ~ "Rural",
      TRUE ~ NA_character_
    ))
}

## ===== NUEVO: construir context_group = SOLO Adult + BMI≥25 y dividir Urban/Rural =====
# 1) Lista de Lanes "Adult"
adult_lanes <- meta %>%
  dplyr::filter(!is.na(Age_group), Age_group == "Adult") %>%
  dplyr::pull(Lane) %>%
  unique()

KAIJU_SUFFIX <- "_kaiju.out"
adult_files  <- paste0(adult_lanes, KAIJU_SUFFIX)

# 2) Columna context_group en kaiju_merged (solo Adult + BMI≥25; Urban/Rural aparte)
kaiju_merged <- kaiju_merged %>%
  dplyr::mutate(
    is_adult_file = file_base %in% adult_files,
    context_group = dplyr::case_when(
      is_adult_file & BMI_group == "BMI≥25" & Group == "Urban" ~ "BMI≥25 Urban",
      is_adult_file & BMI_group == "BMI≥25" & Group == "Rural" ~ "BMI≥25 Rural",
      TRUE ~ NA_character_
    )
  )

## ======================= Helpers genéricos =======================
# Limpieza base (solo Bacteria/Eukaryota a nivel GENUS; excluye Homo)
prep_dat <- function(df_all) {
  df_all %>%
    filter(Domain %in% c("Bacteria","Eukaryota")) %>%
    mutate(Genus = ifelse(is.na(Genus) | Genus=="" | Genus=="Unclassified", NA, str_trim(Genus))) %>%
    filter(!is.na(Genus)) %>%
    filter(!(Domain=="Eukaryota" & Genus=="Homo"))
}

# Métricas por dominio dentro de cada muestra (% dentro del DOMINIO)
domain_norm_tables <- function(dat, group_var = "Group") {
  stopifnot(all(c("file_base", group_var, "Domain", "Genus", "reads") %in% names(dat)))
  totals_domain <- dat %>%
    group_by(.data[[group_var]], file_base, Domain) %>%
    summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")
  genus_reads <- dat %>%
    group_by(.data[[group_var]], file_base, Domain, Genus) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(totals_domain, by = c(group_var, "file_base", "Domain")) %>%
    mutate(pct_domain = if_else(total_reads_domain > 0, 100 * reads / total_reads_domain, 0))
  list(genus_reads = genus_reads, totals_domain = totals_domain)
}

# Elegir géneros “keepers” por dominio (prevalencia y % medio) usando el group_var
pick_keepers <- function(genus_reads, group_var = "Group") {
  prev_mean <- genus_reads %>%
    group_by(.data[[group_var]], Domain, Genus) %>%
    summarise(prevalence = mean(reads > 0),
              mean_pct_dom = mean(pct_domain),
              .groups = "drop")
  
  # NOTA: asumimos 2 niveles en group_var; se ordenan por aparición
  lvls <- sort(unique(prev_mean[[group_var]]))
  
  keepers_bact <- prev_mean %>%
    filter(Domain=="Bacteria") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_BACT,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_BACT),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_BACT) %>% pull(Genus)
  
  keepers_euk <- prev_mean %>%
    filter(Domain=="Eukaryota") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_EUK,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_EUK),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_EUK) %>% pull(Genus)
  
  keepers_bact <- union(keepers_bact,
                        intersect(force_include_bact, unique(genus_reads$Genus[genus_reads$Domain=="Bacteria"])))
  keepers_euk  <- union(keepers_euk,
                        intersect(force_include_euk,  unique(genus_reads$Genus[genus_reads$Domain=="Eukaryota"])))
  
  list(keepers_bact = keepers_bact, keepers_euk = keepers_euk)
}

# Matrices de counts por grupo (según group_var) y CLR
counts_from_long <- function(df_long, group_var, group_value) {
  df_long %>%
    filter(.data[[group_var]] == group_value) %>%
    select(file_base, Genus, reads) %>%
    group_by(file_base, Genus) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Genus, values_from = reads, values_fill = 0) %>%
    column_to_rownames("file_base") %>% as.matrix()
}

clr_from_counts <- function(counts_mat) {
  if (is.null(counts_mat) || nrow(counts_mat) < MIN_SAMPLES_PER_G) return(NULL)
  prop <- sweep(counts_mat, 1, rowSums(counts_mat), "/"); prop[is.na(prop)] <- 0
  logx <- log(prop + PSEUDO); sweep(logx, 1, rowMeans(logx), "-")
}

# Correlaciones BE rápidas con bootstrap
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
  
  # Co-ocurrencia
  pres_abs <- max(min_cooc_abs, ceiling(min_cooc_frac * nrow(counts_mat)))
  pres_mat <- counts_mat > 0
  
  # Spearman masivo = Pearson sobre rangos
  mat_rank <- apply(mat_clr, 2, rank, ties.method = "average")
  R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                  mat_rank[, euk_cols,  drop = FALSE],
                  method = "pearson", use = "pairwise.complete.obs")
  
  # Candidatos: prescreen global
  sel_global <- which(abs(R) >= prescreen_global, arr.ind = TRUE)
  pairs_df <- tibble(
    Genus_1 = bact_cols[sel_global[,"row"]],
    Genus_2 = euk_cols[ sel_global[,"col"]],
    from_whitelist = FALSE
  )
  
  # Candidatos extra: whitelist (más laxo)
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
  
  # Únicos + co-ocurrencia mínima
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

# BH helper
add_q <- function(df) if (nrow(df)) df %>% mutate(q = p.adjust(p, method = "BH")) else df

# Filtro final “muy estricto”
filter_BE <- function(df) {
  if (nrow(df)==0) return(df[0,])
  df %>%
    mutate(abs_rho = abs(rho_mean)) %>%
    mutate(ci_strong = (rho_mean > 0 & CI_lower >= THRESH_BE) |
             (rho_mean < 0 & CI_upper <= -THRESH_BE)) %>%
    filter(
      abs_rho >= THRESH_BE,
      abs(rho) >= THRESH_BE,
      ci_strong,
      sign_consistency >= 0.90,
      q <= ALPHA_Q
    ) %>%
    select(-abs_rho, -ci_strong)
}

make_nodes <- function(edges_df, genus_domain_map) {
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

# --- Mapea etiquetas a tags seguros solo para los nombres de archivo ---
label_to_tag <- function(lbl, group_var = NULL) {
  # Caso especial: BMI_group -> distinción clara en nombres
  if (!is.null(group_var) && group_var == "BMI_group") {
    if (grepl("^\\s*BMI\\s*<\\s*25\\s*$", lbl))  return("BMI24")
    if (grepl("^\\s*BMI\\s*[≥>=]\\s*25\\s*$", lbl)) return("BMI26")
  }
  # Caso especial: context_group (Adult + BMI≥25 Urban/Rural)
  if (!is.null(group_var) && group_var == "context_group") {
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Urban", lbl)) return("BMI26Urban")
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Rural", lbl)) return("BMI26Rural")
  }
  # Fallback general
  out <- gsub("[^A-Za-z0-9]+", "", lbl)
  if (nzchar(out)) out else "grp"
}

## ======================= FUNCIÓN “PIPELINE” GENERAL =======================
`%||%` <- function(a,b) if (!is.null(a)) a else b

run_pipeline_for <- function(df, group_var = "Group", labels = NULL, tag_prefix = NULL) {
  # Asegurar grupos válidos
  df <- df %>% filter(!is.na(.data[[group_var]]))
  stopifnot(length(unique(df[[group_var]])) == 2)
  
  if (is.null(labels)) {
    labels <- sort(unique(df[[group_var]]))
  }
  gA <- labels[1]; gB <- labels[2]
  
  message(">>> Ejecutando para ", group_var, " = {", gA, ", ", gB, "} …")
  
  dat <- prep_dat(df)
  
  # % dentro de dominio
  tabs <- domain_norm_tables(dat, group_var = group_var)
  genus_reads <- tabs$genus_reads
  
  # keepers
  k <- pick_keepers(genus_reads, group_var = group_var)
  keepers_bact <- k$keepers_bact; keepers_euk <- k$keepers_euk
  message("Kept genera — Bacteria=", length(keepers_bact), " | Eukaryota=", length(keepers_euk))
  
  # Filtrado a keepers
  dat_filt <- genus_reads %>%
    filter((Domain=="Bacteria"  & Genus %in% keepers_bact) |
             (Domain=="Eukaryota" & Genus %in% keepers_euk))
  
  # Matrices por grupo
  cnt_A <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gA)
  cnt_B <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gB)
  stopifnot(!is.null(cnt_A), !is.null(cnt_B), nrow(cnt_A) >= MIN_SAMPLES_PER_G, nrow(cnt_B) >= MIN_SAMPLES_PER_G)
  
  mat_A <- clr_from_counts(cnt_A)
  mat_B <- clr_from_counts(cnt_B)
  
  # Alinear columnas
  common_cols <- Reduce(intersect, list(colnames(mat_A), colnames(mat_B)))
  mat_A <- mat_A[, common_cols, drop = FALSE]
  mat_B <- mat_B[, common_cols, drop = FALSE]
  cnt_A <- cnt_A[, common_cols, drop = FALSE]
  cnt_B <- cnt_B[, common_cols, drop = FALSE]
  
  # Genus -> Domain map
  genus_domain_map <- dat_filt %>%
    distinct(Genus, Domain) %>%
    filter(Genus %in% common_cols)
  
  # Correlaciones BE con bootstrap
  message("Computing BE correlations (", gA, ")…")
  edges_A <- fast_pairwise_corr_BE(
    mat_A, cnt_A, genus_domain_map,
    prescreen_global = PRESCREEN_ABS_RHO,
    prescreen_wl     = PRESCREEN_ABS_RHO_WL,
    whitelist = union(force_include_euk, force_include_bact),
    boot_R = BOOT_R
  )
  
  message("Computing BE correlations (", gB, ")…")
  edges_B <- fast_pairwise_corr_BE(
    mat_B, cnt_B, genus_domain_map,
    prescreen_global = PRESCREEN_ABS_RHO,
    prescreen_wl     = PRESCREEN_ABS_RHO_WL,
    whitelist = union(force_include_euk, force_include_bact),
    boot_R = BOOT_R
  )
  
  edges_A <- add_q(edges_A); edges_B <- add_q(edges_B)
  
  # Filtro final muy estricto
  edges_A_f <- filter_BE(edges_A)
  edges_B_f <- filter_BE(edges_B)
  message("BE edges kept | ", gA, ": ", nrow(edges_A_f), " | ", gB, ": ", nrow(edges_B_f))
  
  # Nodos y export
  nodes_A <- make_nodes(edges_A_f, genus_domain_map)
  nodes_B <- make_nodes(edges_B_f, genus_domain_map)
  
  # --- Construcción de tags para los archivos de salida (solo afecta a nombres) ---
  prefix <- (tag_prefix %||% tolower(group_var))
  tagA <- paste0(prefix, "_", label_to_tag(gA, group_var), "_BE_genus_tight")
  tagB <- paste0(prefix, "_", label_to_tag(gB, group_var), "_BE_genus_tight")
  
  write_net(edges_A_f, nodes_A, tagA)
  write_net(edges_B_f, nodes_B, tagB)
  
  invisible(list(
    edges_A = edges_A_f, edges_B = edges_B_f,
    nodes_A = nodes_A, nodes_B = nodes_B,
    genus_domain_map = genus_domain_map
  ))
}

## ======================= CORRER PIPELINE PARA LOS 3 CASOS =======================
# 1) Urban vs Rural
if ("Group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$Group)) && dplyr::n_distinct(na.omit(kaiju_merged$Group)) == 2) {
  res_Group <- run_pipeline_for(
    kaiju_merged %>% filter(!is.na(Group)),
    group_var = "Group",
    labels = c("Rural","Urban"),         # fuerza el orden (opcional)
    tag_prefix = "group"
  )
} else {
  message("Saltando Group: no hay dos niveles en 'Group'.")
}

# 2) BMI≥25 vs BMI<25
if ("BMI_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$BMI_group))) {
  res_BMI <- run_pipeline_for(
    kaiju_merged %>% filter(!is.na(BMI_group)),
    group_var = "BMI_group",
    labels = c("BMI<25","BMI≥25"),       # orden consistente
    tag_prefix = "bmi"
  )
} else {
  message("No hay 'BMI_group' asignado a kaiju_merged (revisa mapeo Lane->file_base).")
}

# 3) context_group: SOLO Adult + BMI≥25, comparando Urban vs Rural (2 salidas más)
if ("context_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$context_group))) {
  # Debe haber exactamente dos niveles: "BMI≥25 Urban" y "BMI≥25 Rural"
  lvls_ctx <- sort(unique(na.omit(kaiju_merged$context_group)))
  if (length(lvls_ctx) == 2) {
    res_CTX <- run_pipeline_for(
      kaiju_merged %>% filter(!is.na(context_group)),
      group_var = "context_group",
      labels = c("BMI≥25 Rural","BMI≥25 Urban"),  # orden explícito
      tag_prefix = "context"
    )
  } else {
    message("Saltando context_group: niveles detectados = ", paste(lvls_ctx, collapse=", "))
  }
} else {
  message("No hay 'context_group' (Adult + BMI≥25) disponible.")
}

message("Listo. Redes exportadas para Group, BMI_group y context_group (si disponibles).")







##############


## =========================================================
## Redes GENUS (Bacteria/Eukaryota) por:
##   1) Group (Urban vs Rural)
##   2) BMI_group (BMI≥25 vs BMI<25)  -> nombres: BMI<25 -> BMI24 ; BMI≥25 -> BMI26
##   3) context_group (SOLO Adult + BMI≥25):
##        - BMI≥25 Urban  -> BMI26Urban (solo en nombres de archivo)
##        - BMI≥25 Rural  -> BMI26Rural (solo en nombres de archivo)
## Además, para cada esquema se generan 3 redes: BE, BB y EE.
## =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(boot); library(readr)
})

set.seed(1234)

## ======================= PARÁMETROS GLOBALES =======================
# Filtros por grupo (prevalencia y % medio, NORMALIZADO POR DOMINIO)
#Prevalencia mínima (MIN_PREV_). Debe aparecer en al menos ese porcentaje de muestras del grupo
MIN_PREV_BACT <- 0.50 #Base .2
MIN_PREV_EUK  <- 0.50 #Base .1

#porcentaje relativo dentro de cada dominio
MIN_MEAN_BACT <- 0.1    #0.06 base # % dentro de Bacteria
MIN_MEAN_EUK  <- 0.1    #0.015 base % dentro de Eukaryota

#el no. máximo de taxones en caso de que sea muy grande
TOP_N_BACT    <- 70
TOP_N_EUK     <- 50

# Correlación / Bootstrap
BOOT_R               <- 300
PSEUDO               <- 1e-6
# Umbrales de prescreen por relación
PRESCREEN_BE         <- 0.5 #0.35
PRESCREEN_BB         <- 0.5 #0.45
PRESCREEN_EE         <- 0.5 #0.30
PRESCREEN_WL         <- 0.5 #0.20   # si toca whitelist

#minimo de muestras para evaluar 
MIN_SAMPLES_PER_G    <- 4

# Co-ocurrencia mínima (muestras donde AMBOS > 0)
MIN_COOC_FRAC <- 0.50          # que la relacion ocurre en al menos ≥35% de las muestras del grupo
MIN_COOC_ABS  <- 4           # y al menos 20 muestras

# Filtro final por relación (umbral de |rho_mean| y FDR)
THRESH_BE <- 0.5
THRESH_BB <- 0.50
THRESH_EE <- 0.5
ALPHA_Q   <- 0.05              # BH FDR

# Whitelist de géneros clave
force_include_euk  <- c("Saccharomyces")
force_include_bact <- character(0)

## ======================= ENTRADAS =======================
# Se asume que ya cargaste 'kaiju_merged' en el ambiente.
# Si no, descomenta y pon tu ruta:
# kaiju_merged <- readr::read_tsv("tu_kaiju_merged.tsv", show_col_types = FALSE)

# Asegura tipos:
kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

# Si no existen columnas taxonómicas separadas pero sí 'taxon_name'
if (!all(c("Domain","Genus") %in% names(kaiju_merged)) && "taxon_name" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    tidyr::separate(
      taxon_name,
      into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
               "Subclass","Order","Family","Genus","Species"),
      sep = ";", fill = "right", extra = "drop"
    )
}

# file_base si hace falta
if (!"file_base" %in% names(kaiju_merged)) {
  if ("file" %in% names(kaiju_merged)) {
    kaiju_merged <- kaiju_merged %>% mutate(file_base = gsub("^.*/", "", file))
  } else {
    stop("No 'file_base' or 'file' column found in kaiju_merged.")
  }
}

## ============= METADATOS: construir BMI_group para BMI≥25 / BMI<25 =============
# Lee tu archivo de metadatos (ajusta la ruta si cambia)
meta <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE, sep = ",", fileEncoding = "latin1",
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

# Limpieza mínima
meta <- meta[, colSums(!is.na(meta)) > 0]
meta <- meta[rowSums(is.na(meta)) < ncol(meta), ]
meta$BMI <- suppressWarnings(as.numeric(meta$BMI))
meta$Age <- suppressWarnings(as.numeric(meta$Age))

# Agrupaciones de IMC (como lo tenías + colapsado a 2 niveles)
meta <- meta %>%
  mutate(
    BMI_group = case_when(
      Percentil_formulas %in% c("Overweight", "Obesity") ~ "Overweight/Obesity",
      Percentil_formulas == "Normal Weight" ~ "Normal Weight",
      Percentil_formulas == "Malnutrition" ~ "Malnutrition",
      TRUE ~ NA_character_
    ),
    BMI_index_group = case_when(
      BMI_group == "Overweight/Obesity" ~ ">25 BMI INDEX",
      BMI_group %in% c("Normal Weight", "Malnutrition") ~ "<25 BMI INDEX",
      TRUE ~ NA_character_
    )
  )

# En tus metadatos, el identificador de secuenciación parece estar en 'Lane'
# y los archivos de Kaiju acaban en "_kaiju.out". Ajusta si tu sufijo es distinto.
files_over25  <- meta %>% filter(BMI_index_group == ">25 BMI INDEX") %>% pull(Lane) %>% unique()
files_under25 <- meta %>% filter(BMI_index_group == "<25 BMI INDEX") %>% pull(Lane) %>% unique()

# Asignar BMI_group a kaiju_merged (etiquetas finales BMI≥25 / BMI<25)
kaiju_merged <- kaiju_merged %>%
  mutate(BMI_group = case_when(
    file_base %in% paste0(files_over25,  "_kaiju.out") ~ "BMI≥25",
    file_base %in% paste0(files_under25, "_kaiju.out") ~ "BMI<25",
    TRUE ~ NA_character_
  ))

## ============= (Opcional) asignar Group (Urban/Rural) si no existe =============
if (!"Group" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    mutate(Group = case_when(
      exists("urban_files") & file_base %in% paste0(urban_files, "_kaiju.out") ~ "Urban",
      exists("rural_files") & file_base %in% paste0(rural_files, "_kaiju.out") ~ "Rural",
      TRUE ~ NA_character_
    ))
}

## ===== NUEVO: construir context_group = SOLO Adult + BMI≥25 y dividir Urban/Rural =====
# 1) Lista de Lanes "Adult"
adult_lanes <- meta %>%
  filter(!is.na(Age_group), Age_group == "Adult",
         !is.na(Gender), Gender == "Male") %>%
  pull(Lane) %>%
  unique()


KAIJU_SUFFIX <- "_kaiju.out"
adult_files  <- paste0(adult_lanes, KAIJU_SUFFIX)

# 2) Columna context_group en kaiju_merged (solo Adult + BMI≥25; Urban/Rural aparte)
kaiju_merged <- kaiju_merged %>%
  dplyr::mutate(
    is_adult_file = file_base %in% adult_files,
    context_group = dplyr::case_when(
      is_adult_file & BMI_group == "BMI≥25" & Group == "Urban" ~ "BMI≥25 Urban",
      is_adult_file & BMI_group == "BMI≥25" & Group == "Rural" ~ "BMI≥25 Rural",
      TRUE ~ NA_character_
    )
  )

## ======================= Helpers genéricos =======================
# Limpieza base (solo Bacteria/Eukaryota a nivel GENUS; excluye Homo)
prep_dat <- function(df_all) {
  df_all %>%
    filter(Domain %in% c("Bacteria","Eukaryota")) %>%
    mutate(Genus = ifelse(is.na(Genus) | Genus=="" | Genus=="Unclassified", NA, str_trim(Genus))) %>%
    filter(!is.na(Genus)) %>%
    filter(!(Domain=="Eukaryota" & Genus=="Homo"))
}

# Métricas por dominio dentro de cada muestra (% dentro del DOMINIO)
domain_norm_tables <- function(dat, group_var = "Group") {
  stopifnot(all(c("file_base", group_var, "Domain", "Genus", "reads") %in% names(dat)))
  totals_domain <- dat %>%
    group_by(.data[[group_var]], file_base, Domain) %>%
    summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")
  genus_reads <- dat %>%
    group_by(.data[[group_var]], file_base, Domain, Genus) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(totals_domain, by = c(group_var, "file_base", "Domain")) %>%
    mutate(pct_domain = if_else(total_reads_domain > 0, 100 * reads / total_reads_domain, 0))
  list(genus_reads = genus_reads, totals_domain = totals_domain)
}

# Elegir géneros “keepers” por dominio (prevalencia y % medio) usando el group_var
pick_keepers <- function(genus_reads, group_var = "Group") {
  prev_mean <- genus_reads %>%
    group_by(.data[[group_var]], Domain, Genus) %>%
    summarise(prevalence = mean(reads > 0),
              mean_pct_dom = mean(pct_domain),
              .groups = "drop")
  
  # NOTA: asumimos 2 niveles en group_var; se ordenan por aparición
  lvls <- sort(unique(prev_mean[[group_var]]))
  
  keepers_bact <- prev_mean %>%
    filter(Domain=="Bacteria") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_BACT,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_BACT),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_BACT) %>% pull(Genus)
  
  keepers_euk <- prev_mean %>%
    filter(Domain=="Eukaryota") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_EUK,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_EUK),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_EUK) %>% pull(Genus)
  
  keepers_bact <- union(keepers_bact,
                        intersect(force_include_bact, unique(genus_reads$Genus[genus_reads$Domain=="Bacteria"])))
  keepers_euk  <- union(keepers_euk,
                        intersect(force_include_euk,  unique(genus_reads$Genus[genus_reads$Domain=="Eukaryota"])))
  
  list(keepers_bact = keepers_bact, keepers_euk = keepers_euk)
}

# Matrices de counts por grupo (según group_var) y CLR
counts_from_long <- function(df_long, group_var, group_value) {
  df_long %>%
    filter(.data[[group_var]] == group_value) %>%
    select(file_base, Genus, reads) %>%
    group_by(file_base, Genus) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Genus, values_from = reads, values_fill = 0) %>%
    column_to_rownames("file_base") %>% as.matrix()
}

clr_from_counts <- function(counts_mat) {
  if (is.null(counts_mat) || nrow(counts_mat) < MIN_SAMPLES_PER_G) return(NULL)
  prop <- sweep(counts_mat, 1, rowSums(counts_mat), "/"); prop[is.na(prop)] <- 0
  logx <- log(prop + PSEUDO); sweep(logx, 1, rowMeans(logx), "-")
}

# -------- Correlaciones GENÉRICAS (BE, BB, EE) con bootstrap --------
fast_pairwise_corr <- function(mat_clr, counts_mat, genus_domain_map,
                               relation = c("BE","BB","EE"),
                               prescreen_global = list(BE = PRESCREEN_BE, BB = PRESCREEN_BB, EE = PRESCREEN_EE),
                               prescreen_wl     = PRESCREEN_WL,
                               min_cooc_frac = MIN_COOC_FRAC,
                               min_cooc_abs  = MIN_COOC_ABS,
                               boot_R = BOOT_R,
                               whitelist = character(0)) {
  
  relation <- match.arg(relation)
  cols <- colnames(mat_clr)
  dom_vec <- setNames(as.character(genus_domain_map$Domain),
                      as.character(genus_domain_map$Genus))
  
  bact_cols <- cols[dom_vec[cols] == "Bacteria"]
  euk_cols  <- cols[dom_vec[cols] == "Eukaryota"]
  
  if (relation == "BE" && (length(bact_cols) == 0 || length(euk_cols) == 0)) return(tibble())
  if (relation == "BB" && length(bact_cols) < 2) return(tibble())
  if (relation == "EE" && length(euk_cols) < 2) return(tibble())
  
  # Co-ocurrencia
  pres_abs <- max(min_cooc_abs, ceiling(min_cooc_frac * nrow(counts_mat)))
  pres_mat <- counts_mat > 0
  
  # Rangos para Spearman masivo
  mat_rank <- apply(mat_clr, 2, rank, ties.method = "average")
  
  if (relation == "BE") {
    R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                    mat_rank[, euk_cols,  drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    sel <- which(abs(R) >= prescreen_global$BE, arr.ind = TRUE)
    pairs_df <- tibble(Genus_1 = bact_cols[sel[, "row"]],
                       Genus_2 = euk_cols[ sel[, "col"]],
                       from_whitelist = FALSE)
    # Whitelist laxo
    if (length(whitelist) > 0) {
      wl_b <- intersect(bact_cols, whitelist)
      wl_e <- intersect(euk_cols,  whitelist)
      if (length(wl_b) > 0) {
        R_wlb <- R[match(wl_b, bact_cols), , drop = FALSE]
        sel_wlb <- which(abs(R_wlb) >= prescreen_wl, arr.ind = TRUE)
        if (length(sel_wlb))
          pairs_df <- bind_rows(pairs_df,
                                tibble(Genus_1 = wl_b[sel_wlb[, "row"]],
                                       Genus_2 = euk_cols[sel_wlb[, "col"]],
                                       from_whitelist = TRUE))
      }
      if (length(wl_e) > 0) {
        R_wle <- R[, match(wl_e, euk_cols), drop = FALSE]
        sel_wle <- which(abs(R_wle) >= prescreen_wl, arr.ind = TRUE)
        if (length(sel_wle))
          pairs_df <- bind_rows(pairs_df,
                                tibble(Genus_1 = bact_cols[sel_wle[, "row"]],
                                       Genus_2 = wl_e[sel_wle[, "col"]],
                                       from_whitelist = TRUE))
      }
    }
  } else if (relation == "BB") {
    R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    upper <- upper.tri(R, diag = FALSE)
    idx <- which(upper & abs(R) >= prescreen_global$BB, arr.ind = TRUE)
    pairs_df <- tibble(Genus_1 = bact_cols[idx[, "row"]],
                       Genus_2 = bact_cols[idx[, "col"]],
                       from_whitelist = FALSE)
  } else { # EE
    R <- stats::cor(mat_rank[, euk_cols, drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    upper <- upper.tri(R, diag = FALSE)
    idx <- which(upper & abs(R) >= prescreen_global$EE, arr.ind = TRUE)
    pairs_df <- tibble(Genus_1 = euk_cols[idx[, "row"]],
                       Genus_2 = euk_cols[idx[, "col"]],
                       from_whitelist = FALSE)
  }
  
  if (nrow(pairs_df) == 0) return(tibble())
  
  # Únicos + co-ocurrencia mínima
  pairs_df <- pairs_df %>%
    distinct(Genus_1, Genus_2, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(cooc = sum(pres_mat[, Genus_1] & pres_mat[, Genus_2])) %>%
    ungroup() %>%
    filter(cooc >= pres_abs)
  
  message("Prescreen(", relation, ") pairs (post co-occur): ", nrow(pairs_df),
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
    scons <- mean(sign(bt$t) == sign(rho0), na.rm = TRUE)
    rel_label <- if (relation=="BE") "Bacteria-Eukaryota" else if (relation=="BB") "Bacteria-Bacteria" else "Eukaryota-Eukaryota"
    tibble(
      Genus_1 = g1, Genus_2 = g2, Relation = rel_label,
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

# BH helper
add_q <- function(df) if (nrow(df)) df %>% mutate(q = p.adjust(p, method = "BH")) else df

# Filtro final por relación
filter_by_relation <- function(df, relation = c("BE","BB","EE")) {
  relation <- match.arg(relation)
  if (nrow(df)==0) return(df[0,])
  
  thr <- switch(relation,
                "BE" = THRESH_BE,
                "BB" = THRESH_BB,
                "EE" = THRESH_EE)
  
  df %>%
    mutate(abs_rho = abs(rho_mean)) %>%
    mutate(ci_strong = (rho_mean > 0 & CI_lower >= thr) |
             (rho_mean < 0 & CI_upper <= -thr)) %>%
    filter(
      abs_rho >= thr,
      abs(rho) >= thr,
      ci_strong,
      sign_consistency >= 0.90,
      q <= ALPHA_Q
    ) %>%
    select(-abs_rho, -ci_strong)
}

make_nodes <- function(edges_df, genus_domain_map) {
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

# --- Mapea etiquetas a tags seguros solo para los nombres de archivo ---
label_to_tag <- function(lbl, group_var = NULL) {
  if (!is.null(group_var) && group_var == "BMI_group") {
    if (grepl("^\\s*BMI\\s*<\\s*25\\s*$", lbl))  return("BMI24")
    if (grepl("^\\s*BMI\\s*[≥>=]\\s*25\\s*$", lbl)) return("BMI26")
  }
  if (!is.null(group_var) && group_var == "context_group") {
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Urban", lbl)) return("BMI26Urban")
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Rural", lbl)) return("BMI26Rural")
  }
  out <- gsub("[^A-Za-z0-9]+", "", lbl)
  if (nzchar(out)) out else "grp"
}

relation_tag <- function(relation_label) {
  if (relation_label == "Bacteria-Eukaryota") return("BE")
  if (relation_label == "Bacteria-Bacteria")  return("BB")
  if (relation_label == "Eukaryota-Eukaryota")return("EE")
  "REL"
}

## ======================= FUNCIÓN “PIPELINE” GENERAL =======================
`%||%` <- function(a,b) if (!is.null(a)) a else b

run_pipeline_for <- function(df, group_var = "Group", labels = NULL, tag_prefix = NULL) {
  # Asegurar grupos válidos
  df <- df %>% filter(!is.na(.data[[group_var]]))
  stopifnot(length(unique(df[[group_var]])) == 2)
  
  if (is.null(labels)) {
    labels <- sort(unique(df[[group_var]]))
  }
  gA <- labels[1]; gB <- labels[2]
  
  message(">>> Ejecutando para ", group_var, " = {", gA, ", ", gB, "} …")
  
  dat <- prep_dat(df)
  
  # % dentro de dominio
  tabs <- domain_norm_tables(dat, group_var = group_var)
  genus_reads <- tabs$genus_reads
  
  # keepers
  k <- pick_keepers(genus_reads, group_var = group_var)
  keepers_bact <- k$keepers_bact; keepers_euk <- k$keepers_euk
  message("Kept genera — Bacteria=", length(keepers_bact), " | Eukaryota=", length(keepers_euk))
  
  # Filtrado a keepers
  dat_filt <- genus_reads %>%
    filter((Domain=="Bacteria"  & Genus %in% keepers_bact) |
             (Domain=="Eukaryota" & Genus %in% keepers_euk))
  
  # Matrices por grupo
  cnt_A <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gA)
  cnt_B <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gB)
  stopifnot(!is.null(cnt_A), !is.null(cnt_B), nrow(cnt_A) >= MIN_SAMPLES_PER_G, nrow(cnt_B) >= MIN_SAMPLES_PER_G)
  
  mat_A <- clr_from_counts(cnt_A)
  mat_B <- clr_from_counts(cnt_B)
  
  # Alinear columnas
  common_cols <- Reduce(intersect, list(colnames(mat_A), colnames(mat_B)))
  mat_A <- mat_A[, common_cols, drop = FALSE]
  mat_B <- mat_B[, common_cols, drop = FALSE]
  cnt_A <- cnt_A[, common_cols, drop = FALSE]
  cnt_B <- cnt_B[, common_cols, drop = FALSE]
  
  # Genus -> Domain map
  genus_domain_map <- dat_filt %>%
    distinct(Genus, Domain) %>%
    filter(Genus %in% common_cols)
  
  # --------- Correlaciones y filtros por RELACIÓN: BE, BB, EE ----------
  relations <- c("BE","BB","EE")
  for (rel in relations) {
    message("Computing ", rel, " correlations (", gA, ")…")
    edges_A <- fast_pairwise_corr(
      mat_A, cnt_A, genus_domain_map,
      relation = rel,
      whitelist = union(force_include_euk, force_include_bact),
      boot_R = BOOT_R
    )
    message("Computing ", rel, " correlations (", gB, ")…")
    edges_B <- fast_pairwise_corr(
      mat_B, cnt_B, genus_domain_map,
      relation = rel,
      whitelist = union(force_include_euk, force_include_bact),
      boot_R = BOOT_R
    )
    
    edges_A <- add_q(edges_A); edges_B <- add_q(edges_B)
    
    edges_A_f <- filter_by_relation(edges_A, relation = rel)
    edges_B_f <- filter_by_relation(edges_B, relation = rel)
    
    message(rel, " edges kept | ", gA, ": ", nrow(edges_A_f), " | ", gB, ": ", nrow(edges_B_f))
    
    nodes_A <- make_nodes(edges_A_f, genus_domain_map)
    nodes_B <- make_nodes(edges_B_f, genus_domain_map)
    
    # Tags de salida
    prefix <- (tag_prefix %||% tolower(group_var))
    rel_tag <- rel
    tagA <- paste0(prefix, "_", label_to_tag(gA, group_var), "_", rel_tag, "_genus_tight")
    tagB <- paste0(prefix, "_", label_to_tag(gB, group_var), "_", rel_tag, "_genus_tight")
    
    write_net(edges_A_f, nodes_A, tagA)
    write_net(edges_B_f, nodes_B, tagB)
  }
  
  invisible(NULL)
}

## ======================= CORRER PIPELINE PARA LOS 3 CASOS =======================
# 1) Urban vs Rural
if ("Group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$Group)) && dplyr::n_distinct(na.omit(kaiju_merged$Group)) == 2) {
  run_pipeline_for(
    kaiju_merged %>% filter(!is.na(Group)),
    group_var = "Group",
    labels = c("Rural","Urban"),         # fuerza el orden (opcional)
    tag_prefix = "group"
  )
} else {
  message("Saltando Group: no hay dos niveles en 'Group'.")
}

# 2) BMI≥25 vs BMI<25
if ("BMI_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$BMI_group))) {
  run_pipeline_for(
    kaiju_merged %>% filter(!is.na(BMI_group)),
    group_var = "BMI_group",
    labels = c("BMI<25","BMI≥25"),       # orden consistente
    tag_prefix = "bmi"
  )
} else {
  message("No hay 'BMI_group' asignado a kaiju_merged (revisa mapeo Lane->file_base).")
}

# 3) context_group: SOLO Adult + BMI≥25, comparando Urban vs Rural
if ("context_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$context_group))) {
  lvls_ctx <- sort(unique(na.omit(kaiju_merged$context_group)))
  if (length(lvls_ctx) == 2) {
    run_pipeline_for(
      kaiju_merged %>% filter(!is.na(context_group)),
      group_var = "context_group",
      labels = c("BMI≥25 Rural","BMI≥25 Urban"),  # orden explícito
      tag_prefix = "context"
    )
  } else {
    message("Saltando context_group: niveles detectados = ", paste(lvls_ctx, collapse=", "))
  }
} else {
  message("No hay 'context_group' (Adult + BMI≥25) disponible.")
}

message("Listo. Redes BE/BB/EE exportadas para Group, BMI_group y context_group (si disponibles).")





















## =========================================================
## Redes GENUS (Bacteria/Eukaryota) por:
##   1) Group (Urban vs Rural)
##   2) BMI_group (BMI≥25 vs BMI<25)  -> nombres: BMI<25 -> BMI24 ; BMI≥25 -> BMI26
##   3) context_group (SOLO Adult + BMI≥25):
##        - BMI≥25 Urban  -> BMI26Urban (solo en nombres de archivo)
##        - BMI≥25 Rural  -> BMI26Rural (solo en nombres de archivo)
## Además, para cada esquema se generan 3 redes: BE, BB y EE.
## =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(boot); library(readr)
})

set.seed(1234)

## ======================= PARÁMETROS GLOBALES =======================
# Filtros por grupo (prevalencia y % medio, NORMALIZADO POR DOMINIO)
#Prevalencia mínima (MIN_PREV_). Debe aparecer en al menos ese porcentaje de muestras del grupo
MIN_PREV_BACT <- 0.30 #Base .2
MIN_PREV_EUK  <- 0.30 #Base .1

#porcentaje relativo dentro de cada dominio
MIN_MEAN_BACT <- 0.1    #0.06 base # % dentro de Bacteria
MIN_MEAN_EUK  <- 0.1    #0.015 base % dentro de Eukaryota

#el no. máximo de taxones en caso de que sea muy grande
TOP_N_BACT    <- 70
TOP_N_EUK     <- 50

# Correlación / Bootstrap
BOOT_R               <- 300
PSEUDO               <- 1e-6
# Umbrales de prescreen por relación
PRESCREEN_BE         <- 0.5 #0.35
PRESCREEN_BB         <- 0.5 #0.45
PRESCREEN_EE         <- 0.5 #0.30
PRESCREEN_WL         <- 0.5 #0.20   # si toca whitelist

#minimo de muestras para evaluar 
MIN_SAMPLES_PER_G    <- 4

# Co-ocurrencia mínima (muestras donde AMBOS > 0)
MIN_COOC_FRAC <- 0.50          # que la relacion ocurre en al menos ≥35% de las muestras del grupo
MIN_COOC_ABS  <- 4           # y al menos 20 muestras

# Filtro final por relación (umbral de |rho_mean| y FDR)
THRESH_BE <- 0.4
THRESH_BB <- 0.55
THRESH_EE <- 0.35
ALPHA_Q   <- 0.05              # BH FDR

# Whitelist de géneros clave
force_include_euk  <- c("Saccharomyces")
force_include_bact <- character(0)

## ======================= ENTRADAS =======================
# Se asume que ya cargaste 'kaiju_merged' en el ambiente.
# Si no, descomenta y pon tu ruta:
# kaiju_merged <- readr::read_tsv("tu_kaiju_merged.tsv", show_col_types = FALSE)

# Asegura tipos:
kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

# Si no existen columnas taxonómicas separadas pero sí 'taxon_name'
if (!all(c("Domain","Genus") %in% names(kaiju_merged)) && "taxon_name" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    tidyr::separate(
      taxon_name,
      into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
               "Subclass","Order","Family","Genus","Species"),
      sep = ";", fill = "right", extra = "drop"
    )
}

# file_base si hace falta
if (!"file_base" %in% names(kaiju_merged)) {
  if ("file" %in% names(kaiju_merged)) {
    kaiju_merged <- kaiju_merged %>% mutate(file_base = gsub("^.*/", "", file))
  } else {
    stop("No 'file_base' or 'file' column found in kaiju_merged.")
  }
}

## ============= METADATOS: construir BMI_group para BMI≥25 / BMI<25 =============
# Lee tu archivo de metadatos (ajusta la ruta si cambia)
meta <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE, sep = ",", fileEncoding = "latin1",
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

# Limpieza mínima
meta <- meta[, colSums(!is.na(meta)) > 0]
meta <- meta[rowSums(is.na(meta)) < ncol(meta), ]
meta$BMI <- suppressWarnings(as.numeric(meta$BMI))
meta$Age <- suppressWarnings(as.numeric(meta$Age))

# Agrupaciones de IMC (como lo tenías + colapsado a 2 niveles)
meta <- meta %>%
  mutate(
    BMI_group = case_when(
      Percentil_formulas %in% c("Overweight", "Obesity") ~ "Overweight/Obesity",
      Percentil_formulas == "Normal Weight" ~ "Normal Weight",
      Percentil_formulas == "Malnutrition" ~ "Malnutrition",
      TRUE ~ NA_character_
    ),
    BMI_index_group = case_when(
      BMI_group == "Overweight/Obesity" ~ ">25 BMI INDEX",
      BMI_group %in% c("Normal Weight", "Malnutrition") ~ "<25 BMI INDEX",
      TRUE ~ NA_character_
    )
  )

# En tus metadatos, el identificador de secuenciación parece estar en 'Lane'
# y los archivos de Kaiju acaban en "_kaiju.out". Ajusta si tu sufijo es distinto.
files_over25  <- meta %>% filter(BMI_index_group == ">25 BMI INDEX") %>% pull(Lane) %>% unique()
files_under25 <- meta %>% filter(BMI_index_group == "<25 BMI INDEX") %>% pull(Lane) %>% unique()

# Asignar BMI_group a kaiju_merged (etiquetas finales BMI≥25 / BMI<25)
kaiju_merged <- kaiju_merged %>%
  mutate(BMI_group = case_when(
    file_base %in% paste0(files_over25,  "_kaiju.out") ~ "BMI≥25",
    file_base %in% paste0(files_under25, "_kaiju.out") ~ "BMI<25",
    TRUE ~ NA_character_
  ))

## ============= (Opcional) asignar Group (Urban/Rural) si no existe =============
if (!"Group" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    mutate(Group = case_when(
      exists("urban_files") & file_base %in% paste0(urban_files, "_kaiju.out") ~ "Urban",
      exists("rural_files") & file_base %in% paste0(rural_files, "_kaiju.out") ~ "Rural",
      TRUE ~ NA_character_
    ))
}

## ===== NUEVO: construir context_group = SOLO Adult + BMI≥25 y dividir Urban/Rural =====
# 1) Lista de Lanes "Adult"
adult_lanes <- meta %>%
  filter(!is.na(Age_group), Age_group == "Adult",
         !is.na(Gender), Gender == "Male") %>%
  pull(Lane) %>%
  unique()


KAIJU_SUFFIX <- "_kaiju.out"
adult_files  <- paste0(adult_lanes, KAIJU_SUFFIX)

# 2) Columna context_group en kaiju_merged (solo Adult + BMI≥25; Urban/Rural aparte)
kaiju_merged <- kaiju_merged %>%
  dplyr::mutate(
    is_adult_file = file_base %in% adult_files,
    context_group = dplyr::case_when(
      is_adult_file & BMI_group == "BMI≥25" & Group == "Urban" ~ "BMI≥25 Urban",
      is_adult_file & BMI_group == "BMI≥25" & Group == "Rural" ~ "BMI≥25 Rural",
      TRUE ~ NA_character_
    )
  )

## ======================= Helpers genéricos =======================
# Limpieza base (solo Bacteria/Eukaryota a nivel GENUS; excluye Homo)
prep_dat <- function(df_all) {
  df_all %>%
    filter(Domain %in% c("Bacteria","Eukaryota")) %>%
    mutate(Genus = ifelse(is.na(Genus) | Genus=="" | Genus=="Unclassified", NA, str_trim(Genus))) %>%
    filter(!is.na(Genus)) %>%
    filter(!(Domain=="Eukaryota" & Genus=="Homo"))
}

# Métricas por dominio dentro de cada muestra (% dentro del DOMINIO)
domain_norm_tables <- function(dat, group_var = "Group") {
  stopifnot(all(c("file_base", group_var, "Domain", "Genus", "reads") %in% names(dat)))
  totals_domain <- dat %>%
    group_by(.data[[group_var]], file_base, Domain) %>%
    summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")
  genus_reads <- dat %>%
    group_by(.data[[group_var]], file_base, Domain, Genus) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(totals_domain, by = c(group_var, "file_base", "Domain")) %>%
    mutate(pct_domain = if_else(total_reads_domain > 0, 100 * reads / total_reads_domain, 0))
  list(genus_reads = genus_reads, totals_domain = totals_domain)
}

# Elegir géneros “keepers” por dominio (prevalencia y % medio) usando el group_var
pick_keepers <- function(genus_reads, group_var = "Group") {
  prev_mean <- genus_reads %>%
    group_by(.data[[group_var]], Domain, Genus) %>%
    summarise(prevalence = mean(reads > 0),
              mean_pct_dom = mean(pct_domain),
              .groups = "drop")
  
  # NOTA: asumimos 2 niveles en group_var; se ordenan por aparición
  lvls <- sort(unique(prev_mean[[group_var]]))
  
  keepers_bact <- prev_mean %>%
    filter(Domain=="Bacteria") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_BACT,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_BACT),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    { if (nrow(.)==0) { message("⚠ Sin métricas para Bacteria en pick_keepers()"); . } else . } %>%
    {
      # Depuración: ¿cuántos fallan por prevalencia vs media?
      tmp <- .
      if (nrow(tmp)) {
        if (all(!tmp$ok)) message("⚠ Ningún género bacteriano pasó los filtros: MIN_PREV_BACT=", MIN_PREV_BACT, 
                                  " y/o MIN_MEAN_BACT=", MIN_MEAN_BACT)
      }
      tmp
    } %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_BACT) %>% pull(Genus)
  
  keepers_euk <- prev_mean %>%
    filter(Domain=="Eukaryota") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_EUK,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_EUK),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    { if (nrow(.)==0) { message("⚠ Sin métricas para Eukaryota en pick_keepers()"); . } else . } %>%
    {
      tmp <- .
      if (nrow(tmp)) {
        if (all(!tmp$ok)) message("⚠ Ningún género eucariota pasó los filtros: MIN_PREV_EUK=", MIN_PREV_EUK, 
                                  " y/o MIN_MEAN_EUK=", MIN_MEAN_EUK)
      }
      tmp
    } %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_EUK) %>% pull(Genus)
  
  keepers_bact <- union(keepers_bact,
                        intersect(force_include_bact, unique(genus_reads$Genus[genus_reads$Domain=="Bacteria"])))
  keepers_euk  <- union(keepers_euk,
                        intersect(force_include_euk,  unique(genus_reads$Genus[genus_reads$Domain=="Eukaryota"])))
  
  if (length(keepers_bact)==0) message("⚠ keepers_bact vacío tras aplicar TOP_N_BACT=", TOP_N_BACT)
  if (length(keepers_euk)==0)  message("⚠ keepers_euk vacío tras aplicar TOP_N_EUK=", TOP_N_EUK)
  
  list(keepers_bact = keepers_bact, keepers_euk = keepers_euk)
}


# Matrices de counts por grupo (según group_var) y CLR
counts_from_long <- function(df_long, group_var, group_value) {
  df_long %>%
    filter(.data[[group_var]] == group_value) %>%
    select(file_base, Genus, reads) %>%
    group_by(file_base, Genus) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Genus, values_from = reads, values_fill = 0) %>%
    column_to_rownames("file_base") %>% as.matrix()
}

clr_from_counts <- function(counts_mat) {
  if (is.null(counts_mat) || nrow(counts_mat) < MIN_SAMPLES_PER_G) return(NULL)
  prop <- sweep(counts_mat, 1, rowSums(counts_mat), "/"); prop[is.na(prop)] <- 0
  logx <- log(prop + PSEUDO); sweep(logx, 1, rowMeans(logx), "-")
}

# -------- Correlaciones GENÉRICAS (BE, BB, EE) con bootstrap --------
fast_pairwise_corr <- function(mat_clr, counts_mat, genus_domain_map,
                               relation = c("BE","BB","EE"),
                               prescreen_global = list(BE = PRESCREEN_BE, BB = PRESCREEN_BB, EE = PRESCREEN_EE),
                               prescreen_wl     = PRESCREEN_WL,
                               min_cooc_frac = MIN_COOC_FRAC,
                               min_cooc_abs  = MIN_COOC_ABS,
                               boot_R = BOOT_R,
                               whitelist = character(0)) {
  
  relation <- match.arg(relation)
  cols <- colnames(mat_clr)
  dom_vec <- setNames(as.character(genus_domain_map$Domain),
                      as.character(genus_domain_map$Genus))
  
  bact_cols <- cols[dom_vec[cols] == "Bacteria"]
  euk_cols  <- cols[dom_vec[cols] == "Eukaryota"]
  
  if (relation == "BE" && (length(bact_cols) == 0 || length(euk_cols) == 0)) {
    message("⚠ Prescreen(", relation, "): no hay columnas suficientes Bacteria/Eukaryota para correlacionar.")
    return(tibble())
  }
  if (relation == "BB" && length(bact_cols) < 2) {
    message("⚠ Prescreen(", relation, "): menos de 2 géneros bacterianos tras keepers.")
    return(tibble())
  }
  if (relation == "EE" && length(euk_cols) < 2) {
    message("⚠ Prescreen(", relation, "): menos de 2 géneros eucariotas tras keepers.")
    return(tibble())
  }
  
  # Co-ocurrencia
  pres_abs <- max(min_cooc_abs, ceiling(min_cooc_frac * nrow(counts_mat)))
  pres_mat <- counts_mat > 0
  
  # Rangos para Spearman masivo
  mat_rank <- apply(mat_clr, 2, rank, ties.method = "average")
  
  if (relation == "BE") {
    R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                    mat_rank[, euk_cols,  drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    sel <- which(abs(R) >= prescreen_global$BE, arr.ind = TRUE)
    if (length(sel)==0) message("⚠ Prescreen(BE): ninguna pareja supera |rho|>=", prescreen_global$BE)
    pairs_df <- tibble(Genus_1 = bact_cols[sel[, "row"]],
                       Genus_2 = euk_cols[ sel[, "col"]],
                       from_whitelist = FALSE)
    # Whitelist laxo
    if (length(whitelist) > 0) {
      wl_b <- intersect(bact_cols, whitelist)
      wl_e <- intersect(euk_cols,  whitelist)
      if (length(wl_b) > 0) {
        R_wlb <- R[match(wl_b, bact_cols), , drop = FALSE]
        sel_wlb <- which(abs(R_wlb) >= prescreen_wl, arr.ind = TRUE)
        if (length(sel_wlb)==0) message("ℹ Prescreen(BE) whitelist-bacteria: ninguna pareja supera |rho|>=", prescreen_wl)
        if (length(sel_wlb))
          pairs_df <- bind_rows(pairs_df,
                                tibble(Genus_1 = wl_b[sel_wlb[, "row"]],
                                       Genus_2 = euk_cols[sel_wlb[, "col"]],
                                       from_whitelist = TRUE))
      }
      if (length(wl_e) > 0) {
        R_wle <- R[, match(wl_e, euk_cols), drop = FALSE]
        sel_wle <- which(abs(R_wle) >= prescreen_wl, arr.ind = TRUE)
        if (length(sel_wle)==0) message("ℹ Prescreen(BE) whitelist-euk: ninguna pareja supera |rho|>=", prescreen_wl)
        if (length(sel_wle))
          pairs_df <- bind_rows(pairs_df,
                                tibble(Genus_1 = bact_cols[sel_wle[, "row"]],
                                       Genus_2 = wl_e[sel_wle[, "col"]],
                                       from_whitelist = TRUE))
      }
    }
  } else if (relation == "BB") {
    R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    upper <- upper.tri(R, diag = FALSE)
    idx <- which(upper & abs(R) >= prescreen_global$BB, arr.ind = TRUE)
    if (length(idx)==0) message("⚠ Prescreen(BB): ninguna pareja supera |rho|>=", prescreen_global$BB)
    pairs_df <- tibble(Genus_1 = bact_cols[idx[, "row"]],
                       Genus_2 = bact_cols[idx[, "col"]],
                       from_whitelist = FALSE)
  } else { # EE
    R <- stats::cor(mat_rank[, euk_cols, drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    upper <- upper.tri(R, diag = FALSE)
    idx <- which(upper & abs(R) >= prescreen_global$EE, arr.ind = TRUE)
    if (length(idx)==0) message("⚠ Prescreen(EE): ninguna pareja supera |rho|>=", prescreen_global$EE)
    pairs_df <- tibble(Genus_1 = euk_cols[idx[, "row"]],
                       Genus_2 = euk_cols[idx[, "col"]],
                       from_whitelist = FALSE)
  }
  
  if (nrow(pairs_df) == 0) return(tibble())
  
  # Únicos + co-ocurrencia mínima
  pairs_df <- pairs_df %>%
    distinct(Genus_1, Genus_2, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(cooc = sum(pres_mat[, Genus_1] & pres_mat[, Genus_2])) %>%
    ungroup()
  
  if (!nrow(pairs_df)) return(tibble())
  n_before <- nrow(pairs_df)
  pairs_df <- pairs_df %>% filter(cooc >= pres_abs)
  if (nrow(pairs_df) == 0) {
    message("⚠ Filtro de co-ocurrencia (", relation, "): 0 de ", n_before,
            " parejas cumplen co-ocurrencia >= ", pres_abs, 
            " (", MIN_COOC_FRAC*100, "% y mínimo absoluto ", MIN_COOC_ABS, ")")
    return(tibble())
  }
  
  message("Prescreen(", relation, ") pairs (post co-occur): ", nrow(pairs_df),
          " (", sum(pairs_df$from_whitelist), " via whitelist)")
  
  if (nrow(pairs_df) == 0) return(tibble())
  
  boot_one <- function(g1, g2) {
    x <- mat_clr[, g1]; y <- mat_clr[, g2]
    if (sd(x) == 0 || sd(y) == 0) {
      message("ℹ Bootstrap(", relation, "): descarto pareja por varianza cero: ", g1, " ~ ", g2)
      return(NULL)
    }
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    rho0 <- unname(ct$estimate); p0 <- ct$p.value
    stat <- function(dd, idx) suppressWarnings(cor(dd[idx,1], dd[idx,2], method = "spearman"))
    bt <- boot(cbind(x, y), statistic = stat, R = boot_R)
    if (all(is.na(bt$t))) {
      message("ℹ Bootstrap(", relation, "): replicados NA para ", g1, " ~ ", g2)
      return(NULL)
    }
    ci <- tryCatch(boot.ci(bt, type = "perc"), error = function(e) NULL)
    scons <- mean(sign(bt$t) == sign(rho0), na.rm = TRUE)
    rel_label <- if (relation=="BE") "Bacteria-Eukaryota" else if (relation=="BB") "Bacteria-Bacteria" else "Eukaryota-Eukaryota"
    tibble(
      Genus_1 = g1, Genus_2 = g2, Relation = rel_label,
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

# BH helper
add_q <- function(df) if (nrow(df)) df %>% mutate(q = p.adjust(p, method = "BH")) else df

# Filtro final por relación
filter_by_relation <- function(df, relation = c("BE","BB","EE")) {
  relation <- match.arg(relation)
  if (nrow(df)==0) return(df[0,])
  
  thr <- switch(relation,
                "BE" = THRESH_BE,
                "BB" = THRESH_BB,
                "EE" = THRESH_EE)
  
  n_in <- nrow(df)
  res <- df %>%
    mutate(abs_rho = abs(rho_mean)) %>%
    mutate(ci_strong = (rho_mean > 0 & CI_lower >= thr) |
             (rho_mean < 0 & CI_upper <= -thr)) %>%
    filter(
      abs_rho >= thr,
      abs(rho) >= thr,
      ci_strong,
      sign_consistency >= 0.90,
      q <= ALPHA_Q
    ) %>%
    select(-abs_rho, -ci_strong)
  
  if (nrow(res)==0) {
    message("⚠ Filtro final (", relation, "): 0 de ", n_in, 
            " pasan: |rho_mean|>=", thr, ", |rho|>=", thr, 
            ", CI fuerte, sign_consistency>=0.90, q<=", ALPHA_Q)
  }
  res
}

make_nodes <- function(edges_df, genus_domain_map) {
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

# --- Mapea etiquetas a tags seguros solo para los nombres de archivo ---
label_to_tag <- function(lbl, group_var = NULL) {
  if (!is.null(group_var) && group_var == "BMI_group") {
    if (grepl("^\\s*BMI\\s*<\\s*25\\s*$", lbl))  return("BMI24")
    if (grepl("^\\s*BMI\\s*[≥>=]\\s*25\\s*$", lbl)) return("BMI26")
  }
  if (!is.null(group_var) && group_var == "context_group") {
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Urban", lbl)) return("BMI26Urban")
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Rural", lbl)) return("BMI26Rural")
  }
  out <- gsub("[^A-Za-z0-9]+", "", lbl)
  if (nzchar(out)) out else "grp"
}

relation_tag <- function(relation_label) {
  if (relation_label == "Bacteria-Eukaryota") return("BE")
  if (relation_label == "Bacteria-Bacteria")  return("BB")
  if (relation_label == "Eukaryota-Eukaryota")return("EE")
  "REL"
}

## ======================= FUNCIÓN “PIPELINE” GENERAL =======================
`%||%` <- function(a,b) if (!is.null(a)) a else b

run_pipeline_for <- function(df, group_var = "Group", labels = NULL, tag_prefix = NULL) {
  # Asegurar grupos válidos
  df <- df %>% filter(!is.na(.data[[group_var]]))
  if (length(unique(df[[group_var]])) != 2) {
    message("⚠ Saltando ", group_var, ": no hay exactamente 2 niveles válidos.")
    return(invisible(NULL))
  }
  
  if (is.null(labels)) {
    labels <- sort(unique(df[[group_var]]))
  }
  gA <- labels[1]; gB <- labels[2]
  
  message(">>> Ejecutando para ", group_var, " = {", gA, ", ", gB, "} …")
  
  dat <- prep_dat(df)
  if (!nrow(dat)) { message("⚠ prep_dat() dejó 0 filas para ", group_var); return(invisible(NULL)) }
  
  # % dentro de dominio
  tabs <- domain_norm_tables(dat, group_var = group_var)
  genus_reads <- tabs$genus_reads
  if (!nrow(genus_reads)) { message("⚠ domain_norm_tables() devolvió 0 filas (", group_var, ")"); return(invisible(NULL)) }
  
  # keepers
  k <- pick_keepers(genus_reads, group_var = group_var)
  keepers_bact <- k$keepers_bact; keepers_euk <- k$keepers_euk
  message("Kept genera — Bacteria=", length(keepers_bact), " | Eukaryota=", length(keepers_euk))
  if (length(keepers_bact)==0 || length(keepers_euk)==0) {
    message("⚠ Filtro keepers: requiere al menos 1 género por dominio. Bact=", length(keepers_bact), ", Euk=", length(keepers_euk))
  }
  
  # Filtrado a keepers
  dat_filt <- genus_reads %>%
    filter((Domain=="Bacteria"  & Genus %in% keepers_bact) |
             (Domain=="Eukaryota" & Genus %in% keepers_euk))
  if (!nrow(dat_filt)) { message("⚠ dat_filt vacío tras keepers para ", group_var); return(invisible(NULL)) }
  
  # Matrices por grupo
  cnt_A <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gA)
  cnt_B <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gB)
  if (is.null(cnt_A) || nrow(cnt_A) < MIN_SAMPLES_PER_G) message("⚠ ", gA, ": no pasa filtro de cantidad de muestras (n=", ifelse(is.null(cnt_A), 0, nrow(cnt_A)), " < ", MIN_SAMPLES_PER_G, ")")
  if (is.null(cnt_B) || nrow(cnt_B) < MIN_SAMPLES_PER_G) message("⚠ ", gB, ": no pasa filtro de cantidad de muestras (n=", ifelse(is.null(cnt_B), 0, nrow(cnt_B)), " < ", MIN_SAMPLES_PER_G, ")")
  stopifnot(!is.null(cnt_A), !is.null(cnt_B), nrow(cnt_A) >= MIN_SAMPLES_PER_G, nrow(cnt_B) >= MIN_SAMPLES_PER_G)
  
  mat_A <- clr_from_counts(cnt_A)
  mat_B <- clr_from_counts(cnt_B)
  if (is.null(mat_A)) message("⚠ ", gA, ": clr_from_counts() devolvió NULL (posible n<mínimo).")
  if (is.null(mat_B)) message("⚠ ", gB, ": clr_from_counts() devolvió NULL (posible n<mínimo).")
  
  # Alinear columnas
  common_cols <- Reduce(intersect, list(colnames(mat_A), colnames(mat_B)))
  if (length(common_cols) == 0) { message("⚠ Sin columnas comunes de géneros entre ", gA, " y ", gB); return(invisible(NULL)) }
  mat_A <- mat_A[, common_cols, drop = FALSE]
  mat_B <- mat_B[, common_cols, drop = FALSE]
  cnt_A <- cnt_A[, common_cols, drop = FALSE]
  cnt_B <- cnt_B[, common_cols, drop = FALSE]
  
  # Genus -> Domain map
  genus_domain_map <- dat_filt %>%
    distinct(Genus, Domain) %>%
    filter(Genus %in% common_cols)
  if (!nrow(genus_domain_map)) { message("⚠ genus_domain_map vacío tras intersección de columnas."); return(invisible(NULL)) }
  
  # --------- Correlaciones y filtros por RELACIÓN: BE, BB, EE ----------
  relations <- c("BE","BB","EE")
  for (rel in relations) {
    message("Computing ", rel, " correlations (", gA, ")…")
    edges_A <- fast_pairwise_corr(
      mat_A, cnt_A, genus_domain_map,
      relation = rel,
      whitelist = union(force_include_euk, force_include_bact),
      boot_R = BOOT_R
    )
    message("Computing ", rel, " correlations (", gB, ")…")
    edges_B <- fast_pairwise_corr(
      mat_B, cnt_B, genus_domain_map,
      relation = rel,
      whitelist = union(force_include_euk, force_include_bact),
      boot_R = BOOT_R
    )
    
    edges_A <- add_q(edges_A); edges_B <- add_q(edges_B)
    
    edges_A_f <- filter_by_relation(edges_A, relation = rel)
    edges_B_f <- filter_by_relation(edges_B, relation = rel)
    
    message(rel, " edges kept | ", gA, ": ", nrow(edges_A_f), " | ", gB, ": ", nrow(edges_B_f))
    
    nodes_A <- make_nodes(edges_A_f, genus_domain_map)
    nodes_B <- make_nodes(edges_B_f, genus_domain_map)
    
    # Tags de salida
    prefix <- (tag_prefix %||% tolower(group_var))
    rel_tag <- rel
    tagA <- paste0(prefix, "_", label_to_tag(gA, group_var), "_", rel_tag, "_genus_tight")
    tagB <- paste0(prefix, "_", label_to_tag(gB, group_var), "_", rel_tag, "_genus_tight")
    
    write_net(edges_A_f, nodes_A, tagA)
    write_net(edges_B_f, nodes_B, tagB)
  }
  
  invisible(NULL)
}

## ======================= CORRER PIPELINE PARA LOS 3 CASOS =======================
# 1) Urban vs Rural
if ("Group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$Group)) && dplyr::n_distinct(na.omit(kaiju_merged$Group)) == 2) {
  run_pipeline_for(
    kaiju_merged %>% filter(!is.na(Group)),
    group_var = "Group",
    labels = c("Rural","Urban"),         # fuerza el orden (opcional)
    tag_prefix = "group"
  )
} else {
  message("Saltando Group: no hay dos niveles en 'Group'.")
}

# 2) BMI≥25 vs BMI<25
if ("BMI_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$BMI_group))) {
  run_pipeline_for(
    kaiju_merged %>% filter(!is.na(BMI_group)),
    group_var = "BMI_group",
    labels = c("BMI<25","BMI≥25"),       # orden consistente
    tag_prefix = "bmi"
  )
} else {
  message("No hay 'BMI_group' asignado a kaiju_merged (revisa mapeo Lane->file_base).")
}

# 3) context_group: SOLO Adult + BMI≥25, comparando Urban vs Rural
if ("context_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$context_group))) {
  lvls_ctx <- sort(unique(na.omit(kaiju_merged$context_group)))
  if (length(lvls_ctx) == 2) {
    run_pipeline_for(
      kaiju_merged %>% filter(!is.na(context_group)),
      group_var = "context_group",
      labels = c("BMI≥25 Rural","BMI≥25 Urban"),  # orden explícito
      tag_prefix = "context"
    )
  } else {
    message("Saltando context_group: niveles detectados = ", paste(lvls_ctx, collapse=", "))
  }
} else {
  message("No hay 'context_group' (Adult + BMI≥25) disponible.")
}

message("Listo. Redes BE/BB/EE exportadas para Group, BMI_group y context_group (si disponibles).")





































## =========================================================
## Redes GENUS (Bacteria/Eukaryota) por:
##   1) Group (Urban vs Rural)
##   2) BMI_group (BMI≥25 vs BMI<25)  -> nombres: BMI<25 -> BMI24 ; BMI≥25 -> BMI26
##   3) context_group (SOLO Adult + BMI≥25):
##        - BMI≥25 Urban  -> BMI26Urban (solo en nombres de archivo)
##        - BMI≥25 Rural  -> BMI26Rural (solo en nombres de archivo)
## Además, para cada esquema se generan 3 redes: BE, BB y EE.
## =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(boot); library(readr)
})

set.seed(1234)

## ======================= PARÁMETROS GLOBALES =======================
# Filtros por grupo (prevalencia y % medio, NORMALIZADO POR DOMINIO)
#Prevalencia mínima (MIN_PREV_). Debe aparecer en al menos ese porcentaje de muestras del grupo
MIN_PREV_BACT <- 0.20 #Base .2
MIN_PREV_EUK  <- 0.10 #Base .1

#porcentaje relativo dentro de cada dominio
MIN_MEAN_BACT <- 0.1    #0.06 base # % dentro de Bacteria
MIN_MEAN_EUK  <- 0.1    #0.015 base % dentro de Eukaryota

#el no. máximo de taxones en caso de que sea muy grande
TOP_N_BACT    <- 70
TOP_N_EUK     <- 50

# Correlación / Bootstrap
BOOT_R               <- 300
PSEUDO               <- 1e-6
# Umbrales de prescreen por relación
PRESCREEN_BE         <- 0.5 #0.35
PRESCREEN_BB         <- 0.5 #0.45
PRESCREEN_EE         <- 0.5 #0.30
PRESCREEN_WL         <- 0.5 #0.20   # si toca whitelist

#minimo de muestras para evaluar 
MIN_SAMPLES_PER_G    <- 4

# Co-ocurrencia mínima (muestras donde AMBOS > 0)
MIN_COOC_FRAC <- 0.50          # que la relacion ocurre en al menos ≥35% de las muestras del grupo
MIN_COOC_ABS  <- 4           # y al menos 20 muestras

# Filtro final por relación (umbral de |rho_mean| y FDR)
THRESH_BE <- 0.4
THRESH_BB <- 0.55
THRESH_EE <- 0.35
ALPHA_Q   <- 0.05              # BH FDR

# Whitelist de géneros clave
force_include_euk  <- c("Saccharomyces")
force_include_bact <- character(0)

## ======================= ENTRADAS =======================
# Se asume que ya cargaste 'kaiju_merged' en el ambiente.
# Si no, descomenta y pon tu ruta:
# kaiju_merged <- readr::read_tsv("tu_kaiju_merged.tsv", show_col_types = FALSE)

# Asegura tipos:
kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

# Si no existen columnas taxonómicas separadas pero sí 'taxon_name'
if (!all(c("Domain","Genus") %in% names(kaiju_merged)) && "taxon_name" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    tidyr::separate(
      taxon_name,
      into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
               "Subclass","Order","Family","Genus","Species"),
      sep = ";", fill = "right", extra = "drop"
    )
}

# file_base si hace falta
if (!"file_base" %in% names(kaiju_merged)) {
  if ("file" %in% names(kaiju_merged)) {
    kaiju_merged <- kaiju_merged %>% mutate(file_base = gsub("^.*/", "", file))
  } else {
    stop("No 'file_base' or 'file' column found in kaiju_merged.")
  }
}

## ============= METADATOS: construir BMI_group para BMI≥25 / BMI<25 =============
# Lee tu archivo de metadatos (ajusta la ruta si cambia)
meta <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE, sep = ",", fileEncoding = "latin1",
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

# Limpieza mínima
meta <- meta[, colSums(!is.na(meta)) > 0]
meta <- meta[rowSums(is.na(meta)) < ncol(meta), ]
meta$BMI <- suppressWarnings(as.numeric(meta$BMI))
meta$Age <- suppressWarnings(as.numeric(meta$Age))

# Agrupaciones de IMC (como lo tenías + colapsado a 2 niveles)
meta <- meta %>%
  mutate(
    BMI_group = case_when(
      Percentil_formulas %in% c("Overweight", "Obesity") ~ "Overweight/Obesity",
      Percentil_formulas == "Normal Weight" ~ "Normal Weight",
      Percentil_formulas == "Malnutrition" ~ "Malnutrition",
      TRUE ~ NA_character_
    ),
    BMI_index_group = case_when(
      BMI_group == "Overweight/Obesity" ~ ">25 BMI INDEX",
      BMI_group %in% c("Normal Weight", "Malnutrition") ~ "<25 BMI INDEX",
      TRUE ~ NA_character_
    )
  )

# En tus metadatos, el identificador de secuenciación parece estar en 'Lane'
# y los archivos de Kaiju acaban en "_kaiju.out". Ajusta si tu sufijo es distinto.
files_over25  <- meta %>% filter(BMI_index_group == ">25 BMI INDEX") %>% pull(Lane) %>% unique()
files_under25 <- meta %>% filter(BMI_index_group == "<25 BMI INDEX") %>% pull(Lane) %>% unique()

# Asignar BMI_group a kaiju_merged (etiquetas finales BMI≥25 / BMI<25)
kaiju_merged <- kaiju_merged %>%
  mutate(BMI_group = case_when(
    file_base %in% paste0(files_over25,  "_kaiju.out") ~ "BMI≥25",
    file_base %in% paste0(files_under25, "_kaiju.out") ~ "BMI<25",
    TRUE ~ NA_character_
  ))

## ============= (Opcional) asignar Group (Urban/Rural) si no existe =============
if (!"Group" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    mutate(Group = case_when(
      exists("urban_files") & file_base %in% paste0(urban_files, "_kaiju.out") ~ "Urban",
      exists("rural_files") & file_base %in% paste0(rural_files, "_kaiju.out") ~ "Rural",
      TRUE ~ NA_character_
    ))
}

## ===== NUEVO: construir context_group = SOLO Adult + BMI≥25 y dividir Urban/Rural =====
# 1) Lista de Lanes "Adult"
adult_lanes <- meta %>%
  filter(!is.na(Age_group), Age_group == "Adult",
         !is.na(Gender), Gender == "Male") %>%
  pull(Lane) %>%
  unique()


KAIJU_SUFFIX <- "_kaiju.out"
adult_files  <- paste0(adult_lanes, KAIJU_SUFFIX)

# 2) Columna context_group en kaiju_merged (solo Adult + BMI≥25; Urban/Rural aparte)
kaiju_merged <- kaiju_merged %>%
  dplyr::mutate(
    is_adult_file = file_base %in% adult_files,
    context_group = dplyr::case_when(
      is_adult_file & BMI_group == "BMI≥25" & Group == "Urban" ~ "BMI≥25 Urban",
      is_adult_file & BMI_group == "BMI≥25" & Group == "Rural" ~ "BMI≥25 Rural",
      TRUE ~ NA_character_
    )
  )

## ======================= Helpers genéricos =======================
# Limpieza base (solo Bacteria/Eukaryota a nivel GENUS; excluye Homo)
prep_dat <- function(df_all) {
  df_all %>%
    filter(Domain %in% c("Bacteria","Eukaryota")) %>%
    mutate(Genus = ifelse(is.na(Genus) | Genus=="" | Genus=="Unclassified", NA, str_trim(Genus))) %>%
    filter(!is.na(Genus)) %>%
    filter(!(Domain=="Eukaryota" & Genus=="Homo"))
}

# Métricas por dominio dentro de cada muestra (% dentro del DOMINIO)
domain_norm_tables <- function(dat, group_var = "Group") {
  stopifnot(all(c("file_base", group_var, "Domain", "Genus", "reads") %in% names(dat)))
  totals_domain <- dat %>%
    group_by(.data[[group_var]], file_base, Domain) %>%
    summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")
  genus_reads <- dat %>%
    group_by(.data[[group_var]], file_base, Domain, Genus) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(totals_domain, by = c(group_var, "file_base", "Domain")) %>%
    mutate(pct_domain = if_else(total_reads_domain > 0, 100 * reads / total_reads_domain, 0))
  list(genus_reads = genus_reads, totals_domain = totals_domain)
}

# Elegir géneros “keepers” por dominio (prevalencia y % medio) usando el group_var
pick_keepers <- function(genus_reads, group_var = "Group") {
  prev_mean <- genus_reads %>%
    group_by(.data[[group_var]], Domain, Genus) %>%
    summarise(prevalence = mean(reads > 0),
              mean_pct_dom = mean(pct_domain),
              .groups = "drop")
  
  # NOTA: asumimos 2 niveles en group_var; se ordenan por aparición
  lvls <- sort(unique(prev_mean[[group_var]]))
  
  keepers_bact <- prev_mean %>%
    filter(Domain=="Bacteria") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_BACT,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_BACT),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_BACT) %>% pull(Genus)
  
  keepers_euk <- prev_mean %>%
    filter(Domain=="Eukaryota") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_EUK,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_EUK),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_EUK) %>% pull(Genus)
  
  keepers_bact <- union(keepers_bact,
                        intersect(force_include_bact, unique(genus_reads$Genus[genus_reads$Domain=="Bacteria"])))
  keepers_euk  <- union(keepers_euk,
                        intersect(force_include_euk,  unique(genus_reads$Genus[genus_reads$Domain=="Eukaryota"])))
  
  list(keepers_bact = keepers_bact, keepers_euk = keepers_euk)
}

# Matrices de counts por grupo (según group_var) y CLR
counts_from_long <- function(df_long, group_var, group_value) {
  df_long %>%
    filter(.data[[group_var]] == group_value) %>%
    select(file_base, Genus, reads) %>%
    group_by(file_base, Genus) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Genus, values_from = reads, values_fill = 0) %>%
    column_to_rownames("file_base") %>% as.matrix()
}

clr_from_counts <- function(counts_mat) {
  if (is.null(counts_mat) || nrow(counts_mat) < MIN_SAMPLES_PER_G) return(NULL)
  prop <- sweep(counts_mat, 1, rowSums(counts_mat), "/"); prop[is.na(prop)] <- 0
  logx <- log(prop + PSEUDO); sweep(logx, 1, rowMeans(logx), "-")
}

# -------- Correlaciones GENÉRICAS (BE, BB, EE) con bootstrap --------
fast_pairwise_corr <- function(mat_clr, counts_mat, genus_domain_map,
                               relation = c("BE","BB","EE"),
                               prescreen_global = list(BE = PRESCREEN_BE, BB = PRESCREEN_BB, EE = PRESCREEN_EE),
                               prescreen_wl     = PRESCREEN_WL,
                               min_cooc_frac = MIN_COOC_FRAC,
                               min_cooc_abs  = MIN_COOC_ABS,
                               boot_R = BOOT_R,
                               whitelist = character(0)) {
  
  relation <- match.arg(relation)
  cols <- colnames(mat_clr)
  dom_vec <- setNames(as.character(genus_domain_map$Domain),
                      as.character(genus_domain_map$Genus))
  
  bact_cols <- cols[dom_vec[cols] == "Bacteria"]
  euk_cols  <- cols[dom_vec[cols] == "Eukaryota"]
  
  if (relation == "BE" && (length(bact_cols) == 0 || length(euk_cols) == 0)) {
    message("⚠ Prescreen(", relation, "): no hay columnas suficientes Bacteria/Eukaryota para correlacionar.")
    return(tibble())
  }
  if (relation == "BB" && length(bact_cols) < 2) {
    message("⚠ Prescreen(", relation, "): menos de 2 géneros bacterianos tras keepers.")
    return(tibble())
  }
  if (relation == "EE" && length(euk_cols) < 2) {
    message("⚠ Prescreen(", relation, "): menos de 2 géneros eucariotas tras keepers.")
    return(tibble())
  }
  
  # Co-ocurrencia
  pres_abs <- max(min_cooc_abs, ceiling(min_cooc_frac * nrow(counts_mat)))
  pres_mat <- counts_mat > 0
  
  # Rangos para Spearman masivo
  mat_rank <- apply(mat_clr, 2, rank, ties.method = "average")
  
  if (relation == "BE") {
    R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                    mat_rank[, euk_cols,  drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    sel <- which(abs(R) >= prescreen_global$BE, arr.ind = TRUE)
    if (length(sel)==0) message("⚠ Prescreen(BE): ninguna pareja supera |rho|>=", prescreen_global$BE)
    pairs_df <- tibble(Genus_1 = bact_cols[sel[, "row"]],
                       Genus_2 = euk_cols[ sel[, "col"]],
                       from_whitelist = FALSE)
    # Whitelist laxo
    if (length(whitelist) > 0) {
      wl_b <- intersect(bact_cols, whitelist)
      wl_e <- intersect(euk_cols,  whitelist)
      if (length(wl_b) > 0) {
        R_wlb <- R[match(wl_b, bact_cols), , drop = FALSE]
        sel_wlb <- which(abs(R_wlb) >= prescreen_wl, arr.ind = TRUE)
        if (length(sel_wlb)==0) message("ℹ Prescreen(BE) whitelist-bacteria: ninguna pareja supera |rho|>=", prescreen_wl)
        if (length(sel_wlb))
          pairs_df <- bind_rows(pairs_df,
                                tibble(Genus_1 = wl_b[sel_wlb[, "row"]],
                                       Genus_2 = euk_cols[sel_wlb[, "col"]],
                                       from_whitelist = TRUE))
      }
      if (length(wl_e) > 0) {
        R_wle <- R[, match(wl_e, euk_cols), drop = FALSE]
        sel_wle <- which(abs(R_wle) >= prescreen_wl, arr.ind = TRUE)
        if (length(sel_wle)==0) message("ℹ Prescreen(BE) whitelist-euk: ninguna pareja supera |rho|>=", prescreen_wl)
        if (length(sel_wle))
          pairs_df <- bind_rows(pairs_df,
                                tibble(Genus_1 = bact_cols[sel_wle[, "row"]],
                                       Genus_2 = wl_e[sel_wle[, "col"]],
                                       from_whitelist = TRUE))
      }
    }
  } else if (relation == "BB") {
    R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    upper <- upper.tri(R, diag = FALSE)
    idx <- which(upper & abs(R) >= prescreen_global$BB, arr.ind = TRUE)
    if (length(idx)==0) message("⚠ Prescreen(BB): ninguna pareja supera |rho|>=", prescreen_global$BB)
    pairs_df <- tibble(Genus_1 = bact_cols[idx[, "row"]],
                       Genus_2 = bact_cols[idx[, "col"]],
                       from_whitelist = FALSE)
  } else { # EE
    R <- stats::cor(mat_rank[, euk_cols, drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    upper <- upper.tri(R, diag = FALSE)
    idx <- which(upper & abs(R) >= prescreen_global$EE, arr.ind = TRUE)
    if (length(idx)==0) message("⚠ Prescreen(EE): ninguna pareja supera |rho|>=", prescreen_global$EE)
    pairs_df <- tibble(Genus_1 = euk_cols[idx[, "row"]],
                       Genus_2 = euk_cols[idx[, "col"]],
                       from_whitelist = FALSE)
  }
  
  if (nrow(pairs_df) == 0) return(tibble())
  
  # Únicos + co-ocurrencia mínima
  pairs_df <- pairs_df %>%
    distinct(Genus_1, Genus_2, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(cooc = sum(pres_mat[, Genus_1] & pres_mat[, Genus_2])) %>%
    ungroup()
  
  if (!nrow(pairs_df)) return(tibble())
  n_before <- nrow(pairs_df)
  pairs_df <- pairs_df %>% filter(cooc >= pres_abs)
  if (nrow(pairs_df) == 0) {
    message("⚠ Filtro de co-ocurrencia (", relation, "): 0 de ", n_before,
            " parejas cumplen co-ocurrencia >= ", pres_abs, 
            " (", MIN_COOC_FRAC*100, "% y mínimo absoluto ", MIN_COOC_ABS, ")")
    return(tibble())
  }
  
  message("Prescreen(", relation, ") pairs (post co-occur): ", nrow(pairs_df),
          " (", sum(pairs_df$from_whitelist), " via whitelist)")
  
  if (nrow(pairs_df) == 0) return(tibble())
  
  boot_one <- function(g1, g2) {
    x <- mat_clr[, g1]; y <- mat_clr[, g2]
    if (sd(x) == 0 || sd(y) == 0) {
      message("ℹ Bootstrap(", relation, "): descarto pareja por varianza cero: ", g1, " ~ ", g2)
      return(NULL)
    }
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    rho0 <- unname(ct$estimate); p0 <- ct$p.value
    stat <- function(dd, idx) suppressWarnings(cor(dd[idx,1], dd[idx,2], method = "spearman"))
    bt <- boot(cbind(x, y), statistic = stat, R = boot_R)
    if (all(is.na(bt$t))) {
      message("ℹ Bootstrap(", relation, "): replicados NA para ", g1, " ~ ", g2)
      return(NULL)
    }
    ci <- tryCatch(boot.ci(bt, type = "perc"), error = function(e) NULL)
    scons <- mean(sign(bt$t) == sign(rho0), na.rm = TRUE)
    rel_label <- if (relation=="BE") "Bacteria-Eukaryota" else if (relation=="BB") "Bacteria-Bacteria" else "Eukaryota-Eukaryota"
    tibble(
      Genus_1 = g1, Genus_2 = g2, Relation = rel_label,
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

# BH helper
add_q <- function(df) if (nrow(df)) df %>% mutate(q = p.adjust(p, method = "BH")) else df

# Filtro final por relación
filter_by_relation <- function(df, relation = c("BE","BB","EE")) {
  relation <- match.arg(relation)
  if (nrow(df)==0) return(df[0,])
  
  thr <- switch(relation,
                "BE" = THRESH_BE,
                "BB" = THRESH_BB,
                "EE" = THRESH_EE)
  
  n_in <- nrow(df)
  res <- df %>%
    mutate(abs_rho = abs(rho_mean)) %>%
    mutate(ci_strong = (rho_mean > 0 & CI_lower >= thr) |
             (rho_mean < 0 & CI_upper <= -thr)) %>%
    filter(
      abs_rho >= thr,
      abs(rho) >= thr,
      ci_strong,
      sign_consistency >= 0.90,
      q <= ALPHA_Q
    ) %>%
    select(-abs_rho, -ci_strong)
  
  if (nrow(res)==0) {
    message("⚠ Filtro final (", relation, "): 0 de ", n_in, 
            " pasan: |rho_mean|>=", thr, ", |rho|>=", thr, 
            ", CI fuerte, sign_consistency>=0.90, q<=", ALPHA_Q)
  }
  res
}

make_nodes <- function(edges_df, genus_domain_map) {
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

# --- Mapea etiquetas a tags seguros solo para los nombres de archivo ---
label_to_tag <- function(lbl, group_var = NULL) {
  if (!is.null(group_var) && group_var == "BMI_group") {
    if (grepl("^\\s*BMI\\s*<\\s*25\\s*$", lbl))  return("BMI24")
    if (grepl("^\\s*BMI\\s*[≥>=]\\s*25\\s*$", lbl)) return("BMI26")
  }
  if (!is.null(group_var) && group_var == "context_group") {
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Urban", lbl)) return("BMI26Urban")
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Rural", lbl)) return("BMI26Rural")
  }
  out <- gsub("[^A-Za-z0-9]+", "", lbl)
  if (nzchar(out)) out else "grp"
}

relation_tag <- function(relation_label) {
  if (relation_label == "Bacteria-Eukaryota") return("BE")
  if (relation_label == "Bacteria-Bacteria")  return("BB")
  if (relation_label == "Eukaryota-Eukaryota")return("EE")
  "REL"
}

## ======================= FUNCIÓN “PIPELINE” GENERAL =======================
`%||%` <- function(a,b) if (!is.null(a)) a else b

run_pipeline_for <- function(df, group_var = "Group", labels = NULL, tag_prefix = NULL) {
  # Asegurar grupos válidos
  df <- df %>% filter(!is.na(.data[[group_var]]))
  if (length(unique(df[[group_var]])) != 2) {
    message("⚠ Saltando ", group_var, ": no hay exactamente 2 niveles válidos.")
    return(invisible(NULL))
  }
  
  if (is.null(labels)) {
    labels <- sort(unique(df[[group_var]]))
  }
  gA <- labels[1]; gB <- labels[2]
  
  message(">>> Ejecutando para ", group_var, " = {", gA, ", ", gB, "} …")
  
  dat <- prep_dat(df)
  if (!nrow(dat)) { message("⚠ prep_dat() dejó 0 filas para ", group_var); return(invisible(NULL)) }
  
  # % dentro de dominio
  tabs <- domain_norm_tables(dat, group_var = group_var)
  genus_reads <- tabs$genus_reads
  if (!nrow(genus_reads)) { message("⚠ domain_norm_tables() devolvió 0 filas (", group_var, ")"); return(invisible(NULL)) }
  
  # keepers
  k <- pick_keepers(genus_reads, group_var = group_var)
  keepers_bact <- k$keepers_bact; keepers_euk <- k$keepers_euk
  message("Kept genera — Bacteria=", length(keepers_bact), " | Eukaryota=", length(keepers_euk))
  if (length(keepers_bact)==0 || length(keepers_euk)==0) {
    message("⚠ Filtro keepers: requiere al menos 1 género por dominio. Bact=", length(keepers_bact), ", Euk=", length(keepers_euk))
  }
  
  # Filtrado a keepers
  dat_filt <- genus_reads %>%
    filter((Domain=="Bacteria"  & Genus %in% keepers_bact) |
             (Domain=="Eukaryota" & Genus %in% keepers_euk))
  if (!nrow(dat_filt)) { message("⚠ dat_filt vacío tras keepers para ", group_var); return(invisible(NULL)) }
  
  # Matrices por grupo
  cnt_A <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gA)
  cnt_B <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gB)
  if (is.null(cnt_A) || nrow(cnt_A) < MIN_SAMPLES_PER_G) message("⚠ ", gA, ": no pasa filtro de cantidad de muestras (n=", ifelse(is.null(cnt_A), 0, nrow(cnt_A)), " < ", MIN_SAMPLES_PER_G, ")")
  if (is.null(cnt_B) || nrow(cnt_B) < MIN_SAMPLES_PER_G) message("⚠ ", gB, ": no pasa filtro de cantidad de muestras (n=", ifelse(is.null(cnt_B), 0, nrow(cnt_B)), " < ", MIN_SAMPLES_PER_G, ")")
  stopifnot(!is.null(cnt_A), !is.null(cnt_B), nrow(cnt_A) >= MIN_SAMPLES_PER_G, nrow(cnt_B) >= MIN_SAMPLES_PER_G)
  
  mat_A <- clr_from_counts(cnt_A)
  mat_B <- clr_from_counts(cnt_B)
  if (is.null(mat_A)) message("⚠ ", gA, ": clr_from_counts() devolvió NULL (posible n<mínimo).")
  if (is.null(mat_B)) message("⚠ ", gB, ": clr_from_counts() devolvió NULL (posible n<mínimo).")
  
  # Alinear columnas
  common_cols <- Reduce(intersect, list(colnames(mat_A), colnames(mat_B)))
  if (length(common_cols) == 0) { message("⚠ Sin columnas comunes de géneros entre ", gA, " y ", gB); return(invisible(NULL)) }
  mat_A <- mat_A[, common_cols, drop = FALSE]
  mat_B <- mat_B[, common_cols, drop = FALSE]
  cnt_A <- cnt_A[, common_cols, drop = FALSE]
  cnt_B <- cnt_B[, common_cols, drop = FALSE]
  
  # Genus -> Domain map
  genus_domain_map <- dat_filt %>%
    distinct(Genus, Domain) %>%
    filter(Genus %in% common_cols)
  if (!nrow(genus_domain_map)) { message("⚠ genus_domain_map vacío tras intersección de columnas."); return(invisible(NULL)) }
  
  # --------- Correlaciones y filtros por RELACIÓN: BE, BB, EE ----------
  relations <- c("BE","BB","EE")
  for (rel in relations) {
    message("Computing ", rel, " correlations (", gA, ")…")
    edges_A <- fast_pairwise_corr(
      mat_A, cnt_A, genus_domain_map,
      relation = rel,
      whitelist = union(force_include_euk, force_include_bact),
      boot_R = BOOT_R
    )
    message("Computing ", rel, " correlations (", gB, ")…")
    edges_B <- fast_pairwise_corr(
      mat_B, cnt_B, genus_domain_map,
      relation = rel,
      whitelist = union(force_include_euk, force_include_bact),
      boot_R = BOOT_R
    )
    
    edges_A <- add_q(edges_A); edges_B <- add_q(edges_B)
    
    edges_A_f <- filter_by_relation(edges_A, relation = rel)
    edges_B_f <- filter_by_relation(edges_B, relation = rel)
    
    message(rel, " edges kept | ", gA, ": ", nrow(edges_A_f), " | ", gB, ": ", nrow(edges_B_f))
    
    nodes_A <- make_nodes(edges_A_f, genus_domain_map)
    nodes_B <- make_nodes(edges_B_f, genus_domain_map)
    
    # Tags de salida
    prefix <- (tag_prefix %||% tolower(group_var))
    rel_tag <- rel
    tagA <- paste0(prefix, "_", label_to_tag(gA, group_var), "_", rel_tag, "_genus_tight")
    tagB <- paste0(prefix, "_", label_to_tag(gB, group_var), "_", rel_tag, "_genus_tight")
    
    write_net(edges_A_f, nodes_A, tagA)
    write_net(edges_B_f, nodes_B, tagB)
  }
  
  invisible(NULL)
}

## ======================= CORRER PIPELINE PARA LOS 3 CASOS =======================
# 1) Urban vs Rural
if ("Group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$Group)) && dplyr::n_distinct(na.omit(kaiju_merged$Group)) == 2) {
  run_pipeline_for(
    kaiju_merged %>% filter(!is.na(Group)),
    group_var = "Group",
    labels = c("Rural","Urban"),         # fuerza el orden (opcional)
    tag_prefix = "group"
  )
} else {
  message("Saltando Group: no hay dos niveles en 'Group'.")
}

# 2) BMI≥25 vs BMI<25
if ("BMI_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$BMI_group))) {
  run_pipeline_for(
    kaiju_merged %>% filter(!is.na(BMI_group)),
    group_var = "BMI_group",
    labels = c("BMI<25","BMI≥25"),       # orden consistente
    tag_prefix = "bmi"
  )
} else {
  message("No hay 'BMI_group' asignado a kaiju_merged (revisa mapeo Lane->file_base).")
}

# 3) context_group: SOLO Adult + BMI≥25, comparando Urban vs Rural
if ("context_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$context_group))) {
  lvls_ctx <- sort(unique(na.omit(kaiju_merged$context_group)))
  if (length(lvls_ctx) == 2) {
    run_pipeline_for(
      kaiju_merged %>% filter(!is.na(context_group)),
      group_var = "context_group",
      labels = c("BMI≥25 Rural","BMI≥25 Urban"),  # orden explícito
      tag_prefix = "context"
    )
  } else {
    message("Saltando context_group: niveles detectados = ", paste(lvls_ctx, collapse=", "))
  }
} else {
  message("No hay 'context_group' (Adult + BMI≥25) disponible.")
}

message("Listo. Redes BE/BB/EE exportadas para Group, BMI_group y context_group (si disponibles).")

## ======================= SECCIÓN ADICIONAL: SOPORTE POR EDGE =======================
## Para cada archivo 'edges_<tag>.txt' generado:
## - Calcula, por grupo y relación, el % de muestras con co-ocurrencia (>0) del par (prevalencia de la relación).
## - Reporta tamaño muestral N.
## - Usa 'sign_consistency' como % de “seguridad” (consistencia bootstrap del signo).
## - Escribe 'edge_support_<tag>.txt' y mensajes-resumen.

try({
  dat_all <- prep_dat(kaiju_merged)
  schemes <- list(
    list(var="Group",        labels=c("Rural","Urban"),             prefix="group"),
    list(var="BMI_group",    labels=c("BMI<25","BMI≥25"),           prefix="bmi"),
    list(var="context_group",labels=c("BMI≥25 Rural","BMI≥25 Urban"), prefix="context")
  )
  relations <- c("BE","BB","EE")
  
  for (sch in schemes) {
    group_var <- sch$var
    if (!group_var %in% names(kaiju_merged)) next
    if (all(is.na(kaiju_merged[[group_var]]))) next
    
    # Tablas por dominio (sin alterar lo anterior)
    tabs <- domain_norm_tables(dat_all, group_var = group_var)
    genus_reads <- tabs$genus_reads
    
    for (lab in sch$labels) {
      # Matriz de counts para este grupo (todas las especies disponibles)
      cnt <- counts_from_long(genus_reads %>% rename(GroupTmp = !!sym(group_var)),
                              "GroupTmp", lab)
      if (is.null(cnt) || nrow(cnt) == 0) {
        message("ℹ Soporte: sin muestras para ", group_var, "=", lab)
        next
      }
      N <- nrow(cnt)
      
      for (rel in relations) {
        tag <- paste0(sch$prefix, "_", label_to_tag(lab, group_var), "_", rel, "_genus_tight")
        f_edges <- paste0("edges_", tag, ".txt")
        if (!file.exists(f_edges)) next
        
        ed <- suppressWarnings(readr::read_tsv(f_edges, show_col_types = FALSE))
        if (!nrow(ed)) next
        
        # Prevalencia de co-ocurrencia (>0 en ambos) por edge
        # (No altera ningún cálculo previo; es informe adicional)
        support <- ed %>%
          rowwise() %>%
          mutate(
            N_samples = N,
            cooc = ifelse(all(c(Source, Target) %in% colnames(cnt)),
                          sum(cnt[, Source] > 0 & cnt[, Target] > 0),
                          NA_integer_),
            prevalence_pct = ifelse(!is.na(cooc), round(100 * cooc / N_samples, 2), NA_real_),
            sign_consistency_pct = round(100 * sign_consistency, 2),
            summary_msg = paste0(
              "En ", lab, " (", group_var, "): '", Source, " ~ ", Target, 
              "' ocurre en ", prevalence_pct, "% de ", N_samples, " muestras; ",
              "seguridad (bootstrap, consistencia de signo) = ", sign_consistency_pct, "%; ",
              "q=", signif(q, 3), ", rho_mean=", signif(rho_mean, 3)
            )
          ) %>%
          ungroup()
        
        # Guardar reporte
        out_file <- paste0("edge_support_", tag, ".txt")
        readr::write_tsv(
          support %>% select(Source, Target, Relation, N_samples, cooc, prevalence_pct,
                             sign_consistency_pct, q, rho_mean, CI_lower, CI_upper, summary_msg),
          out_file
        )
        message("Soporte por edge escrito: ", out_file, 
                " (", nrow(support), " aristas; N=", N, " muestras en ", lab, ")")
        
        # Mensajes sintéticos (una línea por edge puede ser muy verboso;
        # si lo deseas, descomenta la siguiente línea para imprimir todo)
        invisible(lapply(support$summary_msg, message))
      }
    }
  }
}, silent = TRUE)





































### ESTE ES EL FINAL 

## =========================================================
## Redes GENUS (Bacteria/Eukaryota) por:
##   1) Group (Urban vs Rural)
##   2) BMI_group (BMI≥25 vs BMI<25)  -> nombres: BMI<25 -> BMI24 ; BMI≥25 -> BMI26
##   3) context_group (SOLO Adult + BMI≥25):
##        - BMI≥25 Urban  -> BMI26Urban (solo en nombres de archivo)
##        - BMI≥25 Rural  -> BMI26Rural (solo en nombres de archivo)
## Además, para cada esquema se generan 3 redes: BE, BB y EE.
## =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(boot); library(readr)
})

set.seed(1234)

## ======================= PARÁMETROS GLOBALES =======================
# Filtros por grupo (prevalencia y % medio, NORMALIZADO POR DOMINIO)
#Prevalencia mínima (MIN_PREV_). Debe aparecer en al menos ese porcentaje de muestras del grupo
MIN_PREV_BACT <- 0.30 #Base .2
MIN_PREV_EUK  <- 0.30 #Base .1

#porcentaje relativo dentro de cada dominio
MIN_MEAN_BACT <- 0.1    #0.06 base # % dentro de Bacteria
MIN_MEAN_EUK  <- 0.1    #0.015 base % dentro de Eukaryota

#el no. máximo de taxones en caso de que sea muy grande
TOP_N_BACT    <- 70
TOP_N_EUK     <- 50

# Correlación / Bootstrap
BOOT_R               <- 300
PSEUDO               <- 1e-6
# Umbrales de prescreen por relación
PRESCREEN_BE         <- 0.5 #0.35
PRESCREEN_BB         <- 0.5 #0.45
PRESCREEN_EE         <- 0.5 #0.30
PRESCREEN_WL         <- 0.5 #0.20   # si toca whitelist

#minimo de muestras para evaluar 
MIN_SAMPLES_PER_G    <- 4

# Co-ocurrencia mínima (muestras donde AMBOS > 0)
MIN_COOC_FRAC <- 0.50          # que la relacion ocurre en al menos ≥35% de las muestras del grupo
MIN_COOC_ABS  <- 4           # y al menos 20 muestras

# Filtro final por relación (umbral de |rho_mean| y FDR)
THRESH_BE <- 0.35
THRESH_BB <- 0.55
THRESH_EE <- 0.35
ALPHA_Q   <- 0.05              # BH FDR

# Whitelist de géneros clave
force_include_euk  <- c("Saccharomyces")
force_include_bact <- character(0)

## ======================= ENTRADAS =======================
# Se asume que ya cargaste 'kaiju_merged' en el ambiente.
# Si no, descomenta y pon tu ruta:
# kaiju_merged <- readr::read_tsv("tu_kaiju_merged.tsv", show_col_types = FALSE)

# Asegura tipos:
kaiju_merged$reads <- as.numeric(kaiju_merged$reads)

# Si no existen columnas taxonómicas separadas pero sí 'taxon_name'
if (!all(c("Domain","Genus") %in% names(kaiju_merged)) && "taxon_name" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    tidyr::separate(
      taxon_name,
      into = c("Organism","Domain","Supergroup","Kingdom","Phylum","Class",
               "Subclass","Order","Family","Genus","Species"),
      sep = ";", fill = "right", extra = "drop"
    )
}

# file_base si hace falta
if (!"file_base" %in% names(kaiju_merged)) {
  if ("file" %in% names(kaiju_merged)) {
    kaiju_merged <- kaiju_merged %>% mutate(file_base = gsub("^.*/", "", file))
  } else {
    stop("No 'file_base' or 'file' column found in kaiju_merged.")
  }
}

## ============= METADATOS: construir BMI_group para BMI≥25 / BMI<25 =============
# Lee tu archivo de metadatos (ajusta la ruta si cambia)
meta <- read.csv(
  file = "/home/alumno21/axel/files/data_207_3.csv",
  header = TRUE, sep = ",", fileEncoding = "latin1",
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

# Limpieza mínima
meta <- meta[, colSums(!is.na(meta)) > 0]
meta <- meta[rowSums(is.na(meta)) < ncol(meta), ]
meta$BMI <- suppressWarnings(as.numeric(meta$BMI))
meta$Age <- suppressWarnings(as.numeric(meta$Age))

# Agrupaciones de IMC (como lo tenías + colapsado a 2 niveles)
meta <- meta %>%
  mutate(
    BMI_group = case_when(
      Percentil_formulas %in% c("Overweight", "Obesity") ~ "Overweight/Obesity",
      Percentil_formulas == "Normal Weight" ~ "Normal Weight",
      Percentil_formulas == "Malnutrition" ~ "Malnutrition",
      TRUE ~ NA_character_
    ),
    BMI_index_group = case_when(
      BMI_group == "Overweight/Obesity" ~ ">25 BMI INDEX",
      BMI_group %in% c("Normal Weight", "Malnutrition") ~ "<25 BMI INDEX",
      TRUE ~ NA_character_
    )
  )

# En tus metadatos, el identificador de secuenciación parece estar en 'Lane'
# y los archivos de Kaiju acaban en "_kaiju.out". Ajusta si tu sufijo es distinto.
files_over25  <- meta %>% filter(BMI_index_group == ">25 BMI INDEX") %>% pull(Lane) %>% unique()
files_under25 <- meta %>% filter(BMI_index_group == "<25 BMI INDEX") %>% pull(Lane) %>% unique()

# Asignar BMI_group a kaiju_merged (etiquetas finales BMI≥25 / BMI<25)
kaiju_merged <- kaiju_merged %>%
  mutate(BMI_group = case_when(
    file_base %in% paste0(files_over25,  "_kaiju.out") ~ "BMI≥25",
    file_base %in% paste0(files_under25, "_kaiju.out") ~ "BMI<25",
    TRUE ~ NA_character_
  ))

## ============= (Opcional) asignar Group (Urban/Rural) si no existe =============
if (!"Group" %in% names(kaiju_merged)) {
  kaiju_merged <- kaiju_merged %>%
    mutate(Group = case_when(
      exists("urban_files") & file_base %in% paste0(urban_files, "_kaiju.out") ~ "Urban",
      exists("rural_files") & file_base %in% paste0(rural_files, "_kaiju.out") ~ "Rural",
      TRUE ~ NA_character_
    ))
}

## ===== NUEVO: construir context_group = SOLO Adult + BMI≥25 y dividir Urban/Rural =====
# 1) Lista de Lanes "Adult"
adult_lanes <- meta %>%
  filter(!is.na(Age_group), Age_group == "Adult",
         !is.na(Gender), Gender == "Male") %>%
  pull(Lane) %>%
  unique()


KAIJU_SUFFIX <- "_kaiju.out"
adult_files  <- paste0(adult_lanes, KAIJU_SUFFIX)

# 2) Columna context_group en kaiju_merged (solo Adult + BMI≥25; Urban/Rural aparte)
kaiju_merged <- kaiju_merged %>%
  dplyr::mutate(
    is_adult_file = file_base %in% adult_files,
    context_group = dplyr::case_when(
      is_adult_file & BMI_group == "BMI≥25" & Group == "Urban" ~ "BMI≥25 Urban",
      is_adult_file & BMI_group == "BMI≥25" & Group == "Rural" ~ "BMI≥25 Rural",
      TRUE ~ NA_character_
    )
  )

## ======================= Helpers genéricos =======================
# Limpieza base (solo Bacteria/Eukaryota a nivel GENUS; excluye Homo)
prep_dat <- function(df_all) {
  df_all %>%
    filter(Domain %in% c("Bacteria","Eukaryota")) %>%
    mutate(Genus = ifelse(is.na(Genus) | Genus=="" | Genus=="Unclassified", NA, str_trim(Genus))) %>%
    filter(!is.na(Genus)) %>%
    filter(!(Domain=="Eukaryota" & Genus=="Homo"))
}

# Métricas por dominio dentro de cada muestra (% dentro del DOMINIO)
domain_norm_tables <- function(dat, group_var = "Group") {
  stopifnot(all(c("file_base", group_var, "Domain", "Genus", "reads") %in% names(dat)))
  totals_domain <- dat %>%
    group_by(.data[[group_var]], file_base, Domain) %>%
    summarise(total_reads_domain = sum(reads, na.rm = TRUE), .groups = "drop")
  genus_reads <- dat %>%
    group_by(.data[[group_var]], file_base, Domain, Genus) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(totals_domain, by = c(group_var, "file_base", "Domain")) %>%
    mutate(pct_domain = if_else(total_reads_domain > 0, 100 * reads / total_reads_domain, 0))
  list(genus_reads = genus_reads, totals_domain = totals_domain)
}

# Elegir géneros “keepers” por dominio (prevalencia y % medio) usando el group_var
pick_keepers <- function(genus_reads, group_var = "Group") {
  prev_mean <- genus_reads %>%
    group_by(.data[[group_var]], Domain, Genus) %>%
    summarise(prevalence = mean(reads > 0),
              mean_pct_dom = mean(pct_domain),
              .groups = "drop")
  
  # NOTA: asumimos 2 niveles en group_var; se ordenan por aparición
  lvls <- sort(unique(prev_mean[[group_var]]))
  
  keepers_bact <- prev_mean %>%
    filter(Domain=="Bacteria") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_BACT,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_BACT,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_BACT),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_BACT) %>% pull(Genus)
  
  keepers_euk <- prev_mean %>%
    filter(Domain=="Eukaryota") %>%
    group_by(Genus) %>%
    summarise(
      ok = all(prevalence[.data[[group_var]]==lvls[1]] >= MIN_PREV_EUK,
               prevalence[.data[[group_var]]==lvls[2]] >= MIN_PREV_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[1]] >= MIN_MEAN_EUK,
               mean_pct_dom[.data[[group_var]]==lvls[2]] >= MIN_MEAN_EUK),
      overall = mean(mean_pct_dom),
      .groups = "drop"
    ) %>%
    filter(ok) %>% arrange(desc(overall)) %>% slice_head(n = TOP_N_EUK) %>% pull(Genus)
  
  keepers_bact <- union(keepers_bact,
                        intersect(force_include_bact, unique(genus_reads$Genus[genus_reads$Domain=="Bacteria"])))
  keepers_euk  <- union(keepers_euk,
                        intersect(force_include_euk,  unique(genus_reads$Genus[genus_reads$Domain=="Eukaryota"])))
  
  list(keepers_bact = keepers_bact, keepers_euk = keepers_euk)
}

# Matrices de counts por grupo (según group_var) y CLR
counts_from_long <- function(df_long, group_var, group_value) {
  df_long %>%
    filter(.data[[group_var]] == group_value) %>%
    select(file_base, Genus, reads) %>%
    group_by(file_base, Genus) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Genus, values_from = reads, values_fill = 0) %>%
    column_to_rownames("file_base") %>% as.matrix()
}

clr_from_counts <- function(counts_mat) {
  if (is.null(counts_mat) || nrow(counts_mat) < MIN_SAMPLES_PER_G) return(NULL)
  prop <- sweep(counts_mat, 1, rowSums(counts_mat), "/"); prop[is.na(prop)] <- 0
  logx <- log(prop + PSEUDO); sweep(logx, 1, rowMeans(logx), "-")
}

# -------- Correlaciones GENÉRICAS (BE, BB, EE) con bootstrap --------
fast_pairwise_corr <- function(mat_clr, counts_mat, genus_domain_map,
                               relation = c("BE","BB","EE"),
                               prescreen_global = list(BE = PRESCREEN_BE, BB = PRESCREEN_BB, EE = PRESCREEN_EE),
                               prescreen_wl     = PRESCREEN_WL,
                               min_cooc_frac = MIN_COOC_FRAC,
                               min_cooc_abs  = MIN_COOC_ABS,
                               boot_R = BOOT_R,
                               whitelist = character(0)) {
  
  relation <- match.arg(relation)
  cols <- colnames(mat_clr)
  dom_vec <- setNames(as.character(genus_domain_map$Domain),
                      as.character(genus_domain_map$Genus))
  
  bact_cols <- cols[dom_vec[cols] == "Bacteria"]
  euk_cols  <- cols[dom_vec[cols] == "Eukaryota"]
  
  if (relation == "BE" && (length(bact_cols) == 0 || length(euk_cols) == 0)) {
    message("⚠ Prescreen(", relation, "): no hay columnas suficientes Bacteria/Eukaryota para correlacionar.")
    return(tibble())
  }
  if (relation == "BB" && length(bact_cols) < 2) {
    message("⚠ Prescreen(", relation, "): menos de 2 géneros bacterianos tras keepers.")
    return(tibble())
  }
  if (relation == "EE" && length(euk_cols) < 2) {
    message("⚠ Prescreen(", relation, "): menos de 2 géneros eucariotas tras keepers.")
    return(tibble())
  }
  
  # Co-ocurrencia
  pres_abs <- max(min_cooc_abs, ceiling(min_cooc_frac * nrow(counts_mat)))
  pres_mat <- counts_mat > 0
  
  # Rangos para Spearman masivo
  mat_rank <- apply(mat_clr, 2, rank, ties.method = "average")
  
  if (relation == "BE") {
    R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                    mat_rank[, euk_cols,  drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    sel <- which(abs(R) >= prescreen_global$BE, arr.ind = TRUE)
    if (length(sel)==0) message("⚠ Prescreen(BE): ninguna pareja supera |rho|>=", prescreen_global$BE)
    pairs_df <- tibble(Genus_1 = bact_cols[sel[, "row"]],
                       Genus_2 = euk_cols[ sel[, "col"]],
                       from_whitelist = FALSE)
    # Whitelist laxo
    if (length(whitelist) > 0) {
      wl_b <- intersect(bact_cols, whitelist)
      wl_e <- intersect(euk_cols,  whitelist)
      if (length(wl_b) > 0) {
        R_wlb <- R[match(wl_b, bact_cols), , drop = FALSE]
        sel_wlb <- which(abs(R_wlb) >= prescreen_wl, arr.ind = TRUE)
        if (length(sel_wlb)==0) message("ℹ Prescreen(BE) whitelist-bacteria: ninguna pareja supera |rho|>=", prescreen_wl)
        if (length(sel_wlb))
          pairs_df <- bind_rows(pairs_df,
                                tibble(Genus_1 = wl_b[sel_wlb[, "row"]],
                                       Genus_2 = euk_cols[sel_wlb[, "col"]],
                                       from_whitelist = TRUE))
      }
      if (length(wl_e) > 0) {
        R_wle <- R[, match(wl_e, euk_cols), drop = FALSE]
        sel_wle <- which(abs(R_wle) >= prescreen_wl, arr.ind = TRUE)
        if (length(sel_wle)==0) message("ℹ Prescreen(BE) whitelist-euk: ninguna pareja supera |rho|>=", prescreen_wl)
        if (length(sel_wle))
          pairs_df <- bind_rows(pairs_df,
                                tibble(Genus_1 = bact_cols[sel_wle[, "row"]],
                                       Genus_2 = wl_e[sel_wle[, "col"]],
                                       from_whitelist = TRUE))
      }
    }
  } else if (relation == "BB") {
    R <- stats::cor(mat_rank[, bact_cols, drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    upper <- upper.tri(R, diag = FALSE)
    idx <- which(upper & abs(R) >= prescreen_global$BB, arr.ind = TRUE)
    if (length(idx)==0) message("⚠ Prescreen(BB): ninguna pareja supera |rho|>=", prescreen_global$BB)
    pairs_df <- tibble(Genus_1 = bact_cols[idx[, "row"]],
                       Genus_2 = bact_cols[idx[, "col"]],
                       from_whitelist = FALSE)
  } else { # EE
    R <- stats::cor(mat_rank[, euk_cols, drop = FALSE],
                    method = "pearson", use = "pairwise.complete.obs")
    upper <- upper.tri(R, diag = FALSE)
    idx <- which(upper & abs(R) >= prescreen_global$EE, arr.ind = TRUE)
    if (length(idx)==0) message("⚠ Prescreen(EE): ninguna pareja supera |rho|>=", prescreen_global$EE)
    pairs_df <- tibble(Genus_1 = euk_cols[idx[, "row"]],
                       Genus_2 = euk_cols[idx[, "col"]],
                       from_whitelist = FALSE)
  }
  
  if (nrow(pairs_df) == 0) return(tibble())
  
  # Únicos + co-ocurrencia mínima
  pairs_df <- pairs_df %>%
    distinct(Genus_1, Genus_2, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(cooc = sum(pres_mat[, Genus_1] & pres_mat[, Genus_2])) %>%
    ungroup()
  
  if (!nrow(pairs_df)) return(tibble())
  n_before <- nrow(pairs_df)
  pairs_df <- pairs_df %>% filter(cooc >= pres_abs)
  if (nrow(pairs_df) == 0) {
    message("⚠ Filtro de co-ocurrencia (", relation, "): 0 de ", n_before,
            " parejas cumplen co-ocurrencia >= ", pres_abs, 
            " (", MIN_COOC_FRAC*100, "% y mínimo absoluto ", MIN_COOC_ABS, ")")
    return(tibble())
  }
  
  message("Prescreen(", relation, ") pairs (post co-occur): ", nrow(pairs_df),
          " (", sum(pairs_df$from_whitelist), " via whitelist)")
  
  if (nrow(pairs_df) == 0) return(tibble())
  
  boot_one <- function(g1, g2) {
    x <- mat_clr[, g1]; y <- mat_clr[, g2]
    if (sd(x) == 0 || sd(y) == 0) {
      message("ℹ Bootstrap(", relation, "): descarto pareja por varianza cero: ", g1, " ~ ", g2)
      return(NULL)
    }
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    rho0 <- unname(ct$estimate); p0 <- ct$p.value
    stat <- function(dd, idx) suppressWarnings(cor(dd[idx,1], dd[idx,2], method = "spearman"))
    bt <- boot(cbind(x, y), statistic = stat, R = boot_R)
    if (all(is.na(bt$t))) {
      message("ℹ Bootstrap(", relation, "): replicados NA para ", g1, " ~ ", g2)
      return(NULL)
    }
    ci <- tryCatch(boot.ci(bt, type = "perc"), error = function(e) NULL)
    scons <- mean(sign(bt$t) == sign(rho0), na.rm = TRUE)
    rel_label <- if (relation=="BE") "Bacteria-Eukaryota" else if (relation=="BB") "Bacteria-Bacteria" else "Eukaryota-Eukaryota"
    tibble(
      Genus_1 = g1, Genus_2 = g2, Relation = rel_label,
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

# BH helper
add_q <- function(df) if (nrow(df)) df %>% mutate(q = p.adjust(p, method = "BH")) else df

# Filtro final por relación
filter_by_relation <- function(df, relation = c("BE","BB","EE")) {
  relation <- match.arg(relation)
  if (nrow(df)==0) return(df[0,])
  
  thr <- switch(relation,
                "BE" = THRESH_BE,
                "BB" = THRESH_BB,
                "EE" = THRESH_EE)
  
  n_in <- nrow(df)
  res <- df %>%
    mutate(abs_rho = abs(rho_mean)) %>%
    mutate(ci_strong = (rho_mean > 0 & CI_lower >= thr) |
             (rho_mean < 0 & CI_upper <= -thr)) %>%
    filter(
      abs_rho >= thr,
      abs(rho) >= thr,
      ci_strong,
      sign_consistency >= 0.90,
      q <= ALPHA_Q
    ) %>%
    select(-abs_rho, -ci_strong)
  
  if (nrow(res)==0) {
    message("⚠ Filtro final (", relation, "): 0 de ", n_in, 
            " pasan: |rho_mean|>=", thr, ", |rho|>=", thr, 
            ", CI fuerte, sign_consistency>=0.90, q<=", ALPHA_Q)
  }
  res
}

make_nodes <- function(edges_df, genus_domain_map) {
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

# --- Mapea etiquetas a tags seguros solo para los nombres de archivo ---
label_to_tag <- function(lbl, group_var = NULL) {
  if (!is.null(group_var) && group_var == "BMI_group") {
    if (grepl("^\\s*BMI\\s*<\\s*25\\s*$", lbl))  return("BMI24")
    if (grepl("^\\s*BMI\\s*[≥>=]\\s*25\\s*$", lbl)) return("BMI26")
  }
  if (!is.null(group_var) && group_var == "context_group") {
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Urban", lbl)) return("BMI26Urban")
    if (grepl("BMI\\s*[≥>=]\\s*25\\s*Rural", lbl)) return("BMI26Rural")
  }
  out <- gsub("[^A-Za-z0-9]+", "", lbl)
  if (nzchar(out)) out else "grp"
}

relation_tag <- function(relation_label) {
  if (relation_label == "Bacteria-Eukaryota") return("BE")
  if (relation_label == "Bacteria-Bacteria")  return("BB")
  if (relation_label == "Eukaryota-Eukaryota")return("EE")
  "REL"
}

## ======================= FUNCIÓN “PIPELINE” GENERAL =======================
`%||%` <- function(a,b) if (!is.null(a)) a else b

run_pipeline_for <- function(df, group_var = "Group", labels = NULL, tag_prefix = NULL) {
  # Asegurar grupos válidos
  df <- df %>% filter(!is.na(.data[[group_var]]))
  if (length(unique(df[[group_var]])) != 2) {
    message("⚠ Saltando ", group_var, ": no hay exactamente 2 niveles válidos.")
    return(invisible(NULL))
  }
  
  if (is.null(labels)) {
    labels <- sort(unique(df[[group_var]]))
  }
  gA <- labels[1]; gB <- labels[2]
  
  message(">>> Ejecutando para ", group_var, " = {", gA, ", ", gB, "} …")
  
  dat <- prep_dat(df)
  if (!nrow(dat)) { message("⚠ prep_dat() dejó 0 filas para ", group_var); return(invisible(NULL)) }
  
  # % dentro de dominio
  tabs <- domain_norm_tables(dat, group_var = group_var)
  genus_reads <- tabs$genus_reads
  if (!nrow(genus_reads)) { message("⚠ domain_norm_tables() devolvió 0 filas (", group_var, ")"); return(invisible(NULL)) }
  
  # keepers
  k <- pick_keepers(genus_reads, group_var = group_var)
  keepers_bact <- k$keepers_bact; keepers_euk <- k$keepers_euk
  message("Kept genera — Bacteria=", length(keepers_bact), " | Eukaryota=", length(keepers_euk))
  if (length(keepers_bact)==0 || length(keepers_euk)==0) {
    message("⚠ Filtro keepers: requiere al menos 1 género por dominio. Bact=", length(keepers_bact), ", Euk=", length(keepers_euk))
  }
  
  # Filtrado a keepers
  dat_filt <- genus_reads %>%
    filter((Domain=="Bacteria"  & Genus %in% keepers_bact) |
             (Domain=="Eukaryota" & Genus %in% keepers_euk))
  if (!nrow(dat_filt)) { message("⚠ dat_filt vacío tras keepers para ", group_var); return(invisible(NULL)) }
  
  # Matrices por grupo
  cnt_A <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gA)
  cnt_B <- counts_from_long(dat_filt %>% rename(GroupTmp = !!sym(group_var)), "GroupTmp", gB)
  if (is.null(cnt_A) || nrow(cnt_A) < MIN_SAMPLES_PER_G) message("⚠ ", gA, ": no pasa filtro de cantidad de muestras (n=", ifelse(is.null(cnt_A), 0, nrow(cnt_A)), " < ", MIN_SAMPLES_PER_G, ")")
  if (is.null(cnt_B) || nrow(cnt_B) < MIN_SAMPLES_PER_G) message("⚠ ", gB, ": no pasa filtro de cantidad de muestras (n=", ifelse(is.null(cnt_B), 0, nrow(cnt_B)), " < ", MIN_SAMPLES_PER_G, ")")
  stopifnot(!is.null(cnt_A), !is.null(cnt_B), nrow(cnt_A) >= MIN_SAMPLES_PER_G, nrow(cnt_B) >= MIN_SAMPLES_PER_G)
  
  mat_A <- clr_from_counts(cnt_A)
  mat_B <- clr_from_counts(cnt_B)
  if (is.null(mat_A)) message("⚠ ", gA, ": clr_from_counts() devolvió NULL (posible n<mínimo).")
  if (is.null(mat_B)) message("⚠ ", gB, ": clr_from_counts() devolvió NULL (posible n<mínimo).")
  
  # Alinear columnas
  common_cols <- Reduce(intersect, list(colnames(mat_A), colnames(mat_B)))
  if (length(common_cols) == 0) { message("⚠ Sin columnas comunes de géneros entre ", gA, " y ", gB); return(invisible(NULL)) }
  mat_A <- mat_A[, common_cols, drop = FALSE]
  mat_B <- mat_B[, common_cols, drop = FALSE]
  cnt_A <- cnt_A[, common_cols, drop = FALSE]
  cnt_B <- cnt_B[, common_cols, drop = FALSE]
  
  # Genus -> Domain map
  genus_domain_map <- dat_filt %>%
    distinct(Genus, Domain) %>%
    filter(Genus %in% common_cols)
  if (!nrow(genus_domain_map)) { message("⚠ genus_domain_map vacío tras intersección de columnas."); return(invisible(NULL)) }
  
  # --------- Correlaciones y filtros por RELACIÓN: BE, BB, EE ----------
  relations <- c("BE","BB","EE")
  for (rel in relations) {
    message("Computing ", rel, " correlations (", gA, ")…")
    edges_A <- fast_pairwise_corr(
      mat_A, cnt_A, genus_domain_map,
      relation = rel,
      whitelist = union(force_include_euk, force_include_bact),
      boot_R = BOOT_R
    )
    message("Computing ", rel, " correlations (", gB, ")…")
    edges_B <- fast_pairwise_corr(
      mat_B, cnt_B, genus_domain_map,
      relation = rel,
      whitelist = union(force_include_euk, force_include_bact),
      boot_R = BOOT_R
    )
    
    edges_A <- add_q(edges_A); edges_B <- add_q(edges_B)
    
    edges_A_f <- filter_by_relation(edges_A, relation = rel)
    edges_B_f <- filter_by_relation(edges_B, relation = rel)
    
    message(rel, " edges kept | ", gA, ": ", nrow(edges_A_f), " | ", gB, ": ", nrow(edges_B_f))
    
    nodes_A <- make_nodes(edges_A_f, genus_domain_map)
    nodes_B <- make_nodes(edges_B_f, genus_domain_map)
    
    # Tags de salida
    prefix <- (tag_prefix %||% tolower(group_var))
    rel_tag <- rel
    tagA <- paste0(prefix, "_", label_to_tag(gA, group_var), "_", rel_tag, "_genus_tight")
    tagB <- paste0(prefix, "_", label_to_tag(gB, group_var), "_", rel_tag, "_genus_tight")
    
    write_net(edges_A_f, nodes_A, tagA)
    write_net(edges_B_f, nodes_B, tagB)
  }
  
  invisible(NULL)
}

## ======================= CORRER PIPELINE PARA LOS 3 CASOS =======================
# 1) Urban vs Rural
if ("Group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$Group)) && dplyr::n_distinct(na.omit(kaiju_merged$Group)) == 2) {
  run_pipeline_for(
    kaiju_merged %>% filter(!is.na(Group)),
    group_var = "Group",
    labels = c("Rural","Urban"),         # fuerza el orden (opcional)
    tag_prefix = "group"
  )
} else {
  message("Saltando Group: no hay dos niveles en 'Group'.")
}

# 2) BMI≥25 vs BMI<25
if ("BMI_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$BMI_group))) {
  run_pipeline_for(
    kaiju_merged %>% filter(!is.na(BMI_group)),
    group_var = "BMI_group",
    labels = c("BMI<25","BMI≥25"),       # orden consistente
    tag_prefix = "bmi"
  )
} else {
  message("No hay 'BMI_group' asignado a kaiju_merged (revisa mapeo Lane->file_base).")
}

# 3) context_group: SOLO Adult + BMI≥25, comparando Urban vs Rural
if ("context_group" %in% names(kaiju_merged) && any(!is.na(kaiju_merged$context_group))) {
  lvls_ctx <- sort(unique(na.omit(kaiju_merged$context_group)))
  if (length(lvls_ctx) == 2) {
    run_pipeline_for(
      kaiju_merged %>% filter(!is.na(context_group)),
      group_var = "context_group",
      labels = c("BMI≥25 Rural","BMI≥25 Urban"),  # orden explícito
      tag_prefix = "context"
    )
  } else {
    message("Saltando context_group: niveles detectados = ", paste(lvls_ctx, collapse=", "))
  }
} else {
  message("No hay 'context_group' (Adult + BMI≥25) disponible.")
}

message("Listo. Redes BE/BB/EE exportadas para Group, BMI_group y context_group (si disponibles).")

## ======================= SECCIÓN ADICIONAL: SOPORTE POR EDGE (CON LISTA DE FILES) =======================
## No modifica filtros ni funciones anteriores. Solo lee los 'edges_<tag>.txt' ya generados y crea reportes extra.

try({
  dat_all <- prep_dat(kaiju_merged)
  
  # Mapa de info por muestra (file_base) para enriquecer el detalle
  sample_info <- kaiju_merged %>%
    distinct(file_base, Group, BMI_group, context_group)
  
  # Añadir Age desde 'meta' usando Lane -> file_base
  sample_info <- sample_info %>%
    mutate(Lane = sub("_kaiju\\.out$", "", file_base)) %>%
    left_join(meta %>% select(Lane, Age), by = "Lane") %>%
    select(-Lane)
  
  schemes <- list(
    list(var="Group",        labels=c("Rural","Urban"),               prefix="group"),
    list(var="BMI_group",    labels=c("BMI<25","BMI≥25"),             prefix="bmi"),
    list(var="context_group",labels=c("BMI≥25 Rural","BMI≥25 Urban"), prefix="context")
  )
  relations <- c("BE","BB","EE")
  
  for (sch in schemes) {
    group_var <- sch$var
    if (!group_var %in% names(kaiju_merged)) next
    if (all(is.na(kaiju_merged[[group_var]]))) next
    
    # Tablas por dominio (para recomputar presencia >0 sin tocar cálculos)
    tabs <- domain_norm_tables(dat_all, group_var = group_var)
    genus_reads <- tabs$genus_reads
    
    for (lab in sch$labels) {
      # Matriz de counts para este grupo (todas las especies disponibles)
      cnt <- counts_from_long(genus_reads %>% rename(GroupTmp = !!sym(group_var)),
                              "GroupTmp", lab)
      if (is.null(cnt) || nrow(cnt) == 0) {
        message("ℹ Soporte: sin muestras para ", group_var, "=", lab)
        next
      }
      N <- nrow(cnt)
      files_vec <- rownames(cnt)
      
      for (rel in relations) {
        tag <- paste0(sch$prefix, "_", label_to_tag(lab, group_var), "_", rel, "_genus_tight")
        f_edges <- paste0("edges_", tag, ".txt")
        if (!file.exists(f_edges)) next
        
        ed <- suppressWarnings(readr::read_tsv(f_edges, show_col_types = FALSE))
        if (!nrow(ed)) next
        
        # Columnas esperadas en edges: Source, Target, Relation, rho_mean, CI_lower, CI_upper, p, q, sign_consistency
        needed_cols <- c("Source","Target","Relation","rho_mean","CI_lower","CI_upper","p","q","sign_consistency")
        if (!all(needed_cols %in% names(ed))) {
          message("⚠ Formato inesperado en ", f_edges, " (faltan columnas). Me salto este archivo.")
          next
        }
        
        # --- Resumen por edge ---
        support <- ed %>%
          rowwise() %>%
          mutate(
            N_samples = N,
            cooc = ifelse(all(c(Source, Target) %in% colnames(cnt)),
                          sum(cnt[, Source] > 0 & cnt[, Target] > 0),
                          NA_integer_),
            prevalence_pct = ifelse(!is.na(cooc), round(100 * cooc / N_samples, 2), NA_real_),
            sign_consistency_pct = round(100 * sign_consistency, 2), # en %
            summary_msg = paste0(
              "En ", lab, " (", group_var, "): '", Source, " ~ ", Target,
              "' ocurre en ", prevalence_pct, "% de ", N_samples, " muestras; ",
              "seguridad (bootstrap, consistencia de signo) = ", sign_consistency_pct, "%; ",
              "q=", signif(q, 3), ", rho_mean=", signif(rho_mean, 3)
            )
          ) %>%
          ungroup()
        
        # Guardar resumen
        out_file <- paste0("edge_support_", tag, ".txt")
        readr::write_tsv(
          support %>% select(Source, Target, Relation, N_samples, cooc, prevalence_pct,
                             sign_consistency_pct, q, rho_mean, CI_lower, CI_upper, summary_msg),
          out_file
        )
        
        # --- Detalle por muestra (lista de files que soportan cada edge) ---
        
        
        # Para cada edge, extraer file_base donde ambos géneros > 0
        details_list <- lapply(seq_len(nrow(ed)), function(i) {
          g1 <- ed$Source[i]; g2 <- ed$Target[i]
          if (!all(c(g1,g2) %in% colnames(cnt))) return(NULL)
          
          mask <- (cnt[, g1] > 0) & (cnt[, g2] > 0)
          if (!any(mask)) return(NULL)
          
          files_hit <- files_vec[mask]
          tibble(
            Source = g1,
            Target = g2,
            file_base = files_hit
          )
        })
        details <- bind_rows(details_list)
        
        if (!is.null(details) && nrow(details)) {
          details <- details %>%
            left_join(sample_info, by = "file_base") %>%
            mutate(
              Group_info         = Group,
              Age_info           = Age,
              BMI_group_info     = BMI_group,
              context_group_info = context_group,
              group_var          = group_var,
              group_label        = lab,
              tag                = tag
            ) %>%
            select(Source, Target, file_base,
                   Group_info, Age_info, BMI_group_info, context_group_info,
                   group_var, group_label, tag)
          
          out_details <- paste0("edge_support_details_", tag, ".txt")
          readr::write_tsv(details, out_details)
        }
        
        message("Soporte por edge escrito: ", out_file,
                if (exists("out_details")) paste0(" + ", out_details) else "",
                " (", nrow(support), " aristas; N=", N, " muestras en ", lab, ")")
      }
    }
  }
}, silent = TRUE)

