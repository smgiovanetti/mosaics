#setwd("/vf/users/SBGE/simone/fastqs_111725")
#in fastqs111725, files have been processed by 1bcountvalidbc.py and 2classify_final.py
#output of these steps are in annotated_tables_barcode_new
#each file is a sample (wt or mut for 403), replicate (1,2), bin (1,2,3,4), concentration (7,9)
#also reading in at same time the read counts for each fastq so can compute ratio 
library(tidyverse)

input_dir <- "annotated_tables_barcode_new"
files <- list.files(input_dir,
                    pattern = "_frequency\\.txt$",
                    full.names = TRUE)

variant_cols_raw <- c(
  "439","440","441","443","445","449","459",
  "478","483","484","486","490","493","494",
  "498","501","505"
)

read_counts <- read.delim("5-7-26readcount.txt",
                          header = TRUE,
                          stringsAsFactors = FALSE)
read_one_file <- function(path) {
  fn <- basename(path)
  
  sample_info <- read_counts %>%
    filter(name == fn)
  
  if (nrow(sample_info) == 0) {
    stop(paste("No metadata found for", fn))
  }
  
  read_count <- sample_info$read_count
  cell <- sample_info$cell
  
  #concentration 0, no bin
  m0 <- str_match(
    fn,
    "^([mw])(1|2)-0_S\\d+_R1_001\\.fastq\\.gz_frequency\\.txt$"
  )
  
  #concentration 7 or 9, with bin
  m1 <- str_match(
    fn,
    "^([mw])(1|2)-(\\d)-([79])_S\\d+_R1_001\\.fastq\\.gz_frequency\\.txt$"
  )
  
  if (!all(is.na(m0[1, ]))) {
    sample_type <- m0[1, 2]
    replicate <- as.integer(m0[1, 3])
    bin <- NA_integer_
    concentration <- 0
  } else if (!all(is.na(m1[1, ]))) {
    sample_type <- m1[1, 2]
    replicate <- as.integer(m1[1, 3])
    bin <- as.integer(m1[1, 4])
    concentration <- as.integer(m1[1, 5])
  } else {
    stop("Filename did not match expected pattern: ", fn)
  }
  
  # m = S at 403, w = R at 403
  first_genotype_pos <- ifelse(sample_type == "m", "S", "R")
  
  df <- read.delim(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    encoding = "latin1"
  )
  
  freq_col <- grep("_freq$", names(df), value = TRUE)
  if (length(freq_col) != 1) {
    stop("Expected exactly one *_freq column in: ", fn)
  }
  
  variant_cols <- if (all(variant_cols_raw %in% names(df))) {
    variant_cols_raw
  } else if (all(paste0("X", variant_cols_raw) %in% names(df))) {
    paste0("X", variant_cols_raw)
  } else {
    stop("Missing variant columns in: ", fn)
  }
  
  df$barcode <- toupper(trimws(as.character(df$barcode)))
  
  # genotype = first position (403) from filename + positions 439-505 from the table
  df$genotype <- apply(df[, variant_cols, drop = FALSE], 1, paste, collapse = "_")
  df$genotype <- paste(first_genotype_pos, df$genotype, sep = "_")
  
  df$frequency <- as.numeric(df[[freq_col]])
  df$replicate <- replicate
  df$bin <- bin
  df$concentration <- concentration
  df$read_count <- read_count
  df$cell <- cell
  df %>%
    transmute(
      genotype,
      barcode,
      replicate,
      bin,
      concentration,
      cell,
      read_count,
      frequency
    )}

combined_df <- map_dfr(files, read_one_file) %>%
  group_by(genotype, barcode, replicate, bin, concentration, cell, read_count) %>%
  summarise(frequency = sum(frequency), .groups = "drop")

###combined_df has columns: genotype barcode replicate bin concentration frequency
write.table(
  combined_df,
  "combined_genotype_barcode_table.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#now reading in whitelist dfs and then will filter for genotypes
setwd("/vf/users/SBGE/simone")
WT_whitelist <- read.delim("WT_whitelist.txt",
                           header = TRUE,
                           stringsAsFactors = FALSE)

MUT_whitelist <- read.delim("MUT_whitelist.txt",
                            header = TRUE,
                            stringsAsFactors = FALSE)

WT_whitelist$pattern <- paste0("R_", WT_whitelist$pattern)
MUT_whitelist$pattern <- paste0("S_", MUT_whitelist$pattern)

combined_whitelist <- bind_rows(WT_whitelist, MUT_whitelist)

###filtering 
combined_df_whitelisted <- combined_df %>%
  semi_join(combined_whitelist, by = c("barcode", "genotype" = "pattern"))

##calcualte ratio (frequency/read_count)
combined_df_whitelisted <- combined_df_whitelisted %>%
  mutate(
    read_count = as.numeric(gsub(",", "", read_count)),
    cell = as.numeric(gsub(",", "", cell)),
    ratio = frequency / read_count,
    ratio_cell = cell * (frequency / read_count)
  )




#now pivot wider. matching rows with the same barcode genotype replicate concentration 
#but first need to separate bin0 (since not tied to concentration)

# separate rows with no bin
df_na_bin <- combined_df_whitelisted %>%
  filter(is.na(bin) | tolower(as.character(bin)) == "N/A") %>%
  mutate(bin = 0)

# keep only bin 1-4 rows
df_bins <- combined_df_whitelisted %>%
  filter(!is.na(bin), tolower(as.character(bin)) != "n/a", bin %in% 1:4)

# reshape to wide, keep ratio and frequency(for filtering)
df_final <- df_bins %>%
  select(genotype, barcode, replicate, concentration, bin, ratio, ratio_cell,frequency) %>%
  mutate(bin = paste0("bin", bin)) %>%
  pivot_wider(
    names_from = bin,
    values_from = c(ratio, ratio_cell, frequency),
    names_glue = "{bin}_{.value}",
    values_fill = 0
  )

# make sure the columns are in the order you want
df_final <- df_final %>%
  select(
    genotype,
    barcode,
    replicate,
    concentration,
    bin1_frequency,
    bin2_frequency,
    bin3_frequency,
    bin4_frequency,
    bin1_ratio,
    bin2_ratio,
    bin3_ratio,
    bin4_ratio,
    bin1_ratio_cell,
    bin2_ratio_cell,
    bin3_ratio_cell,
    bin4_ratio_cell
  )

## now adding bin0 values
df_na_bin <- df_na_bin %>%
  group_by(genotype, barcode, replicate) %>%
  summarise(
    bin0_ratio = first(ratio),
    bin0_ratio_cell = first(ratio_cell),
    bin0_frequency = sum(frequency, na.rm = TRUE),
    .groups = "drop"
  )

df_final_combined <- df_final %>%
  left_join(df_na_bin,
            by = c("genotype", "barcode", "replicate")) %>%
  mutate(
    bin0_ratio = replace_na(bin0_ratio, 0),
    bin0_ratio_cell = replace_na(bin0_ratio_cell, 0),
    bin0_frequency = replace_na(bin0_frequency, 0)
  )


###saving this
write.table(df_final_combined,
            "7-7-26collapsed_ratio.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

df_final_combined <- read.delim(
  "5-7-26collapsed_ratio.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

df_final_combined <- read.delim(
  "7-7-26collapsed_ratio.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

###calculate mean_bin (weighted)
df_final_combined <- df_final_combined %>%
  mutate(
    mean_bin = (
      bin1_ratio * 4 +
        bin2_ratio * 3 +
        bin3_ratio * 2 +
        bin4_ratio * 1
    ) / (
      bin1_ratio +
        bin2_ratio +
        bin3_ratio +
        bin4_ratio
    ),
    
    mean_bin_cell = (
      bin1_ratio_cell * 4 +
        bin2_ratio_cell * 3 +
        bin3_ratio_cell * 2 +
        bin4_ratio_cell * 1
    ) / (
      bin1_ratio_cell +
        bin2_ratio_cell +
        bin3_ratio_cell +
        bin4_ratio_cell
    )
  )

positions <- c(
  "403","439","440","441","443","445","449","459",
  "478","483","484","486","490","493","494","498","501","505"
)

pos_cols  <- c(
  "403","439","440","441","443","445","449","459",
  "478","483","484","486","490","493","494","498","501","505"
)

df_final_combined <- df_final_combined %>%
  separate(genotype,
           into = positions,
           sep = "_",
           remove = FALSE)

##add column to ask if 'window' (contiguous S) or 'mixed'
is_window <- function(x) {
  s_idx <- which(x == "S")
  
  # no S or all S counts as a window
  if (length(s_idx) <= 1) return(TRUE)
  
  # contiguous if every gap is exactly 1
  all(diff(s_idx) == 1)
}

df_final_combined <- df_final_combined %>%
  rowwise() %>%
  mutate(
    window = if_else(
      is_window(c_across(all_of(pos_cols))),
      "window",
      "mixed"
    )
  ) %>%
  ungroup()

#18 VIOLINS, RESIDUE EFFECT (for window chimeras, add #window filter)
long_positions <- df_final_combined %>%
  pivot_longer(
    cols = all_of(positions),
    names_to = "position",
    values_to = "pos_value"
  ) %>%
  filter(!is.na(mean_bin))
long_positions %>%
  filter(
    concentration == 7, #window == "window",
    (bin1_frequency +
       bin2_frequency +
       bin3_frequency +
       bin4_frequency) > 10,
    (bin1_frequency > 0) +
      (bin2_frequency > 0) +
      (bin3_frequency > 0) +
      (bin4_frequency > 0) >= 2
  )  %>%
  ggplot(aes(x = pos_value, y = mean_bin_cell)) +
  geom_violin() +
  stat_summary(fun = mean, geom = "point", color = "#d7df23", size = 1.5) +
  labs(
    x = "background",
    y = "mean bin",
    title = "-7 ace2, pacbio whitelist, bin1-4 > 10, observed in >1 bin, cell ratio"
  ) +
  facet_wrap(~position, scales = "free_x", ncol = 6) +
  theme_minimal() +
  theme(legend.position = "none") 
  
ggsave("9-8-26allconc7_pbwhitelist_1234_10_2bins.pdf", width = 6, height = 6)


#R/S DISTRIBUTION
position_props <- df_final_combined %>%
  filter(
    concentration == 7,
    (bin1_frequency +
       bin2_frequency +
       bin3_frequency +
       bin4_frequency) > 10, (bin1_frequency > 0) + (bin2_frequency > 0) + (bin3_frequency > 0) + (bin4_frequency > 0) >= 2
  ) %>%
  select(genotype, all_of(pos_cols)) %>%
  distinct() %>%
  pivot_longer(
    cols = all_of(pos_cols),
    names_to = "position",
    values_to = "allele"
  ) %>%
  count(position, allele) %>%
  group_by(position) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
ggplot(position_props,
       aes(x = factor(position, levels = pos_cols),
           y = proportion,
           fill = allele)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "R" = "#d7df23",
      "S" = "#2b23df"
    )
  ) +
  theme_minimal() +
  labs(
    x = "position",
    y = "proportion",
    fill = "genotype",
    title = "-7 ace2, whitelist, bin1-4>10, observed in >1 bin, R/S distribution"
  )
ggsave("9-8-26proportionposition.pdf", height = 6, width = 6)


##FREQUENCY OF PAIRS
#deviation of pairs (25:25:25:25)
pos_cols_17  <- c(
  "439","440","441","443","445","449","459",
  "478","483","484","486","490","493","494","498","501","505"
)
geno_unique <- df_final_combined %>%
  filter(concentration == 7,  (bin1_frequency +
                                 bin2_frequency +
                                 bin3_frequency +
                                 bin4_frequency) > 10,
         (bin1_frequency > 0) +
           (bin2_frequency > 0) +
           (bin3_frequency > 0) +
           (bin4_frequency > 0) >= 2
  ) %>%
  select(genotype, all_of(pos_cols_17)) %>%
  distinct()

pairs <- combn(pos_cols_17, 2, simplify = FALSE)

pair_freq <- map_dfr(pairs, function(x) {
  
  p1 <- x[1]
  p2 <- x[2]
  
  tab <- geno_unique %>%
    count(
      allele1 = .data[[p1]],
      allele2 = .data[[p2]]
    ) %>%
    complete(
      allele1 = c("R", "S"),
      allele2 = c("R", "S"),
      fill = list(n = 0)
    ) %>%
    mutate(
      prop = n / sum(n),
      combo = paste0(allele1, allele2)
    )
  
  tab %>%
    select(combo, prop) %>%
    mutate(
      position1 = p1,
      position2 = p2
    )
})


pair_plot <- pair_freq %>%
  mutate(
    position1 = factor(position1, levels = pos_cols_17),
    position2 = factor(position2, levels = rev(pos_cols_17)),
    combo = factor(combo, levels = c("RR", "RS", "SR", "SS"))
  )

ggplot(
  pair_plot,
  aes(
    x = position1,
    y = position2,
    fill = prop
  )
) +
  geom_tile(color = "white") +
  #geom_text(
  #aes(label = sprintf("%d", round(prop * 100))),
  #color = "black",
  #size = 3
  #) +
  facet_wrap(~combo, nrow = 2) +
  scale_fill_gradient2(
    midpoint = 0.25,
    low = "#878ec6",
    mid = "white",
    high = "#d7df23",
    limits = c(0, 0.5),
    name = "frequency"
  ) +
  theme_minimal() +
  labs(
    x = "position",
    y = "position",
    title = "-7 ace2, whitelist, bins1-4 > 10, observed in > 1 bin, white = 25"
  )
ggsave("9-28-26variantheatmap.pdf", height = 8, width = 8)

#deviation of pairs, calculate expected based on single variant distribution
single_freq <- geno_unique %>%
  pivot_longer(
    cols = all_of(pos_cols_17),
    names_to = "position",
    values_to = "allele"
  ) %>%
  count(position, allele) %>%
  group_by(position) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()

# observed vs expected for each pair and each combination
pair_combo_dev <- map_dfr(pairs, function(x) {
  
  p1 <- x[1]
  p2 <- x[2]
  
  obs <- geno_unique %>%
    count(
      allele1 = .data[[p1]],
      allele2 = .data[[p2]]
    ) %>%
    complete(
      allele1 = c("R", "S"),
      allele2 = c("R", "S"),
      fill = list(n = 0)
    ) %>%
    mutate(
      observed = n / sum(n),
      combo = paste0(allele1, allele2)
    ) %>%
    select(allele1, allele2, combo, observed)
  
  p1_freq <- single_freq %>%
    filter(position == p1) %>%
    select(allele, prop) %>%
    rename(allele1 = allele, prop1 = prop)
  
  p2_freq <- single_freq %>%
    filter(position == p2) %>%
    select(allele, prop) %>%
    rename(allele2 = allele, prop2 = prop)
  
  expected <- expand_grid(
    allele1 = c("R", "S"),
    allele2 = c("R", "S")
  ) %>%
    left_join(p1_freq, by = "allele1") %>%
    left_join(p2_freq, by = "allele2") %>%
    mutate(expected = prop1 * prop2) %>%
    select(allele1, allele2, expected)
  
  obs %>%
    left_join(expected, by = c("allele1", "allele2")) %>%
    mutate(
      deviation = observed - expected,
      position1 = p1,
      position2 = p2
    )
})

pair_combo_dev <- pair_combo_dev %>%
  mutate(
    position1 = factor(position1, levels = pos_cols_17),
    position2 = factor(position2, levels = rev(pos_cols_17)),
    combo = factor(combo, levels = c("RR", "RS", "SR", "SS"))
  )

ggplot(pair_combo_dev, aes(x = position1, y = position2, fill = deviation)) +
  geom_tile(color = "white") +
  facet_wrap(~combo, nrow = 2) +
  #geom_text(aes(label = sprintf("%.1f", deviation)), size = 2) +
  scale_fill_gradient2(
    midpoint = 0,
    low = "#878ec6",
    mid = "white",
    high = "#d7df23",
    limits = c(-0.25, 0.25),
    name = "obs - exp"
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    x = "position",
    y = "position",
    title = "deviation from expected (observed - expected)"
  )
ggsave("9-28-26deviationfromexpecteddist.pdf", height = 8, width = 8)

pair_combo_dev %>%
  summarise(
    median_abs_dev = median(abs(deviation)),
    mean_abs_dev = mean(abs(deviation)),
    max_abs_dev = max(abs(deviation))
  )

##VARIANT DISTRIBUTION
#all R and all S mean mean_bin
df_final_combined %>%
  filter(
    concentration == 7,
    (bin1_frequency +
       bin2_frequency +
       bin3_frequency +
       bin4_frequency) > 10,
    (bin1_frequency > 0) +
      (bin2_frequency > 0) +
      (bin3_frequency > 0) +
      (bin4_frequency > 0) >= 2
  ) %>%
  mutate(
    all_R = if_all(c(`403`,`439`,`440`,`441`,`443`,`445`,`449`,`459`,
                     `478`,`483`,`484`,`486`,`490`,`493`,`494`,
                     `498`,`501`,`505`),
                   ~ . == "R"),
    all_S = if_all(c(`403`,`439`,`440`,`441`,`443`,`445`,`449`,`459`,
                     `478`,`483`,`484`,`486`,`490`,`493`,`494`,
                     `498`,`501`,`505`),
                   ~ . == "S")
  ) %>%
  filter(all_R | all_S) %>%
  mutate(category = ifelse(all_R, "all R", "all S")) %>%
  group_by(category) %>%
  summarise(
    mean_mean_bin = mean(mean_bin_cell),
    n = n(),
    .groups = "drop"
  )

#all R = 1.13, all S = 3.03

df_final_combined %>%
  filter(concentration == 7,
         (bin1_frequency +
            bin2_frequency +
            bin3_frequency +
            bin4_frequency) > 10,
         (bin1_frequency > 0) +
           (bin2_frequency > 0) +
           (bin3_frequency > 0) +
           (bin4_frequency > 0) >= 2
  ) %>%
  ggplot(aes(x = "", y = mean_bin_cell)) +
  #geom_histogram(binwidth = 0.05) +
  #geom_vline(xintercept = 1.13, color = "#d7df23", linewidth = .5) +
  #geom_vline(xintercept = 3.03, color = "#d7df23", linewidth = .5) +
  geom_violin()+
  geom_hline(
    yintercept = 1.13,
    color = "#d7df23",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 3.03,
    color = "#d7df23",
    linewidth = 0.5
  ) +
  labs(
    y = "mean bin",
    title = "-7 ace2, bins1-4>10, observed in >1 bin") +
  theme_minimal()
ggsave("9-8-26vdis_conc7_pbwhitelist_1-4_10_cell.pdf", width = 6, height = 6)
ggsave("9-6-26vdis_conc7_pbwhitelist_1-4_10_cell_VIOLIN.pdf", width = 6, height = 6)



##EFFECT OF NUMBER OF S ON MEAN BIN
df_final_combined <- df_final_combined %>%
  rowwise() %>%
  mutate(num_S = sum(c_across(all_of(positions)) == "S", na.rm = TRUE)) %>%
  ungroup()

df_final_combined %>%
  filter(concentration == 7, 
         (bin1_frequency +
            bin2_frequency +
            bin3_frequency +
            bin4_frequency) > 10,
         (bin1_frequency > 0) +
           (bin2_frequency > 0) +
           (bin3_frequency > 0) +
           (bin4_frequency > 0) >= 2
  ) %>%
  ggplot(aes(x = factor(num_S), y = mean_bin_cell)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = "number of S",
    y = "mean_bin",
    title = "concentration 7, sum bin1-4 >10, observed in >1 bin"
  )
ggsave("9-8-26numS_conc7_pbwhitelist_1-4_10_cell.pdf", width = 6, height = 6)

library(patchwork)

#MODELS
library(broom)
df_model <- df_final_combined %>%
  filter(concentration == 7, 
         (bin1_frequency +
            bin2_frequency +
            bin3_frequency +
            bin4_frequency) > 10,
         (bin1_frequency > 0) +
           (bin2_frequency > 0) +
           (bin3_frequency > 0) +
           (bin4_frequency > 0) >= 2
         ) %>%
  select(mean_bin_cell, all_of(positions)) %>%
  rename_with(~ paste0("pos_", .x), all_of(positions)) %>%
  mutate(across(starts_with("pos_"), as.factor))

#main effect
m1_18 <- lm(
  mean_bin_cell ~ pos_403 + pos_439 + pos_440 + pos_441 + pos_443 +
    pos_445 + pos_449 + pos_459 + pos_478 + pos_483 +
    pos_484 + pos_486 + pos_490 + pos_493 + pos_494 +
    pos_498 + pos_501 + pos_505,
  data = df_model
)

significant_terms <- tidy(m1_18) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    p_adj = p.adjust(p.value, method = "fdr")
  ) %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)

lm_terms <- broom::tidy(
  m1_18,
  conf.int = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    position = stringr::str_extract(term, "\\d+"),
    p_adj = p.adjust(p.value, method = "fdr")
  )

ggplot(
  lm_terms,
  aes(
    x = estimate,
    y = position
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    width = 0.15
  ) +
  geom_point(size = 2) +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  labs(
    x = " effect of S",
    y = "position"
  )
ggsave("9-17-26positioneffect.pdf", width = 6, height = 6)


#lm on window chimeras
df_model_window <- df_final_combined %>%
  filter(
    concentration == 7,
    window == "window",
    (bin1_frequency +
       bin2_frequency +
       bin3_frequency +
       bin4_frequency) > 10,
    (
      (bin1_frequency > 0) +
        (bin2_frequency > 0) +
        (bin3_frequency > 0) +
        (bin4_frequency > 0)
    ) >= 2
  ) %>%
  select(mean_bin_cell, all_of(positions)) %>%
  rename_with(
    ~ paste0("pos_", .x),
    all_of(positions)
  ) %>%
  mutate(
    across(starts_with("pos_"), as.factor)
  )
m1_18_window <- lm(
  mean_bin_cell ~ pos_403 + pos_439 + pos_440 + pos_441 + pos_443 +
    pos_445 + pos_449 + pos_459 + pos_478 + pos_483 +
    pos_484 + pos_486 + pos_490 + pos_493 + pos_494 +
    pos_498 + pos_501 + pos_505,
  data = df_model_window
)
lm_terms_window <- broom::tidy(
  m1_18_window,
  conf.int = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    position = stringr::str_extract(term, "\\d+"),
    p_adj = p.adjust(p.value, method = "fdr")
  ) %>%
  arrange(p_adj)

lm_terms_window
#end window chimeras


m1_6 <- lm(
  mean_bin_cell ~ pos_449 + pos_486 + pos_493 + pos_498 + pos_501 + pos_505,
  data = df_model
)

m1_7 <- lm(
  mean_bin_cell ~ pos_403 + pos_449 + pos_486 + pos_493 + pos_498 + pos_501 + pos_505,
  data = df_model
)

m1_8 <- lm(
  mean_bin_cell ~ pos_403 + pos_439 + pos_449 + pos_486 + pos_493 + pos_498 + pos_501 + pos_505,
  data = df_model
)

# pairwise interactions
m2_18 <- lm(
  mean_bin_cell ~ (pos_403 + pos_439 + pos_440 + pos_441 + pos_443 +
                     pos_445 + pos_449 + pos_459 + pos_478 + pos_483 +
                     pos_484 + pos_486 + pos_490 + pos_493 + pos_494 +
                     pos_498 + pos_501 + pos_505)^2,
  data = df_model
)

m2_6 <- lm(
  mean_bin_cell ~ (pos_449 + pos_486 + pos_493 + pos_498 + pos_501 + pos_505)^2,
  data = df_model
)

m2_7 <- lm(
  mean_bin_cell ~ (pos_403 + pos_449 + pos_486 + pos_493 + pos_498 + pos_501 + pos_505)^2,
  data = df_model
)

m2_8 <- lm(
  mean_bin_cell ~ (pos_403 + pos_439 + pos_449 + pos_486 + pos_493 + pos_498 + pos_501 + pos_505)^2,
  data = df_model
)


# three-way interactions
m3_18 <- lm(
  mean_bin_cell ~ (pos_403 + pos_439 + pos_440 + pos_441 + pos_443 +
                     pos_445 + pos_449 + pos_459 + pos_478 + pos_483 +
                     pos_484 + pos_486 + pos_490 + pos_493 + pos_494 +
                     pos_498 + pos_501 + pos_505)^3,
  data = df_model
)

m3_6 <- lm(
  mean_bin_cell ~ (pos_449 + pos_486 + pos_493 + pos_498 + pos_501 + pos_505)^3,
  data = df_model
)

m3_7 <- lm(
  mean_bin_cell ~ (pos_403 + pos_449 + pos_486 + pos_493 + pos_498 + pos_501 + pos_505)^3,
  data = df_model
)

m3_8 <- lm(
  mean_bin_cell ~ (pos_403 + pos_439 + pos_449 + pos_486 + pos_493 + pos_498 + pos_501 + pos_505)^3,
  data = df_model
)

mods <- list(
  m1_6, m1_7, m1_8, m1_18,
  m2_6, m2_7, m2_8, m2_18,
  m3_6, m3_7, m3_8, m3_18
)

data.frame(
  model = c(
    "6 main", "7 main", "8 main", "18 main",
    "6 pairwise", "7 pairwise", "8 pairwise", "18 pairwise",
    "6 3-way", "7 3-way", "8 3-way", "18 3-way"
  ),
  adj_r2 = sapply(mods, function(x) summary(x)$adj.r.squared)
)


plot_df <- data.frame(
  residues = factor(c(6, 7, 8, 18), levels = c(6, 7, 8, 18)),
  main    = c(0.4109732, 0.5377056, 0.5564632, 0.5614255),
  pairwise = c(0.4352176, 0.5674043, 0.5931325, 0.6038801),
  three_way = c(0.4680235, 0.6123675, 0.6461569, 0.6639046)
) |>
  tidyr::pivot_longer(
    cols = c(main, pairwise, three_way),
    names_to = "model",
    values_to = "adj_r2"
  ) |>
  mutate(
    model = factor(
      model,
      levels = c("main", "pairwise", "three_way"),
      labels = c("main effects", "pairwise interactions", "3 way interactions")
    )
  )

ggplot(plot_df, aes(x = model, y = adj_r2, fill = residues)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c(
    "6"  = "#2b23df",
    "7"  = "#8993FF",
    "8"  = "#74a9b6",
    "18" = "#d7df23"
  )) +
  theme_minimal() +
  labs(
    x = NULL,
    y = "adjusted r2",
    fill = NULL
  ) +
  ylim(0, 1)
ggsave("9-8-26modelcomparison.pdf", width = 6, height = 6)

#pairwise p adj
pairwise_results <- broom::tidy(
  m2_18
) %>%
  filter(
    grepl(":", term)
  ) %>%
  mutate(
    p_adj = p.adjust(
      p.value,
      method = "BH"
    ),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

pairwise_results %>%
  filter(
    grepl("pos_493", term) &
      grepl("pos_483|pos_484|pos_494", term)
  )

pairwise_results %>%
  filter(
    grepl("pos_498", term) &
      grepl("pos_483|pos_484|pos_494", term)
  )



###looking at 1 position vs all
pos_cols <- c(
  "403","439","440","441","443","445","449","459",
  "478","483","484","486","490","493","494","498","501","505"
)

anchor_pos <- "498"
other_positions <- setdiff(pos_cols, anchor_pos)

df_long_interactions <- df_final_combined %>%
  filter(concentration == 7, 
         (bin1_frequency +
            bin2_frequency +
            bin3_frequency +
            bin4_frequency) > 10,
         (bin1_frequency > 0) +
           (bin2_frequency > 0) +
           (bin3_frequency > 0) +
           (bin4_frequency > 0) >= 2
  ) %>%
  pivot_longer(
    cols = all_of(other_positions),
    names_to = "other_position",
    values_to = "other_value"
  ) %>%
  mutate(
    interaction = paste(.data[[anchor_pos]], other_value, sep = "_")
  )

ggplot(df_long_interactions, aes(x = interaction, y = mean_bin_cell)) +
  geom_violin(alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "#d7df23", size = 1.5) +
  facet_wrap(~other_position, scales = "free_x", ncol = 6) +
  theme_minimal() +
  labs(
    x = paste(anchor_pos, "vs position"),
    y = "mean_bin",
    title = paste(anchor_pos, "vs all positions, bins1-4 > 10, observed in >1")
  ) +
  theme(legend.position = "none")
ggsave("9-8-26post498vsall_7_ratio.pdf", width = 12, height = 6)

#484, 493, 494
df_3pos_493 <- df_final_combined %>%
  filter(
    concentration == 7,
    (bin1_frequency +
       bin2_frequency +
       bin3_frequency +
       bin4_frequency) > 10,
    (bin1_frequency > 0) +
      (bin2_frequency > 0) +
      (bin3_frequency > 0) +
      (bin4_frequency > 0) >= 2
  ) %>%
  mutate(
    genotype_3 = paste0(
      `484`,
      `493`,
      `494`
    )
  )

order_df_3_493 <- df_3pos_493 %>%
  group_by(genotype_3) %>%
  summarise(
    mean_bin_cell = mean(mean_bin_cell, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_bin_cell)

df_3pos_493 %>%
  mutate(
    genotype_3 = factor(
      genotype_3,
      levels = order_df_3_493$genotype_3
    )
  ) %>%
  ggplot(aes(x = genotype_3, y = mean_bin_cell)) +
  geom_violin(alpha = 0.9) +
  #geom_jitter(width = 0.15, alpha = 0.15, size = 0.5) +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "#d7df23",
    size = 2
  ) +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "483 484 493",
    y = "mean_bin_cell"
  )

genotype_levels_493 <- order_df_3_493$genotype_3

genotype_grid_493 <- order_df_3_493 %>%
  select(genotype_3) %>%
  mutate(
    genotype_3 = factor(
      genotype_3,
      levels = genotype_levels_493
    ),
    `484` = substr(as.character(genotype_3), 1, 1),
    `493` = substr(as.character(genotype_3), 2, 2),
    `494` = substr(as.character(genotype_3), 3, 3)
  ) %>%
  pivot_longer(
    cols = c(`484`, `493`, `494`),
    names_to = "position",
    values_to = "allele"
  ) %>%
  mutate(
    position = factor(
      position,
      levels = c("484", "493", "494")
    )
  )

p_violin_493 <- df_3pos_493 %>%
  mutate(
    genotype_3 = factor(
      genotype_3,
      levels = genotype_levels_493
    )
  ) %>%
  ggplot(
    aes(
      x = genotype_3,
      y = mean_bin_cell
    )
  ) +
  geom_violin(alpha = 0.9) +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "#d7df23",
    size = 2
  ) +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "484 493 494",
    y = "mean_bin_cell"
  )

p_grid_493 <- ggplot(
  genotype_grid_493,
  aes(
    x = position,
    y = genotype_3,
    fill = allele
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "R" = "#d7df23",
      "S" = "#8993FF"
    )
  ) +
  theme_minimal() +
  labs(
    x = "position",
    y = NULL,
    fill = "allele"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p_violin_493 + p_grid_493 +
  plot_layout(
    widths = c(3, 1)
  )

ggsave("10-8-26_484_493_494violins_grid.pdf", width = 6, height = 6)

#484, 494, 498
df_4pos_498 <- df_final_combined %>%
  filter(
    concentration == 7,
    (bin1_frequency +
       bin2_frequency +
       bin3_frequency +
       bin4_frequency) > 10,
    (bin1_frequency > 0) +
      (bin2_frequency > 0) +
      (bin3_frequency > 0) +
      (bin4_frequency > 0) >= 2
  ) %>%
  mutate(
    genotype_4 = paste0(
      `484`,
      `494`,
      `498`
    )
  )

order_df_4_498 <- df_4pos_498 %>%
  group_by(genotype_4) %>%
  summarise(
    mean_bin_cell = mean(mean_bin_cell, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_bin_cell)

df_4pos_498 %>%
  mutate(
    genotype_4 = factor(
      genotype_4,
      levels = order_df_4_498$genotype_4
    )
  ) %>%
  ggplot(aes(x = genotype_4, y = mean_bin_cell)) +
  geom_violin(alpha = 0.9) +
  #geom_jitter(width = 0.15, alpha = 0.15, size = 0.5) +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "#d7df23",
    size = 2
  ) +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "483 484 494 498",
    y = "mean_bin_cell"
  )
ggsave("9-15-26_483_484_494_498violins.pdf", width = 6, height = 6)



genotype_levels_498 <- order_df_4_498$genotype_4

genotype_grid_498 <- order_df_4_498 %>%
  select(genotype_4) %>%
  mutate(
    genotype_4 = factor(
      genotype_4,
      levels = genotype_levels_498
    ),
    `484` = substr(as.character(genotype_4), 1, 1),
    `494` = substr(as.character(genotype_4), 2, 2),
    `498` = substr(as.character(genotype_4), 3, 3)
  ) %>%
  pivot_longer(
    cols = c(`484`, `494`, `498`),
    names_to = "position",
    values_to = "allele"
  ) %>%
  mutate(
    position = factor(
      position,
      levels = c("484", "494", "498")
    )
  )

p_violin_498 <- df_4pos_498 %>%
  mutate(
    genotype_4 = factor(
      genotype_4,
      levels = genotype_levels_498
    )
  ) %>%
  ggplot(
    aes(
      x = genotype_4,
      y = mean_bin_cell
    )
  ) +
  geom_violin(alpha = 0.9) +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "#d7df23",
    size = 2
  ) +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "484 494 498",
    y = "mean_bin_cell"
  )

p_grid_498 <- ggplot(
  genotype_grid_498,
  aes(
    x = position,
    y = genotype_4,
    fill = allele
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "R" = "#d7df23",
      "S" = "#8993FF"
    )
  ) +
  theme_minimal() +
  labs(
    x = "position",
    y = NULL,
    fill = "allele"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_violin_498 + p_grid_498 +
  plot_layout(
    widths = c(3, 1)
  )

ggsave("10-8-26_484_494_498violins_grid.pdf", width = 6, height = 6)

#downsampling attempt

full_effects <- map_dfr(positions, function(pos) {
  
  means <- df_model %>%
    group_by(state = .data[[paste0("pos_", pos)]]) %>%
    summarise(
      mean_pheno = mean(mean_bin_cell, na.rm = TRUE),
      .groups = "drop"
    )
  
  tibble(
    position = pos,
    full_effect = means$mean_pheno[means$state == "S"] -
      means$mean_pheno[means$state == "R"]
  )
})

sample_sizes <- c(
  100, 250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000
)

set.seed(123)

n_reps <- 1000

estimate_effects <- function(dat) {
  
  map_dfr(positions, function(pos) {
    
    means <- dat %>%
      group_by(state = .data[[paste0("pos_", pos)]]) %>%
      summarise(
        mean_pheno = mean(mean_bin_cell, na.rm = TRUE),
        .groups = "drop"
      )
    
    tibble(
      position = pos,
      effect = means$mean_pheno[means$state == "S"] -
        means$mean_pheno[means$state == "R"]
    )
  })
}

sampling_results <- map_dfr(
  sample_sizes,
  function(n) {
    
    map_dfr(
      seq_len(n_reps),
      function(rep) {
        
        sampled <- df_model %>%
          slice_sample(
            n = min(n, nrow(df_model))
          )
        
        sampled_effects <- estimate_effects(sampled)
        
        sampled_effects %>%
          left_join(
            full_effects,
            by = "position"
          ) %>%
          summarise(
            rmsd = sqrt(
              mean(
                (effect - full_effect)^2,
                na.rm = TRUE
              )
            )
          ) %>%
          mutate(
            n_genotypes = n,
            replicate = rep
          )
      }
    )
  }
)

sampling_summary <- sampling_results %>%
  group_by(n_genotypes) %>%
  summarise(
    median_rmsd = median(rmsd),
    q25 = quantile(rmsd, 0.25),
    q75 = quantile(rmsd, 0.75),
    .groups = "drop"
  )

ggplot(
  sampling_summary,
  aes(x = n_genotypes, y = median_rmsd)
) +
  geom_ribbon(
    aes(ymin = q25, ymax = q75),
    alpha = 0.2
  ) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  theme_minimal() +
  labs(
    x = "number of genotypes",
    y = "rmsd from library main effects"
  )
ggsave('9-30-26downsample_genotypes_rmsd_100reps.pdf', width = 6, height = 6)

sample_sizes <- c(
  10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500
)

set.seed(123)

n_reps <- 1000


downsample_effects <- map_dfr(
  sample_sizes,
  function(n) {
    
    map_dfr(
      seq_len(n_reps),
      function(rep) {
        
        sampled <- df_model %>%
          slice_sample(
            n = min(n, nrow(df_model))
          )
        
        estimate_effects(sampled) %>%
          mutate(
            n_genotypes = n,
            replicate = rep
          )
      }
    )
  }
)

positions8 <- c(
  "403", "439", "449", "486",
  "493", "498", "501", "505"
)

direction_results <- downsample_effects %>%
  filter(position %in% positions8) %>%
  left_join(
    full_effects %>%
      select(position, full_effect),
    by = "position"
  ) %>%
  mutate(
    same_direction =
      sign(effect) == sign(full_effect)
  )

direction_by_site <- direction_results %>%
  group_by(position, n_genotypes) %>%
  summarise(
    fraction_correct_direction =
      mean(same_direction, na.rm = TRUE),
    .groups = "drop"
  )
residue_cols <- c(
  "403" = "#440154",
  "439" = "#692A99",
  "449" = "#A3319F",
  "486" = "#CA3C93",
  "493" = "#F2637F",
  "498" = "#F78C79",
  "501" = "#F6A97A",
  "505" = "#EDD9A3"
)

ggplot(
  direction_by_site,
  aes(
    x = n_genotypes,
    y = fraction_correct_direction,
    group = position,
    color = position
  )
) +
  geom_line() +
  geom_point() +
  scale_color_manual(
    values = residue_cols
  ) +
  scale_x_continuous(
    limits = c(0, 500)
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  theme_minimal() +
  labs(
    x = "number of genotypes",
    y = "fraction with correct effect direction",
    color = "position"
  )
ggsave("9-30-26downsample_genotype.pdf", width = 6, height = 6)



#downsampling 'interaction' detection: 493
df_493_context <- df_model %>%
  mutate(
    genotype_3 = paste0(
      pos_484,
      pos_493,
      pos_494
    )
  )

test_493_interaction <- function(data, n_sample) {
  
  # 1. Randomly sample from the ENTIRE library
  sampled <- data %>%
    slice_sample(
      n = min(n_sample, nrow(data)),
      replace = FALSE
    )
  
  # 2. Only now retain the four genotype classes
  focal <- sampled %>%
    filter(
      genotype_3 %in% c(
        "RRR",
        "RRS",
        "SSR",
        "SSS"
      )
    ) %>%
    mutate(
      background = case_when(
        genotype_3 %in% c("RRR", "RRS") ~ "RR",
        genotype_3 %in% c("SSR", "SSS") ~ "SS"
      ),
      allele_493 = case_when(
        genotype_3 %in% c("RRR", "SSR") ~ "R",
        genotype_3 %in% c("RRS", "SSS") ~ "S"
      ),
      background = factor(
        background,
        levels = c("RR", "SS")
      ),
      allele_493 = factor(
        allele_493,
        levels = c("R", "S")
      )
    )
  
  # Need observations in all four classes
  counts <- focal %>%
    count(genotype_3)
  
  if (!all(
    c("RRR", "RRS", "SSR", "SSS") %in%
    counts$genotype_3
  )) {
    return(
      tibble(
        n_sample = n_sample,
        n_focal = nrow(focal),
        interaction_estimate = NA_real_,
        interaction_p = NA_real_
      )
    )
  }
  
  # 3. Pairwise model
  model <- lm(
    mean_bin_cell ~ background * allele_493,
    data = focal
  )
  
  # 4. Pull out interaction term
  interaction <- broom::tidy(model) %>%
    filter(term == "backgroundSS:allele_493S")
  
  tibble(
    n_sample = n_sample,
    n_focal = nrow(focal),
    interaction_estimate = interaction$estimate,
    interaction_p = interaction$p.value
  )
}

sample_sizes_493 <- c(
  100, 250, 500,
  1000, 2000, 3000, 4000, 5000,
  7500, 10000,
  15000, 20000, 30000,
  40000, 50000
)

set.seed(123)

n_reps <- 1000

interaction_downsample <- map_dfr(
  sample_sizes_493,
  function(n) {
    
    map_dfr(
      seq_len(n_reps),
      function(rep) {
        
        test_493_interaction(
          df_493_context,
          n_sample = n
        ) %>%
          mutate(replicate = rep)
      }
    )
  }
)

interaction_downsample <- interaction_downsample %>%
  mutate(
    significant = interaction_p < 0.05
  )

interaction_power <- interaction_downsample %>%
  group_by(n_sample) %>%
  summarise(
    n_valid = sum(!is.na(interaction_p)),
    fraction_valid = mean(!is.na(interaction_p)),
    fraction_significant = mean(
      significant,
      na.rm = TRUE
    ),
    median_n_focal = median(n_focal),
    .groups = "drop"
  )

ggplot(
  interaction_power,
  aes(
    x = n_sample,
    y = fraction_significant
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  theme_minimal() +
  labs(
    x = "number of genotypes",
    y = "fraction with 493 × (483/484)"
  )
ggsave("10-8-26downsample_493x.pdf", width = 6, height = 6)

#downsampling 'interaction' detection: 498
df_498_context <- df_model %>%
  mutate(
    genotype_3 = paste0(
      pos_484,
      pos_494,
      pos_498
    )
  )

test_498_interaction <- function(data, n_sample) {
  
  # 1. Randomly sample from the ENTIRE library
  sampled <- data %>%
    slice_sample(
      n = min(n_sample, nrow(data)),
      replace = FALSE
    )
  
  # 2. Only now retain the four genotype classes
  focal <- sampled %>%
    filter(
      genotype_3 %in% c(
        "RRR",
        "RRS",
        "SSR",
        "SSS"
      )
    ) %>%
    mutate(
      background = case_when(
        genotype_3 %in% c("RRR", "RRS") ~ "RR",
        genotype_3 %in% c("SSR", "SSS") ~ "SS"
      ),
      allele_498 = case_when(
        genotype_3 %in% c("RRR", "SSR") ~ "R",
        genotype_3 %in% c("RRS", "SSS") ~ "S"
      ),
      background = factor(
        background,
        levels = c("RR", "SS")
      ),
      allele_498 = factor(
        allele_498,
        levels = c("R", "S")
      )
    )
  
  # Need observations in all four classes
  counts <- focal %>%
    count(genotype_3)
  
  if (!all(
    c("RRR", "RRS", "SSR", "SSS") %in%
    counts$genotype_3
  )) {
    return(
      tibble(
        n_sample = n_sample,
        n_focal = nrow(focal),
        interaction_estimate = NA_real_,
        interaction_p = NA_real_
      )
    )
  }
  
  # 3. Interaction model
  model <- lm(
    mean_bin_cell ~ background * allele_498,
    data = focal
  )
  
  # 4. Pull out interaction term
  interaction <- broom::tidy(model) %>%
    filter(term == "backgroundSS:allele_498S")
  
  tibble(
    n_sample = n_sample,
    n_focal = nrow(focal),
    interaction_estimate = interaction$estimate,
    interaction_p = interaction$p.value
  )
}

sample_sizes_498 <- c(
  100, 250, 500,
  1000, 2000, 3000, 4000, 5000,
  7500, 10000,
  15000, 20000, 30000,
  40000, 50000
)

set.seed(123)

n_reps <- 1000

interaction_downsample <- map_dfr(
  sample_sizes_498,
  function(n) {
    
    map_dfr(
      seq_len(n_reps),
      function(rep) {
        
        test_498_interaction(
          df_498_context,
          n_sample = n
        ) %>%
          mutate(replicate = rep)
      }
    )
  }
)

interaction_downsample <- interaction_downsample %>%
  mutate(
    significant = interaction_p < 0.05
  )

interaction_power <- interaction_downsample %>%
  group_by(n_sample) %>%
  summarise(
    n_valid = sum(!is.na(interaction_p)),
    fraction_valid = mean(!is.na(interaction_p)),
    fraction_significant = mean(
      significant,
      na.rm = TRUE
    ),
    median_n_focal = median(n_focal),
    .groups = "drop"
  )

ggplot(
  interaction_power,
  aes(
    x = n_sample,
    y = fraction_significant
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  theme_minimal() +
  labs(
    x = "number of genotypes",
    y = "fraction with 498 × (483/494)"
  )
ggsave("10-8-26downsample_498x.pdf", width = 6, height = 6)


