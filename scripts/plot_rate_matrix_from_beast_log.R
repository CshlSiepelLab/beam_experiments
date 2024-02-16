#Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# Specify your log file and output file
log_file <- commandArgs(trailingOnly = TRUE)[1]
primary_tissue <- commandArgs(trailingOnly = TRUE)[2]

# log_file <- "/Users/staklins/projects/crispr-barcode-cancer-metastasis/stephen_data/longer_10million_mcmc_unsymmetrical_machina_m8_sims_compare_beast_machina_fixedtreeanalysis_default_2_12_24/T_seed3_unlabeled_true_tree_final_input_xml.log"
# primary_tissue <- "P"

# burnin <- 0.1   ### Seems like Tracer values are closer to averages with burnin kept in, so I turned this off for now.

# Read in log data
log_df <- read.table(log_file, header=TRUE, comment.char = "#", sep="\t")

# # Discard burnin samples
# total_rows <- nrow(log_df)
# burnin <- round(total_rows * burnin)
# log_df <- read.table(log_file, header=TRUE, comment.char = "#", sep="\t", skip=burnin)

# Subset columns to only those relevant to the rates
log_df <- log_df[grep("geoSubstModelLogger.relGeoRate_", names(log_df))]

# Rename columns to only have source and recipient tissue labels
names(log_df) <- gsub("geoSubstModelLogger.relGeoRate_", "", names(log_df))

# Melt the DataFrame to long format
melted_df <- reshape2::melt(log_df)
melted_df <- separate(melted_df, variable, into = c("Source", "Recipient"), sep = "_")
grouped_df <- melted_df %>%
  # filter(value > 0) %>%
  group_by(Source, Recipient) %>%
  summarise(mean_rate = mean(value, na.rm = TRUE))


order_source <- c(grouped_df$Source)
order_recipient <- c(grouped_df$Recipient)

order <- unique(c(order_source, order_recipient))
order <- c(primary_tissue, sort(unique(order[order != primary_tissue])))

add_rows <- setdiff(order, order_source)
add_cols <- setdiff(order, order_recipient)

heatmap_df <- pivot_wider(grouped_df, names_from = Recipient, values_from = mean_rate)
heatmap_df <- column_to_rownames(heatmap_df, var = "Source")

# Add rows for source tissues missing
heatmap_df <- rbind(heatmap_df, setNames(data.frame(matrix(NA, ncol = ncol(heatmap_df), nrow = length(add_rows))), colnames(heatmap_df)))
rownames(heatmap_df)[(nrow(heatmap_df) - length(add_rows) + 1):nrow(heatmap_df)] <- add_rows

# Add cols for recipient tissues missing
heatmap_df[,add_cols] <- NA

heatmap_df <- heatmap_df %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

# Create a ggplot2 heatmap
heatmap <- ggplot(heatmap_df, aes(x = factor(colname, levels = order), y = factor(rowname, levels = order), fill=value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
  theme_minimal() +
  labs(x="Recipient tissue", y="Source tissue", fill="Rate") +
  theme(axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=18, color="black"),
        axis.title=element_text(size=20, color="black"),
        legend.text=element_text(size=18, color="black"),
        legend.title=element_text(size=20, color="black"),
        panel.grid = element_blank())

output_file <- sub("\\.log$", "_log.pdf", log_file)
ggsave(output_file, heatmap)
