#Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(reshape2)

# # Specify your log file and output file
log_file <- commandArgs(trailingOnly = TRUE)[1]
primary_tissue <- commandArgs(trailingOnly = TRUE)[2]

# log_file <- "beast_gundem_2015_2_21_24/A29_asym/A29_unlabeled_tree_final_input_xml.log"
# primary_tissue <- "prostate"

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
melted_df <- melt(log_df)
melted_df <- separate(melted_df, variable, into = c("Source", "Recipient"), sep = "_")
grouped_df <- melted_df %>%
  # filter(value > 0) %>%
  group_by(Source, Recipient) %>%
  summarise(mean_rate = mean(value, na.rm = TRUE))


order_source <- c(grouped_df$Source)
order_recipient <- c(grouped_df$Recipient)

order <- unique(c(order_source, order_recipient))

order <- unique(order[order != primary_tissue])
order_numeric <- as.numeric(sub("M", "", order))
if (!any(is.na(order_numeric))) {
order <- paste("M", sort(order_numeric), sep = "")
}

order <- c(primary_tissue, order)

add_rows <- setdiff(order, order_source)
add_cols <- setdiff(order, order_recipient)

heatmap_df <- pivot_wider(grouped_df, names_from = Recipient, values_from = mean_rate)


if (length(add_rows) != 0) {
# Add rows for source tissues missing
heatmap_df <- rbind(heatmap_df, setNames(data.frame(matrix(NA, ncol = ncol(heatmap_df), nrow = length(add_rows))), colnames(heatmap_df)))
rownames(heatmap_df)[(nrow(heatmap_df) - length(add_rows) + 1):nrow(heatmap_df)] <- add_rows
}

if (length(add_cols) != 0) {
# Add cols for recipient tissues missing
heatmap_df[,add_cols] <- NA
}

### NEED TO FIX THIS AND DECIDE HOW TO PLOT THE RATES
subst <- 0
num_states <- nrow(heatmap_df)
state_freq <- rep(1 / num_states, num_states)
for (i in 1:nrow(heatmap_df)) {
  row_sum <- sum(heatmap_df[i, ], na.rm = TRUE) - heatmap_df[i, i]
  heatmap_df[i, i] <- -row_sum
  subst <- subst + (row_sum * state_freq[i])
}

# Normalize all rate matrix entries by dividing them by the negative sum of the off-diagonal elements in the row
for (i in 1:nrow(heatmap_df)) {
  row_sum <- sum(heatmap_df[i, ], na.rm = TRUE) - heatmap_df[i, i]
  heatmap_df[i, ] <- heatmap_df[i, ] / subst
}

heatmap_df <- rownames_to_column(heatmap_df, var = "Recipient")
heatmap_df <- column_to_rownames(heatmap_df, var = "Source")

heatmap_df <- heatmap_df %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

# Increase color scale when rates go above 1 to the nearest 0.5 above the highest rate
max_limit=1
max_value = max(na.omit(heatmap_df$value))
if ( max_value > 1) {
  max_limit <- ceiling(max_value * 2) / 2
}

# optional to set a permanent scale bar limit
max_limit <- 1

# Create a ggplot2 heatmap
heatmap <- ggplot(heatmap_df, aes(x = factor(colname, levels = order), y = factor(rowname, levels = order), fill=value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, max_limit)) +
  theme_minimal() +
  labs(x="Recipient tissue", y="Source tissue", fill="Rate") +
  theme(axis.text.x=element_text(size=18, color="black", angle = 90, hjust = 1),
        axis.text.y=element_text(size=18, color="black"),
        axis.title=element_text(size=20, color="black"),
        legend.text=element_text(size=18, color="black"),
        legend.title=element_text(size=20, color="black"),
        panel.grid = element_blank())

output_file <- sub("\\.log$", "_log.pdf", log_file)
ggsave(output_file, heatmap)
