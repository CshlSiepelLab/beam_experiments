
library(EvoTraceR)
args = commandArgs(trailingOnly=TRUE)

input_dir <- args[1]
output_dir <- args[2]
asv_threshold <- as.numeric(args[3])
barcode_version <- args[4]
cores <- as.numeric(args[5])

# Set parameters based on barcode version
if (barcode_version == "BC10v2") {
  ref_name <- "BC10v2"
  ref_seq <- "TCTACACGCGCGTTCAACCGAGGAAAACTACACACGCATTCAACCACGGTTTATTACACGCACATTCAACCGTGGACTGCTACACACGCGCTCAACCACGGATATTTACGCACACGTTCAACCGCGGATTGTTACACCCGCATTCAACCGAGGTCACCTACACCCGCACTCAACCGGGGTACGCGACACGTGCGATCAACCGAGGCTTACTACCCGCACGTTCAACTGGGGAACACTGCACGCGAGTTCGACCGGGGATC"
  ref_flank_right <- "ACCGGGGATC"
} else if (barcode_version == "BC10v0") {
  ref_name <- "BC10v0"
  ref_seq <- "TCTACACGCGCGTTCAACCGAGGAAAACTACACACACGTTCAACCACGGTTTTTTACACACGCATTCAACCACGGACTGCTACACACGCACTCAACCGTGGATATTTACATACTCGTTCAACCGTGGATTGTTACACCCGCGTTCAACCAGGGTCAGATACACCCACGTTCAACCGTGGTACTATACTCGGGCATTCAACCGCGGCTTTCTGCACACGCCTACAACCGCGGAACTATACACGTGCATTCACCCGTGGATC"
  ref_flank_right <- "CCCGTGGATC"
} else {
  stop("Invalid barcode version. Must be either 'BC10v0' or 'BC10v2'")
}

# Common parameters for both versions
ref_flank_left <- "TCTAC"
ref_cut_sites <- c(17, 43, 69, 95, 121, 147, 173, 199, 225, 251)
ref_border_sites <- c(1, 26, 52, 78, 104, 130, 156, 182, 208, 234)

# unzip files
# List all files in the directory
files <- list.files(input_dir)
# Filter files with ".zip" extension
zip_files <- files[grepl("\\.zip$", files)]
# Check if there are any zip files
if (length(zip_files) > 0) {
  # Loop through each zip file and unzip its contents
  for (zip_file in zip_files) {
    # Specify the full path to the zip file
    zip_file_path <- file.path(input_dir, zip_file)
    # Unzip the file
    unzip(zip_file_path, exdir = input_dir)
    cat("Unzipped:", zip_file, "\n")
  }
} else {
  cat("No zip files found in the directory.\n")
}

trimmomatic_path <- Sys.getenv("TRIMMOMATIC_PATH")
flash_path <- Sys.which("flash")

EvoTraceR_object <-
  initialize_EvoTraceR(
    input_dir = input_dir,
    output_dir = output_dir,
    trimmomatic_path = trimmomatic_path,
    flash_path = flash_path)

EvoTraceR_object <-
  asv_analysis(EvoTraceR_object = EvoTraceR_object,
               ref_name = ref_name,
               ref_seq = ref_seq,
               ref_flank_left = paste0("^", ref_flank_left),
               ref_flank_right = paste0(ref_flank_right, "$"),
               ref_cut_sites = ref_cut_sites,
               ref_border_sites = ref_border_sites,
               flanking_filtering = "right",
               output_figures = TRUE,
               asv_count_cutoff = asv_threshold, # minimum number of ASVs to be counted; decided as 3 on: 03/25/22
               # pair-wise alignment parameters between un-edited barcode and edited barcode (ASV)
               pwa_type = "global", # based on AmpliCan (global = Needleman-Wunsch)
               pwa_gapOpening = -25, # based on AmpliCan: -25
               pwa_gapExtension = 0, # based on AmpliCan: 0
               pwa_match = 15, # based on AmpliCan: 15
               pwa_mismatch = -4, # based on AmpliCan: -4
               cleaning_window = c(10, 10), # cleaning window +/- from Cas9 editing size (nucleotide 17 in guide) is considered as an edit 
               batch_size = 100,
               cores = cores               
               )

EvoTraceR_object <-
  analyse_mutations(EvoTraceR_object = EvoTraceR_object)

EvoTraceR_object <-
  infer_phylogeny(EvoTraceR_object = EvoTraceR_object, mutations_use = "del_ins")

EvoTraceR_object <-
  create_df_summary(EvoTraceR_object)

save.image(paste0(output_dir, "/", "evotracer.RData"))
