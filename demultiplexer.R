#' Searches a fastq file for reads containing specified barcodes
#'
#' @author Joshua Jackson <jjack29_uwo.ca>
#' @examples

extract_reads_with_barcodes <- function(path_to_input_fastq_r1, path_to_input_fastq_r2, path_to_barcodes_txt, path_to_output_fastq_r1, path_to_output_fastq_r2, length_UMI = 0, char_trail = "T", length_trail = 0, provide_updates = TRUE) {
  # Loads package with enhanced version of data.frame structure (speeding up data manipulation)
  library(data.table)
  
  # Reads the barcode input file
  barcodes <- readLines(path_to_barcodes_txt)
  
  # Prepends periods to ignore the UMI in the read (number of periods prepended is specified by the function's length_UMI parameter) 
  mod_barcodes <- paste0(strrep(".", length_UMI), barcodes)
  
  # Appends a trailing character if one occurs after the UMI and barcode (character type and amount appended is specified by char_trail and length_trail parameters)
  mod_barcodes <- paste0(mod_barcodes, strrep(char_trail, length_trail))
  
  # Prepends a caret to ensure only the beginning of the read is being searched
  mod_barcodes <- paste0("^", mod_barcodes)
  
  # Reads the input fastq file (putting it in the data structure loaded from the data.table library)
  input_fastq_r1 <- fread(path_to_input_fastq_r1, sep="\n", header=FALSE)
  input_fastq_r2 <- fread(path_to_input_fastq_r2, sep="\n", header=FALSE)
  
  # Creates a list storing all sequences containing specified barcodes (updating progress along the way)
  seq_containing_barcodes <- list()
  for (i in 1:length(barcodes)) {
    bc <- mod_barcodes[i]
    matched_rows <- grep(bc, input_fastq_r1$V1)
    seq_containing_barcodes[[i]] <- matched_rows
    if (provide_updates == TRUE) {
      cat("Barcode:", barcodes[i], "- Matched sequences:", length(matched_rows), "\n")
      cat(paste0(round(i/length(barcodes)*100, 1), "% finished"), "\n")
    }
  }
  
  # Converts list to vector (to allow for further data manipulation)
  vec_seq_containing_barcodes <- unlist(seq_containing_barcodes)
  
  # Creates a vector to store the full 4-line fastq version of each sequence
  fastq_containing_barcodes <- c()
  for (i in 0:3) {
    fastq_containing_barcodes <- c(fastq_containing_barcodes, vec_seq_containing_barcodes - 1 + i)
  }
  
  # order lines correctly (so consecutive 4 lines correspond to the same fastq entry)
  fastq_containing_barcodes <- sort(fastq_containing_barcodes)
  
  # Output the correct entries
  output_data_r1 <- input_fastq_r1[fastq_containing_barcodes,]
  fwrite(output_data_r1, path_to_output_fastq_r1, col.names=FALSE, row.names=FALSE, sep="\n", quote=FALSE)
  # Keep the corresponding reads in read2
  output_data_r2 <- input_fastq_r2[fastq_containing_barcodes,]
  fwrite(output_data_r2, path_to_output_fastq_r2, col.names=FALSE, row.names=FALSE, sep="\n", quote=FALSE)
  
  if (provide_updates == TRUE) {
    cat(paste0("All done, find fastq file containing all reads with the specified barcodes at this location: ", path_to_output_fastq_r1, " and ", path_to_input_fastq_r2))
  }
}




demultiplex_reads <- function(read1_file, read2_file, output_read2_file, UMI_start, UMI_end, barcode_start, barcode_end) {
  library(data.table)
  # Read in the FASTQ files and create data.tables
  cat(paste0("Reading ", read1_file, " and ", read2_file, "\n"))
  read1_dt <- fread(read1_file, sep = "\n", header = FALSE, skip = 0, nrows = Inf)
  read2_dt <- fread(read2_file, sep = "\n", header = FALSE, skip = 0, nrows = Inf)
  
  # Extract the UMI and barcode
  cat("Extracting UMIs and barcodes from read1 sequence lines\n")
  UMIs <- substr(read1_dt$V1, UMI_start, UMI_end)[seq(2, nrow(read1_dt), by = 4)]
  barcodes <- substr(read1_dt$V1, barcode_start, barcode_end)[seq(2, nrow(read1_dt), by = 4)]
  
  # Update the read2 headers
  cat("Prepending UMIs and barcodes to read2 header lines\n")
  header_indices <- seq(1, nrow(read2_dt), by = 4)
  read2_dt[header_indices, V1 := paste0("@UMI_", UMIs, ":BC_", barcodes, ":", sub("^@", "", V1))]
  
  # Write the output lines to the new read2 file
  fwrite(read2_dt, output_read2_file, sep = "\n", col.names = FALSE, row.names = FALSE)
  cat(paste0("All done, find fastq file containing demultiplexed reads at this location: ", output_read2_file))
}



extract_reads_with_barcodes(path_to_input_fastq_r1="~/demultiplexData/01_R1.fastq", path_to_input_fastq_r2="~/demultiplexData/01_R2.fastq", path_to_barcodes_txt="~/demultiplexData/barcodes.txt", path_to_output_fastq_r1="~/demultiplexData/02_contain_bar_R1.fastq", path_to_output_fastq_r2="~/demultiplexData/02_contain_bar_R2.fastq", length_UMI = 6, char_trail = "T", length_trail = 8, provide_updates = TRUE)
demultiplex_reads("~/demultiplexData/02_contain_bar_R1.fastq", "~/demultiplexData/02_contain_bar_R2.fastq", "~/demultiplexData/03_demultiplexed.fastq", UMI_start=0, UMI_end=5, barcode_start=6, barcode_end=11)

