#' Searches a fastq file for reads containing specified barcodes
#'
#' @param path_to_input_fastq The path to the fastq file being searched
#' @param path_to_barcodes_txt The path to the text file containing the barcodes you want to search for
#' @param path_to_output_fastq The path to the file you want to store the outputted fastq file containing reads with the barcodes
#' @param length_UMI The length of the UMI at the beginning of the each read (default = 0 for reads containing no UMI)
#' @param char_trail The character trailing the barcode in each read (default = 'T'), used to avoid false positives
#' @param length_trail The minimum number of occurrences of the character trailing the barcode in each read (default = 0 for reads containing no trailing character), used to avoid false positives
#' @param provide_updates Leave as TRUE if you want to be updated on the progress of the function along with the number of reads containing each barcode (default = TRUE)
#' @return Does not return anything. Outputs a fastq file (located at param path_to_output_fastq) containing all reads with specified barcode
#' @examples
#' extract_reads_with_barcodes("~/rawData/read1.fastq", "~rawData/barcodes.txt", "~/rawData/reads_containing_barcodes.fastq", length_UMI = 6, length_trail = 6)

extract_reads_with_barcodes <- function(path_to_input_fastq, path_to_barcodes_txt, path_to_output_fastq, length_UMI = 0, char_trail = "T", length_trail = 0, provide_updates = TRUE) {
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
  input_fastq <- fread(path_to_input_fastq, sep="\n", header=FALSE)
  
  # Creates a list storing all sequences containing specified barcodes (updating progress along the way)
  seq_containing_barcodes <- list()
  for (i in 1:length(barcodes)) {
    bc <- mod_barcodes[i]
    matched_rows <- grep(bc, input_fastq$V1)
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
  output_data <- input_fastq[fastq_containing_barcodes,]
  fwrite(output_data, path_to_output_fastq, col.names=FALSE, row.names=FALSE, sep="\n", quote=FALSE)
  if (provide_updates == TRUE) {
    cat(paste0("All done, find fastq file containing all reads with the specified barcodes at this location: ", path_to_output_fastq))
  }
}
