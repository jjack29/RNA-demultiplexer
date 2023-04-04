# RNA-demultiplexer
R package that demultiplexes RNA reads in the format &lt;UMI>&lt;sample barcode>&lt;polyN>&lt;cDNA>

Function:
extract_reads_with_barcodes(path_to_input_fastq, path_to_barcodes_txt, path_to_output_fastq, length_UMI = 0, char_trail = "T", length_trail = 0, provide_updates = TRUE)
