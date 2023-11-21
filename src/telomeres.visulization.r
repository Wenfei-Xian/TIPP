args<-commandArgs(T)

library(ggplot2)
library(Biostrings)
library(stringdist)
library(stats)

read_fasta <- function(fasta_file_path) {
  dna_seq <- readDNAStringSet(fasta_file_path)
  lengths <- width(dna_seq)
  names <- names(dna_seq)
  sequences <- as.character(dna_seq)
  return(data.frame(name=names, length=lengths, sequence=sequences))
}

read_bed <- function(bed_file_path) {
  bed_df <- read.delim(bed_file_path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(bed_df) < 3) {
    stop("BED file must contain at least three columns.")
  }
  names(bed_df) <- c("chrom", "start", "end")
  bed_df$start <- as.numeric(bed_df$start)
  bed_df$end <- as.numeric(bed_df$end)
  bed_df$color <- "orange" 
  return(bed_df)
}

calculate_distance_matrix <- function(fasta_df) {
  distance_matrix <- stringdist::stringdistmatrix(fasta_df$sequence, fasta_df$sequence, method = "cosine")
  rownames(distance_matrix) <- fasta_df$name
  colnames(distance_matrix) <- fasta_df$name
  return(as.dist(distance_matrix)) 
}

cluster_sequences <- function(fasta_df) {
  distance_matrix <- calculate_distance_matrix(fasta_df)
  hc <- hclust(distance_matrix) 
  fasta_df$ordered_names <- fasta_df$name[hc$order] 
  fasta_df <- fasta_df[order(fasta_df$ordered_names),] 
  return(fasta_df)
}

plot_sequences <- function(fasta_df, bed_df) {
  fasta_df$name_factor <- factor(fasta_df$name, levels = fasta_df$ordered_names)
  
  bed_df$ymin <- match(bed_df$chrom, levels(fasta_df$name_factor)) - 0.15
  bed_df$ymax <- match(bed_df$chrom, levels(fasta_df$name_factor)) + 0.15
  
  max_length <- max(fasta_df$length)
  
  plot_height <- length(levels(fasta_df$name_factor)) * 0.4 + 1  
  
  p <- ggplot(data = fasta_df) +
    geom_rect(aes(xmin = 0, xmax = length, ymin = as.numeric(name_factor) - 0.15, ymax = as.numeric(name_factor) + 0.15), fill = "grey80") +
    geom_rect(data = bed_df, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = color), color = "black", size = 0.2) +
    geom_text(aes(x = length, y = as.numeric(name_factor) - 0.3, label = name), angle = 0, hjust = 1, size = 3, color = "black") +
    scale_fill_identity() +
    scale_x_continuous(breaks=seq(0, max_length, by=5000), labels=seq(0, max_length, by=5000)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle=45, hjust=1),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color="grey", linetype="solid"),
      panel.background = element_blank()
    ) +
    labs(y = "Telomere") +
    ylim(0.5, length(levels(fasta_df$name_factor)) + 0.5) 
  
  return(list(plot = p, height = plot_height))
}

fasta_file_path <- args[1]
bed_file_path <- args[2]

fasta_df <- read_fasta(fasta_file_path)
fasta_df <- cluster_sequences(fasta_df) 
bed_df <- read_bed(bed_file_path)

plot_data <- plot_sequences(fasta_df, bed_df)

ggsave(args[3], plot_data$plot, width = 6.69, height = plot_data$height, units = "in", limitsize = FALSE)
