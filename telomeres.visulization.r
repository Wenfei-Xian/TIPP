args<-commandArgs(T)

library(ggplot2)
library(Biostrings)
library(stringdist)
library(stats)

# 函数：读取FASTA文件并包含序列
read_fasta <- function(fasta_file_path) {
  dna_seq <- readDNAStringSet(fasta_file_path)
  lengths <- width(dna_seq)
  names <- names(dna_seq)
  sequences <- as.character(dna_seq) # 提取序列作为字符向量
  return(data.frame(name=names, length=lengths, sequence=sequences))
}

# 函数：读取BED文件
read_bed <- function(bed_file_path) {
  bed_df <- read.delim(bed_file_path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(bed_df) < 3) {
    stop("BED file must contain at least three columns.")
  }
  names(bed_df) <- c("chrom", "start", "end")
  bed_df$start <- as.numeric(bed_df$start)
  bed_df$end <- as.numeric(bed_df$end)
  bed_df$color <- "orange" # 设置颜色
  return(bed_df)
}

# 新函数：计算FASTA序列间的距离矩阵
calculate_distance_matrix <- function(fasta_df) {
  # 使用stringdist包计算基于字符串的距离
  distance_matrix <- stringdist::stringdistmatrix(fasta_df$sequence, fasta_df$sequence, method = "cosine")
  rownames(distance_matrix) <- fasta_df$name
  colnames(distance_matrix) <- fasta_df$name
  return(as.dist(distance_matrix)) # 转换为距离对象
}

# 新函数：对FASTA序列进行聚类排序
cluster_sequences <- function(fasta_df) {
  distance_matrix <- calculate_distance_matrix(fasta_df)
  hc <- hclust(distance_matrix) # 层次聚类
  fasta_df$ordered_names <- fasta_df$name[hc$order] # 根据聚类结果重新排序名称
  fasta_df <- fasta_df[order(fasta_df$ordered_names),] # 重新排序数据框
  return(fasta_df)
}

# 修改后的绘图函数
plot_sequences <- function(fasta_df, bed_df) {
  # 创建序列名的因子，确保顺序
  fasta_df$name_factor <- factor(fasta_df$name, levels = fasta_df$ordered_names)
  
  # 添加ymin和ymax到bed数据框，用于定位rect
  bed_df$ymin <- match(bed_df$chrom, levels(fasta_df$name_factor)) - 0.15
  bed_df$ymax <- match(bed_df$chrom, levels(fasta_df$name_factor)) + 0.15
  
  # 计算绘图所需的最大长度
  max_length <- max(fasta_df$length)
  
  # 根据序列的数量动态设置PDF文件的高度
  plot_height <- length(levels(fasta_df$name_factor)) * 0.4 + 1  # 增加额外的空间用于放置文本
  
  # 创建ggplot对象
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
    ylim(0.5, length(levels(fasta_df$name_factor)) + 0.5) # 为了确保每个矩形有足够的空间
  
  # 返回绘图对象和绘图高度
  return(list(plot = p, height = plot_height))
}

# 主流程保持不变，只是调用了修改后的绘图函数
fasta_file_path <- args[1]
bed_file_path <- args[2]

fasta_df <- read_fasta(fasta_file_path)
fasta_df <- cluster_sequences(fasta_df) # 使用聚类排序
bed_df <- read_bed(bed_file_path)

# 绘制图形
plot_data <- plot_sequences(fasta_df, bed_df)

# 输出图形到PDF，文件宽度设为6.69 = 170mm，高度根据plot_data的height决定
ggsave(args[3], plot_data$plot, width = 6.69, height = plot_data$height, units = "in", limitsize = FALSE)

