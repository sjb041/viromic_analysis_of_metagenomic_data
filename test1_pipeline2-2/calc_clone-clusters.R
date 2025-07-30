# 基于聚类结果表 cluster.tsv 计算克隆数，
# 首先将实验组序列和对照组序列合并为一个序列文件，
# 对这个合并的序列文件99.4聚类，
# 取出那些同时含有实验组序列和对照组序列的cluster，
# 所有cluster中有多少条不同的实验组序列，就是实验组的克隆数，对照组同理。

# 输入：
# control_path（对照组序列文件的路径）
# treatment_path 的路径（实验组序列文件的路径）
# clusters_path（聚类结果表的路径）
# 输出：
# control_clone（对照组中的克隆数）
# treatment_clone（实验组中的克隆数）

library(seqinr)
library(dplyr)
library(optparse)

# =============== 定义命令行参数 ===============

option_list <- list(
  make_option(c("-c", "--clusters_path"), type = "character", 
              help = "Path to the clusters file (e.g., clusters.tsv)"),
  
  make_option(c("-t", "--treatment_path"), type = "character", 
              help = "Path to the treatment FASTA file"),
  
  make_option(c("-r", "--control_path"), type = "character", 
              help = "Path to the control FASTA file")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# 使用参数
clusters_path <- args$clusters_path
treatment_path <- args$treatment_path
control_path <- args$control_path

# =============== 定义 get_nameslist 函数 ===============
# 将contig_id保存为列表元素，
# 列表中每一个元素是一个contig_id

get_nameslist <- function(name_character){
  
  # 将字符串按 '' 分割，得到每一组
  groups <- unlist(strsplit(name_character, "''"))
  
  # 创建一个空的列表
  nameslist <- list()
  
  # 循环处理每一组数据
  for (group in groups) {
    # 将每一组按 ',' 分割
    nodes <- unlist(strsplit(group, ","))
    
    # 将分割后的节点添加到结果列表中
    nameslist <- c(nameslist, nodes)
  }
  
  return (nameslist)
}

  
# =============== 主要计算流程 ===============

# 读文件
clusters <- read.delim(file = clusters_path, header = F, stringsAsFactors = F, check.names = F)
control_seq <- read.fasta(file = control_path,
                          as.string = TRUE,
                          forceDNAtolower = FALSE,
                          whole.header = FALSE)   # 注意序列名不要取整个描述行
treatment_seq <- read.fasta(file = treatment_path,
                            as.string = TRUE,
                            forceDNAtolower = FALSE,
                            whole.header = FALSE)
  
# 把数据框 clusters 中第二列含有逗号的行提取出来得到新的数据框（建议检查第二列是否同时含有实验组和对照组序列）
clusters_filtered <- dplyr::filter(clusters, grepl(",", clusters[[2]]))
  
# 把 cluster_filtered 的第二列提取出来转化为列表 list
clusters_filtered_names <- get_nameslist(clusters_filtered[[2]])
# 把 序列文件的 names 转化为 list
control_seq_names <- get_nameslist(names(control_seq))
treatment_seq_names <- get_nameslist(names(treatment_seq))
  
# 计算 clone 个数
control_clone <- length(intersect(control_seq_names, clusters_filtered_names))
treatment_clone <- length(intersect(treatment_seq_names, clusters_filtered_names))
  
message("Clones of control: ", control_clone)
message("Clones of treatment: ", treatment_clone)

#message("Clones of ", tools::file_path_sans_ext(basename(control_path)), ": ", control_clone)
#message("Clones of ", tools::file_path_sans_ext(basename(treatment_path)), ": ", treatment_clone)

