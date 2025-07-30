# 基于比对结果表 ani.tsv 计算克隆数
# 首先筛选比对结果中 pid >= 99.4, tcov >= 85
# 删除 qname 和 tname 同时属于 control_seq 或 treatment_seq 的行（去除内部匹配）
# 取数据框的 qname 和 tname，合并然后去重得到 unique_names
# unique_names 中有几个属于实验组，那么实验组就有几条克隆，对照组同理

# 输入：
# control_path（对照组序列文件的路径）
# treatment_path 的路径（实验组序列文件的路径）
# ani_path（比对结果表的路径）
# 输出：
# control_clone（对照组中的克隆数）
# treatment_clone（实验组中的克隆数）

library(seqinr)
library(dplyr)
library(optparse)

# =============== 定义命令行参数 ===============

option_list <- list(
  make_option(c("-a", "--ani_path"), type = "character", 
              help = "Path to the clusters file (e.g., ani.tsv)"),
  
  make_option(c("-t", "--treatment_path"), type = "character", 
              help = "Path to the treatment FASTA file"),
  
  make_option(c("-r", "--control_path"), type = "character", 
              help = "Path to the control FASTA file")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# 使用参数
ani_path <- args$ani_path
treatment_path <- args$treatment_path
control_path <- args$control_path

# =============== 主要计算流程 ===============
  
# 读取文件  
control_seq <- read.fasta(file = control_path,
                          as.string = TRUE,
                          forceDNAtolower = FALSE,
                          whole.header = FALSE)     # 注意序列名不要取整个描述行
treatment_seq <- read.fasta(file = treatment_path,
                            as.string = TRUE,
                            forceDNAtolower = FALSE,
                            whole.header = FALSE)
ani <- read.delim(file = ani_path, header = T, stringsAsFactors = F, check.names = F)
  
# 相似度筛选
ani <- dplyr::filter(ani, pid >= 99.4, tcov >= 85 )
# 删除 qname 和 tname 同时属于 control_seq 或 treatment_seq 的行（去除内部匹配）
ani_filtered <- ani[!((ani$qname %in% names(control_seq) & ani$tname %in% names(control_seq)) |
                      (ani$qname %in% names(treatment_seq) & ani$tname %in% names(treatment_seq))), ]
  
# 取数据框的两列合并然后去重，再计算数量
unique_names <- unique(c(ani_filtered$qname, ani_filtered$tname))
control_clone <- length(intersect(names(control_seq), unique_names))
treatment_clone <- length(intersect(names(treatment_seq), unique_names))
  
  
message("Clones of control: ", control_clone)
message("Clones of treatment: ", treatment_clone)
