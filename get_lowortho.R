library(ggplot2)
library(cowplot)
# read file. Each line one orthogroup
ortho = scan("Orthogroups.txt", what=character(), sep="\n")
ortho = strsplit(ortho, " ") # split genes
# exclude orthogroup name and keep species name only
ortho_spp = lapply(ortho, function(x) sapply(x[2:length(x) ], function(y) strsplit(y, "\\.")[[1]][1] ))
# vector of unique species names
spp = unique(unlist(ortho_spp))
# subsample orthogroups that include all species
ortho_spp_complete = ortho_spp[ sapply(ortho_spp, function(x) all(spp %in% x) ) ]
ortho_spp_complete_len = sapply(ortho_spp_complete, length) # get lengths from complete data set
# ggplot2
# histogram of orthogroup sizes
theme_set(theme_gray(base_size = 18))
ggplot(data=NULL, aes(ortho_spp_complete_len[-1])) + # exclude first (biggest)
  geom_histogram(bins=100) + xlab("number of gene copies") +
  geom_label(data=NULL, aes(x=1000, y=1000, label="3278 orthogroups"))
dev.copy2pdf(file="histogram_almost_all_orthogroups.pdf") # save to file
# group size vs number of orthogroups
df.len = data.frame(size=60:80, number=sapply(60:80, function(x) sum(ortho_spp_complete_len == x) ))
ggplot(df.len, aes(size, number)) +
  geom_point() +geom_line() +
  scale_x_continuous(breaks=60:80) +
  ylab("number of orthogroups") +
  xlab("number of gene copies")
dev.copy2pdf(file="scatterplot_ortho_vs_copies.60-80.pdf") # save to file
df.len = data.frame(size=60:200, number=sapply(60:200, function(x) sum(ortho_spp_complete_len == x) ))
ggplot(df.len, aes(size, number)) +
  geom_point() +geom_line() +
  scale_x_continuous(breaks=seq(60, 200, 10)) +
  ylab("number of orthogroups") +
  xlab("number of gene copies")
dev.copy2pdf(file="scatterplot_ortho_vs_copies.60-200.pdf") # save to file
