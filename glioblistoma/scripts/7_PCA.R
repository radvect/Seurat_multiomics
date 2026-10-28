library(tidyverse)
library(tibble)

#here I took the top first HVG discovered previously, took 4 cells where they are nonzero and computed PCa with 2 ways to do that





#gene_names <- top_genes[1:2]
#nonzero <- Matrix::colSums(
#  rna_norm[gene_names, , drop = FALSE] > 0
#) == 2
#cell_indices <- which(nonzero)[1:4]
#cells <- head(colnames(rna_norm)[nonzero], 4)

#cells

cells <- colnames(rna_norm)[c(14, 44, 62, 65)]

gene_names <- top_genes[1:2]

Gene_A <- as.numeric(rna_norm[gene_names[1], cells])
Gene_B <- as.numeric(rna_norm[gene_names[2], cells])

expression_tibble <- tibble(
  cells,
  !!gene_names[1] := as.numeric(rna_norm[gene_names[1], cells]),
  !!gene_names[2] := as.numeric(rna_norm[gene_names[2], cells])
)


# ggplot(expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
#   geom_point(color = "cornflowerblue") +
#   geom_text(hjust = 0, vjust = -1) +
#   xlim(0, 2) +
#   ylim(0, 2) +
#   theme_bw() +
#   xlab("Gene A") +
#   ylab("Gene B") +
#   ggtitle("Example Expression Values from Four Cells") +
#   theme(plot.title = element_text(hjust = 0.5))


# Gene_A_mean <- mean(expression_tibble[[gene_names[1]]])
# Gene_B_mean <- mean(expression_tibble[[gene_names[2]]])
# 
# center_of_data <- c(Gene_A_mean, Gene_B_mean)
# names(center_of_data) <- gene_names

# center_of_data

# 
# p <- ggplot(expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
#   geom_point(color = "cornflowerblue") +
#   annotate("point", x = Gene_A_mean, y = Gene_B_mean, color = "red", size = 3) +
#   geom_text(hjust = 0, vjust = -1) +
#   annotate("text", x = Gene_A_mean, y = Gene_B_mean, color = "red", label="New Center", hjust = 0, vjust = -1) +
#   xlim(0, 2) +
#   ylim(0, 2) +
#   theme_bw() +
#   xlab("Gene A") +
#   ylab("Gene B") +
#   ggtitle("Example Expression Values from Four Cells") +
#   theme(plot.title = element_text(hjust = 0.5))

# recentered_expression_tibble <- expression_tibble %>%
#   mutate(
#     !!gene_names[1] := .data[[gene_names[1]]] - Gene_A_mean,
#     !!gene_names[2] := .data[[gene_names[2]]] - Gene_B_mean
#   )
 
# recentered_expression_tibble



# p1 <- ggplot(recentered_expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
#   geom_point(color = "cornflowerblue") +
#   annotate("point", x = 0, y = 0, color = "red", size = 3) +
#   geom_text(hjust = 0, vjust = -1) +
#   annotate("text", x = 0, y = 0, color = "red", label="New Center", hjust = 0, vjust = -1) +
#   xlim(-2, 2) +
#   ylim(-2, 2) +
#   theme_bw() +
#   xlab("Gene A") +
#   ylab("Gene B") +
#   ggtitle("Example Re-centered Expression Values from Four Cells") +
#   theme(plot.title = element_text(hjust = 0.5))
# 


# recentered_expression_matrix <- recentered_expression_tibble %>% 
#   column_to_rownames("cells") %>% 
#   as.matrix()


# cov_matrix <- cov(recentered_expression_matrix)
# 
# 
# print(cov_matrix)
# covariance <- cov(recentered_expression_matrix[,  gene_names[1]], recentered_expression_matrix[,  gene_names[2]])
# eig <- eigen(cov_matrix)
# print(eig)
# 
# pc_scores <- recentered_expression_matrix %*% eig$vectors
# colnames(pc_scores) <- c("PC_1", "PC_2")

# print(pc_scores)


# pct_var_explained <- (eig$values / sum(eig$values)) * 100
# 
# names(pct_var_explained) <- c("PC_1", "PC_2")
# print(pct_var_explained)



prcomp_PCA <- expression_tibble %>% 
  column_to_rownames("cells") %>% 
  prcomp()

print(prcomp_PCA$center)

prcomp_eigenvalues <- prcomp_PCA$sdev ** 2
names(prcomp_eigenvalues) <- c("PC_1", "PC_2")

print(prcomp_eigenvalues)
print(prcomp_PCA$rotation)
print(prcomp_PCA$x)


pct_var_explained <- (prcomp_eigenvalues / sum(prcomp_eigenvalues)) * 100

names(pct_var_explained) <- c("PC_1", "PC_2")
print(pct_var_explained)





