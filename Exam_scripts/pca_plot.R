read_pca_output <- function(path){
  pca_df = read.table(path)
  names(pca_df) = c('Population', 'Individual', 'PC1', 'PC2', 'PC3')
  # enforce the order of populations
  pca_df$Individual <- factor(pca_df$Individual , levels =c('CEU','YRI', 'CHB', 'MXL'))
  pca_df = pca_df[order(pca_df$Individual),]
  return(pca_df)
}

plot_pca_basic <- function(pca_df, pc, lab, title){
  plot(pca_df$PC1 , pc, col = pca_df$Individual, pch = 16, main = title, xlab = "PC1", ylab = lab)
  legend(x="topleft", legend = levels(pca_df$Individual), fill = palette()[1:4])
}


pca_clean = read_pca_output('filter.eigenvec')
head(pca_clean)

# Make some plots and save them to a file

png("PC1v2pruned.png")
plot_pca_basic(pca_clean, pc = pca_clean$PC2, lab = "PC2", title = 'Principal Component 1 vs Principal Component 2')
dev.off()

png("PC1v2pruned.png")
plot_pca_basic(pca_clean, pc = pca_clean$PC3, lab = "PC3", title = 'Principal Component 1 vs Principal Component 3')
dev.off()

print("Done plotting PCAs")