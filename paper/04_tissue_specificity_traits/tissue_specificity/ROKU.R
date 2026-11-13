library(TCC)
library(tidyr)
library(dplyr)
library(stringr)

### tissue specific==============================================================
t <- read.table("v11_raw_fit_fit_results.txt",sep="\t", header=T)

# Convert from long to wide format
tt <- t %>%
  pivot_wider(
    id_cols = promoter,
    names_from = tissue,
    values_from = overall_mean
  )

df <- as.data.frame(tt)
row.names(df) <- df$promoter
d <- df[,-1]

# ROKU
result <- ROKU(d, upper.limit = 0.25, sort = TRUE)

# Z_test
roku_ztest <- function(d, result, alpha = 0.05) {
  data <- as.matrix(d)
  outliers <- result$outlier
  outliers <- outliers[rownames(data), colnames(data)]
  
  true_outliers <- outliers
  
  # Create a matrix to hold the p values
  p_values <- matrix(NA, nrow = nrow(outliers), ncol = ncol(outliers))
  rownames(p_values) <- rownames(outliers)
  colnames(p_values) <- colnames(outliers)
  
  # Process each promoter (row)
  for(i in 1:nrow(data)) {
    print(paste0("Processing promoter ", i, " of ", nrow(data)))
    # Get the expression values of this promoter in all tissues
    promoter_data <- data[i, ]
    
    # Get the ROKU outlier results of this promoter
    promoter_outliers <- outliers[i, ]
    
    # Compute the mean and SD of the non-outlier tissues (those marked 0)
    normal_indices <- which(promoter_outliers == 0 & !is.na(promoter_data))
    outlier_indices <- which(!is.na(promoter_outliers) & promoter_outliers != 0)

    if(length(outlier_indices) == 0) {next}
    
    if(length(normal_indices) >= 2) {  # Make sure there are enough samples
      # Use the mean and SD
      normal_mean <- mean(promoter_data, na.rm = TRUE)
      # Add a small constant to the SD to avoid division by 0
      normal_sd <- max(sd(promoter_data, na.rm = TRUE), 1e-10)
      
      # Run the z test only on data points already marked as outliers
      for(j in outlier_indices) {
        # Handle only data points marked as outlier (-1 or 1) and not NA
        if(!is.na(promoter_outliers[j]) && promoter_outliers[j] != 0 && !is.na(promoter_data[j])) {
          # Compute the z score
          z_score <- (promoter_data[j] - normal_mean) / normal_sd
          
          # Compute the two-tailed p value
          p_value <- 2 * pnorm(-abs(z_score))
          p_values[i, j] <- p_value
          
          # If the p value is larger than alpha, it is not a significant outlier; mark it 0
          if(p_value > alpha) {
            print(paste0("False Positive for ", names(promoter_outliers)[j]))
            true_outliers[i, j] <- 0  # Reset to non-outlier
          }
        }
      }
    } else {
      # If there are not enough non-outlier tissues to serve as reference, mark all outliers of this promoter 0
        true_outliers[i, outlier_indices] <- 0
    }
  }
  
  # Return the results
  return(list(
    true_outliers = true_outliers,
    p_values = p_values
  ))
}


# Run the z test and get the results
true_outlier_result <- roku_ztest(d, result)

result$true_outliers <- true_outlier_result$true_outliers
result$p_values <- true_outlier_result$p_values

saveRDS(result, file="v11_tissue_specific_ROKU_results.rds")



### tissues plots==============================================================================
library(data.table)
load("/home/users/nus/e1327989/YuhaoDong/GTEx/v11_filter_recal/Promoter_info_filter_5.RData")
head(info)

load("/home/users/nus/e1327989/YuhaoDong/GTEx/v11_filter_recal/filter_5_recal_att_filter.RData")
#recal_att_f <- fread("tissue_specific/v10_raw_fit_fit_results.txt")
recal_att_f <- as.data.frame(recal_att_f)
head(recal_att_f)[,1:8]


t <- readRDS("v11_tissue_specific_ROKU_results.rds")

out <- as.data.frame(t$true_outliers)
rank <- as.data.frame(t$rank)
colnames(rank) <- "rank"
rank$test <- 1
head(rank,20)
head(rank[row.names(rank) %in% row.names(subset(out, Brain_Cortex != 0)), ])

###
head(out)
tt <- as.data.frame(t(out))
head(tt)[,1:8]

tt$rowsum <- 0
for(i in 1:nrow(tt)){
  s <- sum(abs(tt[i,]), na.rm=T)
  tt[i,]$rowsum <- s
}
tt$tissue <- row.names(tt)
ttt <- tt[c('tissue', 'rowsum')]
ttt <- ttt[order(ttt$rowsum, decreasing = T),]
ttt$tissue <- factor(ttt$tissue, levels = ttt$tissue[order(ttt$rowsum, decreasing = F)])
head(ttt)


### tissues plots
library(ggplot2)

### top 20 bar
library(dplyr)
library(ggplot2)
num <- 20
# 1. Sort and select the top 20
top <- ttt %>%
  arrange(desc(rowsum)) %>%
  slice(1:num)

# 2. Merge the remaining tissues into "Others"
others_sum <- ttt %>%
  arrange(desc(rowsum)) %>%
  slice((num + 1):n()) %>%
  summarise(rowsum = sum(rowsum)) %>%
  mutate(tissue = "Others")

# 3. Combine the top 20 and Others
plot_data <- rbind(top, others_sum)
  # Reorder the factor so that "Others" comes last
  plot_data$tissue <- factor(plot_data$tissue, levels = rev(plot_data$tissue))

# 4. Plot
plots <- ggplot(plot_data, aes(x = tissue, y = rowsum)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = paste0("Tissue Specific Count of ",sum(plot_data$rowsum)," promoters"),
       x = NULL, y = "Count") +
  theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
  coord_flip()

pdf("v11_top20_bar.pdf")
print(plots)
dev.off()



###single promoter====================================================================
library(ggplot2)
load("/home/users/nus/e1327989/YuhaoDong/GTEx/v11_filter_recal/Promoter_info_filter_5.RData")
load("/home/users/nus/e1327989/YuhaoDong/GTEx/v11_filter_recal/Promoter_info_all_v11.RData")


subset(info, geneId == "ENSG00000114331.16")
subset(info_all, geneId == "ENSG00000114331.16" & internalPromoter == FALSE)
# TMEM38B
info_all[grep("ENSG00000095209", info_all$geneId),]
# TMEM161B
info_all[grep("ENSG00000164180", info_all$geneId),]
# TMEM161B-DT
info_all[grep("ENSG00000247828", info_all$geneId),]



raw <- read.table("v11_raw_fit_fit_results.txt",sep="\t", header=T)
head(raw)

result <- readRDS("v11_tissue_specific_ROKU_results.rds")
out <- as.data.frame(result$true_outliers)

for(i in c(24,25,26,27,28,29,30,31,32)){
    p <- paste0("PromoterID_417", i);
    o <- as.data.frame(t(out[p,]));
    print(paste0(p, ": ",sum(o[,1], na.rm=T)))
    }

promoterID <- "PromoterID_41728"

o <- as.data.frame(t(out["PromoterID_41728",]))
oo <- subset(o, o[1] > 0 | o[1] < 0)


p <- subset(raw, promoter == "PromoterID_41728")

p$sig <- "none"
p[p$tissue %in% row.names(oo),]$sig <- "sig"


g <- ggplot(p, aes(x = tissue, y = overall_mean, group = 1, color = sig)) + 
  geom_point(size = 3) +
  geom_line(aes(group = 1), color = "black", linetype = "dashed") +  # Add connecting lines
  theme_minimal() +
  labs(title = promoterID) +
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("sig" = "red", "none" = "grey"))


pdf("AKAP6_PromoterID_41728.pdf", width = 12, height = 6)
print(g)
dev.off()


info_all["PromoterID_17115",]
info[grep("ENSG00000290418", info$geneId),c(1:7)]
subset(info_all[grep("ENSG00000290418", info_all$geneId),c(1:7)], internalPromoter == FALSE)

p1 <- subset(raw, promoter == "PromoterID_134459")
p2 <- subset(raw, promoter == "PromoterID_134460")

p3 <- subset(raw, promoter == "PromoterID_134462")

p4 <- subset(raw, promoter == "PromoterID_134466")

p5 <- subset(raw, promoter == "PromoterID_101380")


##
pp <- ggplot() +
  ## ---------- p1 ----------
geom_point(
  data = p1,
  aes(x = tissue, y = overall_mean, color = "p1"),  # Put the colour in aes, with the value set as a character string
  size = 3
) +
  geom_line(
    data = p1,
    aes(x = tissue, y = overall_mean, group = 1, color = "p1"),
    linetype = "dashed"
  ) +
  
  ## ---------- p2 ----------
geom_point(
  data = p2,
  aes(x = tissue, y = overall_mean, color = "p2"),
  size = 3
) +
  geom_line(
    data = p2,
    aes(x = tissue, y = overall_mean, group = 1, color = "p2"),
    linetype = "dashed"
  ) +
  
  ## ---------- Common theme ----------
theme_minimal() +
  labs(title = "PromoterID_17115-PromoterID_17117",
       y = "mean value",
       color = "Dataset") +          # Legend title
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  ## ---------- Set colours manually ----------
scale_color_manual(
  values = c(p1 = "red", p2 = "blue"),
  labels = c(p1 = "PromoterID_17115",
             p2 = "PromoterID_17117")
)

pdf("ENSG00000108306.13.pdf", width = 12, height = 6)
print(pp)
dev.off()




pp <- ggplot() +
  ## ---------- p1 ----------
geom_point(
  data = p1,
  aes(x = tissue, y = overall_mean, color = "p1"),  # Put the colour in aes, with the value set as a character string
  size = 3
) +
  geom_line(
    data = p1,
    aes(x = tissue, y = overall_mean, group = 1, color = "p1"),
    linetype = "dashed"
  ) +
  
  ## ---------- p2 ----------
geom_point(
  data = p2,
  aes(x = tissue, y = overall_mean, color = "p2"),
  size = 3
) +
  geom_line(
    data = p2,
    aes(x = tissue, y = overall_mean, group = 1, color = "p2"),
    linetype = "dashed"
  ) +
  
  ## ---------- p3 ----------
geom_point(
  data = p3,
  aes(x = tissue, y = overall_mean, color = "p3"),
  size = 3
) +
  geom_line(
    data = p3,
    aes(x = tissue, y = overall_mean, group = 1, color = "p3"),
    linetype = "dashed"
  ) +
  
  ## ---------- Common theme ----------
theme_minimal() +
  labs(title = "PromoterID_134459-60-62",
       y = "mean value",
       color = "Dataset") +          # Legend title
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  ## ---------- Set colours manually ----------
scale_color_manual(
  values = c(p1 = "red", p2 = "blue", p3 = "green"),
  labels = c(p1 = "PromoterID_134459",
             p2 = "PromoterID_134460",
             p3 = "PromoterID_134462")
)


pdf("ENSG00000290418.2.pdf", width = 12, height = 6)
print(pp)
dev.off()



##
pp <- ggplot() +
  ## ---------- p1 ----------
geom_point(
  data = p1,
  aes(x = tissue, y = overall_mean, color = "p1"),  # Put the colour in aes, with the value set as a character string
  size = 3
) +
  geom_line(
    data = p1,
    aes(x = tissue, y = overall_mean, group = 1, color = "p1"),
    linetype = "dashed"
  ) +
  
  ## ---------- p2 ----------
geom_point(
  data = p2,
  aes(x = tissue, y = overall_mean, color = "p2"),
  size = 3
) +
  geom_line(
    data = p2,
    aes(x = tissue, y = overall_mean, group = 1, color = "p2"),
    linetype = "dashed"
  ) +
  
  ## ---------- p3 ----------
geom_point(
  data = p3,
  aes(x = tissue, y = overall_mean, color = "p3"),
  size = 3
) +
  geom_line(
    data = p3,
    aes(x = tissue, y = overall_mean, group = 1, color = "p3"),
    linetype = "dashed"
  ) +
  
  ## ---------- p4 ----------
geom_point(
  data = p4,
  aes(x = tissue, y = overall_mean, color = "p4"),
  size = 3
) +
  geom_line(
    data = p4,
    aes(x = tissue, y = overall_mean, group = 1, color = "p4"),
    linetype = "dashed"
  ) +

  ## ---------- Common theme ----------
theme_minimal() +
  labs(title = "PromoterID_109388-109389-109392-109393",
       y = "mean value",
       color = "Dataset") +          # Legend title
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  ## ---------- Set colours manually ----------
scale_color_manual(
  values = c(p1 = "red", p2 = "blue", p3 = "green"),
  labels = c(p1 = "PromoterID_109388",
             p2 = "PromoterID_109389",
             p3 = "PromoterID_109392",
             p4 = "PromoterID_109393")
)



pdf("ENSG00000257261.8.pdf", width = 12, height = 6)
print(pp)
dev.off()




##
pp <- ggplot() +
  ## ---------- p1 ----------
geom_point(
  data = p1,
  aes(x = tissue, y = overall_mean, color = "p1"),  # Put the colour in aes, with the value set as a character string
  size = 3
) +
  geom_line(
    data = p1,
    aes(x = tissue, y = overall_mean, group = 1, color = "p1"),
    linetype = "dashed"
  ) +
  
  ## ---------- p2 ----------
geom_point(
  data = p2,
  aes(x = tissue, y = overall_mean, color = "p2"),
  size = 3
) +
  geom_line(
    data = p2,
    aes(x = tissue, y = overall_mean, group = 1, color = "p2"),
    linetype = "dashed"
  ) +
  
  ## ---------- p3 ----------
geom_point(
  data = p3,
  aes(x = tissue, y = overall_mean, color = "p3"),
  size = 3
) +
  geom_line(
    data = p3,
    aes(x = tissue, y = overall_mean, group = 1, color = "p3"),
    linetype = "dashed"
  ) +
  
  ## ---------- p4 ----------
geom_point(
  data = p4,
  aes(x = tissue, y = overall_mean, color = "p4"),
  size = 3
) +
  geom_line(
    data = p4,
    aes(x = tissue, y = overall_mean, group = 1, color = "p4"),
    linetype = "dashed"
  ) +

  ## ---------- p5 ----------
geom_point(
  data = p5,
  aes(x = tissue, y = overall_mean, color = "p5"),
  size = 3
) +
  geom_line(
    data = p5,
    aes(x = tissue, y = overall_mean, group = 1, color = "p5"),
    linetype = "dashed"
  ) +

  ## ---------- Common theme ----------
theme_minimal() +
  labs(title = "PromoterID_5026-PromoterID_5035",
       y = "mean value",
       color = "Dataset") +          # Legend title
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  ## ---------- Set colours manually ----------
scale_color_manual(
  values = c(p1 = "red", p2 = "blue", p3 = "green", p4 = "orange", p5 = "pink"),
  labels = c(p1 = "PromoterID_5026",
             p2 = "PromoterID_5027",
             p3 = "PromoterID_5033",
             p4 = "PromoterID_5034",
             p5 = "PromoterID_5035")
)



pdf("ENSG00000065809.14.pdf", width = 12, height = 6)
print(pp)
dev.off()







r <- cbind(recal_att_f[1:4], recal_att_f[promoterID])
r <- na.omit(r)
r$sig <- "none"
r[r$SMTSD %in% row.names(oo),]$sig <- "sig"

head(r)


ggplot(r, aes(x = SMTSD, y = r[[promoterID]], fill = sig)) + 
  geom_violin(width=3, alpha=0.7) + theme_minimal() +
  labs(title = promoterID) +
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("sig" = "red", "none" = "grey"))



info["PromoterID_20342",]
subset(info, geneId == "ENSG00000114331.16")




























