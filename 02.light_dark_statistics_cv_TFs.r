bb <- read.csv("count.tables/All.TFs.in.light.exp.csv", header = TRUE)
head(bb)
dim(bb)

colData <- data.frame(
  Sample = colnames(bb)[-1],
  Group = factor(
    c(rep("constant", 8), rep("natural", 8)),
    levels = c("constant", "natural")
  )
)


library(dplyr)
library(tidyr)

## filter out genes at least 3 smaples  >0 
expressed_genes <- norm_data_long %>%
  group_by(Group, Gene) %>%
  summarise(non_zero_count = sum(Expression > 0), .groups = "drop") %>%
  group_by(Gene) %>%
  filter(all(non_zero_count >= 3)) %>%  
  distinct(Gene)

## trans to long 


norm_data_long <- bb %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(colData, by = "Sample")

##### 
norm_data_long_filtered <- norm_data_long %>%
  filter(Gene %in% expressed_genes$Gene)



############# batch calculate the statistics result of each gene

stats_summary <- norm_data_long_filtered %>%
  group_by(Group, Gene) %>%
  summarise(
    Mean = mean(Expression),
    SD = sd(Expression),
    CV = ifelse(Mean > 0, SD / Mean, NA), 
    Min = min(Expression),
    Max = max(Expression),
    Range = Max - Min,
    .groups = "drop"
  )


head(stats_summary)

write.csv(stats_summary, file="All.TFs.in_light_dark.statistics.csv")

# How many genes CV >1 in each group
cv_above_0.5_counts <- stats_summary %>%
  filter(CV > 0.5 & !is.na(CV)) %>%
  group_by(Group) %>%
  summarise(
    High_Var_Gene_Count = n(),
    .groups = "drop"
  )

#### transfer to wide table 
stats_summary_wide <- stats_summary %>%
  pivot_wider(
    names_from = Group,
    values_from = c(Mean, SD, CV, Min, Max, Range),
    names_glue = "{.value}_{Group}"
  )


## output the result
write.csv(stats_summary_wide, "count.tables/All.TFs_in.light_dark_stats_summary.csv", row.names = FALSE)

head(stats_summary_wide)


### select the various gene 
high_var_genes_with_group <- stats_summary %>%
  filter(CV > 0.5 & !is.na(CV)) %>%
  select(Gene, Group, CV)

print(high_var_genes_with_group)
write.csv(high_var_genes_with_group, "count.tables/lifecycle_high_var_genes_with_group.csv", row.names = FALSE)



## CV distribution 
library(ggplot2)
library(dplyr)


ggplot(stats_summary %>% filter(!is.na(CV)), aes(x = Group, y = CV, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of CV Across Genes (DESeq2 Normalized)",
       x = "Group", y = "Coefficient of Variation (CV)") +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")) +
  theme(legend.position = "none")
ggsave("cv_distribution_light_dark.pdf", width = 8, height = 6)


cv_summary <- stats_summary %>%
  filter(!is.na(CV)) %>%
  group_by(Group) %>%
  summarise(
    Median_CV = median(CV),
    Q1_CV = quantile(CV, 0.25),
    Q3_CV = quantile(CV, 0.75),
    Max_CV = max(CV)
  )
print(cv_summary)
