library(tidyverse)

# Start with the df combined_df_edit 
gene <- combined_df_subset %>%
  filter(
    grepl("ferritin", Matched_Product, ignore.case = TRUE) |
      grepl("ferritin", upstream_Product, ignore.case = TRUE) |
      grepl("ferritin", downstream_Product, ignore.case = TRUE)
  )
gene <- gene[!is.na(gene$percentage),]
gene <- gene %>%
  mutate(
    week       = as.numeric(week),
    percentage = as.numeric(percentage)
  )
slopes <- gene %>%
  group_by(Sample, Position) %>%
  summarise(
    slope = coef(lm(percentage ~ week))[2],
    .groups = "drop"
  )
gene2 <- gene %>%
  left_join(slopes, by = c("Sample","Position")) %>%
  mutate(pct_adj = ifelse(slope < 0, 1 - percentage, percentage)) %>%
  filter(OLS_slope > abs(0.005)) %>%
  filter(OLS_pvalue <= abs(0.1)) %>%
  group_by(Sample, Species, week, Mutation_Type) %>%
  summarise(
    n        = n(),
    mean_pct = mean(pct_adj, na.rm = TRUE),
    se       = if_else(n > 1,
                       sd(pct_adj, na.rm = TRUE) / sqrt(n),
                       0),
    .groups = "drop"
  )
colnames(iHMP_metadata_disease_summarized)[1] <- "Sample"
gene2 <- left_join(gene2, iHMP_metadata_disease_summarized)
gene2$Mutation_Type[is.na(gene2$Mutation_Type) | gene2$Mutation_Type == ""] <- "Intergenic"


my_colors <- c(#"H4015"="#e0d7ca", "H4017"="#e0caab", "M2028"="#d4a96e","P6009"="#9e7741","C3001"="#e0d7ca", "H4006"="#d4a96e", "M2014"="#e0d7ca", "M2068"="#69491d", "M2034"="#e0d7ca",
               #"C3002"="#e0d7ca", "C3030"="#e0caab", "H4006"="#f0cb99","H4015"="#d4a96e","H4017"="#e39c39", "H4020"="#9e7741", "M2008"="#805821", "M2085"="#69491d",
               #"H4007"="#e0d7ca", "H4015"="#f0cb99", "H4017"="#d4a96e","H4038"="#9e7741","M2028"="#69491d",
               "H4007"="#e0d7ca", "H4015"="#d4a96e","M2028"="#69491d",               
               #"M2069"="#9f8fbd", "H4027"="#9f8fbd", "C3005"="#d5cce6", "H4044"="#9f8fbd", "M2083"="#634a91", "P6013" ="#9f8fbd", "P6012" ="#9f8fbd",
               "M2064"="#d5cce6", "M2069"="#9f8fbd", "P6012"="#634a91", "M2103"="#392957",
               "P6018"="#81b58c", "H4024"="#b2d6ba", "M2039"="#81b58c", "P6014"="#81b58c", "C3022"="#b2d6ba", "M2072"="#81b58c",
               "H4023"="#c3e0ca", "M2042"="#8abf97", "M2072"="#5db36f", "P6018"="#3c7548")

table(gene2$Sample, gene2$diagnosis)

ggplot(gene2, aes(x = week, y = mean_pct*100, color = Sample, shape = Species)) +
  geom_line(se=F) +
  geom_point(size = 1, stroke = 1, fill = NA, alpha=1) +
  scale_color_manual(values=my_colors) +
  labs(x = "Week", y = "Mutation Frequency (%)") +
  facet_wrap(~Mutation_Type, nrow=4) + 
  custom_theme +
  theme(legend.position = "right")

ggsave("../Figures/mutations/feoAB_ihmp_mutations_s005_p01.png", units = "in", width = 7, height = 6)

