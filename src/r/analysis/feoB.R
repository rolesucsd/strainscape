library(data.table)
feoB <- snvs[snvs$gene %like% "feoB",]
feoB <- left_join(feoB, meta)

#table(left_join(snvs[!duplicated(snvs[,c("patient_id","bin")]),],meta)$group)
#CD     nonIBD     UC 
#1015    693    608 

feoB_map <- feoB[!duplicated(feoB[,c("patient_id","bin","mutation_type","is_sweep","chromosome")]),]
feoB_table <- as.data.frame(table(feoB_map$group, feoB_map$mutation_type, feoB_map$is_sweep))
feoB_table[feoB_table$Var1 %in% "nonIBD", "Freq"] <- feoB_table[feoB_table$Var1 %in% "nonIBD", "Freq"] / 693 * 100
feoB_table[feoB_table$Var1 %in% "UC", "Freq"] <- feoB_table[feoB_table$Var1 %in% "UC", "Freq"] / 608 * 100
feoB_table[feoB_table$Var1 %in% "CD", "Freq"] <- feoB_table[feoB_table$Var1 %in% "CD", "Freq"] / 1015 * 100

feoB_table <- feoB_table %>%
  mutate(Var1 = factor(Var1, levels = c("nonIBD", "UC", "CD")))
ggplot(feoB_table, aes(x = Var2, y = Freq, fill = Var1)) +
  facet_wrap(~Var3, scales = "free_y", labeller = labeller(Var3 = c(`TRUE` = "Sweeping mutation", `FALSE` = "Non-sweeping mutation"))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("UC" = "#bb9df5",
                               "nonIBD" = "#71d987",
                               "CD" = "#ffc169")) +
  labs(x = NULL, y = "% of isolates") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 11)
  )

ggsave(file.path(OUT, "feoB.png"), width = 4.75, height = 3.5, dpi = 320)

# For each mutation that occurd intragenic, find the gene length, and the position of the mutation
feoB_mutation <- feoB %>%
  filter(coding_status == "Coding") %>%
  select(c(1:5, 22:39, 79:85)) %>%
  mutate(
    gene_length = abs(as.numeric(stop) - as.numeric(start)),
    gene_position = if_else(
      strand == "+",
      as.numeric(position) - as.numeric(start),
      as.numeric(stop) - as.numeric(position)
    ),
    mutation_type = case_when(
      mutation_type %in% c("Nonsense", "Missense") ~ "Nonsynonymous",
      TRUE ~ mutation_type
    )
  ) %>%
  filter(gene_length < 5000) %>%
  filter(gene_length > 1)

ggplot(feoB_mutation, aes(x=gene_length))+
  geom_histogram() +
  scale_x_continuous(breaks = pretty_breaks(n = 20))

ggplot(feoB_mutation, aes(x=gene_position/3, fill=is_sweep)) +
  geom_histogram(bins=70) +
#  geom_vline(xintercept=1132/3, linetype = "dashed") +
#  geom_vline(xintercept=1325/3, linetype = "dashed") +
#  geom_vline(xintercept=1766/3, linetype = "dashed") +
#  geom_vline(xintercept=1989/3, linetype = "dashed") +
#  geom_vline(xintercept=2092/3, linetype = "dashed") +
  scale_fill_manual(values=c("grey","red")) +
  labs(x = NULL, y = "No. of variants") +
  facet_wrap(~mutation_type, ncol = 1, scales="free") +
  scale_x_continuous(breaks = pretty_breaks(n = 8))
ggsave(file.path(OUT, "alr_length.png"), width = 6, height = 3.5, dpi = 320)


# Looking at the trend over time
feoB_map <- feoB_map %>%
  pivot_longer(
    cols = all_of(names(feoB_map)[suppressWarnings(!is.na(as.numeric(names(feoB_map))))]),
    names_to = "week",
    values_to = "proportion"
  )
feoB_map <- feoB_map[!is.na(feoB_map$proportion),]
feoB_map$week <- as.numeric(feoB_map$week)

ggplot(feoB_map[feoB_map$is_sweep,], aes(
  x = week,
  y = proportion,
  group = interaction(patient_id, chromosome, position),
  color = group,
  shape = mutation_type
)) +
  geom_line()+
  geom_point()+
  scale_color_manual(values = c(
    "UC"     = "#bb9df5",
    "nonIBD" = "#71d987",
    "CD"     = "#ffc169"
  )) +
  facet_wrap(~patient_id+bin)
ggsave(file.path(OUT, "feoB_sweeping.png"), width = 30, height = 30, dpi = 320)


