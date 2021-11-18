# If run by itself, source 'main_script' first
# source('main_script.R')

# Correlation between all measures and covariates. ----

corr.matrix <- complete.frame %>% 
  select(JSD, NJSD, Number.Of.Authors, Number.Of.Publications, gini) %>% 
  mutate(Number.Of.Authors = log(Number.Of.Authors), Number.Of.Publications = log(Number.Of.Publications)) %>% 
  corr.test()
gt(as_tibble(corr.matrix$stars)) %>%
  cols_width(everything() ~ pct(20)) %>% 
  tab_header("p-Value of variable correlation")
gt(as_tibble(corr.matrix$t)) %>% 
  cols_width(everything() ~ pct(20)) %>% 
  fmt_number(everything(), decimals = 3) %>% 
  tab_header("t-Value of variable correlation")

png(filename = "../supp_figures/corr_matrix.png", width = 3000,
    height = 1500, units = "px", res = 300)
complete.frame %>% 
  select(JSD, NJSD,
         Number.Of.Authors, Number.Of.Publications, gini) %>%
  mutate(Number.Of.Authors = log(Number.Of.Authors), Number.Of.Publications = log(Number.Of.Publications)) %>% 
  pairs.panels(method = "pearson",
               density = TRUE,
               ellipses = TRUE,
               ci = TRUE,
               lm = TRUE)
dev.off()


# Other measures ----

# Non normalized JSD
jsd.reg <- lm(data = complete.frame, JSD ~ Label + gini + 
                Number.Of.Publications.Log + Number.Of.Authors.Log) 
jsd.frame <- jsd.reg %>% 
  emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>% 
  rename(Journal = Label) %>% 
  left_join(journal.metadata)
jsd.plot <- jsd.frame %>% 
  ggplot(aes(x = Journal, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = CogSciSoc)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "JSD") +
  coord_flip() + 
  theme_cowplot() + 
  background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(breaks = seq(0.0, 1, 0.05))

jsd.cluster.reg <- lm(data = complete.frame, cluster.JSD ~ Label + gini + 
                        Number.Of.Publications.Log + Number.Of.Authors.Log) 
jsd.cluster.frame <- jsd.cluster.reg %>% 
  emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>% 
  rename(Journal = Label) %>% 
  left_join(journal.metadata)
jsd.cluster.plot <- jsd.cluster.frame %>% 
  ggplot(aes(x = Journal, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = CogSciSoc)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "JSD (Cluster)") +
  coord_flip() + 
  theme_cowplot() + 
  background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(breaks = seq(0.0, 1, 0.05)) +
  theme(axis.text.y = element_blank())

# Normalized JSD based on reduced matrix
njsd.cluster.reg <- lm(data = complete.frame, NJSD.cluster ~ Label + gini + 
                         Number.Of.Publications.Log + Number.Of.Authors.Log) 
njsd.cluster.frame <- njsd.cluster.reg %>% 
  emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>% 
  rename(Journal = Label) %>% 
  left_join(journal.metadata)

njsd.cluster.plot <- njsd.cluster.frame %>% 
  ggplot(aes(x = Journal, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = CogSciSoc)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "NJSD (CLUSTER)") +
  coord_flip() + 
  theme_cowplot() + 
  background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(breaks = seq(0.0, 0.5, 0.01)) +
  theme(axis.text.y = element_blank())

# Mean pairwise JSD

pairwise.reg <- lm(data = complete.frame, 
                   pairwise.JSD ~ Label + gini + 
                     Number.Of.Authors.Log + Number.Of.Publications.Log) 
pairwise.frame <- pairwise.reg %>% 
  emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>%
  rename(Journal = Label) %>% 
  left_join(journal.metadata)
pairwise.plot <- pairwise.frame %>% 
  ggplot(aes(x = Journal, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = CogSciSoc)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "Pairwise JSD") +
  coord_flip() + 
  theme_cowplot() + 
  background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(breaks = seq(0.0, 0.5, 0.01))

pairwise.cluster.reg <- lm(data = complete.frame, 
                           pairwise.cluster.JSD ~ Label + gini + 
                             Number.Of.Authors.Log + Number.Of.Publications.Log) 
pairwise.cluster.frame <- pairwise.cluster.reg %>% 
  emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>%
  rename(Journal = Label) %>% 
  left_join(journal.metadata)

pairwise.cluster.plot <- pairwise.cluster.frame %>% 
  ggplot(aes(x = Journal, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = CogSciSoc)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "Pairwise JSD (CLUSTER)") +
  coord_flip() + 
  theme_cowplot() + 
  background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(breaks = seq(0.0, 0.5, 0.01)) +
  theme(axis.text.y = element_blank())

# Mean pairwise euclidean distance
euclidean.reg <- lm(data = complete.frame, 
                    euclidean ~ Label + gini + 
                      Number.Of.Authors.Log + Number.Of.Publications.Log) 
euclidean.frame <- euclidean.reg %>% 
  emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>%
  rename(Journal = Label) %>% 
  left_join(journal.metadata)

euclidean.plot <- euclidean.frame %>% 
  ggplot(aes(x = Journal, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = CogSciSoc)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "Mean Euclidean") +
  coord_flip() + 
  theme_cowplot() + 
  background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(breaks = seq(10, 25, 1))

euclidean.cluster.reg <- lm(data = complete.frame, 
                            euclidean.cluster ~ Label + gini + 
                              Number.Of.Authors.Log + Number.Of.Publications.Log) 
euclidean.cluster.frame <- euclidean.cluster.reg %>% 
  emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>%
  rename(Journal = Label) %>% 
  left_join(journal.metadata)

euclidean.cluster.plot <- euclidean.cluster.frame %>% 
  ggplot(aes(x = Journal, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = CogSciSoc)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "Mean Euclidean (CLUSTER)") +
  coord_flip() + 
  theme_cowplot() + 
  background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(breaks = seq(20, 30, 1)) +
  theme(axis.text.y = element_blank())

# Make a grid with all 8 measures

measure.grid <- plot_grid(NULL, njsd.plot, njsd.cluster.plot, 
                          NULL, jsd.plot, jsd.cluster.plot,
                          NULL, pairwise.plot, pairwise.cluster.plot,
                          NULL, euclidean.plot, euclidean.cluster.plot,
                          ncol = 3, rel_widths = c(0.1, 1.1, 1), align = "h",
                          labels = c("(A)", "", "", 
                                     "(B)", "", "", 
                                     "(C)", "", "", 
                                     "(D)", "", ""),
                          label_size = 30)
ggsave(plot = measure.grid, filename = "../supp_figures/measure_grid.png", width = 15,
       height = 20, dpi = 300, bg = "white")



# Plot journals -----

theory.levels <- c("bayesian", "connectionism", "distributed", "dynamical",
                   "ecological", "embodied", "enactive", "symbolic")
names(theory.levels) <- theory.levels

# Bootstrap confidence intervals for each theory for the journals
journal.boot <- map_dfc(theory.levels, function(x){
  theory.boot <- authors.and.abstracts.cogsci %>% 
    group_by(Journal) %>% 
    do(data.frame(rbind(Hmisc::smean.cl.boot(unlist(.[,x]), B = 10000))))
  colnames(theory.boot) <- c("Journal", paste0(x, "Mean"), paste0(x, "Lower"), paste0(x, "Upper"))
  return(theory.boot)
}) %>% 
  select(Journal = `Journal...1`, contains("Mean"), 
         contains("Lower"), contains("Upper")) %>% 
  group_by() %>% 
  left_join(journal.metadata) %>% 
  select(-Journal, Journal = Label)

# Use minimum and maximum for the plots
min.journal <- journal.boot %>% 
  select(contains("Upper"), -Journal) %>% 
  summarize_all(min) %>% 
  mutate_all(floor)
max.journal <- journal.boot %>% 
  select(contains("Upper"), -Journal) %>% 
  summarize_all(max)
max.min.journal <- bind_rows(max.journal, min.journal)
colnames(max.min.journal) <- theory.levels
png("../main_figures/journal_chart.png", width = 10, height = 5, units = 'in', res = 300)
par(mfrow=c(2,4))
for(journal in journal.boot$Journal){
  plot.journal(journal)
}
dev.off()

# Plot all other authors if you want ----

# Get all authors with more than 2 publications
# all.authors <- all.author.stats$author %>% unique()
# 
# for(author in all.authors){
#   print(author)
#   plot.author(author, save = TRUE)
# }
