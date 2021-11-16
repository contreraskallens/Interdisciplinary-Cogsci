library(tidyverse)
library(cowplot)
library(fmsb)
library(psych)

# When re-running everything, uncomment source() lines and comment read_rds() lines
# If not, just load the RDS objects that result from them.

source('functions.R')
# source('JSD.R')
# source('theories.R')

# This is the result of JSD.R
# Gini column includes nan because of the formula
complete.frame <- read_rds('../saved_objects/complete_frame_ids.rds') %>%   
  filter(!(is.nan(gini))) %>% 
  mutate(label = Journal)

# Make a table to fix labels for neater visualization



new.labels <- tibble(Journal = sort(unique(complete.frame$Journal)),
                     Label = c("Annu. Rev. Neurosci.", "Annu. Rev. Psychol.", "Behav. Brain Sci.",
                             "Brain Cognition", "Cognition", "Cogn. Psychol.", "Cogn. Sci.",
                             "J. Cognitive Neurosci.", "J. Mem. Lang.", "Nat. Neurosci",
                             "Neuron", "Neuropsychologia", "Neuropsychology", "Neuropsychol. Rev.", 
                             "Neurosci. Biobehav. R.", "Psychol. Bull.", "Psychol. Rev.",
                             "Top. Cogn. Sci.", "Trends Cogn. Sci", "Trends Neurosci."))
complete.frame <- complete.frame %>% 
  left_join(new.labels)

complete.frame <- complete.frame %>% 
  mutate(cogsci = ifelse(Label %in% c("Cogn. Sci.", "Top. Cogn. Sci."),TRUE, FALSE))

# Bootstrapping Confidence Intervals of measures ---------

# Get the bootstrapped CI of each measure for each journal

# Used for plotting
color.frame <- select(complete.frame, Label, cogsci) %>% distinct

jsd.boot.plot <- complete.frame %>% 
  group_by(Label) %>% 
  do(data.frame(rbind(Hmisc::smean.cl.boot(.$JSD, B = 1000)))) %>% 
  left_join(color.frame) %>% 
  rename(Journal = Label)
njsd.boot.plot <- complete.frame %>% 
  group_by(Label) %>% 
  do(data.frame(rbind(Hmisc::smean.cl.boot(.$NJSD, B = 1000)))) %>% 
  left_join(color.frame) %>% 
  rename(Journal = Label)

# Plot with CI

jsd.boot.plot <- jsd.boot.plot %>% 
  ggplot(aes(x = Journal, y = Mean,
             ymin = Lower, ymax = Upper, color, color = cogsci)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "NJSD") +
  coord_flip() + 
  cowplot::theme_cowplot() + 
  cowplot::background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  labs(title = "Mean JSD by Journal", subtitle = "Bootstrapped Confidence Intervals") +
  scale_color_manual(values = c("black", "red"))
jsd.boot.plot
ggsave(plot = jsd.boot.plot, filename = "../extra_figures/jsd_boot.png", width = 15,
       height = 6, dpi = 300)

njsd.boot.plot <- njsd.boot.plot %>% 
  ggplot(aes(x = Journal, y = Mean,
             ymin = Lower, ymax = Upper, color, color = cogsci)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "NJSD") +
  coord_flip() + 
  cowplot::theme_cowplot() + 
  cowplot::background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  labs(title = "Mean NJSD by Journal", subtitle = "Bootstrapped Confidence Intervals") +
  scale_color_manual(values = c("black", "red"))
njsd.boot.plot
ggsave(plot = njsd.boot.plot, filename = "../extra_figures/njsd_boot.png", width = 15,
       height = 6, dpi = 300)


# Correlation between all measures and covariates.

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

# Regressing out covariates --------

# Covariates: Number of Publications and Gini 
# (Number of authors is used in JSD as it's known they are highly correlated)

color.frame.journal <- select(complete.frame, Journal = Label, cogsci) %>% 
  distinct()

# Because of the insanely long right tail of number of publications nad number of authors, log them
complete.frame <- complete.frame %>%
  mutate(Number.Of.Publications.Log = log(Number.Of.Publications),
         Number.Of.Authors.Log = log(Number.Of.Authors))

jsd.reg <- lm(data = complete.frame, 
              JSD ~ Label + gini + Number.Of.Authors.Log + Number.Of.Publications.Log) 
jsd.reg.2 <- lm(data = complete.frame, 
              JSD ~ Label + gini + Number.Of.Authors.Log) 
jsd.frame <- jsd.reg %>% 
  emmeans::emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>%
  rename(Journal = Label) %>% 
  left_join(color.frame.journal)

jsd.frame.2 <- jsd.reg.2 %>% 
  emmeans::emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>%
  rename(Journal = Label) %>% 
  left_join(color.frame.journal)

njsd.reg <- lm(data = complete.frame, NJSD ~ Label + gini + Number.Of.Publications.Log) 
njsd.reg.2 <- lm(data = complete.frame, NJSD ~ Label + gini) 
summary(njsd.reg)
summary(njsd.reg.2)
# Get semi partial R^2 of Journal
njsd.reg.reduced <- lm(data = complete.frame, NJSD ~ Number.Of.Publications.Log + gini) 
summary(njsd.reg)$adj.r.squared - summary(njsd.reg.reduced)$adj.r.squared
anova(njsd.reg.reduced, njsd.reg)


njsd.frame <- njsd.reg.2 %>% 
  emmeans::emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>% 
  rename(Journal = Label) %>% 
  left_join(color.frame.journal)
njsd.frame.2 <- njsd.reg.2 %>% 
  emmeans::emmeans("Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>% 
  rename(Journal = Label) %>% 
  left_join(color.frame.journal)

# Plot
jsd.plot <- jsd.frame %>% 
  ggplot(aes(x = Journal, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = cogsci)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "JSD") +
  coord_flip() + 
  cowplot::theme_cowplot() + 
  cowplot::background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("black", "red"))
jsd.plot
ggsave(plot = jsd.plot, filename = "../extra_figures/jsd_reg_plot.png", width = 15,
       height = 6, dpi = 300)

njsd.plot <- njsd.frame %>% 
  ggplot(aes(x = Journal, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = cogsci)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "NJSD") +
  coord_flip() + 
  cowplot::theme_cowplot() + 
  cowplot::background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("black", "red"))

njsd.plot
ggsave(plot = njsd.plot, filename = "../extra_figures/njsd_reg_plot.png", width = 15,
       height = 6, dpi = 300)

# Theories per author and journals ---------

# Load result of theories.R
authors.and.abstracts.cogsci <- read_rds("../saved_objects/cogsci_theories.rds")

# For each author with more than 2 publications, gets the mean score for each theory

get.author.stats <- function(author){
  print(author)
  author.df <- authors.and.abstracts.cogsci %>%
    filter_all(any_vars(str_detect(., author))) %>%
    select(bayesian, connectionism, ecological, distributed, embodied, dynamical, enactive, symbolic) %>%
    sweep(1, rowSums(.), FUN = "/")
  if(nrow(author.df) <= 2){return(NA)}

  author.df <- author.df %>%
    colMeans() %>%
    t() %>%
    as_tibble %>%
    add_column(author = author)
  return(author.df)
}

# all.authors <- authors.and.abstracts.cogsci %>%
#   select(contains("name")) %>%
#   unlist() %>%
#   unique() %>%
#   sort()
# all.author.stats <- map(all.authors, get.author.stats)
# all.author.stats %>%
#   write_rds("../saved_objects/allauthorstats.rds")

all.author.stats <- read_rds("../saved_objects/allauthorstats.rds")
all.author.stats <- all.author.stats[!(is.na(all.author.stats))]
all.author.stats <- bind_rows(all.author.stats)

# Get the maximum and minimum score for each theory for plotting
theory.max <- all.author.stats  %>% 
  select(-author) %>% 
  summarize_all(max)
theory.min <- all.author.stats %>% 
  select(-author) %>% 
  summarize_all(min)
# Also get the mean to use as baseline
theory.mean <- all.author.stats %>%
  select(-author) %>%
  summarize_all(mean)

# Normalize as proportion of maximum
all.author.normalized <- all.author.stats %>% 
  mutate(bayesian = bayesian / theory.max$bayesian,
         connectionism = connectionism / theory.max$connectionism,
         ecological = ecological / theory.max$ecological,
         distributed = distributed / theory.max$distributed,
         embodied = embodied / theory.max$embodied,
         dynamical = dynamical / theory.max$dynamical,
         enactive = enactive / theory.max$enactive,
         symbolic = symbolic / theory.max$symbolic)

plot.author <- function(this.author, save = FALSE){
  # Function to plot authors as "power charts" or "radar charts"
  # The vertices of the radar are the maximum of each theory and the
  # Base is the minimum of each theory. A grey polygon is drawn to show the mean
  # score for each theory.
  # If you want to save the file instead of viewing it, save = TRUE
  
  author.wide <- all.author.stats %>% 
    filter(author == this.author) %>% 
    select(-author)
  author.initial <- str_extract(this.author, "[:alpha:]")
  if(save){
    author.dir <- paste0("../author_figures/", author.initial, "/")
    if(!(dir.exists(author.dir))){
      dir.create(author.dir)
    }
    png(paste0("../author_figures/", author.initial, "/", this.author, ".png"), width = 1500, height = 1500,
        res = 300, units = "px")}

  max.min <- bind_rows(theory.max, theory.min)
  author.wide <- bind_rows(max.min, theory.mean, author.wide) 
  author.plot <- fmsb::radarchart(df = author.wide,
             pcol = c("black", "#dc3091"), # Color for border
             pty = c(32, 32), # Symbol of points
             plwd = c(1, 3), # Linewidth of border
             plty = c(3, 1), # Linetype
             pfcol = c(rgb(0.3, 0.3, 0.3, 0.25), # Fill color
                       rgb(0.62, 0.02, 0.36, 0.5)),
             cglty=1, # Linetype for grid
             cglwd = 1, # Line width for grid
             cglcol= "#3c0d76", # Line color for grid
             axistype = 0, # No axis labels
             seg = 4, # Number of segments of the radar
             axislabcol = "#3c0d76", # Color of axis label
             vlcex = 1, # Font size for labels
             title = paste(this.author) # Title
  )
  if(save){
    author.plot
    dev.off()}
}

# Get all authors with more than 2 publications
all.authors <- all.author.stats$author %>% unique()
# Loop through authors and save all the charts
for(author in all.authors){
  print(author)
  plot.author(author, save = TRUE)
}


# Plot journals
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
  select(Journal, contains("Mean"), contains("Lower"), contains("Upper")) %>% 
  group_by() %>% 
  left_join(new.labels) %>% 
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

plot.journal <- function(this.journal){
  journal.df <- journal.boot %>% 
    filter(Journal == this.journal)
  upper <- journal.df %>% 
    select(contains("Upper"))
  colnames(upper) <- theory.levels
  lower <- journal.df %>% 
    select(contains("Lower"))
  colnames(lower) <- theory.levels
  mean <- journal.df %>% 
    select(contains("Mean"))
  colnames(mean) <- theory.levels
  
  all.journal.data <- bind_rows(max.min.journal,
                                upper, mean)
  print(all.journal.data)
  fmsb::radarchart(df = all.journal.data,
             pcol = c("black", "black"),
             pty = c(32, 32),
             plwd = c(1, 3),
             plty = c(3, 1),
             pfcol = c(rgb(0.3, 0.3, 0.3, 0.25),
                       rgb(0.26, 0.45, 0.29, 0.5)),
             cglty=1,
             cglwd = 1,
             cglcol= "#3c0d76",
             axistype = 0,
             seg = 4,
             axislabcol = "#3c0d76",
             vlcex = 1,
             title = paste(this.journal)
  )
}

png("../main_figures/journal_chart.png", width = 10, height = 5, units = 'in', res = 300)
par(mfrow=c(2,4))
for(journal in journal.boot$Journal){
  plot.journal(journal)
}
dev.off()
