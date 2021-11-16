# Load libraries, scripts and resources ----
library(tidyverse)
library(cowplot)
library(psych)
library(emmeans)
source('functions.R')

# This is the result of JSD.R
# Gini column includes nan because of the formula
complete.frame <- read_rds('../saved_objects/complete_frame_ids_cluster.rds') %>%   
  filter(!(is.nan(gini))) %>% 
  mutate(label = Journal)
# Get neater names and categories for Journals
journal.metadata <- read_csv("journal_info.csv")
complete.frame <- left_join(complete.frame, journal.metadata)


# Get covariate-controlled njsd ----

# Because of the insanely long right tail of number of publications and number of authors, log them
complete.frame <- complete.frame %>%
  mutate(Number.Of.Publications.Log = log(Number.Of.Publications),
         Number.Of.Authors.Log = log(Number.Of.Authors))
njsd.reg <- lm(data = complete.frame, NJSD ~ Label + gini + Number.Of.Publications.Log + Number.Of.Authors.Log) 
summary(njsd.reg)

# Get semi partial R^2 of Journal
njsd.reg.reduced <- lm(data = complete.frame, NJSD ~ gini + Number.Of.Publications.Log + Number.Of.Authors.Log) 
summary(njsd.reg)$adj.r.squared - summary(njsd.reg.reduced)$adj.r.squared

# Get model comparison for including Journal factor
anova(njsd.reg.reduced, njsd.reg)
  
# Get estimated mean NJSD that's controlled for variables
njsd.frame <- emmeans::emmeans(njsd.reg, "Label") %>% 
  summary() %>% 
  arrange(desc(emmean)) %>% 
  left_join(journal.metadata)

# Plot njsd 
njsd.plot <- njsd.frame %>% 
  ggplot(aes(x = Label, y = emmean,
             ymin = lower.CL, ymax = upper.CL, color = CogSciSoc)) +
  geom_errorbar(width = 0.5, size = 1) +
  geom_linerange(size = 1) +
  geom_point(size = 2) + 
  labs(y = "NJSD", x = "Journal") +
  coord_flip() + 
  theme_cowplot() + 
  background_grid() + 
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("red", "black")) +
  scale_y_continuous(breaks = seq(0.2, 0.5, 0.01))

# Make hypothetical distributions to illustrate plot

high.dists <- ggplot() +
  geom_area(data = data.frame(values = dbeta(seq(0, 1, length = 1000), 2, 7),
                              index = c(1:1000)), aes(y = values, x = index),
            fill = "seagreen", color = "black", alpha = 0.75) +
  geom_area(data = data.frame(values = dbeta(seq(0, 1, length = 1000), 7, 2),
                              index = c(1:1000)), aes(y = values, x = index),
            fill = "skyblue", color = "black", alpha = 0.75) +
  theme_nothing()

low.dists <- ggplot() +
  geom_area(data = data.frame(values = dbeta(seq(0, 1, length = 1000), 5, 6),
                              index = c(1:1000)), aes(y = values, x = index),
            fill = "seagreen", color = "black", alpha = 0.75) +
  geom_area(data = data.frame(values = dbeta(seq(0, 1, length = 1000), 6, 5),
                              index = c(1:1000)), aes(y = values, x = index),
            fill = "skyblue", color = "black", alpha = 0.75) +
  theme_nothing()
# Make arrows for plot
arrow.right <- ggplot()+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), size = 2,
               arrow = arrow(length = unit(0.2, "inches")), lineend = "round", linejoin = "round") +
  theme_void()
arrow.left <- ggplot()+
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 0), size = 2,
               arrow = arrow(length = unit(0.2, "inches")), lineend = "round", linejoin = "round") +
  theme_void()

njsd.with.dists <- plot_grid(plot_grid(NULL, low.dists, arrow.left, arrow.right, high.dists, rel_widths = c(0.11, 0.22, 0.22, 0.22, 0.22), nrow = 1),
                             njsd.plot, nrow = 2,
                             rel_heights = c(0.1, 0.9))
njsd.with.dists
ggsave(plot = njsd.with.dists, filename = "../main_figures/njsd_with_dists.png", width = 15,
       height = 7, dpi = 300, bg = "white")


# Time analyses ----

complete.frame <- mutate(complete.frame, YearNum = as.numeric(Year))

# Journal by journal
year.reg.journal <- lm(data = complete.frame, NJSD ~ Label + YearNum + Label * YearNum + gini + Number.Of.Authors.Log + Number.Of.Publications.Log)
summary(year.reg)
year.journal.frame <- summary(ref_grid(year.reg.journal, at = list(YearNum = (2009:2018)))) %>% 
  left_join(journal.metadata)
journal.year.plot <- ggplot() + 
  geom_point(data = filter(year.journal.frame, CogSciSoc != "CogSciSoc"), 
             aes(x = YearNum, y = prediction), color = "black", size = 0.5) +
  geom_point(data = filter(year.journal.frame, CogSciSoc == "CogSciSoc"), 
             aes(x = YearNum, y = prediction), color = "red", size = 2) +
  theme_cowplot() +
  labs(y = "NJSD") +
  scale_x_continuous(breaks = c(2009:2018), labels = c(2009:2018))
journal.year.plot

# Plot CogSci time series versus other time series
year.reg.topic <- lm(data = complete.frame, NJSD ~ YearNum + Topic + YearNum * Topic + 
                       Number.Of.Publications.Log + Number.Of.Authors.Log + gini)
year.topic.frame <- summary(ref_grid(year.reg.topic, at = list(YearNum = c(2009:2018))))

topic.year.plot <- ggplot(data = filter(year.topic.frame, Topic != "Other"), 
                          aes(x = YearNum, y = prediction, group = Topic)) +
  geom_line(aes(color = Topic), size = 1) +
  geom_point(aes(fill = Topic), shape = 23, size = 3) +
  theme_cowplot() +
  scale_x_continuous(breaks = c(2009:2018), labels = c(2009:2018)) +
  labs(y = "NJSD")
  
year.grid <- plot_grid(journal.year.plot, topic.year.plot, nrow = 1)
year.grid
ggsave(plot = year.grid, filename = "../main_figures/year_grid.png", width = 15,
       height = 7, dpi = 300, bg = "white")

# Interaction between cogsci-noncogsci and year
type.reg <- lm(formula = NJSD ~ CogSciSoc + gini + Number.Of.Authors.Log +
                 Number.Of.Publications.Log +  YearNum + CogSciSoc * YearNum, 
               data = complete.frame)
summary(type.reg)

# Author plots ----

authors.and.abstracts.cogsci <- read_rds("../saved_objects/cogsci_theories.rds")

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
theory.max <- select(all.author.stats, -author) %>% 
  summarize_all(max)
theory.min <- select(all.author.stats, -author) %>% 
  summarize_all(min)
# Also get the mean to use as baseline
theory.mean <- select(all.author.stats, -author) %>%
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

max.min <- bind_rows(theory.max, theory.min)

png("../main_figures/authors_chart.png", width = 10, height = 5, units = 'in', res = 300, )
par(mfrow = c(1, 2))
author1.wide <- bind_rows(max.min, theory.mean, select(filter(all.author.stats, author == "Anderson, JR"), -author))
fmsb::radarchart(df = author1.wide, pcol = c("black", rgb(0, 0.5, 0.5, 0.9)), # Color for border
                                pty = c(32, 32), # Symbol of points
                                plwd = c(1, 3), # Linewidth of border
                                plty = c(3, 1), # Linetype
                                pfcol = c("grey", rgb(0, 0.5, 0.5, 0.5)),
                                cglty=1, # Linetype for grid
                                cglwd = 1, # Line width for grid
                                cglcol= "#3c0d76", # Line color for grid
                                axistype = 0, # No axis labels
                                seg = 4, # Number of segments of the radar
                                axislabcol = "#3c0d76", # Color of axis label
                                vlcex = 1, # Font size for labels
                                title = "Author 1") # Title
author2.wide <- bind_rows(max.min, theory.mean, select(filter(all.author.stats, author == "Turvey, MT"), -author))
fmsb::radarchart(df = author2.wide, pcol = c("black", rgb(0.9, 0, 0.1, 0.9)), # Color for border
                 pty = c(32, 32), # Symbol of points
                 plwd = c(1, 3), # Linewidth of border
                 plty = c(3, 1), # Linetype
                 pfcol = c("grey", rgb(0.9, 0, 0.1, 0.5)),
                 cglty=1, # Linetype for grid
                 cglwd = 1, # Line width for grid
                 cglcol= "#3c0d76", # Line color for grid
                 axistype = 0, # No axis labels
                 seg = 4, # Number of segments of the radar
                 axislabcol = "#3c0d76", # Color of axis label
                 vlcex = 1, # Font size for labels
                 title = "Author 2") # Title

dev.off()


# Supp figures ----

# Correlation between all measures and covariates. ----

png(filename = "../supp_figures/corr_matrix.png", width = 3000,
    height = 1500, units = "px", res = 300)
complete.frame %>% 
  select(JSD, NJSD, 
         Number.Of.Authors, Number.Of.Publications, gini) %>%
  mutate(Number.Of.Authors = log(Number.Of.Authors)) %>% 
  pairs.panels(method = "pearson",
               density = TRUE,
               ellipses = TRUE,
               ci = TRUE,
               lm = TRUE)
dev.off()


# Other measures ----






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
