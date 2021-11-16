library(tidyverse)
library(viridis)
library(parallelDist)
library(dendextend)
library(far)
library(ggwordcloud)
source('functions.R')
set.seed(123)

# Load and pre-process data -------


# Read document by term matrix created using python. This are the documents used in the Scientometrics paper
theories.dbt <- read_csv('../theories_data/dbt_matrix_scientometrics.csv')
theories.dbt[is.na(theories.dbt)] <- 0

# Get rid of words that appear 5 or fewer times
theories.dbt <- theories.dbt[, which(colSums(theories.dbt) > 5)] 
# Load the catalog that includes information on every paper
theories.catalog <- read_csv('../theories_data/catalog_scientometrics.csv') %>%
  select(-id)

# To select the journals in new database to be used, look at the journals that appeared in the scientometrics dataset
# Only use journals that have more than 1% of the papers with which the semantic space was constructed
journals.both <- theories.catalog$journal %>% 
  table() %>% 
  as_tibble %>%
  mutate(proportion = n / sum(n)) %>% 
  filter(proportion > .01)
journals.both <- journals.both$.

# Weight the terms using the entropy
theories.dbt <- entropy.weight(theories.dbt)

# Replace NA with 0.
theories.dbt[is.na(theories.dbt)] <- 0 

# Singular Value Decomposition
reduced.theories.dbt <- svd(theories.dbt)
theory.list <- unique(theories.catalog$theory)
names(theory.list) <- theory.list

# D  is the number of dimensions used for the analysis
D = 100
  
# Create document loadings by multiplying the U matrix of the SVD 
# by a diagonilized version of the singular values
document.loadings <- reduced.theories.dbt$u[,1:D] %*% diag(reduced.theories.dbt$d)[1:D, 1:D] %>% 
  scale() # Scaling has shown good results
term.loadings <- reduced.theories.dbt$v[,1:D] %*% diag(reduced.theories.dbt$d)[1:D, 1:D]
row.names(term.loadings) <- colnames(theories.dbt)
row.names(document.loadings) <- row.names(theories.dbt)

# Add the theory of each document to the document loadings
document.loadings <- document.loadings %>%
  as_tibble() %>%
  mutate(theory = as.factor(theories.catalog$theory))



# Observe theory structure --------

# Get cosine distance matrix for quick theory structure
distance.matrix <- document.loadings %>%
  select(1:100, -theory) %>% # Reduce space as we're not gonna use all 100 dimensions for classification
  as.matrix %>%
  parDist(method = "cosine") %>%
  as.matrix()

# Dendrogram and heatmap are big picture of theory structure

# Get mean distances to every theory
distances <- map_dfc(theory.list, function(theory){
  # Returns a DF where each row is a document and each column is a theory.
  # Each cell is the mean distance of document x to documents of theory y
  theory.mask <- theories.catalog$theory == theory
  theory.distances <- distance.matrix[,theory.mask]
  theory.distances <- rowMeans(theory.distances, na.rm = T)
  return(theory.distances)
})

distances$theory <- theories.catalog$theory
mean.distances <- distances %>%
  group_by(theory) %>%
  summarise_all(mean, na.rm = TRUE)
mean.distances %>% 
  gather("theory2", "value", bayesian:symbolic) %>% 
  ggplot(aes(x = theory, y = theory2, fill = value, label = round(value, digits = 3))) +
  geom_tile() +
  scale_fill_viridis() +
  geom_text()

dendrogram <- hclust(dist(mean.distances, upper = T), method = "ward.D2")
dendrogram$labels <- theory.list
par(cex = 1.5, srt = 30)
png(filename = "../main_figures/dendrogram.png", res = 300, width = 4000,
    height = 1500, units = "px")
dendrogram %>% 
  as.dendrogram(hang = 0.1) %>% 
  set("labels_cex", 1.5) %>% 
  as.ggdend() %>% 
  ggplot(horiz = TRUE) +
  theme(plot.margin = margin(0, 10, 0, 0, "pt"))
dev.off()

# Train GLMs --------


# Train one glm for theory with dimensions that best predict that theory.
# Same methodology as scientometrics paper.

theories <- unique(document.loadings$theory)
names(theories) <- theories

all.dims <- map(theories, function(this.theory){
  print(this.theory)
  # Turn theory into a binary field
  theory.papers <- document.loadings %>%
    mutate(is.this.theory = as.factor((theory == this.theory))) %>%
    select(-theory)
  # Extract coefficients of glm and sort them in decreasing absolute value
  glm <- glm(data = theory.papers, is.this.theory ~ ., family = "binomial")
  # Return the highest 30 coefficients per theory among the 100 dimensions
  dims <- glm$coefficients %>% 
    .[2:length(.)] %>% 
    abs %>% 
    sort(decreasing = TRUE) %>% 
    .[1:30] %>% 
    names()
  return(dims)
})

all.glm <- map(theories, function(this.theory){
  # Train a new set of glm using only the most predictive 30 dimensions 
  # picked in the previous step
  print(this.theory)
  theory.papers <- document.loadings %>%
    mutate(is.this.theory = as.factor((theory == this.theory))) %>%
    select(-theory)
  relevant.dimensions <- all.dims[[this.theory]]
  theory.papers <- theory.papers[,c(relevant.dimensions, "is.this.theory")]
  glm <- glm(data = theory.papers, is.this.theory ~ ., family = "binomial")
})

# Load the document by term matrix of the new papers
# Takes quite some time to read, so you might want to uncompress and then load the .csv
# topics.dbt <- read_csv('../theories_data/all_dbt.zip') %>%
#   select(-1) %>%
#   as.matrix()
# 
# Determine the shared terms between the Scientometrics dbt and the new dbt
# shared.terms <- colnames(topics.dbt)[colnames(topics.dbt) %in% colnames(theories.dbt)]
# topics.dbt.reduced <- topics.dbt[,shared.terms] %>%
#   as.matrix()
# rm(topics.dbt)

# Save reduced dbt to avoid loading entire thing every time
# topics.dbt.reduced %>% 
#   as_data_frame() %>% 
#   write_csv("all_dbt_reduced.csv")

# Load product of previous steps
topics.dbt.reduced <- read_csv("../theories_data/all_dbt_reduced.zip") %>% 
  as.matrix()

# Fold-in new semantic space -------

# Add the terms that are in Contreras Kallens & Dale (2018) 
# matrix and not in the Topics matrix
other.terms <- matrix(0,
                      ncol = sum(!(colnames(theories.dbt) %in% colnames(topics.dbt.reduced))),
                      nrow = nrow(topics.dbt.reduced))
colnames(other.terms) <- colnames(theories.dbt)[!(colnames(theories.dbt) %in% colnames(topics.dbt.reduced))]
topics.dbt.reduced <- cbind(topics.dbt.reduced, other.terms)

# Get them in the same order
topics.dbt.reduced <- topics.dbt.reduced[,colnames(theories.dbt)]

# Fold-in new journal dbt onto LSA space
projected.topics <- topics.dbt.reduced %*% term.loadings
colnames(projected.topics) <- colnames(document.loadings)[1:D]

# Load new abstracts -----

# Get frame with abstracts. 
# Filter out papers without abstracts. 
# Return a dataframe with authors, abstracts, title and year

years <- read_rds('../saved_objects/journalData.rds')

authors.and.abstracts <- years %>% 
  map(function(y){
    map(y, function(x){
      x <- x %>% 
        select(contains("wos_standard"), contains("abstract"), contains("title"), contains("journal"), contains("year")) %>% 
        filter(`static_data.summary.pub_info.@attributes.has_abstract` == "Y")
      return(x)}) %>% 
      reduce(bind_rows)
  }) %>% 
  reduce(bind_rows) %>% 
  mutate(abstract = ifelse(!(is.na(static_data.fullrecord_metadata.abstracts.abstract.abstract_text.p)),
                                                      static_data.fullrecord_metadata.abstracts.abstract.abstract_text.p, 
                           static_data.fullrecord_metadata.abstracts.abstract.abstract_text.p.0))

# Continue with journal list from line 25

authors.and.abstracts.cogsci <- filter(authors.and.abstracts, 
                                       static_data.summary.titles.title.0 %in% journals.both)

# Get only those papers from the folded-in matrix
projected.topics <- projected.topics %>% 
  .[which(authors.and.abstracts$static_data.summary.titles.title.0 %in% journals.both),] %>% 
  scale() # Scaling seems to work well, and the GLMs were trained on scaled data


# Predict theories using GLM -------

# Same column names as document.loadings to use relevant dims
colnames(projected.topics) <- colnames(document.loadings)[1:100] 

predictions <- map_dfc(theories, function(theory){
  relevant.dimensions <- all.dims[[theory]]
  relevant.lsa <- projected.topics[,relevant.dimensions] %>% 
    as_tibble()
  model <- all.glm[[theory]]
  predicted <- predict(model, newdata = relevant.lsa,
                       type = "response")
})

authors.and.abstracts.cogsci <- bind_cols(authors.and.abstracts.cogsci,
                                          predictions)
authors.and.abstracts.cogsci <- authors.and.abstracts.cogsci %>%
  rename(Journal = static_data.summary.titles.title.0,
         Year = `static_data.summary.pub_info.@attributes.pubyear`)

authors.and.abstracts.cogsci %>% 
  write_rds("../saved_objects/cogsci_theories.rds", compress = "xz")

# Generate the wordclouds ---------

# Get a varimax rotation of the term loadings.
# Technical details can be gotten from Latent semantic analysis - Nicholas E. Evangelopoulos (2013)
terms.varimax <- varimax(term.loadings)
varimax.terms <- unclass(terms.varimax$loadings)

# Rotato document loadings by the same rotation matrix
varimax.documents <- as.matrix(document.loadings[,1:100]) %*% terms.varimax$rotmat %>% 
  as_tibble()
colnames(varimax.documents) <- colnames(document.loadings)[1:100]
varimax.documents$theory <- document.loadings$theory
colnames(varimax.terms) <- colnames(document.loadings)[1:100]

# Get dimensions for each theory in the rotated space
# Divide them into positive or negative depending on coefficient
# Get the 30 higher ones
dims.for.wc <- map(theories, function(this.theory){
  print(this.theory)
  theory.papers <- varimax.documents %>%
    mutate(is.this.theory = as.factor((theory == this.theory))) %>%
    select(-theory)
  glm <- glm(data = theory.papers, is.this.theory ~ ., family = "binomial")
  order.of.dims <- glm$coefficients %>% 
    .[2:length(.)] %>% 
    abs %>% 
    sort(decreasing = TRUE) %>% 
    .[1:10] %>% 
    names()
  dims <- glm$coefficients[order.of.dims]
  positive <- names(dims)[dims > 0]
  negative <- names(dims)[dims < 0]
  return(list(positive = positive, negative = negative))
})


# For each theory, compute a positive and a negative matrix with
# Dimension x Term. Then, take the (scaled) mean of those terms in the dimensions
words.valences <- map(theories, function(theory){
  positive.valences <- dims.for.wc[[theory]][["positive"]]
  negative.valences <- dims.for.wc[[theory]][["negative"]]
  
  positive.matrix <- varimax.terms[,positive.valences] %>% 
    rowMeans %>% 
    sort(decreasing = TRUE)
  negative.matrix <- varimax.terms[,negative.valences] %>% 
    rowMeans %>% 
    sort(decreasing = TRUE)
  positive.terms <-  names(positive.matrix)
  positive.values <-  positive.matrix
  
  negative.terms <- names(negative.matrix)
  negative.values <- negative.matrix
  
  return(list(positive = tibble(terms = positive.terms, values = positive.values),
              negative = tibble(terms = negative.terms, values = negative.values)))
})


# For each theory, generate a positive and a negative word cloud with the 
# 20 words with the highest mean scores in their respective dimensions.
wordclouds <- map(theories, function(theory){
  positive <- words.valences[[theory]][["positive"]]
  negative <- words.valences[[theory]][["negative"]]
  positive <- ggplot(positive[1:25,], aes(label = terms, size = values)) + 
    geom_text_wordcloud(rm_outside = TRUE, shape = "square") +
    scale_size_area(aes(color = values), max_size = 14) + 
    scale_color_viridis() + 
    theme_cowplot() + 
    theme(axis.line  = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  negative <- ggplot(negative[1:20,], aes(label = terms, size = values)) + 
    geom_text_wordcloud(rm_outside = TRUE, shape = "square") +
    scale_size_area(aes(color = values), max_size = 14) + 
    scale_color_viridis() + 
    theme_cowplot() + 
    theme(axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  ggsave(positive, filename = paste0("../supp_figures/word_clouds/positive_wordcloud_", theory, ".png"), scale = 1.5, device = "png")
  ggsave(negative, filename = paste0("../supp_figures/word_clouds/negative_wordcloud_", theory, ".png"), scale = 1.5, device = "png")
  return(list(positive = positive, negative = negative))
})


# Save them as a grid
plot_grid(wordclouds$bayesian$positive, wordclouds$bayesian$negative, 
          wordclouds$connectionism$positive, wordclouds$connectionism$negative,
          wordclouds$distributed$positive, wordclouds$distributed$negative,
          wordclouds$dynamical$positive, wordclouds$dynamical$negative,
          wordclouds$ecological$positive, wordclouds$ecological$negative,
          wordclouds$embodied$positive, wordclouds$embodied$negative,
          wordclouds$enactive$positive, wordclouds$enactive$negative,
          wordclouds$symbolic$positive, wordclouds$symbolic$negative,
          ncol = 2, align = "h", labels = c("Bayesian - Positive", "Bayesian - Negative", "Connectionism - Positive", "Connectionism - Negative",
                                            "Distributed - Positive", "Distributed - Negative", "Dynamical - Positive", "Dynamical - Negative", 
                                            "Ecological - Positive", "Ecological - Negative", "Embodied - Positive", "Embodied - Negative",
                                            "Enactive - Positive", "Enactive - Negative", "Symbolic - Positive", "Symbolic - Negative"), 
          label_size = 30)
ggsave("../supp_figures/wordcloud_panel.png", dpi = 600, width = 20, height = 25, units = "in")
