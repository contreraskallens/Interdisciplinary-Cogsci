library(Matrix)
library(tidyverse)
library(fastcluster)
library(cluster)
library(cowplot)
source('functions.R')

# Load results from data_processing.R ----
years <- read_rds('../saved_objects/journalData.rds') # Saved product from previous
publication.count.matrix <- read_rds('../saved_objects/authorByJournal_id.rds')

# JSD ----

# Loop through every journal to prepare data for JSD computation.
# Years is a list where every element is a year.
# Each year contains a list for each of the journals in the data.
# Each journal list is cleaned to keep only information relevant to compute jsd.
# Loop returns a list where every element is an article.
data.for.jsd <- years %>%
  map(.x = ., .f = function(yearlist){
    map(.x = yearlist, .f = function(mylist){
      mylist <- clean.for.JSD(mylist)}) %>%
      unlist(recursive = F)
  }) %>%
  unlist(recursive = F)

names(data.for.jsd) <- c()

# For every article in the database, get the authors ids. Returns as list.
authors.per.article <- map(data.for.jsd, get.authors) 

# Eliminate from data.for.jsd and authors.per.article 
# all the articles that have only 1 author
n.authors <- map_dbl(authors.per.article, length)
data.for.jsd <- data.for.jsd[n.authors > 1]
authors.per.article <- authors.per.article[n.authors > 1] 
unique.authors <- authors.per.article %>% 
  unlist() %>% 
  unique()
authors.with.data <- unique.authors %in% row.names(publication.count.matrix)
unique.authors <- unique.authors[authors.with.data] # only authors with data
publication.count.matrix <- publication.count.matrix[unique.authors, ] # Keeps only co-authors with data

# Mean pubs per author and unique authors per journal before cleanup
rowSums(publication.count.matrix) %>% mean()
rowSums(publication.count.matrix) %>% sd()
(publication.count.matrix > 0) %>% colSums() %>% mean()
(publication.count.matrix > 0) %>% colSums() %>% sd()

author.threshold <- 2 # Define a threshold of number of unique authors to keep a journal in the matrix.
journal.count.matrix <- publication.count.matrix > 0
journals.over.threshold <- colSums(journal.count.matrix) > author.threshold
considered.journals <- colnames(publication.count.matrix)[journals.over.threshold]
publication.count.matrix <- publication.count.matrix[, considered.journals]

pub.threshold <- 1 # Define a threshold of number of publications per author to be included.
pubs.per.author <- rowSums(publication.count.matrix)
considered.authors <- unique.authors[pubs.per.author > pub.threshold]
publication.count.matrix <- publication.count.matrix[considered.authors, ]

# Do cleanup again to remove stragglers. This is one is enough.
journal.count.matrix <- publication.count.matrix > 0
journals.over.threshold <- colSums(journal.count.matrix) > author.threshold
considered.journals <- colnames(publication.count.matrix)[journals.over.threshold]
publication.count.matrix <- publication.count.matrix[, considered.journals]

pubs.per.author <- rowSums(publication.count.matrix)
considered.authors <- considered.authors[pubs.per.author > pub.threshold]
publication.count.matrix <- publication.count.matrix[considered.authors, ]

# Now do a second pass through all the articles to see which ones 
# remain after previous cleanup.
co.authored.with.data <- map(data.for.jsd, function(authors){
  article.authors <- get.authors(authors)
  authors.with.data <- article.authors[article.authors %in% considered.authors]
  if(length(authors.with.data) <= 1){
    return(NA)
  }else{
    these.authors.matrix <- publication.count.matrix[authors.with.data, considered.journals]
    these.authors.journals <- colSums(these.authors.matrix)
    these.authors.journals <- colnames(these.authors.matrix)[these.authors.journals > 0]
    return(list("authors" = authors.with.data, "journals" = these.authors.journals))
  }
})
data.for.jsd <- data.for.jsd[!(is.na(co.authored.with.data))]
co.authored.with.data <- co.authored.with.data[!(is.na(co.authored.with.data))]
co.authors.per.article <- map(co.authored.with.data, function(x){return(x$authors)})
considered.authors <- co.authors.per.article %>% 
  unlist() %>% 
  unique()
considered.journals <- map(co.authored.with.data, function(x){return(x$journals)}) %>% 
  unlist() %>% 
  unique()
# This deletes the authors that are co-authors exclusively with people eliminated in the cleanup
publication.count.matrix <- publication.count.matrix[considered.authors, considered.journals]

# Mean pubs per author and unique authors per journal after cleanup
rowSums(publication.count.matrix) %>% mean()
rowSums(publication.count.matrix) %>% sd()
(publication.count.matrix > 0) %>% colSums() %>% mean()
(publication.count.matrix > 0) %>% colSums() %>% sd()


# Build a journal correlation distance matrix for clustering ----
# Beware. The code on this section takes LONG to run. 
# journal.count.matrix <- as.matrix(publication.count.matrix)
# distance.matrix <- bigcor(journal.count.matrix, size = 2000)
# distance.matrix <- 1 - distance.matrix
# dimnames(distance.matrix) <- list(colnames(publication.count.matrix),
#                                   colnames(publication.count.matrix))
# dist.matrix <- as.dist(distance.matrix)
# write_rds(dist.matrix, "../saved_objects/dist_matrix.Rds", compress = "xz")

dist.matrix <- read_rds("../saved_objects/dist_matrix.Rds")
cluster.solution <- fastcluster::hclust(dist.matrix, method = "ward.D2")

silhouette_k = function(k, clusters, distances){
  print(k)
  ## Coerce clusters to hclust, cut with k clusters, and calculate silhouette index
  clusters %>% 
    as.hclust() %>% 
    cutree(k = k) %>% 
    cluster::silhouette(distances) %>% 
    as('matrix') %>% 
    as_tibble()
}

all.silhouettes <- seq(50, 6000, 50) %>% 
  set_names() %>% 
  map_dfr(silhouette_k, cluster.solution, dist.matrix, 
          .id = 'k') %>% 
  mutate(k = as.integer(k))

sil.plot = ggplot(all.silhouettes, 
                  aes(x = k, y = sil_width, 
                      label = k, 
                      group = k)) +
  stat_summary(fun = mean, 
               geom = 'line', 
               group = 1L) +
  stat_summary(fun = mean, geom = 'point', size = 2) +
  labs(x = "Number of clusters", y = "Silhouette") +
  scale_x_continuous(breaks = seq(0, 6100, 500)) +
  scale_y_continuous(breaks = seq(0, 0.08, 0.01)) +
  theme_cowplot() +
  background_grid() +
  expand_limits(x = c(-100, 6100),
                y = 0.08)
sil.plot
ggsave(plot = sil.plot, filename = "../supp_figures/sil_plot.png", width = 15,
       height = 6, dpi = 300, bg = "white")

# Global max seems to be at k = 3000
journal.clusters <- cutree(cluster.solution, k = 3000)
journal.clusters <- data.frame(journal = colnames(publication.count.matrix),
                               cluster = journal.clusters)
cluster.columns <- map(unique(journal.clusters$cluster), function(this.cluster){
  print(this.cluster)
  this.cluster <- filter(journal.clusters, cluster == this.cluster)$journal
  if(length(this.cluster) == 1){
    return(publication.count.matrix[,this.cluster])
  }else{
    cluster.matrix <- publication.count.matrix[,this.cluster]
    return(rowSums(cluster.matrix))}
})
new.columns <- matrix(0, nrow = nrow(publication.count.matrix), ncol = length(unique(journal.clusters$cluster)))
for(i in 1:length(cluster.columns)){
  new.columns[,i] <- cluster.columns[[i]]
}
row.names(new.columns) <- row.names(publication.count.matrix)
complete.publication.count <- publication.count.matrix
publication.count.matrix <- new.columns

# Transform the sparse matrix into a probability matrix now.
# First, transpose so that each column is now the journal count of each author
publication.probability.matrix <- sweep(publication.count.matrix, 1, FUN =  "/", rowSums(publication.count.matrix))
complete.pub.probability <- sweep(complete.publication.count, 1, FUN =  "/", rowSums(complete.publication.count))

# Calculate distance measures ---------

# Calculate the JSD of every article 

# Loop  through list of authors and list of articles 
# and measures various distances between the publication histories of the authors of an article
# returns a list. NA indicates articles with undefined JSD.
all.results <- map2(co.authors.per.article, data.for.jsd, function(authors, paper){
  # print(authors %in% considered.authors)
  authors <- authors[authors %in% row.names(publication.count.matrix)]
  print(authors)
  # Save unique authors of paper (names) for reference
  all.authors <- paper %>% 
    unlist()
  all.authors <- all.authors[str_detect(names(all.authors), "wos_standard")] %>% 
    unique()
  
  # Unlist the authors, which is now a vector where every entry is an ID in the form of a string.
  authors <- unlist(authors)
  # Slice both id x journal matrices using the author ids.
  author.probability.matrix <- publication.probability.matrix[authors, ]
  author.count.matrix <- publication.count.matrix[authors, ]
  complete.author.probability <- complete.pub.probability[authors, ]
  complete.author.count <- complete.publication.count[authors, ]  
  
  # Compute: Generalized JSD between all authors at the same time,
  JSD.cluster <- get.JSD(author.probability.matrix)
  JSD.complete <- get.JSD(complete.author.probability)
  
  # Other measures require pairwise combinations
  unique.combos <- expand.grid(1:length(authors), 1:length(authors)) %>% 
    apply(1, sort) %>% 
    t() %>% 
    unique()
  unique.combos <- unique.combos[unique.combos[,1] != unique.combos[,2],]
  # Next line saves papers with 2 authors
  n.unique <- ifelse(is.null(nrow(unique.combos)), 1, nrow(unique.combos))
  # Mean pairwise JSD between authors of a paper
  if(n.unique == 1){
    pairwise.jsd.cluster <- JSD.cluster
    pairwise.jsd.complete <- JSD.complete
  } else {
    pairwise.jsd.cluster <- map_dbl(1:nrow(unique.combos), function(x){
      these.authors <- unique.combos[x,]
      small.matrix <- author.probability.matrix[these.authors,]
      return(get.JSD(small.matrix))
    })
    pairwise.jsd.complete <- map_dbl(1:nrow(unique.combos), function(x){
      these.authors <- unique.combos[x,]
      small.matrix <- complete.author.probability[these.authors,]
      return(get.JSD(small.matrix))
    })
  }
  
  # Mean pairwise Euclidean distance between the publication count
  # histories of authors of an article
  
  if(n.unique == 1){
    pairwise.euclidean.cluster <- as.numeric(dist(author.count.matrix))
    pairwise.euclidean.complete <- as.numeric(dist(complete.author.count))
  } else {
    pairwise.euclidean.cluster <- map_dbl(1:nrow(unique.combos), function(x){
      these.authors <- unique.combos[x,]
      small.matrix <- author.count.matrix[these.authors,]
      return(as.numeric(dist(small.matrix)))
    })
    pairwise.euclidean.complete <- map_dbl(1:nrow(unique.combos), function(x){
      these.authors <- unique.combos[x,]
      small.matrix <- complete.author.count[these.authors,]
      return(as.numeric(dist(small.matrix)))
    })
  }
  # Mean number of articles published previously by the authors of the article
  # Also, Standard Deviation of number of publications
  # And Gini of the number of publications of authors
  mean.number.articles <- mean(rowSums(author.count.matrix))
  sd.number.articles <- sd(rowSums(author.count.matrix))
  
  # Total number of previous publications for the authors of the article
  number.of.publications <- sum(author.count.matrix)
  
  # Gini index of previous publications
  gini.index <- get.gini(author.count.matrix)
  
  # For each article, return all information:
  # Generalized JSD, mean pairwise JSD, mean euclidean distance, mean number of articles,
  # SD of number of articles, number of authors considered, name of journal, names of the considered authors,
  # the total (considered and not considered) authors of the paper, year of publication, title of the paper,
  # total number of articles of the authors
  this.JSD <- tibble(JSD = JSD.complete,
                     cluster.JSD = JSD.cluster,
                     pairwise.JSD = mean(pairwise.jsd.complete),
                     pairwise.cluster.JSD = mean(pairwise.jsd.cluster),
                     euclidean = mean(pairwise.euclidean.complete),
                     euclidean.cluster = mean(pairwise.euclidean.cluster),
                     gini = gini.index,
                     Number.Of.Authors = length(authors),
                     Number.Of.Publications = number.of.publications,
                     Mean.Articles = mean.number.articles, SD.Articles = sd.number.articles,
                     Journal = paper$static_data.summary.titles.title.0,
                     Considered.Authors = paste("--", paste(authors, collapse = "--"), sep = ""),
                     Authors = paste("--", paste(all.authors, collapse = "--"), sep = ""),
                     Year = paper$`static_data.summary.pub_info.@attributes.pubyear`,
                     Title = paper$static_data.summary.titles.title.5)
  return(this.JSD)
})

# Remove NA entries (single authors)
complete.frame <- all.results %>% 
  bind_rows()
complete.frame <- complete.frame %>% 
  # Add boundary to calculate the normalized JSD
  mutate(bound = log2(Number.Of.Authors), 
         NJSD = JSD / bound, 
         NJSD.cluster = cluster.JSD / bound,
         Journal = factor(Journal)) %>% 
  # Filter articles belonging to an unplanned journal
  filter(Journal != "NEXT-GENERATION OPTICAL NETWORKS FOR DATA CENTERS AND SHORT-REACH LINKS III") %>% 
  mutate(Journal = droplevels(Journal))

# Fix the levels, as some of them include volumes.
levels(complete.frame$Journal) <- c(rep("ANNUAL REVIEW OF NEUROSCIENCE", 10),
                                    rep("ANNUAL REVIEW OF PSYCHOLOGY", 9),
                                    "BEHAVIORAL AND BRAIN SCIENCES", "BRAIN AND COGNITION",
                                    "COGNITION", "COGNITIVE PSYCHOLOGY", "COGNITIVE SCIENCE",
                                    "JOURNAL OF COGNITIVE NEUROSCIENCE", "JOURNAL OF MEMORY AND LANGUAGE",
                                    "NATURE NEUROSCIENCE", "NATURE NEUROSCIENCE", "NEURON",
                                    "NEUROPSYCHOLOGIA", "NEUROPSYCHOLOGY", "NEUROPSYCHOLOGY REVIEW",
                                    "NEUROSCIENCE AND BIOBEHAVIORAL REVIEWS",
                                    "PSYCHOLOGICAL BULLETIN", "PSYCHOLOGICAL REVIEW",
                                    "TOPICS IN COGNITIVE SCIENCE", "TRENDS IN COGNITIVE SCIENCES",
                                    "TRENDS IN NEUROSCIENCES")

write_rds(complete.frame, "../saved_objects/complete_frame.rds")