library(Matrix)
library(tidyverse)
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
n.authors <- map_dbl(authors.per.article, length)
data.for.jsd <- data.for.jsd[n.authors > 1]
co.authored.articles <- authors.per.article[n.authors > 1] 
short.authors.per.article <- co.authored.articles %>% 
  unlist() %>% 
  unique()
authors.with.data <- short.authors.per.article %in% row.names(publication.count.matrix)
short.authors.per.article <- short.authors.per.article[authors.with.data]
pub.threshold <- 1
pub.per.author <- rowSums(publication.count.matrix[short.authors.per.article,])
considered.authors <- short.authors.per.article[pub.per.author > pub.threshold]

author.threshold <- 2
journal.count.matrix <- publication.count.matrix > 0
author.per.journal.count <- colSums(journal.count.matrix) > author.threshold
considered.journals <- colnames(publication.count.matrix)[author.per.journal.count]


# Now go over the articles to see what authors/journals are still considered in the study
# Eliminates authors that co-author only with the ones that were eliminated in previous steps
co.authored.with.data <- map(co.authored.articles, function(authors){
  authors.with.data <- authors[authors %in% considered.authors]
  if(length(authors.with.data) <= 1){
    print(authors.with.data)
    return(NA)
  }else{
    these.authors.matrix <- publication.count.matrix[authors.with.data,considered.journals]
    authors.with.data <- authors.with.data[rowSums(these.authors.matrix) > pub.threshold]
    if(length(authors.with.data) <= 1){
      print(authors.with.data)
      return(NA)
    }
    these.authors.journals <- colSums(these.authors.matrix)
    these.authors.journals <- colnames(these.authors.matrix)[these.authors.journals > 0]
    return(list("authors" = authors.with.data, "journals" = these.authors.journals))
  }
})
data.for.jsd <- data.for.jsd[!(is.na(co.authored.with.data))]
co.authored.with.data <- co.authored.with.data[!(is.na(co.authored.with.data))]

co.authors <- map(co.authored.with.data, function(x){return(x$authors)}) %>% 
  unlist() %>% 
  unique()
journals <- map(co.authored.with.data, function(x){return(x$journals)}) %>% 
  unlist() %>% 
  unique()

publication.count.matrix <- publication.count.matrix[co.authors, journals]
publication.count.matrix <- publication.count.matrix[,colSums(publication.count.matrix) > 0]
journal.count.matrix <- t(publication.count.matrix) %>% 
  as.matrix()

# journal.pca <- irlba::prcomp_irlba(journal.count.matrix, n = 300)

# plot(1:300, journal.pca$sdev)
# journal.count.matrix <- journal.pca$x[,1:150]

# Euclidean distance and dot product: 
# sum(a_i^2) + sum(a_i^2) - 2(sum(a_i * b_i)) = sum(a - b)^2

# lengths <- rowSums(journal.count.matrix^2)
# dot.product <- tcrossprod(journal.count.matrix)
# dot.product <- dot.product * -2
# eu.dist <- sweep(dot.product, 1, lengths, "+") # plus length 1
# eu.dist <- sweep(eu.dist, 2, lengths, "+") # plus length 2
# eu.dist <- sqrt(eu.dist)
# diag(eu.dist) <- 0
# 
distance.matrix <- matrix(data = NA,
                          ncol = ncol(publication.count.matrix),
                          nrow = ncol(publication.count.matrix))
buffer.size <- 500
t.journal.count <- t(journal.count.matrix)
for(i in 1:ceiling(nrow(journal.count.matrix) / buffer.size)){
  start <- (i * buffer.size - buffer.size) + 1
  end <- i * buffer.size
  if(end > nrow(journal.count.matrix)){end <- nrow(journal.count.matrix)}
  this.cor <- cor(t.journal.count[,start:end], t.journal.count, method = "spearman")
  distance.matrix[start:end,] <- this.cor
  print(i)
  gc()
}

distance.matrix <- 1 - distance.matrix
dimnames(distance.matrix) <- list(colnames(publication.count.matrix),
                                  colnames(publication.count.matrix))

# write.csv(distance.matrix, "distance_matrix.csv")
# distance.matrix <- read.csv("distance_matrix.csv")[,-1]

dist.matrix <- as.dist(distance.matrix)
write_rds(dist.matrix, "dist_matrix.Rds", compress = "xz")

dist.matrix <- read_rds("dist_matrix.Rds")

cluster.solution <- cluster::agnes(dist.matrix, 
                                   diss = TRUE, 
                                   keep.diss = F, 
                                   keep.data = F,
                                   method = "ward")

write_rds(cluster.solution, "cluster_solution.Rds")
cluster.solution <- read_rds("cluster_solution.Rds")

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

all.silhouettes <- c(10, seq(50, 6000, 50)) %>% 
  set_names() %>% 
  map_dfr(silhouette_k, cluster.solution, dist.matrix, 
          .id = 'k') %>% 
  mutate(k = as.integer(k))

sil.plot <- ggplot(all.silhouettes, 
                  aes(x = k, y = sil_width, 
                      label = k, 
                      group = k)) +
  stat_summary(fun = mean, 
               geom = 'line', 
               group = 1L) +
  stat_summary(fun = mean, geom = 'label', size = 5) +
  # ggbeeswarm::geom_beeswarm(alpha = .1)
  labs(x = "Number of clusters", y = "Silhouette")
sil.plot

write_rds(all.silhouettes, "all_silhouettes.Rds")

journal.clusters <- cutree(cluster.solution, k = 3500)
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
old.publication.count <- publication.count.matrix
publication.count.matrix <- new.columns

# Transform the sparse matrix into a probability matrix now.
# First, transpose so that each column is now the journal count of each author
publication.probability.matrix <- sweep(publication.count.matrix, 1, FUN =  "/", rowSums(publication.count.matrix))
old.pub.probability <- sweep(old.publication.count, 1, FUN =  "/", rowSums(old.publication.count))
# Several articles of Neuropsychologia 2018 have no authors with ids. returned as NULL. 
# Also Anonymous papers. Both these groups are stored in papers.no.ids for reference.
papers.no.ids <- data.for.jsd[which(sapply(authors.per.article, is.null) == TRUE)]
data.for.jsd <- data.for.jsd[setdiff(1:length(data.for.jsd), which(sapply(authors.per.article, is.null) == TRUE))]
# Calculate distance measures ---------

# Calculate the JSD of every article 

authors.per.article <- map(data.for.jsd, get.authors) 

data.for.jsd <- data.for.jsd[setdiff(1:length(data.for.jsd), which(sapply(authors.per.article, is.null) == TRUE))]

# Loop  through list of authors and list of articles 
# and measures various distances between the publication histories of the authors of an article
# returns a list. NA indicates articles with undefined JSD.
all.results <- map2(authors.per.article, data.for.jsd, function(authors, paper){
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
  old.author.probability <- old.pub.probability[authors, ]
  old.author.count <- old.publication.count[authors, ]  # Compute: Generalized JSD between all authors at the same time,
  # Mean pairwise JSD between authors of a paper
  # Mean pairwise Euclidean distance between the publication count
  # histories of authors of an article
  # Mean number of articles published previously by the authors of the article
  # Also, Standard Deviation of number of publications
  # And Gini of the number of publications of authors

  mean.number.articles <- mean(rowSums(author.count.matrix))
  sd.number.articles <- sd(rowSums(author.count.matrix))

  generalized.JSD <- get.JSD(author.probability.matrix)
  old.jsd <- get.JSD(old.author.probability)
  # Total number of previous publications for the authors of the article
  number.of.publications <- sum(author.count.matrix)

  # Gini index of previous publications
  gini.index <- get.gini(author.count.matrix)

  # For each article, return all information:
  # Generalized JSD, mean pairwise JSD, mean euclidean distance, mean number of articles,
  # SD of number of articles, number of authors considered, name of journal, names of the considered authors,
  # the total (considered and not considered) authors of the paper, year of publication, title of the paper,
  # total number of articles of the authors
  this.JSD <- tibble(JSD = generalized.JSD,
                     old.JSD = old.jsd,
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
all.results <- all.results[which(!(is.na(all.results)))]
complete.frame <- all.results %>% 
  bind_rows()

complete.frame <- complete.frame %>% 
  # Add boundary to calculate the normalized JSD
  mutate(bound = log2(Number.Of.Authors), 
         NJSD = JSD / bound, 
         old.NJSD = old.JSD / bound,
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

write_rds(complete.frame, "../saved_objects/complete_frame_ids_cluster.rds")
