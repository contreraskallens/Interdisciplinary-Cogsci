require(Matrix)
require(tidyverse)
source('functions.R')

# Load results from data_processing.R ----
years <- read_rds('../saved_objects/journalData.rds') # Saved product from previous
publication.count.matrix <- read_rds('../saved_objects/authorByJournal_id.rds')

# Transform the sparse matrix into a probability matrix now.
# First, transpose so that each column is now the journal count of each author
publication.probability.matrix <- publication.count.matrix %>% t()

# mat <- as.matrix(publication.probability.matrix)

# publication.probability.matrix@x holds the non-zero entries of the matrix.
# publication.probability.matrix@p holds the number of non-zero elements in each column
# Thus, we can make a vector where the sum of each column is repeated x times,
# where x is the number of non-zero entries in that column, and element-wise
# divide the non-zero entries by this vector to get the within-author proportion
# instead of the raw count.

publication.probability.matrix@x <- publication.probability.matrix@x / rep.int(colSums(publication.probability.matrix), diff(publication.probability.matrix@p))
publication.probability.matrix <- t(publication.probability.matrix)

publication.probability.matrix <- irlba::prcomp_irlba(publication.probability.matrix, n = 200)
publication.count.matrix <- irlba::prcomp_irlba(publication.count.matrix, n = 200)

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

# Several articles of Neuropsychologia 2018 have no authors with ids. returned as NULL. 
# Also Anonymous papers. Both these groups are stored in papers.no.ids for reference.
papers.no.ids <- data.for.jsd[which(sapply(authors.per.article, is.null) == TRUE)]


# Calculate distance measures ---------

# Calculate the JSD of every article 
considered.authors <- list() # Save authors whose JSD was calculated

# Loop  through list of authors and list of articles 
# and measures various distances between the publication histories of the authors of an article
# returns a list. NA indicates articles with undefined JSD.
all.results <- map2(authors.per.article, data.for.jsd, function(authors, paper){
  print(authors)
  # Save unique authors of paper (names) for reference
  all.authors <- paper %>% 
    unlist()
  all.authors <- all.authors[str_detect(names(all.authors), "wos_standard")] %>% 
    unique()
  
  # Unlist the authors, which is now a vector where every entry is an ID in the form of a string.
  authors <- unlist(authors)
  
  # If there's only 1 or less than 1 author left after cleanup, skip (undefined JSD)
  if(is.null(authors)){ 
    return(NA)
  }
  if(length(authors) <= 1){ 
    return(NA)
  }
  
  # Filter authors that are not included in the data. 
  # These are the authors that authored papers in our database, 
  # but don't return a personal publication history
  if(FALSE %in% (authors %in% row.names(publication.probability.matrix))){
    remaining.authors <- authors[authors %in% row.names(publication.probability.matrix)]
    missing.authors <- authors[!(authors %in% row.names(publication.probability.matrix))]
    authors <- remaining.authors
  }
  
  # Slice both id x journal matrices using the author ids.
  author.probability.matrix <- publication.probability.matrix[authors, ]
  author.count.matrix <- publication.count.matrix[authors, ]
  
  # If after removing 0s there's only 1 author left, 
  # distribution becomes "double", so filter those out.
  if(typeof(author.probability.matrix) == "double"){ 
    return(NA)
  }

  # Filter authors that have no publications in the previously constructued
  # matrices of past publication history
  author.probability.matrix <- author.probability.matrix[rowSums(author.probability.matrix) > 0, ] 
  author.count.matrix <- author.count.matrix[rowSums(author.count.matrix) > 0, ]
  
  # Filter single-author papers again
  if(typeof(author.probability.matrix) == "double"){ 
    return(NA)
  }
  
  # Compute: Generalized JSD between all authors at the same time,
  # Mean pairwise JSD between authors of a paper
  # Mean pairwise Euclidean distance between the publication count 
  # histories of authors of an article
  # Mean number of articles published previously by the authors of the article
  # Also, Standard Deviation of number of publications
  # And Gini of the number of publications of authors
  
  mean.number.articles <- mean(rowSums(author.count.matrix))
  sd.number.articles <- sd(rowSums(author.count.matrix))
  
  generalized.JSD <- get.JSD(author.probability.matrix)
  
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
                     gini = gini.index,
                    Number.Of.Authors = length(authors),
                     Number.Of.Publications = number.of.publications,
                     Mean.Articles = mean.number.articles, SD.Articles = sd.number.articles,
                     Journal = paper$static_data.summary.titles.title.0, 
                     Considered.Authors = paste("--", paste(authors, collapse = "--"), sep = ""), 
                     Authors = paste("--", paste(all.authors, collapse = "--"), sep = ""),
                     Year = paper$`static_data.summary.pub_info.@attributes.pubyear`,
                     Title = paper$static_data.summary.titles.title.5)
  })

# Remove NA entries (single authors)
all.results <- all.results[which(!(is.na(all.results)))]
complete.frame <- all.results %>% 
  bind_rows()


complete.frame <- complete.frame %>% 
  # Add boundary to calculate the normalized JSD
  mutate(bound = log2(Number.Of.Authors), NJSD = JSD / bound, Journal = factor(Journal)) %>% 
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
write_rds(complete.frame, "../saved_objects/complete_frame_ids.rds")