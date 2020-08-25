clean.for.JSD <- function(df){
  # Cleans the dataframe of useless columns to loop more efficiently when calculating JSD.
  df <- df %>%
    # Keep the names in wos_standard format
    select(contains("wos_standard"),
           # Keep the role of each author
           matches("static_data\\.summary\\.names\\.name\\.[[:digit:]]\\.@attributes\\.role"),
           # Keep the year of publication
           contains("static_data.summary.pub_info.@attributes.pubyear"), 
           # Keep the journal name
           contains("static_data.summary.titles.title.0"),
           # Keep the title of the article
           contains("static_data.summary.titles.title.5"),
           # Keep the ids of the authors
           contains("daisng"))
  # Turn into a list for more efficient looping
  df <- pmap(as.list(df), list)
  # Remove NA
  df <- lapply(df, function(x){
    x <- x[!is.na(x)]
  })
  return(df)
}


get.entropy <- function(probability.distributions) {
  # Get the entropy of a probability distribution.
  # Input must be a numeric vector that sums 1
  log.probability.distribution <- probability.distributions * log(probability.distributions)
  log.probability.distribution[is.nan(log.probability.distribution)] <- 0
  entropy <- sum(log.probability.distribution) * -1
  return(entropy)
}


get.JSD <- function(matrix.of.distributions) {
  # Get the Generalized Jensen-Shannon Divergence between the rows of a matrix.
  # The rows of the matrix must each be a probability distribution, 
  # that is, the rows of the matrix must each sum 1.

  # Outputs the Jensen-Shannon Divergence between the rows as a scalar.
  # For simplicity, weights are uniform (1 / length)
  
  weight <- 1 / nrow(matrix.of.distributions)
  joint.distribution <- matrix.of.distributions * weight
  # This sums all the rows of the distribution matrix to get a weighted joint distribution
  joint.distribution.sum <- colSums(joint.distribution)
  # Entropy of the joint distribution
  joint.entropy <- get.entropy(joint.distribution.sum) 
  # Sum of the entropy of each distribution
  individual.entropies <- weight * (apply(matrix.of.distributions, 1, get.entropy))
  entropy.of.individuals <- sum(individual.entropies)  
  JSD <- (joint.entropy - entropy.of.individuals)
  return(JSD)
}

get.authors <- function(article){
  # Receives list with WOS attributes. Returns the unique IDs of all the authors.
  authors <- article %>%
    t() %>%
    as_tibble() %>% 
    select(contains("daisng")) %>% 
    unlist()
  # Safeguard so editors are not included in each article's data
  roles <- article %>% 
    t() %>% 
    as_tibble() %>% 
    select(contains("role")) %>% 
    unlist()
  is.author <- map_lgl(roles, function(x){x == "author"})
  author.places <- str_extract(string = names(roles)[is.author], pattern = "[[:digit:]]")
  # Extract the IDs of only author places
  authors <- authors[str_detect(names(authors), paste(author.places, collapse = "|"))] %>% 
    .[!(is.na(.))] %>% 
    unique()
  return(authors)
}

get.gini <- function(publication.count.matrix){
  # Calculate gini for each article.
  # using https://www.statsdirect.com/help/nonparametric_methods/gini_coefficient.htm
  
  # Get total number of publications per author, excluding those with 0 publications.
  number.of.publications <- rowSums(publication.count.matrix) %>% 
    sort()
  number.of.publications <- number.of.publications[number.of.publications > 0]
  
  # Get the ordered difference from the mean number of papers
  mean.number.publications <- mean(number.of.publications)
  number.authors <- length(number.of.publications)
  differences <- number.of.publications - mean.number.publications
  
  # Multiply differences from the mean by their rank
  differences <- differences * c(1:number.authors)
  differences <- sum(differences)
  numerator <- differences * 2
  denominator <- (number.authors ^ 2) * mean.number.publications
  return(numerator / denominator)
  
}  


entropy.weight <- function(matrix) {
  # Weights a document by term matrix by the entropy of the terms for LSA
  
  # Transposes for division
  temp.matrix <- t(matrix) 
  # Generates a probability matrix
  a <- temp.matrix / rowSums(temp.matrix) 
  # Calculate entropy of terms 
  # (1 + (sum of probability times log2 probability, divided by the total number of documents)).
  b <- 1 + ((rowSums(a * log2(a), na.rm = T)) / log2(dim(matrix)[1])) 
  # Log normalizes frequency matrix and deletes 0 entropy ones.
  temp.matrix <- log(temp.matrix[which(b > 0, arr.ind = T), ] + 1) 
  # Weight log normalized matrix by multiplying terms with 
  # entropy higher than 0 by its entropy.
  temp.matrix <- temp.matrix * b[b > 0] 
  # Returns non-transposed matrix.
  return(t(temp.matrix)) 
}

