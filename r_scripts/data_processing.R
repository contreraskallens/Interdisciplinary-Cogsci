# IMPORTANT: This script requires uncompressing the data.7z archive in the data folder.
# All writing of objects is commented out. Uncomment if you want to save again.
# Packages ------------

library(ndjson)
library(tidyverse)
library(Matrix)
library(data.table)



# Load journal data -------

# Load all years of 20 journals and put them in a list with year as
# name of every element. Each journal is a dataframe with all the articles.
# Because data is saved as JSON, flatten before.
year.paths <- list.files('../data/journals/', recursive = T, full.names = T, pattern = ".json")
years <- list()

for(file in year.paths){
  journal.name <- str_extract(file, "(?<=[:digit:]\\.).+(?=.json)")
  year <- str_extract(file, "[:digit:]{4}")
  journal.df <- stream_in(file) %>%
    unlist() %>%
    ndjson::flatten() %>%
    as_tibble()
  years[[year]][[journal.name]] <- journal.df
}

# write_rds(years, '../saved_objects/journalData.rds', compress = "xz", compression = 9)


# Load author data --------


# Each author's JSON was previously transformed into a CSV. 
# Load all of them and put them in an author list.
# The resulting list has author with reduced number of columns.

author.paths <- list.files('../data/authors/', full.names = T, pattern = ".csv")
authors <- list()

for (file in author.paths){
  # Skip empty files
  if(file.info(file)$size == 0){
    next
  }  
  author.name <- str_extract(file, "(?<=authors/).+(?=.csv)")
  
  # This line allows to resume process if for some reason 
  # it was interrupted
  if(author.name %in% names(authors)){
    next
  }
  print(author.name)
  # Remove useless columns: full name, citation count, separate versions of the name, dates, 
  # contact information
  authors[[author.name]] <- fread(file) %>%
    select(-contains("full_name"), -contains("count"), -contains("last_name"), -contains("first_name"),
           -contains("display_name"), -contains("reprint"), -contains("seq_no"), -contains("role"),
           -contains("date"), -contains("pubmonth"), -contains("pubyear"), -contains("email"), -contains("page"),
           -contains("addr"), -contains("suffix"))
}

write_rds(authors, '../saved_objects/allAuthors.rds', compress = "xz", compression = 9)


# Obtain IDs of each author ---------


# Loop through list of authors to fetch id and journals
authors.ids <- list()
for(author in new.authors){
  # Print author to monitor
  print(author)
  author.data <- authors[[author]]
  # Match by last name to avoid problems with differing/missing middle names
  author.name <- str_extract(author, ".*(?=,)") 
  # For every column of the author's dataframe,
  # get the position the author has in that paper (integer)
  # Then get the WoS id of that position.
  author.journals <- map(.x = c(1:nrow(author.data)),
                       .f = function(x){
                         this.row <- author.data[x,]
                         # Reduce df to have only names
                         only.names <- select(this.row, contains("wos_standard"))
                         # Get name of journal
                         this.journal <- this.row$static_data.summary.titles.title.0
                         # Get columns that contain the author's last name
                         these.places <- which(str_detect(only.names, author.name)) 
                         # Get the places that the author appear in each row by extracting
                         # the digit of columns that match the name
                         these.places <- names(only.names)[these.places] %>% 
                           str_extract("[[:digit:]]")
                         # For each of the places where the author appears, extract the id by pasting
                         # the extracted digit in the column namefor each id in the paper
                         these.ids <- map_chr(these.places, function(place){
                           # When author is the only author, "these.places" contains an NA,
                           # And the author's id is in the static_data.summary.names.name.@attributes.daisng_id
                           # column instead of a numbered column.
                           if(!(is.na(place))){
                             id.variable <- paste0("static_data.summary.names.name.", place, ".@attributes.daisng_id")
                           } else{id.variable <- "static_data.summary.names.name.@attributes.daisng_id"}
                           return(unlist(author.data[x, ..id.variable]))
                         })
                         # Return a list where [[ids]] contains the ids in the papers, and [[journal]] contains
                         # the journals where those ids appear.
                         ids.and.journals <- list(ids = these.ids, journal = this.journal)
                         return(ids.and.journals)
                          })
  # Author.journals now contains a list where every element is the
  # ids.and.journals variable of each row in the author's dataframe
  for(id.journal in author.journals){
    # Loop through each id and each associated journals. If the id is not in the
    # aggregated list, Add an entry to the overall list with id as a name and a
    # vector of journal names as the content. If it's already there, append this
    # vector of journal names to the vector already contained.
    for(author.id in id.journal$ids){
      if(author.id %in% names(authors.ids)){
        authors.ids[[author.id]] <- c(authors.ids[[author.id]], id.journal$journal)
      } else{authors.ids[[author.id]] <- id.journal$journal}
    }
    
  }
}

# write_rds(author.ids, "../saved_objects/authors_ids.rds")


# Contruct ID by Journal matrix ----------

# Get all the unique journals that appear
all.journals <- unique(unlist(authors.ids)) %>% sort()

# Fill the id x journal matrix with 0 and then loop through authors to
# fill their corresponding row. IDs are rows, journals are columns.
# Each column is the count of past publication in that journal.
get.journal.count <- function(author, author.journal.data = author.journals, journal.data = all.journals){
  print(author)
  author.journal.data <- author.journal.data[[author]]
  all.pubs <- vector("numeric", length = length(journal.data))
  author.pubs <- match(author.journal.data, journal.data) %>%
    table()
  all.pubs[as.numeric(names(author.pubs))] <- author.pubs
  return(all.pubs)
}

# Saved as sparse matrix as it's pretty big
journal.matrix <- lapply(names(authors.ids), get.journal.count, author.journal.data = authors.ids)
journal.matrix <- journal.matrix %>% 
  unlist() %>%
  Matrix(sparse = T, data = .,
         ncol = length(authors.ids))
colnames(journal.matrix) <- names(authors.ids)
row.names(journal.matrix) <- sort(all.journals)

# Transpose to get author x journal instead.
journal.matrix <- journal.matrix %>%
  t()

write_rds(journal.matrix, "../saved_objects/authorByJournal_id.rds", compress = "xz")