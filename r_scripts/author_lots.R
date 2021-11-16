# Get all authors with more than 2 publications
all.authors <- all.author.stats$author %>% unique()
# Loop through authors and save all the charts
for(author in all.authors){
  print(author)
  plot.author(author, save = TRUE)
}
