# this script takes the clean database produced by clean_database and produces a two different versions of the
# database: a stemmed version and a lemmatized version, as lemmatized_catalog.txt and stemmed_catalog.txt#

import nltk
from nltk.tokenize import word_tokenize
from nltk.stem import WordNetLemmatizer

import pandas as pd
lemmat = WordNetLemmatizer()

tag_dictionary = {'CC': 0, 'CD': 1, 'DT': 2, 'EX': 3, 'IN': 4, 'LS': 5, 'MD': 6, 'PDT': 7, 'POS': 8, 'PRP': 9,
                  'PRP$': 10, 'RB': 11, 'RBS': 12, 'RBR': 13, 'RP': 14, 'SYM': 15, 'TO': 16, 'UH': 17, 'WDT': 18,
                  'WP': 19, 'WP$': 20, 'WRB': 21}  # removes everything but adjectives,
whole_text = pd.read_csv("clean_all_abstracts.csv")
abstracts = whole_text["value"]

abstracts = [word_tokenize(abstract) if type(abstract) is str else " " for abstract in abstracts] # every abstract a list of words
abstracts = [nltk.pos_tag(abstract) for abstract in abstracts] # every word a tuple of word, tag
abstracts = [[word_tag[0] for word_tag in abstract if word_tag[1] not in tag_dictionary and len(word_tag[0]) > 2] for abstract in abstracts]
abstracts = [[lemmat.lemmatize(word) for word in abstract] for abstract in abstracts]
abstracts = [' '.join(abstract) for abstract in abstracts]

whole_text["value"] = abstracts
whole_text.to_csv("clean_all_abstracts.csv", index=False)



