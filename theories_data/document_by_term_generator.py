# this script creates a document by (unique) term of the catalog specified as database (default:
# 'lemmatized_catalog.txt'). Writes it down as 'document_by_term.txt'#

import re
import pandas as pd
import numpy as np

def flatten(mylist):  # flatten one level deep list
    new_list = []
    for item in mylist:
        n = 0
        while n < len(item):
            new_list.append(item[n])
            n = n + 1
    return new_list


abstracts = pd.read_csv('clean_all_abstracts.csv')
word_bag = [abstract.split(" ") if type(abstract) is str else " " for abstract in abstracts["value"]]
unique_terms = [word for abstract in word_bag for word in abstract]
unique_terms = sorted(set(unique_terms))
pd.DataFrame(unique_terms).T.to_csv("all_dbt.csv", header=False, index=False)

# Do it in 100 batches
limit = len(abstracts) / 100
i = 0
j = 0
abstract_dbt = []

for abstract in word_bag:
    abstract_row = []
    for word in unique_terms:
        abstract_row.append(abstract.count(word))
    abstract_dbt.append(abstract_row)
    i = i + 1
    if i >= limit:
        i = 0
        print(j)
        j = j + 1
        abstract_dbt = pd.DataFrame.from_records(abstract_dbt)
        abstract_dbt.to_csv("all_dbt.csv", mode="a", header=False, index=False)
        abstract_dbt = []

abstract_dbt
abstract_dbt = pd.DataFrame.from_records(abstract_dbt) # the last two ones are for the remaining papers after the loop ends
abstract_dbt.to_csv("all_dbt.csv", mode="a", header=False, index=False)
# manually removed columns that start with numbers in csv of rep
