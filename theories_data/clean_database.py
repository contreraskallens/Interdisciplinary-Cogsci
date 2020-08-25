# this script produces a version of the catalog of papers without the keywords used for searching each theory,
# in lower case, missing symbols and grammatical markers.

import re
import pandas as pd

text = pd.read_csv("new_abstracts.csv")
abstracts = text["value"]

abstracts = [abstract.lower() if type(abstract) is str else "" for abstract in abstracts]
abstracts = [re.sub(r'\([Cc]\).+\.', '', abstract) for abstract in abstracts] # removes copyright
abstracts = [re.sub(r'[^\w\s]s', '', abstract) for abstract in abstracts] # removes possessives
abstracts = [re.sub(r'[^\w\s]', '', abstract) for abstract in abstracts]  # removes punctuantion and symbols

text["value"] = abstracts
text.to_csv("clean_all_abstracts.csv", index=False)
