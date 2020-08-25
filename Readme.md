Script and data for assessing interdisciplinarity in cognitive science

# Data

## The *data* folder
All the data used to produce the results is contained in this repository.
The *data* folder includes the data used in the JSD - Gini part of the paper.
It contains two subfolders:
* The *journals* subfolder includes all papers published in 20 journals from 2009 to 2018 as obtained from Web of Science. This folder is compressed in the *journals.7z* archive and must be decompressed before running scripts to re-generate the data.
* The *authors* subfolder includes one .csv for each of the authors that coauthored any of the papers. Each .csv has the publications made by that authors in the period before and during the data contained in *journals*. The authors were divided into 4 .7z archives that must be decompressed before running scripts to re-generate the data.

## The *theories_data* folder

This folder includes the data **and scripts** for the LSA theory analysis part of the paper.
The *.py* scripts take a catalog and transform it into a lemmatized document by term matrix.
The scripts have to be run in the following order:

1. clean_database.py
2. lemmatizer_stemmer.py
3. document_by_term_generator.py

To produce the data for the study, all 3 scripts must be run on *new_abstracts.csv* and *catalog.scientometrics.csv*.
The output of these scripts generate *clean_database.csv* (Scientometrics paper data) and *clean_all_abstracts.csv* (this paper's new data) through scripts *1* and *2*, and *dbt_matrix_scientometrics.csv* (Scientometrics data) and *all_dbt.csv* (New data, compressed into *all_dbt.zip*). *all_dbt_reduced.zip* contains a compressed version of the reduced document by term matrix produced in the *theories.R* script.

# R scripts

The *r_scripts* folder contains the main analysis scripts for the paper. The *saved_objects* folder contains *.Rds* files that save the intermediate objects generated in these scripts to save time.
*data_processing.R* takes the files in the *data* folder and produces a usable version for the rest of the analysis. This script takes a very long time to run as each author file has to be incorporated. The output of this script is saved in *journalData.Rds* and *authorByJournal.rds*.
*JSD.R* calculates the Jensen-Shannon divergence and Gini Index for all papers in the *journals* data folder. The output of this script is contained in *complete_frame.rds*.
*theories.R* builds semantic space with each theory included in the Scientometrics paper, and applies them to the *journal* papers in Cognitive Science journals through a random forest. The predictions and scores for each theory in each article of the CogSci journals are saved in *cogsci_theories.rds*.
*main_script.R* takes the objects produced by the other scripts and generates the figures and models reported in the paper. Optionally, it can be commented/uncommented to run all scripts and not rely on the *saved_objects* folder.
*functions.R* contains most of the general purpose functions used in the rest of the scripts.

# Figures

* Folder *author_figures* includes the theory distribution for each author who (co-)authored a paper in the three cognitive science journals being considered.
* Folder *main_figures* includes the figures used in the main text of the paper.
* Folder *supp_figures* includes the figures used as supplementary material in the paper.
