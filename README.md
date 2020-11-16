# pubmed_textmining()
Term-based text-Mining of the 'PubMed' repository.
Requires internet access!

The function takes vectors of string as term sets for fix-terms (research topic focus) and pub-terms (terms which can pivot around your research focus).
Install the function with `devtools::install_github("Moshkante/pubmed_mining")` or simply via CRAN.
Mining results are given in text format with pointwise mutual information ('PMI'-score) and related article titles and publishing years.
More features will be added in the next version.

Example:  
`fixterms = "bike" #or multiple terms    
pubterms = c("dangerous", "extreme")  
output = getwd() #or "YOUR/DESIRED/PATHWAY"  
pubmed_textmining(fixterms, pubterms, output)`  
