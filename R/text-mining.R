#' Pubmed-textmining
#'
#'  Function for text-mining the PubMed repository based on defined sets of terms.
#' The relationship between fixterms (related to your research topic) and pubterms (terms which pivot around your research focus)
#' is calculated using the pointwise mutual information algorithm. A text file is generated with the pmi-scores for each fixterm.
#' Then for each collocation pairs (a fixterm + a pubterm), a text file is generated with related article titles and publishing years.
#' Additional Author section will follow in the next version updates.
#'
#' @param fixterms The input argument fixterms refers to any term that strongly relates to your focus of research. Either a string or a vector of strings.
#' @param pubterms The input argument pubterms recapitulates any terms you wish to pivot around your research focus. Either a string or a vector of strings.
#' @param output The output argument defines the directory you wish to store the generated text files. Default is your current location.
#' @return Returns text files with pmi-scores for each fixterm and text files for every collocation pairs bearing the article titles and publishing years. Authors will follow in the next version.
#' @export
#' @examples
#' \dontrun{
#' fixterms = c("bike", "downhill")
#' pubterms = c("dangerous", "extreme", "injuries")
#' output = getwd() #or "YOUR/DESIRED/OUTPUT/PATHWAY"
#' pubmed_textmining(fixterms, pubterms, output)
#' }

pubmed_textmining <- function(fixterms, pubterms, output){

  #set working directory
  if (missing(output)) {
  output = getwd()
  cat(paste("Attention, you did not enter any output path.", "Default pathway is set to your current location.", "\n\n"))
  }

  if (missing(pubterms) | missing(fixterms))
    stop("You are missing a character vector for the pubterms and/or fixterms argument")

  if (!is.character(fixterms) | !is.character(pubterms))
    stop("The pubterms and/or fixterms you have entered are not character vectors. Please review the sets of terms.")

  #################################################
  ### Some useful functions from previous years ###
  #################################################

  # Counting term frequency in pubmed (requires internet access!)
  pubmed_count <- function(query)
  {
    res = easyPubMed::get_pubmed_ids(query)

    count <- res$Count
    return(as.numeric(count))
  }

  # pointwise mutual information between two PubMed search terms
  pmi  <- function(term1, term2)
  {
    #get current size of Pubmed (for normalization)
    cur_size <- pubmed_count("1800:2100[dp]")

    count1 <- pubmed_count(term1)
    count2 <- pubmed_count(term2)
    count12 <- pubmed_count(paste(term1,'+AND+',term2,sep=""))

    return( log( count12 * cur_size / (count1 * count2)) )
  }

  #################################################
  ### term related text-mining functions ##########
  #################################################

  #pointwise mutual information definition
  pmi_definition = "Definition of Pointwise Mutual Information (PMI) scoring from Wikipedia:\n\nGood collocation pairs have high PMI because the probability of co-occurrence\nis only slightly lower than the probabilities of occurrence of each word. Conversely,\na pair of words whose probabilities of occurrence are considerably higher than their\nprobability of co-occurrence gets a small PMI score.\n"

  #pmi of one fixterm and several pubterms

  ######################################################################
  ### fixterms = terms that define the focus of your research ##########
  ### pubterms = terms whish can pivot around your focus as you wish ###
  ######################################################################

  #calculate and extract pmi
  for (j in 1:length(fixterms))
  {
    #Initialize pubterms score dataframe for each fixterm
    pubterms_score = {}
    for(i in 1:length(pubterms))
    {
      cat(paste("Collocation pairs:", fixterms[j], "~", pubterms[i]))
      pubterms_score[i] = pmi(fixterms[j], pubterms[i])
      cat(paste(" - PMI-score:",pubterms_score[i], "\n"))
    }
    #Save pmi score plus definition into txt file
    filename = paste(fixterms[j], format(Sys.time(), "%H-%M-%S") ,"pmi-scores.txt", sep="-")
    utils::write.table(pmi_definition, file.path(output, filename), col.names = FALSE, row.names = FALSE, quote = FALSE)
    utils::write.table(paste(fixterms[j], "~", pubterms, "PMI-score:", pubterms_score, sep = " "), file.path(output, filename), append=TRUE, col.names=FALSE,row.names=FALSE,sep="\t", quote=FALSE)
    cat(paste("\n"))
  }

  ##########################################################
  ### Extract relevant article titles ######################
  ##########################################################

  #loop to extract article titles which match the given query
  for(k in 1:length(fixterms))
  {
    for(t in 1:length(pubterms))
    {
      #reset collections
      my_query = {}
      my_titles = {}
      my_year = {}
      my_years = {}
      #load query
      my_query[t] <- paste(fixterms[k],'+AND+',pubterms[t],sep="")
      my_entrez_id <- easyPubMed::get_pubmed_ids(my_query[t])
      if (isTRUE(my_entrez_id$Count == 0)){
        next
      } else {
        cat(paste("Writing text file for", fixterms[k], "~", pubterms[t], ":"))
        #extract xml file of related articles
        my_xml <- easyPubMed::fetch_pubmed_data(pubmed_id_list = my_entrez_id)
        #extract title and publication year
        my_titles <- easyPubMed::custom_grep(my_xml, "ArticleTitle", "char")
        my_year <- easyPubMed::custom_grep(my_xml, "PubDate", "char")
        #extract years vector for each title for this pair
        for (n in 1:length(my_year))
        {
          #some publication years are displayed with "medlinedate" instead of "pubdate"
          if (is.na(my_year[n]) == TRUE | is.null(my_year[n]) == TRUE) {
            my_years[n] <- NA
          }
          else if (stringr::str_detect(my_year[n], "MedlineDate") == TRUE) {
            my_years[n] <- easyPubMed::custom_grep(my_year[n],"MedlineDate", "char")
          } else {
            my_years[n] <- easyPubMed::custom_grep(my_year[n],"Year", "char")
          }
        }
        #save year and title in text files separated for collocation pairs
        filename <- paste(my_query[t], ".txt", sep="")
        utils::write.table(paste(my_years, my_titles, sep = " "), file.path(output, filename), col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
        cat(paste(" done", "\n"))
      }
    }
    cat(paste("\n"))
  }
  cat(paste("Results saved in", output, "\n", "Please note that pairs with pmi-score of -Inf have been skipped.", "\n"))
}
###########################################################
### END TEXT MINING #######################################
###########################################################
