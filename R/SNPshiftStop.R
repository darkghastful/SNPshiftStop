
# if aa.mut.pos not porvided results may not be accurate
# if ref aa not porvided results may not be accurate
# if alt aa not porvided results may not be accurate

#' @importFrom magrittr %>%
#' @importFrom stringr str_split str_length str_remove
#' @importFrom Biostrings translate DNAString
#' @importFrom bqutils content.from.endpoint
#' SNPshiftStop
#'
#' @param gene.symbol Gene symbol with frameshift mutation
#' @param taxon provide taxon character string or number; found using NCBI (default human (9606))
#' @param mutation.from.start Location of mutation relative to start
#' @param ref.bp reference base pair (default is '"-"')
#' @param alt.bp alternate base pair (default is '"-"')
#' @param ref.aa reference amino acid (default is NA)
#' @param alt.aa alternate amino acid (default is NA)
#' @param aa.mut.pos amino acid mutation position (default is NA)
#' @param start.loss does this mutation contain a start.loss (default is TRUE)
#'
#' @return (gene.symbol, aa.seq, mutated.aa.seq, shortened.of.total.aa, dna.seq, mutated.dna.seq)
#' @export
#'
#' @examples
#' \donttest{
#' SNPshiftStop(gene.symbol="CFTR", taxon="human", mutation.from.start=1213, ref.bp="T", ref.aa="F", alt.aa="L", aa.mut.pos=405)
#' }
SNPshiftStop <- function(ncbi.api.key=Sys.getenv("ncbi.api.key"), gene.symbol, taxon="human", mutation.from.start, ref.bp="-", alt.bp="-", ref.aa=NA, alt.aa=NA, aa.mut.pos=NA, start.loss=FALSE){

  # Stop conditions
  if(!nzchar(ncbi.api.key)){
    stop("Please provide ncbi.api.key")
  }

  if((ref.bp=="-" & alt.bp=="-")){
    stop("Please provide either ref.bp or alt.bp.")
  }

  # Warning conditions

  if(is.na(ref.aa) | is.na(alt.aa) | is.na(aa.mut.pos)){
    if(is.na(ref.aa) | is.na(aa.mut.pos)){
      warning("The reference amino acid can not be checked against the sequence pulled from NCBI.")
    }else if(is.na(alt.aa) | is.na(aa.mut.pos)){
      warning("The alternate amino acid can not be checked against the generated mutant protein.")
    }

    if(is.na(ref.aa)){
      warning("To ensure accuracy please provide the reference amino acid (ref.aa).")
    }
    if(is.na(alt.aa)){
      warning("To ensure accuracy please provide the alternate amino acid (alt.aa).")
    }
    if(is.na(aa.mut.pos)){
      warning("To ensure accuracy please provide the location of the amino acid mutation (aa.mut.pos).")
    }
  }



  if(ref.bp=="-"){
    ref.bp <- NA
    ins.del <- "insertion"
    ref.alt.bp <- alt.bp
  }else if(alt.bp=="-"){
    alt.bp <- NA
    ins.del <- "deletion"
    ref.alt.bp <- ref.bp
  }

  ##### Pull gene info
  url <- paste0("https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/", gene.symbol, "/taxon/", taxon, "?api_key=", ncbi.api.key)
  gene.info <- bqutils::content.from.endpoint(url, content_type="")$reports[1,1]

  ######################## Ensure that the provided position and bp mutation are consistant with the sequence from entrez

  ##### Pull protein coding information from the product report of the NCBI page for the gene symbol
  url <- paste0("https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/", gene.symbol, "/taxon/", taxon, "/product_report?types=PROTEIN_CODING&api_key=", ncbi.api.key)
  proteins <- bqutils::content.from.endpoint(url, content_type="")$reports[1,1]

  # Identify the protein of interest (based on MANE column)
  protein <- proteins[,"transcripts"][[1]] %>%
    .[!is.na(.[,"select_category"]),] # Matched ensemble dna and protein sequence

  protein.aa.length <- protein[,"protein"][,"length"]

  ##### Pull the coding dna sequence (cds) accession and the bp range of the protein within the gene from protein.cds.accession/protein.dna.accession
  protein.cds.accession <- protein[,"cds"][,"accession_version"]
  protein.coding <- protein[,"cds"][,"range"][[1]]


  #################### may or may not be correct should maybe rely on the dna string pulled as I know that the pulled dna sliced using "protein.coding"
  protein.dna.length <- protein[,"length"]

  ##### Pull the accession number provided for the protein of interest
  protein.dna.accession <- protein[,"accession_version"]


  ########## Ensure accession number provided for the protein of interest matches that given for the "cds" of the protein of interest
  if(protein.cds.accession!=protein.dna.accession){
    warning("Protein accession values (both from the same NCBI Dataset) do not match. Do not trust further results.")
  }

  # if(length(protein.dna.accession)>1){
  #   protein.dna.accession <- protein[protein[,"select_category"]=="MANE_SELECT", "accession_version"]
  # }
  if(length(protein.dna.accession)>1){
    protein.dna.accession <- protein %>%
      .[.[,"name"]=="transcript variant 3", "accession_version"]
  }


  # DNA
  ##### Pull dna seq from accession
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=", protein.dna.accession)
  dna.seq <- bqutils::content.from.endpoint(url, accept="text/plain", fromJSON=FALSE)
  dna.seq <- str_split(dna.seq, "\n")[[1]] %>%
    .[2:length(.)]
  dna.seq <- paste(dna.seq, collapse="")

  ##########
  if(protein.dna.length!=str_length(dna.seq)){ # Ensure the protein lengths are the same
    warning("Protein length from NCBI Datasets does not match that from Entrez. Results may not be accurate.")
  }

  #################### protein.dna.seq <- protein[,"protein"][,"length"]



  # AA
  # Pull aa seq from accession
  protein.aa.accession <- protein[,"protein"][,"accession_version"]

  if(length(protein.aa.accession)>1){
    protein.aa.accession <- protein[,"protein"] %>%
      .[.[,"isoform_name"]=="isoform 3", "accession_version"]
  }


  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id=", protein.aa.accession)
  aa.seq <- bqutils::content.from.endpoint(url, accept="text/plain", fromJSON=FALSE)
  aa.seq <- str_split(aa.seq, "\n")[[1]] %>%
    .[2:length(.)]
  aa.seq <- paste(aa.seq, collapse="")

  ##########
  if(protein.aa.length!=str_length(aa.seq)){ # Ensure the protein lengths are the same
    warning("Protein length from NCBI Datasets does not match that from Entrez. Results may not be accurate.")
  }

  ## Reassign the proein aa length just in case the previous warning has a true statement (can be removed when made into a function)
  protein.aa.length <- str_length(aa.seq)



  # Double check that both the accessed amino acid sequence and translated dna sequence match before proceeding
  ### Pull the DNA sequences for each row in protein.coding
  for(row in 1:nrow(protein.coding)){
    begin <- protein.coding[row, "begin"]
    end <- protein.coding[row, "end"]
    protein.coding[row,"sequence"] <- paste(str_split(dna.seq, "")[[1]][begin:end], collapse="")
  }

  #################### there may never be more than one row for the protein coding but make sure that this is the case
  ### Collapse the protein coding columns in to one usable sequence
  coding.seq <- paste0(protein.coding[,"sequence"], collapse="")

  ### Translate the "coding.seq" to aa
  translated.aa.seq <- translate(DNAString(coding.seq))

  # Remove stop from "translated.aa.seq" (for comparison with the accessed "aa.seq")
  translated.aa.seq <- as.character(translated.aa.seq) %>%
    str_remove(., "\\*")

  ########## Double check that the "aa.seq" of the protein of interest matches the "dna.seq" translation ("translated.aa.seq") of the protein of interest
  if(aa.seq!=translated.aa.seq){ #
    warning("The amino acid sequence generated from the dna accession does not match that from the protein accession. Do not trust further results.")
  }



  # Analyze the mutated sequence
  ####################
  # if(ceiling(mutation.from.start/3)!=aa.mut.pos){
  #   warning("The provided mutation location from the start of the coding region does not match the provided location for the aa mutation.")
  # }

  if(ins.del=="deletion"|ins.del=="del"){
    mutated.dna.seq <- paste0(str_split(coding.seq, "")[[1]][-c(mutation.from.start)], collapse="")
  }else if(ins.del=="insertion"|ins.del=="ins"){
    mutation.split <- str_split(coding.seq, "")[[1]]
    section.1 <- paste0(mutation.split[1:min(mutation.from.start)-1], collapse="")
    section.2 <- paste0(mutation.split[min(mutation.from.start):length(mutation.split)], collapse="")
    mutated.dna.seq <- paste0(section.1, alt.bp, section.2, collapse="")
    # alt.bp==paste(str_split(mutated.dna.seq, "")[[1]][mutation.from.start], collapse="")
  }

  # Add check to see if the mutation is actually a frameshift
  if(!is.na(ref.bp)){
    if(length(str_split(ref.bp, "")[[1]])%%3==0){
      warning("The provided mutation is not a proper frameshift. The deletion is contains full codons.")
    }
  }else if(!is.na(alt.bp)){
    if(length(str_split(alt.bp, "")[[1]])%%3==0){
      warning("The provided mutation is not a proper frameshift. The insertion is contains full codons.")
    }
  }

  if(start.loss){
    mutated.dna.seq <- paste0("ATG", str_split(mutated.dna.seq, "ATG", 2)[[1]][2])
  }
  mutated.aa.seq <- suppressWarnings(as.character(translate(DNAString(mutated.dna.seq))))
  # str_split(mutated.aa.seq, "")[[1]][aa.mut.pos]


  if(ins.del=="deletion"|ins.del=="del"){
    if(ref.bp!=paste0(str_split(coding.seq, "")[[1]][mutation.from.start], collapse="")&!start.loss){
      stop(paste0("The deleted base pair (at the provided muation site) does not match the reference base pair in the sequence aquired from NCBI.",
                  " Provided:", ref.bp, " NCBI:", paste0(str_split(coding.seq, "")[[1]][mutation.from.start], collapse="")))
      # stop()
    }
  }else if(ins.del=="insertion"|ins.del=="ins"){
    if(alt.bp!=paste0(str_split(mutated.dna.seq, "")[[1]][mutation.from.start], collapse="")&!start.loss){
      stop(paste0("The inserted base pair at the provided muation site does not match the alternate base pair in the sequence aquired from NCBI.",
                  " Provided:", alt.bp, " NCBI:", paste0(str_split(mutated.dna.seq, "")[[1]][mutation.from.start], collapse="")))
    }
  }

  if(!is.na(ref.aa) & !is.na(aa.mut.pos)){
    if(ref.aa!=str_split(aa.seq, "")[[1]][aa.mut.pos]&!is.na(ref.aa)&!start.loss){
      stop(paste0("The provided reference amino acid does not match the amino acid in the mutation position of the protein sequence aquired from NCBI.",
                  " Provided:", ref.aa, " NCBI:", str_split(aa.seq, "")[[1]][aa.mut.pos]))
    }
  }

  if(!is.na(alt.aa) & !is.na(aa.mut.pos)){
    if(alt.aa!=str_split(mutated.aa.seq, "")[[1]][aa.mut.pos]&!is.na(alt.aa)&!start.loss){
      stop(paste0("The mutation introduced to the protein sequence aquired from NCBI does not match the amino acid mutation provided
                  (at the position provided).",
                  " Provided:", alt.aa, " NCBI:", str_split(mutated.aa.seq, "")[[1]][aa.mut.pos]))
    }
  }


  mutated.protein.length <- str_length(str_split(mutated.aa.seq, "\\*")[[1]][1])

  missing.aa <- protein.aa.length-(mutated.protein.length+1)
  shortened.of.total.aa <- paste0(mutated.protein.length+1, "/", protein.aa.length)

  # print(aa.seq)
  # print(mutated.aa.seq)
  # print(gene.symbol)
  # print(ref.alt.bp)
  # print(ins.del)
  # print(mutation.from.start)
  # print(aa.mut.pos)
  # print(str_split(aa.seq, "")[[1]][aa.mut.pos])
  # print(paste0(str_split(mutated.aa.seq, "")[[1]][aa.mut.pos]))

  cat(paste(
    "\n\033[1mReference:\033[22m\n", aa.seq,
    "\n\033[1mAlternate:\033[22m\n", mutated.aa.seq,
    "\n\n\033[1mIn the gene", gene.symbol, "there was an", ref.alt.bp, ins.del, "at position", paste0(mutation.from.start, "."),
    "The frameshift began at", aa.mut.pos, "with a change of", str_split(aa.seq, "")[[1]][aa.mut.pos], "to", paste0(str_split(mutated.aa.seq, "")[[1]][aa.mut.pos], "."),
    "This introduced a premature stop at position", shortened.of.total.aa, "(shortening the protein by", missing.aa, "amino acids).\033[22m\n\n"
  ))
  # cat(paste("\nThe frameshift introduced a premature stop at position", shortened.of.total.aa, "(shortening the protein by", missing.aa, "amino acids)"))

  return(list(gene.symbol=gene.symbol, aa.seq=aa.seq, mutated.aa.seq=mutated.aa.seq, shortened.of.total.aa=shortened.of.total.aa, dna.seq=dna.seq, mutated.dna.seq=mutated.dna.seq))
}


