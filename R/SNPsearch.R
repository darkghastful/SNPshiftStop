#' SNPsearch
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_split str_length str_remove
#' @importFrom Biostrings translate DNAString
#' @importFrom bqutils content.from.endpoint subset.object
#' @importFrom methods is
#'
#' @param rsid rsid of SNP
#'
#' @return frame of mutations asssociated with rsid
#' @export
#'
#' @examples
#' SNPsearch("rs333")
SNPsearch <- function(rsid){
  url <- paste0("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/", str_remove(rsid, "rs"))
  snp.info <- content.from.endpoint(url)

  results <- list()
  allele_annotations <- snp.info$primary_snapshot_data$allele_annotations
  for(ann.number in 1:nrow(allele_annotations)){
    ann <- allele_annotations[ann.number,]
    asms <- ann$assembly_annotation[[1]]
    for(asm.number in 1:nrow(asms)) {
      asm <- asms[asm.number,]
      if(asm$annotation_release != "Homo sapiens Annotation Release 110") next
      genes <- asm$genes[[1]]
      for(gene.number in 1:length(genes)) {
        gene <- genes[gene.number,]$locus
        rnas <- genes[gene.number, "rnas"][[1]]
        for(rna.number in 1:length(rnas)) {
          rna <- rnas[rna.number,]
          prot <- rna$protein
          if(is.null(prot)){
            next
          }
          var <- prot$variant
          if(all(is.na(var))){
            next
          }
          frame.list <- list(
            sequence_ontology="prot$sequence_ontology[[1]]$name",
            ref.dna="rna$id",
            # bp.pos="rna$codon_aligned_transcript_change$position",
            ref.bp="rna$codon_aligned_transcript_change$deleted_sequence",
            alt.bp="rna$codon_aligned_transcript_change$inserted_sequence",
            ref.pro="rna$product_id",
            # aa.pos="rna$protein$variant$spdi$position",
            # ref.aa="rna$protein$variant$spdi$deleted_sequence",
            # alt.aa="rna$protein$variant$spdi$inserted_sequence",
            hgvs="rna$hgvs"
          )
          append.list <- list(gene=gene)
          for(element in 1:length(frame.list)){
            if(is(try(eval(parse(text=frame.list[element]))), "try-error")| is.null(try(eval(parse(text=frame.list[element]))))){
              append.list[[element+1]] <- NA
            }else{
              append.list[[element+1]] <- eval(parse(text=frame.list[element]))
            }
            names(append.list)[element+1] <- c(names(frame.list)[element])
          }
          results <- append(results, list(append.list))
        }
      }
    }
  }

  # Convert results to a data frame
  frame <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
  return(frame)
}
