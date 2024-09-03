#' run_rsubread
#'
#' @param dir path to directory containing all the bam files
#' @param gff_file path to gff file for reference genome
#'
#' @return a count matrix
#' @import magrittr
#' @import Rsubread
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
#'
#' @examples
run_rsubread <- function(dir, gff_file){


  bam_files <- tibble::as_tibble(dir)

  count_tibble <- bam_files %>%
    dplyr::mutate(read_count = purrr::map(value, function(ii){
      cd <- Rsubread::featureCounts(files = ii,
                                    annot.ext = gff_file,
                                    isGTFAnnotationFile = TRUE,
                                    isPairedEnd = TRUE, nthreads = 8,
                                    strandSpecific = 2) # ref: https://github.com/igordot/genomics/blob/master/notes/rna-seq-strand.md
      cc <- cd$counts %>% as.data.frame()%>%
        tibble::rownames_to_column("gene_name")
    }))

  count_mat <- count_tibble %>%
    dplyr::select(2) %>%
    tidyr::unnest(cols = c(read_count)) %>%
    tidyr::pivot_longer(cols = -c(gene_name)) %>%
    tidyr::drop_na() %>%
    tidyr::pivot_wider(names_from = name,values_from = value) %>%
    dplyr::rename_if(is.numeric,
                     dplyr::funs(stringr::str_replace(., "_star_alignAligned.sortedByCoord.out.bam", "")))

  return(count_mat)

}
