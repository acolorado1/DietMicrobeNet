# load necessary packages 
packages <- c("readr", "dplyr", "argparse")
for (package in packages) {
  if(!require(package, character.only = T)){
    install.packages(package)
    library(package)
  }
}

# create parser
parser <- ArgumentParser()
parser$add_argument('--diet_file', help = "CSV file containing food items")
parser$add_argument('--content_file', help = 'CSV file containing content info on food')
parser$add_argument('--ExDes_file', help = 'CSV file containing external descriptors on compounds')
parser$add_argument('--output_file', help = 'file path for output')
args <- parser$parse_args()

#' From foods found in FooDB get list of KEGG compounds found in their collective metabolomes
#'
#' @param diet_df CSV generated through included website 
#' @param content_df content CSV  from FooDB
#' @param external_descriptor_df CompoundExternalDescriptor CSV from FooDB
#' @param output_path file path for KEGG compounds
#'
#' @return file with a list of KEGG Compounds 
#'
#' @examples
#' 
get_diet_fooDB_compounds <- function(diet_df, 
                                     content_df, 
                                     external_descriptor_df, 
                                     output_path){ 
  
  # load in data 
  foods <- read_csv(diet_df)
  Content <- read_csv(content_df)
  CExtDes <- read_csv(external_descriptor_df)
  
  # get food IDs
  ids <- foods$id

  # Get compound IDs, I believe these are labeled as the Source IDs in the contents file
  food_content <- Content[Content$food_id %in% ids, ]
  comp_ids <- unique(food_content$source_id)

  # get external descriptors
  extDes <- CExtDes[CExtDes$compound_id %in% comp_ids, ]
  extDes <- extDes %>%
    select(external_id, compound_id) %>%
    na.omit()

  # get only external descriptors that are KEGG related
  comp_KEGG <- unique(subset(extDes, grepl("^C[0-9]{5}$", external_id)))

  write.table(comp_KEGG$external_id, output_path, row.names = F, sep = '\n', col.names = F, quote = F)
}

# call function
get_diet_fooDB_compounds(diet_df = args$diet_file,
                        content_df = args$content_file,
                        external_descriptor_df = args$ExDes_file,
                        output_path = args$output_file)
