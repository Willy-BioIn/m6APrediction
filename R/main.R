
#' Nucleotide-level Encoding for DNA String Vectors
#'
#' Splits an **equal-length DNA string vector** into a factorized data frame at the nucleotide level.
#' Each column corresponds to a nucleotide position, and factor levels for all columns are standardized to "A", "T", "C", "G"
#' for subsequent feature processing in machine learning models.
#'
#' @param dna_strings A character vector where each element is a pure DNA sequence (containing only "A", "T", "C", "G").
#' All sequences **must have the same length** (the function uses the length of the first sequence as the standard).
#' @return A factorized data frame with column names formatted as "nt_pos1", "nt_pos2", ..., "nt_posN" (N = length of DNA sequences).
#' Each column has fixed factor levels: c("A", "T", "C", "G"). Rows correspond to each DNA string in the input.
#' @examples
#' # Example:
#' example_data<-read.csv(system.file("extdata", "m6A_input_example.csv",package="m6APrediction"))
#' dna_5mer_vec <- example_data$DNA_5mer
#' dna_encoded_df <- dna_encoding(dna_5mer_vec)
#' head(dna_encoded_df, 3)
#' @export

dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Batch Prediction of m6A Modification Status
#'
#' Uses a pre-trained machine learning model and a batch feature data frame to predict m6A modification probability
#' and classification status for each sample. Automatically standardizes factor levels for RNA type and RNA region,
#' processes 5-mer sequences via the dna_encoding function, and returns the complete data frame with prediction results.
#'
#' @param ml_fit A pre-trained machine learning model (e.g. randomForest).
#' Support the `predict()` method with `type = "prob"` to output a probability matrix (containing a "Positive" column).
#' @param feature_df A data frame **required to contain the following 7 columns** (order does not matter):
#' gc_content (numeric, GC content), RNA_type (character/factor, RNA type),
#' RNA_region (character/factor, RNA region), exon_length (numeric, exon length),
#' distance_to_junction (numeric, distance to junction), evolutionary_conservation (numeric, evolutionary conservation score),
#' DNA_5mer (character, 5-mer DNA sequence, containing only "A", "T", "C", "G").
#' @param positive_threshold A numeric value defining the probability threshold for positive ("Positive") predictions. Defaults to 0.5.
#' @return A data frame with 2 additional columns appended to the input `feature_df`:
#' - predicted_m6A_prob: Numeric, probability of being predicted as m6A-positive.
#' - predicted_m6A_status: Factor, predicted status with levels: "Negative" and "Positive".
#' @seealso \code{\link{dna_encoding}} (dependency function for DNA_5mer processing),
#' \code{\link{prediction_single}} (single-sample prediction function)
#' @examples
#' rf_fit<-readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' example_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' pred_result <- prediction_multiple(ml_fit = rf_fit, feature_df = example_data)
#' head(pred_result, 5)
#'
#' @import randomForest
#' @export




prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  feature_df_t <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))

  predicted_m6A_prob_t<-predict(ml_fit,feature_df_t,type = "prob")
  predicted_m6A_prob<-predicted_m6A_prob_t[,"Positive"]

  predicted_m6A_status_t<-ifelse(predicted_m6A_prob>positive_threshold,"Positive","Negative")
  predicted_m6A_status<-factor(predicted_m6A_status_t,levels = c("Negative", "Positive"))

  feature_df<-cbind(feature_df,predicted_m6A_prob,predicted_m6A_status)
  return(feature_df)
}



#' Single-Sample Prediction of m6A Modification Status
#'
#' Takes feature parameters for a single sample (e.g., GC content, RNA type) and calls the prediction_multiple function
#' to predict m6A status. Eliminates the need to manually construct a data frame and directly returns a named vector
#' containing "predicted probability" and "predicted status" for easy quick calls.
#'
#' @param ml_fit A pre-trained machine learning model, identical to the `ml_fit` parameter in prediction_multiple (must support probability output).
#' @param gc_content A numeric value(0-1), GC content of the single sample (e.g., 0.45 = 45%).
#' @param RNA_type A character string, RNA type of the single sample. **Allowed values only**: "mRNA", "lincRNA", "lncRNA", "pseudogene".
#' @param RNA_region A character string, RNA region type of the single sample. **Allowed values only**: "CDS", "intron", "3'UTR", "5'UTR".
#' @param exon_length A numeric value, exon length of the single sample.
#' @param distance_to_junction A numeric value, distance of the single sample to the junction.
#' @param evolutionary_conservation A numeric value, evolutionary conservation score of the single sample.
#' @param DNA_5mer A character string, 5-mer DNA sequence of the single sample. **Must be 5 characters long** (containing only "A", "T", "C", "G").
#' @param positive_threshold A numeric value defining the probability threshold for positive predictions. Defaults to 0.5 (same as prediction_multiple).
#' @return A named character vector with 2 elements:
#' - "predicted_m6A_prob": Character-formatted numeric, probability of being predicted as m6A-positive (convert to numeric with as.numeric()).
#' - "predicted_m6A_status": Character, predicted status ("Positive" or "Negative").
#' @seealso \code{\link{prediction_multiple}} (underlying batch prediction function called by this function)
#' @examples
#' rf_fit<-readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' result<-prediction_single(
#'   ml_fit = rf_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "3'UTR",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.7,
#'   DNA_5mer = "ATCGT",
#'   positive_threshold = 0.5
#' )
#' print(result)
#' @export



prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){

  data_offer<-data.frame(
    gc_content=gc_content,
    RNA_type=RNA_type,
    RNA_region=RNA_region,
    exon_length=exon_length,
    distance_to_junction=distance_to_junction,
    evolutionary_conservation=evolutionary_conservation,
    DNA_5mer=DNA_5mer
  )

  res<-prediction_multiple(ml_fit,data_offer,positive_threshold)

  predicted_m6A_prob<-res[,"predicted_m6A_prob"]
  predicted_m6A_status<-as.character(res[,"predicted_m6A_status"])


  returned_vector<-c(predicted_m6A_prob,predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector)
}

