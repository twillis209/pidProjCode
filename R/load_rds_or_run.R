#' @title Check for RDS file and if absent, make it
#'
#' @param file_path Path to file
#' @param function_to_run Function to generate object
#'
#' @return value of function_to_run
#'
#' @export
load_rds_or_run <- function(file_path, function_to_run) {
  if(file.exists(file_path)) {
    readRDS(file_path)
  } else {
    result <- function_to_run()

    saveRDS(result, file = file_path)

    result
  }
}
