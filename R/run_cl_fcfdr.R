#' @title Run fcFDR program with command-line arguments
#'
#' @param cl_args Vector of command-line arguments
#' @param PID_ROOT Path to root of PID project directory tree
#' @import data.table
#' @import argparse
#' @importFrom fcfdr flexible_cfdr
#' @export
run_cl_fcfdr <- function(cl_args, PID_ROOT = Sys.getenv('pidRoot')) {
  parser <- ArgumentParser(description = 'Run the fcFDR procedure on a set of GWAS results')
  parser$add_argument('-if', '--input_file', type = 'character', help = 'Path to input file', default = file.path(PID_ROOT, 'joinedData/pid/pid_nine_imd_two_contr.tsv.gz'))
  parser$add_argument('-chr', '--chromosome', type = 'character', help = 'Label of chromosome column in GWAS file', default = 'CHR38')
  parser$add_argument('-bp', '--basepair', type = 'character', help = 'Label of BP column in GWAS file', default = 'BP38')
  parser$add_argument('-m', '--maf', type = 'character', help = 'Label of minor allele frequency column in GWAS file', default = 'ALT_FREQ')
  parser$add_argument('-ma', '--no_match', action = 'store_true', help = 'Do not match MAF distribution between whole set and independent subset', default = FALSE)
  parser$add_argument('-p', '--principal', type = 'character', help = 'Label of principal p-value column in GWAS file', required = T)
  parser$add_argument('-q', '--auxiliary', action = 'store', dest = 'auxiliary', type = 'character', help = 'Comma-delimited list of labels of auxiliary trait columns in GWAS file.', required  =  T)
  parser$add_argument('-w', '--weights', action = 'store', dest = 'weights', type = 'character', help = 'Comma-delimited list of labels of weight columns. Order should match that of auxiliary trait labels', required  =  T)
  parser$add_argument('-v', '--v_values', action = 'store', dest = 'v_values', type = 'character', help = 'Comma-delimited list of labels of v-value columns. Order should match that of auxiliary trait labels', required  =  T)
  parser$add_argument('-ac', '--add_columns', nargs = '*', type = 'character', help = 'Additional columns to include in the output to identify SNPs', default = c('SNPID', 'REF', 'ALT'))
  parser$add_argument('-at', '--aux_transform', nargs = '+', type = 'character', help = 'Comma-delimited list of transformations to apply to auxiliary values. Order should match that of auxiliary trait labels. Current valid values: \'identity\', \'log\', \'z\'', required = T)
  parser$add_argument('-op', '--outputPath', type = 'character', help = 'Path to output file', required = T)
  parser$add_argument('-nt', '--noOfThreads', type = 'integer', help = 'Number of threads to use', default = 1)

  trans <- list(log = log, identity = identity, z = function(x) qnorm(x/2))

  args <- parser$parse_args(cl_args)

  args$auxiliary <- unlist(strsplit(args$auxiliary, ','))
  args$weights <- unlist(strsplit(args$weights, ','))
  args$v_values <- unlist(strsplit(args$v_values, ','))
  args$aux_transform <- unlist(strsplit(args$aux_transform, ','))

  if(length(args$aux_transform) != length(args$auxiliary)) {
    stop("No. of specified auxiliary transformations does not match number of auxiliary covariates")
  }

  if(any(!is.element(args$aux_transform, c('identity', 'log', 'z')))) {
    stop("Specified auxiliary transformations contain values not matching one of following valid transformation identifiers: \'identity\', \'log\', \'z\'")
  }

  setDTthreads(threads=args$noOfThreads)

  gwasCols <- c(args$chromosome, args$basepair, args$add_columns, args$maf, args$principal, args$auxiliary, args$weights)

  dat <- fread(args$input_file, sep = '\t', header = T, select = gwasCols)

  if(any(!is.element(gwasCols, names(dat)))) {
    stop(sprintf("data.table object is missing the following columns specified in the command-line arguments: %s\n", paste(gwasCols[!is.element(gwasCols, names(dat))], collapse = ', ')))
  }

  dat <- dat[!(get(args$chromosome) == 6 & get(args$basepair) %between% c(24e6, 45e6))]

  # If ur-principal p-value is absent, can't do anything
  dat <- dat[!is.na(get(args$principal)) & !is.na(get(args$maf))]

  # Transforming p-values of 0 to Z-scores yields Inf values
  for(x in c(args$principal, args$auxiliary)) {
    dat[ get(x) < 1e-300, (x) := 1e-300]
  }

  prin_pvalues <- c(args$principal, args$v_values)

  for(i in seq_along(args$auxiliary)) {
    # prin_pvalues[i] should never contain NA values as we drop the ur_prin_pvalues NA rows then carry over
    # If we only want to include values for which we have LDAK weight values, include !is.na(get(aux_weights[i]))
    print(sprintf('Iteration %d', i))

    p <- dat[!is.na(get(args$auxiliary[i])), get(prin_pvalues[i])]

    q <- dat[!is.na(get(args$auxiliary[i])), get(args$auxiliary[i])]

    q_trans <- trans[[args$aux_transform[i]]](q)

    # TODO check for NA in transformed values
    if(any(is.na(q_trans))) {
      stop(sprintf("Iteration %d: NA value(s) in transformed auxiliary values", i))
    }

    ind_indices <- dat[!is.na(get(args$auxiliary[i]))][get(args$weights[i]) > 0, which = T]

    maf <- dat[!is.na(get(args$auxiliary[i])), get(args$maf)]

    if(args$no_match) {
      fcfdr_result <- fcfdr::flexible_cfdr(p = p,
                                  q = q_trans,
                                  indep_index = ind_indices,
                                  maf = NULL,
                                  check_indep_cor = F,
                                  enforce_p_q_cor = F)
    } else {
      fcfdr_result <- fcfdr::flexible_cfdr(p = p,
                                    q = q_trans,
                                    indep_index = ind_indices,
                                    maf = maf,
                                    check_indep_cor = F,
                                    enforce_p_q_cor = F)
    }

    dat[, prin_pvalues[i+1] := Inf]

    dat[!is.na(get(args$auxiliary[i])), prin_pvalues[i+1] := fcfdr_result[[1]]$v]

    dat[is.infinite(get(prin_pvalues[i+1])), (prin_pvalues[i+1]) := get(prin_pvalues[i])]

    if(nrow(dat[is.infinite(get(prin_pvalues[i+1]))]) > 0) {
      stop(sprintf("%d infinite-valued rows left in %s", nrow(dat[is.infinite(get(prin_pvalues[i+1]))]), prin_pvalues[i+1]))
    }

    fwrite(dat, file = args$outputPath, sep  = '\t', row.names = F, col.names = T, quote = F)
  }
}
