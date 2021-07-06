#' @title Run cFDR program with command-line arguments
#'
#' @param cl_args Vector of command-line arguments
#' @param PID_ROOT Path to root of PID project directory tree
#' @import data.table
#' @import argparse
#' @importFrom parallel mclapply mcmapply
#' @importFrom cfdr vl il fit.2g
#' @export
run_cl_cfdr <- function(cl_args, PID_ROOT = Sys.getenv('pidRoot')) {
  parser <- ArgumentParser(description = 'Run the fcFDR procedure on a set of GWAS results')
  parser$add_argument('-if', '--input_file', type = 'character', help = 'Path to input file', default = file.path(PID_ROOT, 'joinedData/pid/pid_nine_imd_two_contr.tsv.gz'))
  parser$add_argument('-chr', '--chromosome', type = 'character', help = 'Label of chromosome column in GWAS file', default = 'CHR38')
  parser$add_argument('-bp', '--basepair', type = 'character', help = 'Label of BP column in GWAS file', default = 'BP38')
  parser$add_argument('-p', '--principal', type = 'character', help = 'Label of principal p-value column in GWAS file', required = T)
  parser$add_argument('-q', '--auxiliary', action = 'store', dest = 'auxiliary', type = 'character', help = 'Comma-delimited list of labels of auxiliary trait columns in GWAS file.', required  =  T)
  parser$add_argument('-pt', '--p_threshold', type = 'double', help = 'Threshold for principal p-values', default = 1e-6)
  parser$add_argument('-w', '--weights', action = 'store', dest = 'weights', type = 'character', help = 'Comma-delimited list of labels of weight columns. Order should match that of auxiliary trait labels', required  =  T)
  parser$add_argument('-v', '--v_values', action = 'store', dest = 'v_values', type = 'character', help = 'Comma-delimited list of labels of v-value columns. Order should match that of auxiliary trait labels', required  =  T)
  parser$add_argument('-ac', '--add_columns', nargs = '*', type = 'character', help = 'Additional columns to include in the output to identify SNPs', default = c('SNPID', 'REF', 'ALT'))
  parser$add_argument('-o', '--outputPath', type = 'character', help = 'Path to output file', required = T)
  parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)
  parser$add_argument('-si', '--start_index', type = 'integer', help = 'Index of auxiliary trait to take as first auxiliary covariate. Used for restarting jobs where results file contains partial results.', default = 1)

  args <- parser$parse_args(cl_args)

  args$auxiliary <- unlist(strsplit(args$auxiliary, ','))
  args$weights <- unlist(strsplit(args$weights, ','))
  args$v_values <- unlist(strsplit(args$v_values, ','))

  if(args$start_index < 1 || args$start_index > length(args$auxiliary)) {
    stop("Invalid start index. Must be a positive integer and no greater than the number of auxiliary covariates supplied")
  }

  setDTthreads(threads=args$no_of_threads)

  gwasCols <- c(args$chromosome, args$basepair, args$add_columns, args$principal, args$auxiliary, args$weights)

  if(args$start_index > 1) {
    gwasCols <- c(gwasCols, args$v_values[1:(args$start_index-1)])
  }

  dat <- fread(args$input_file, sep = '\t', header = T, select = gwasCols)

  if(any(!is.element(gwasCols, names(dat)))) {
    stop(sprintf("data.table object is missing the following columns specified in the command-line arguments: %s\n", paste(gwasCols[!is.element(gwasCols, names(dat))], collapse = ', ')))
  }

#  if(args$start_index > 1 & any(!is.element(args$v_values[1:args$start_index], names(dat)))) {
#    missing_v_values <- paste((args$v_values[1:args$start_index])[!is.element(args$v_values[1:args$start_index], names(dat))], collapse = ', ')
#    stop(sprintf("v-values columns %s missing from input file, cannot restart at index %d", missing_v_values, args$start_index))
#  }

  dat <- dat[!(get(args$chromosome) == 6 & get(args$basepair) %between% c(24e6, 45e6))]

  # If ur-principal p-value is absent, can't do anything
  dat <- dat[!is.na(get(args$principal))]

  # Transforming p-values of 0 to Z-scores yields Inf values
  for(x in c(args$principal, args$auxiliary)) {
    dat[ get(x) < 1e-300, (x) := 1e-300]
  }

  prin_pvalues <- c(args$principal, args$v_values)

  for(i in args$start_index:length(args$auxiliary)) {
    # prin_pvalues[i] should never contain NA values as we drop the ur_prin_pvalues NA rows then carry over
    print(sprintf('Iteration %d', i))

    if(i > 1) {
      cols <- c(args$chromosome, args$basepair, args$principal, prin_pvalues[i], args$auxiliary[i], args$weights[i])
    } else {
      cols <- c(args$chromosome, args$basepair, args$principal, args$auxiliary[i], args$weights[i])
    }

    sub_dat <- dat[!is.na(get(args$auxiliary[i])), ..cols]

    cols <- c(args$principal, args$weights[i], args$auxiliary[i])

    q0_dat <- dat[ get(args$principal) > 0.5 & get(args$weights[i]) > 0, ..cols]

    est_q0_pars <- fit.2g(P = q0_dat[[args$auxiliary[i]]], weights = q0_dat[[args$weights[i]]])$pars

    rm(q0_dat)

    folds <- mclapply(unique(sub_dat[[args$chromosome]]), function(x) which(sub_dat[[args$chromosome]]==x), mc.cores = args$no_of_threads)

    candidate_indices <- sub_dat[get(args$principal) < args$p_threshold, which = T]

    # Organise the candidate indices by fold
    ind <- mclapply(folds, function(x) intersect(candidate_indices, x), mc.cores = args$no_of_threads)

    # Exclude ind, folds with no data points
    non_empty_indices <- ind[sapply(ind, function(x) length(x) > 0)]

    folds_with_indices <- folds[sapply(ind, function(x) length(x) > 0)]

    # Compute L-regions
    v <- mcmapply(function(x,y) vl(sub_dat[[prin_pvalues[i]]], sub_dat[[args$auxiliary[i]]], indices = x, mode = 2, fold = y), non_empty_indices, folds_with_indices, mc.cores = args$no_of_threads, SIMPLIFY = F)

    # il calls are fast enough not to justify their being parallelised
    # Integrate over L-regions to obtain v-values
    for(j in 1:length(non_empty_indices)) {
        sub_dat[non_empty_indices[[j]], (prin_pvalues[i+1]) := il(v[[j]], pi0_null = est_q0_pars[1], sigma_null = est_q0_pars[2], distx = "norm")];
    }

    # still assuming one SNP at most per basepair, so we can join on chrom and bp
    cols <- c(args$chromosome, args$basepair, prin_pvalues[i+1])

    dat <- merge(dat, sub_dat[, ..cols], all.x = T, suffixes = c('', ''))

    # Impute as previous iteration's p-/v-value and replace NA values produced by il (see notes.org)
    dat[is.na(get(prin_pvalues[i+1])), (prin_pvalues[i+1]) := get(prin_pvalues[i])]

    fwrite(dat, file = args$outputPath, sep  = '\t', row.names = F, col.names = T, quote = F)
  }
}
