################ IMPORT PARALLEL LIBRARIES #####################
library('doParallel')
library('yaml')

############### SOURCE LOCAL FILES NECCESARY ###################
source('SL_Metrics/Metrics.R')
source('SL_Metrics/Properties.R')

############### LOAD PARAMETER FILE #################
load.search.grid <- function(parameter.filename) {
  # TODO: IMPLEMENT
  parsed.parameters <- read_yaml(parameter.filename)
  search.grid <- expand.grid(parsed.parameters$alpha_values,
                             parsed.parameters$complexity,
                             parsed.parameters$average_levels,
                             parsed.parameters$structure_types,
                             parsed.parameters$trial_seeds,
                             parsed.parameters$number_of_samples,
                             parsed.parameters$number_of_variables
                             )
  workers <- parsed.parameters$workers
  return(list('search.grid' = search.grid, 'workers' = workers))
}

algorithm.results <- function(algorithm, true.dag, true.fit, clean.data, sink.node, algorithm.params) {
  # Get start time to measure delta after algorithm 
  start.time <- Sys.time()
  # Call the passed algorithm, with the passed algorithm parameters
  learnt.dag <- do.call(algorithm, algorithm.params) 
  # Calculate the algorithm time taken
  time.taken <- difftime(Sys.time(), start.time, units = "secs")[[1]]
  # Calculate Skeleton Results
  skeleton.results <- Calculate.F1(true = true.dag, learnt = learnt.dag, skeleton = TRUE)
  
  # Calculate DAG Results
  dag.results <- Calculate.F1(true = true.dag, learnt = learnt.dag, skeleton = FALSE)
  
  return(list('Precision' = precision, 'Recall' = recall, 'V.Precision' = v.precision, 
              'V.Recall' = v.recall, 'OR.Correct' = OR.Correct, 'Time.Taken' = time.taken))
}

try.fit <- function(structure, seed, alpha, max.levels) {
  fit <- tryCatch({
    fit <- withTimeout({
      generate.parameters(structure, seed = seed, alpha = alpha, max.levels = max.levels)
    }, timeout= 300); 
  }, TimeoutException = function(ex) {
    return('Gen Timeout')
  })
  return(fit)
}


experiment <- function(algorithm, nodes, alpha, struct.type, trial, samples, complexity, method.args, levels, timeout = 600) {
  
  # --------------------------------------- GROUND TRUTH GENERATION --------------------------------------- #
  
  structure <- create.structure(algorithm, nodes, seed = trial, method = struct.type, method.args = method.args)
  sink.node <- identify.sink(structure$structure)
  
  parent.num <- c()
  for (node in structure$structure$nodes) {
    parent.num <- c(parent.num, length(node$parents))
  }
  median.edges <-  median(parent.num)
  
  fit <- try.fit(structure, seed = trial, alpha = alpha, max.levels = levels)
  
  if(fit == 'Gen Timeout') {
    return(list('Structure Type' = struct.type, 'Complexity' = complexity, 'Median Edges' = median.edges, 'Alpha' = alpha, 
                'Number of Nodes' = nodes, 'Number of Samples' = samples, 'Max.Levels' = levels, 'N.Params' = 'Gen Timeout',
                'degree.mean' = get.degree.mean(structure$structure), 'degree.variance' = get.degree.variance(structure$structure),
                'Trial Number' = trial, 'PCTABU Precision' = 'Gen Timeout', 'PCTABU Recall' = 'Gen Timeout', 
                'PCTABU V Precision' = 'Gen Timeout', 'PCTABU V Recall' = 'Gen Timeout', 
                'PCTABU OR Correct' = 'Gen Timeout', 'PCTABU Time Taken' = 'Gen Timeout'))
  }
  
  clean.data <- generate.synthetic.dataset(fit, sample.size = samples, seed = trial)
  
  # --------------------------------------- END GROUND TRUTH GENERATION --------------------------------------- #
  
  alg.results <- tryCatch({
    alg.results <- withTimeout({
      algorithm.results(algorithm, structure, fit, clean.data, sink.node, alpha = NULL)
    }, timeout= timeout); 
  }, TimeoutException = function(ex) {
    list('Precision' = 'Timeout', 'Recall' = 'Timeout', 'V.Precision' = 'Timeout', 'V.Recall' = 'Timeout', 'OR.Correct' = 'Timeout', 'Time.Taken' = 'Timeout')
  })
  
  return(list('Structure Type' = struct.type, 'Complexity' = complexity, 'Median Edges' = median.edges, 'Alpha' = alpha, 
              'Number of Nodes' = nodes, 'Number of Samples' = samples, 'Max.Levels' = levels, 'N.Params' = bnlearn::nparams(fit),
              'degree.mean' = get.degree.mean(structure$structure), 'degree.variance' = get.degree.variance(structure$structure),
              'Trial Number' = trial, alg.results))
  
}


perform.experiments <- function(algorithm, parameter.filename) {
  parameters <- load.search.grid(parameter.filename)
  search.grid <- parameters$search.grid
  
  print(paste0('Performing ', nrow(search.grid), ' experiments.'))
  
  #Multi-Processing stuff
  cores <- detectCores()
  cluster <- makeCluster(parameters$workers, outfile = "")
  registerDoParallel(cluster)
  options(width = 250, cores = parameters$workers)
  
  alpha.df <- foreach (i = 1:nrow(search.grid),
                       .combine = rbind,
                       .packages=c('SynBN', 'bnlearn', 'pcalg', 'R.utils'),
                       .export = c('experiment', 'try.fit', 'algorithm.results')) %dopar% {
    alpha <- search.grid[i, 1]
    compl <- search.grid[i, 2]
    levels <- search.grid[i, 3]
    struct.type <- search.grid[i, 4]
    trial <- search.grid[i, 5]
    samples <- search.grid[i, 6]
    nodes <- search.grid[i, 7]
    exp.results <- experiment(algorithm, nodes, alpha, struct.type, trial, samples, 1, NULL, levels)
    exp.results
  }
  colnames(alpha.df) <- colnames(df)
  df <- rbind(df, alpha.df)
  rownames(df) <- 1:nrow(df)
  print(tail(df, 10))
  write.csv(df, file = output.filename, row.names = F)
  stopCluster(cluster)
}