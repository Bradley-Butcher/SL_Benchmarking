library(doParallel)

v.structs.match <- function(comparison.structure, true.structure) {
  #True V structure in a dataframe
  true.v.structs <- data.frame(vstructs(cpdag(true.structure, moral = TRUE), moral = TRUE))
  
  #V structure to compare in a dataframe
  comp.v.structs <- data.frame(vstructs(cpdag(comparison.structure)))
  
  n.v.struct.true <- nrow(true.v.structs)
  
  n.v.struct.comp.true <- nrow(merge(comp.v.structs, true.v.structs))
  
  n.v.struct.comp <- nrow(comp.v.structs)
  
  precision <- if (n.v.struct.comp > 0) n.v.struct.comp.true / n.v.struct.comp else 0
  recall <- if (n.v.struct.true > 0) n.v.struct.comp.true / n.v.struct.true else 0
  
  return(list('precision' = precision,
              'recall' = recall))
}

get.results.wvstructs <- function(true, learnt, label) {
    labels <- c(paste0(label, ' Precision'), paste0(label, ' Recall'), paste0(label, ' F1'), 
                paste0(label, ' V Precision'), paste0(label, ' V Recall'), paste0(label, ' V F1'),
                paste0(label, ' DAG Precision'), paste0(label, ' DAG Recall'), paste0(label, ' DAG F1'))
    
    v.results <-  v.structs.match(learnt, true)
    skeleton.results <- compare(skeleton(true),
                                         skeleton(learnt))
    dag.results <- compare(true, learnt)
    
    skeleton.precision <- skeleton.results$tp / (skeleton.results$tp + skeleton.results$fp)
    skeleton.recall <- skeleton.results$tp / (skeleton.results$tp + skeleton.results$fn)
    skeleton.f1 <- 2 * (skeleton.precision * skeleton.recall) / (skeleton.precision + skeleton.recall)
    
    v.f1 <- 2 * (v.results$precision * v.results$recall) / (v.results$precision + v.results$recall)
    
    dag.precision <- dag.results$tp / (dag.results$tp + dag.results$fp)
    dag.recall <- dag.results$tp / (dag.results$tp + dag.results$fn)
    dag.f1 <- 2 * (dag.precision * dag.recall) / (dag.precision + dag.recall)
    
    res <- c(skeleton.precision, skeleton.recall, skeleton.f1,
      v.results$precision, v.results$recall, v.f1, 
      dag.precision, dag.recall, dag.f1)
    
  
  names(res) <- labels
  
  return(res)
}

get.ensemble.results <- function(true, learnt.models, label) {
  results <- c()
  n <- length(learnt.models)
  for(learnt.model in learnt.models) {
    single.res <- get.results.wvstructs(true, model2network(learnt.model), label)
    results <- c(results, unlist(single.res))
  }
  summarized.results <- c('mean.skeleton.precision' = mean(results[(1:n * 9) - 8]),
              'std.skeleton.precision' = sd(results[(1:n * 9) - 8]),
              'mean.skeleton.recall' = mean(results[(1:n * 9) - 7]),
              'std.skeleton.recall' = sd(results[(1:n * 9) - 7]),
              'mean.skeleton.f1' = mean(results[(1:n * 9) - 6]),
              'std.skeleton.f1' = sd(results[(1:n * 9) - 6]),
              'mean.v.precision' = mean(results[(1:n * 9) - 5]),
              'std.v.precision' = sd(results[(1:n * 9) - 5]),
              'mean.v.recall' = mean(results[(1:n * 9) - 4]),
              'std.v.recall' = sd(results[(1:n * 9) - 4]),
              'mean.v.f1' = mean(results[(1:n * 9) - 3]),
              'std.v.f1' = sd(results[(1:n * 9) - 3]),
              'mean.dag.precision' = mean(results[(1:n * 9) - 2]),
              'std.dag.precision' = sd(results[(1:n * 9) - 2]),
              'mean.dag.recall' = mean(results[(1:n * 9) - 1]),
              'std.dag.recall' = sd(results[(1:n * 9) - 1]),
              'mean.dag.f1' = mean(results[(1:n * 9)]),
              'std.dag.f1' = sd(results[(1:n * 9)])
              )
  names(summarized.results) <- paste0(names(summarized.results), ' ', label)
  return(summarized.results)
}

gen.results <- function(results.filename, subsample.directory, modelstring.directory, 
                        algorithm, algorithm.name, algorithm.parameters, algorithm.source,
                        result.savedir, n.workers = 5, ensemble = FALSE) {
  #Load Results
  results.df <- read.csv(results.filename)
  
  struct.conv <- list('forest fire' = 'ForestFire', 'ic-dag' = 'icDAG', 'preferential attachment' = 'PrefAttach')
  
  results <- NULL
  experiments.todo <- NULL
  
  if(file_test("-f", result.savedir)) {
    results <- read.csv(result.savedir)
    performed.experiments <- results$Experiment
    experiments.todo <- setdiff(1:nrow(results.df), performed.experiments)
    print('Resuming Experiment')
  } else {
    experiments.todo <- 1:nrow(results.df)
  }
  
  cluster <<- makeCluster(n.workers, outfile='')
  registerDoParallel(cluster)
  on.exit(print('Stopping Cluster...'))
  on.exit(stopCluster(cluster))
  
  #Iterate over rows
  final.results <- foreach(i=experiments.todo,
                           .export=c('get.results.wvstructs', 'get.ensemble.results', 'v.structs.match'),
                           .packages=c('bnlearn', 'readr'), 
                           .combine=rbind) %dopar% {
    source(algorithm.source, local = TRUE)
    print(paste0('Starting experiment ', i))
    row <- results.df[i, ]
    #Load True Structure
    true.string <- row$truth_modelstring
    true.model <- model2network(as.character(true.string))
    
    #Get subsample data
    datafile <- row$input_datafile
    datafile <- tail(strsplit(as.character(datafile), split="/", fixed=TRUE)[[1]], 1)
    datafile <- paste0(subsample.directory, '/', row$input_structuretype, '/', datafile)
    datafile <- paste0(strsplit(datafile, '\\.')[[1]][1], '_sussex.csv')
    data <- read.csv(datafile)
    data[] <- lapply(data, as.factor)
    
    #Get Number of Parameters
    n.params <- nparams(true.model, data)
    
    #Get Refs learnt model
    modelstring.file.start <- paste0(modelstring.directory, '/REFS_', row$input_structuretype,
                                     '_node', row$input_nodes, '_sample', row$input_samples, '_*')
    
    modelstring.file <- Sys.glob(modelstring.file.start)[1]
    refs.learnt <- model2network(strsplit(read_file(modelstring.file), '[\\\\]|[^[:print:]]')[[1]][1])
    
    #Calculate REFS results
    refs.results <- get.results.wvstructs(true.model, refs.learnt, 'REFS')
    
    #Get Algorithm Learnt Model
    algorithm.learnt <- do.call(algorithm, append(list(data = data), algorithm.parameters))
    
    algorithm.results <- NULL
    ensemble.results <- NULL
    #Get algorithm result
    if(!ensemble) {
      algorithm.results <- get.results.wvstructs(true.model, algorithm.learnt, algorithm.name)
    } else {
      algorithm.results <- get.results.wvstructs(true.model, algorithm.learnt$consensus.dag, algorithm.name)
      ensemble.results <- get.ensemble.results(true.model, algorithm.learnt$ensemble.models, algorithm.name)
    }
    
    #Store results in row
    results.row <- c('Experiment' = i, 'Structure Type' = as.character(row$input_structuretype), 'Number of Samples' = row$input_samples,
                     'Number of Nodes' = row$input_nodes, 'N.Params' = n.params, refs.results, algorithm.results)
    
    if(ensemble) {
      results.row <- c(results.row, ensemble.results)
    }

    if(!file_test("-f", result.savedir)) {
      print('Saving New CSV')
      cat.conn <- file(result.savedir)
      cat(names(results.row), file = cat.conn, sep = ",")
      close(cat.conn)
      
      cat.conn <- file(result.savedir, open = "a")
      cat('\n', file = cat.conn)
      cat(unname(results.row), file = cat.conn, sep = ",")
      close(cat.conn)
      
    } else {
      print('Appending to CSV')
      cat.conn <- file(result.savedir, open = "a")
      cat('\n', file = cat.conn)
      cat(unname(results.row), file = cat.conn, sep = ",")
      close(cat.conn)
    }
    results.row
  }
  write.csv(final.results, result.savedir, row.names = FALSE)
  stopCluster(cluster)
  
  return(results)
}

####################### EXAMPLE USAGE ############################

setwd('/home/surgo')
source('SL_Methods/Cluster_PCTABU.R')

results <- gen.results('SL_Benchmarking/REFS_Results/notes_20181121_01.csv',
                       'SL_Benchmarking/REFS_Results/data/',
                       'SL_Benchmarking/REFS_Results/',
                       cluster.pctabu.ensemble,
                       'Clustered.Fragments.PCTABU',
                       list(n = 20),
                       'SL_Methods/Cluster_PCTABU.R',
                       'Clustered_Fragments_PCTABU_ENS.csv',
                       n.workers = 10,
                       ensemble = TRUE)