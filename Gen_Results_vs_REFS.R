library(readr)

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
  res <- tryCatch({
    labels <- c(paste0(label, ' Precision'), paste0(label, ' Recall'), paste0(label, ' F1'), 
                paste0(label, ' V Precision'), paste0(label, ' V Recall'), paste0(label, ' V F1'),
                paste0(label, ' DAG Precision'), paste0(label, ' DAG Recall'), paste0(label, ' DAG F1'))
    
    v.results <-  v.structs.match(learnt, true)
    skeleton.results <- bnlearn::compare(bnlearn::skeleton(true),
                                         bnlearn::skeleton(learnt))
    dag.results <- bnlearn::compare(true, learnt)
    
    skeleton.precision <- skeleton.results$tp / (skeleton.results$tp + skeleton.results$fp)
    skeleton.recall <- skeleton.results$tp / (skeleton.results$tp + skeleton.results$fn)
    skeleton.f1 <- 2 * (skeleton.precision * skeleton.recall) / (skeleton.precision + skeleton.recall)
    
    v.f1 <- 2 * (v.results$precision * v.results$recall) / (v.results$precision + v.results$recall)
    
    dag.precision <- dag.results$tp / (dag.results$tp + dag.results$fp)
    dag.recall <- dag.results$tp / (dag.results$tp + dag.results$fn)
    dag.f1 <- 2 * (dag.precision * dag.recall) / (dag.precision + dag.recall)
    
    c(skeleton.precision, skeleton.recall, skeleton.f1,
      v.results$precision, v.results$recall, v.f1, dag.precision, dag.recall, dag.f1)
    
  }, error = function(e){
    return(rep('ERROR', 12))
  })
  
  names(res) <- labels
  
  return(res)
}

gen.results <- function(results.filename, subsample.directory, modelstring.directory, 
                        algorithm, algorithm.name, algorithm.parameters, result.savedir) {
  #Load Results
  results.df <- read.csv(results.filename)
  
  struct.conv <- list('forest fire' = 'ForestFire', 'ic-dag' = 'icDAG', 'preferential attachment' = 'PrefAttach')
  
  results <- NULL
  
  start.idx <- 1
  
  if(file_test("-f", result.savedir)) {
    results <- read.csv(result.savedir)
    start.idx <- nrow(results) + 1
    print(paste0('Resuming from experiment #', start.idx))
  }
  
  #Iterate over rows
  for(i in start.idx:nrow(results.df)) {
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
    n.params <- bnlearn::nparams(true.model, data)
    
    #Get Refs learnt model
    modelstring.file.start <- paste0(modelstring.directory, '/REFS_', row$input_structuretype,
                                     '_node', row$input_nodes, '_sample', row$input_samples, '_*')
    modelstring.file <- Sys.glob(modelstring.file.start)[1]
    refs.learnt <- model2network(strsplit(read_file(modelstring.file), '[\\\\]|[^[:print:]]')[[1]][1])
    
    #Calculate REFS results
    refs.results <- get.results.wvstructs(true.model, refs.learnt, 'REFS')
    
    #Get Algorithm Learnt Model
    algorithm.learnt <- do.call(algorithm, append(list(data = data), algorithm.parameters))
    
    #Get algorithm results
    algorithm.results <- get.results.wvstructs(true.model, algorithm.learnt, algorithm.name)
    
    #Store results in row
    results.row <- c('Structure Type' = as.character(names(struct.conv)[which(struct.conv %in% row$input_structure)]), 'Number of Samples' = row$input_samples,
                     'Number of Nodes' = row$input_nodes, 'N.Params' = n.params, refs.results, algorithm.results)
    
    results <- rbind(results, results.row)
    write.csv(results, result.savedir, row.names = NULL)
  }
  return(results)
}

results <- gen.results('~/Desktop/REFS_Results/notes_20181121_01.csv',
                       '~/Desktop/REFS_Results/data/',
                       '~/Desktop/REFS_Results/',
                       dag.threshold.ensemble,
                       'Clustered.Fragments.PCTABU',
                       list(n = 20),
                       'Clustered_Fragments_PCTABU.csv')