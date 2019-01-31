Calculate.F1 <- function(true, learnt, skeleton = FALSE) {
  if(skeleton) {
    true <- skeleton(true)
    learnt <- skeleton(learnt)
  }
  comparison <- bnlearn::compare(true, learnt)
  precision <- comparison$tp / (comparison$tp + comparison$fp)
  recall <- comparison$tp / (comparison$tp + comparison$fn)
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  return(list('F1' = f1, 'Precision' = precision, 'Recall' = recall))
}

sensitivity.analysis <- function(bn.structure, bn.fit, algorithm, skeleton, algorithm.parameters, samples, tests) {
  scores <- c()
  for (i in 1:tests) {
    print(paste0('Test: ', i, ', Samples: ', samples))
    data.sample <- SynBN::generate.synthetic.dataset(bn.fit, sample.size = samples, seed = i)
    learnt.structure <- do.call(algorithm, append(list(data = data.sample), algorithm.parameters))
    scores <- c(scores, Calculate.F1(bn.structure, learnt.structure, skeleton = skeleton)$F1)
  }
  return(list(scores = scores, samples = rep(samples, tests)))
}

sensitivity.plot <- function(bn.structure, bn.fit, algorithm, skeleton = FALSE, algorithm.parameters = NULL,
                             samples=c(1000, 5000, 10000), tests = 10, algorithm.name = 'TEST') {
  library(ggplot2)
  scores <- c()
  sample.labels <- c()
  for(sample.size in samples) {
    results <- sensitivity.analysis(bn.structure, bn.fit, 
                                    algorithm, skeleton, 
                                    algorithm.parameters, 
                                    sample.size, tests)
    scores <- append(scores, results$scores)
    sample.labels <- append(sample.labels, results$samples)
  }
  df <- data.frame(sample.labels, scores)
  
  p <- ggplot(df, aes(x = sample.labels, y = scores, group = sample.labels)) +
    geom_boxplot(colour = "black", fill = "#56B4E9") +
    scale_y_continuous(name = "DAG F1 score") +
    xlab('Samples') +
    ggtitle(paste0("Sensitivity analysis of "), algorithm.name)
  print(p)
}
