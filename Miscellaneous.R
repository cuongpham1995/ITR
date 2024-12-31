summarize_data <- function(data, filter_condition) {
  data %>%
    filter(!!filter_condition) %>%
    select(
      gender, age, bmi, time_from_diagnosis, levodopa_dose, A, updrs_1, updrs_2, updrs_3,
      dys, fluc, nausea, sleep, ortho, Time_since_first_Levo, outcome_updrs3
    ) %>%
    summarise(
      Gender = n.pct(gender),
      BMI = mean.sd(bmi),
      Age = mean.sd(age),
      time_from_diagnosis = mean.sd(time_from_diagnosis),
      Levodopa_dosage = mean.sd(levodopa_dose),
      Levodopa_duration = mean.sd(Time_since_first_Levo),
      DRA = n.pct(A),
      updrs_1 = mean.sd(updrs_1),
      updrs_2 = mean.sd(updrs_2),
      updrs_3 = mean.sd(updrs_3),
      updrs_3_outcome = mean.sd(outcome_updrs3),
      Dyskinesia = n.pct(dys),
      Fluctuation = n.pct(fluc),
      Nausea = n.pct(nausea),
      insomnia = n.pct(sleep),
      orthostasis = n.pct(ortho)
    )
}


n.pct = function(x, dec.pt = 1){
  res = paste( format(round(sum(x == 1),dec.pt), nsmall = 1) , " (", format( round( mean(x == 1)*100 ,dec.pt), nsmall = 1) , "%",")", sep = "")
  return(res)
}

mean.sd = function(x, dec.pt = 1){
  res = paste( format( round(mean(x),dec.pt), nsmall = 1) , " (", format( round(sd(x),dec.pt), nsmall = 1 ) , ")", sep = "")
  return(res)
}

range.func = function(x, dec.pt = 1){
  x.range = round(range(x), dec.pt)
  res = paste( "(", format(x.range[1], nsmall = 1) , ", ", format(x.range[2], nsmall = 1) , ")", sep = "" )
  return(res)
}


calculate.statistics = function(all_results, FUN){
  all_results = all_results[sapply(all_results, FUN = function(x) length(x) == length(all_results[[1]]))]
  average_return = lapply(2:length(all_results[[1]]), function(j){
  component_values <- lapply(all_results, function(result) result[[j]])
  mat_comb = do.call(rbind, component_values)
    return(apply(mat_comb,2, FUN, na.rm = T))
    
  })
  return(average_return)
}

