#install.packages("doParallel")
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)


ptime <- system.time({
  foreach(i=1:10000) %do% sqrt(i)
})

ptime

ptime <- system.time({
foreach(i=1:10000) %dopar% sqrt(i)
})

ptime

