###############################
## 01 Download PlanetScope
for(numSite in 1){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 run_script_01.sh ',numSite,sep=''))
}

