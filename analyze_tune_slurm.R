#!/usr/bin/Rscript

#===============================================================================
#   This script launches a series of jobs on a Slurm cluster such as UPPMAX
#-------------------------------------------------------------------------------

# Specify which folds to run in each job
folds <- list(1:9, 10:17, 18:25)


#-------------------------------o
#   Launch!

dir.create("runcontrol", showWarnings=FALSE)
for(i in seq_along(folds)){
    batch.script <- sprintf(
"#! /bin/bash -l
#SBATCH -A b2010028
#SBATCH -p node -N 1
#SBATCH -t 9:00:00
#SBATCH --qos=b2010028_4nodes
#SBATCH -J subtypes
#SBATCH --output=runcontrol/analyze_tune_%i.out
#SBATCH --error=runcontrol/analyze_tune_%i.err

R -f analyze_tune.R --args %s",
    i, i, paste(folds[[i]], collapse=" "))
    
    batch.file <- sprintf("runcontrol/analyze_tune_%i.sh", i)
    cat(batch.script, file=batch.file)
    system(paste("sbatch", batch.file))
}

