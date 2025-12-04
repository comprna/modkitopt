
-- dependencies
  -- If in an HPC environment, these can be loaded via modules

-- installation
-- If pipeline is interrupted, run again with -resume (Nextflow automatically checkpoints intermediate results and supports resuming interrupted runs)

-- How to run
  -- Include profile (local, pbs, slurm) via -profile flag

-- Explain how threads work (for pbs, slurm then 8 threads allocated per modkit pileup call)
  -- Automatically handled by nextflow