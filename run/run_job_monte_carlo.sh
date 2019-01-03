for i in {1..50}
  do msub -N true_covariance -q singlenode -l nodes=1:ppn=16 -l walltime=03:00:00:00 job_monte_carlo.sh
done
