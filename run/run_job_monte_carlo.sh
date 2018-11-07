for i in {1..10}
  do msub -N true_covariance -q singlenode -l nodes=1:ppn=16 -l walltime=00:12:00:00 job_monte_carlo.sh
done
