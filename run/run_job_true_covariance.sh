#msub -N true_covariance -q fat -l nodes=1:ppn=16 -l walltime=02:00:00:00 -l pmem=64000mb job_true_covariance.sh
msub -N true_covariance -l nodes=1:ppn=4 -l walltime=02:00:00:00 -l pmem=10000mb job_true_covariance.sh
