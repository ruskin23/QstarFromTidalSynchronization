rm convergence_test.txt
touch convergence_test.txt
echo $'system_kic   quantile    r   quantile_value  burn-in total_steps thin' >convergence_test.txt
for S in *.h5; do python3 h5_analysis.py -s $S -c convergence_test.txt; done

