


for S in 6312521 6359798 8746310 4947726 3241344 3973504 6525196 5393558 11200773 4579321 6949550 6927629 3838496 4839180 10965963 3834364 8580438; do
nohup python3 h5_analysis.py -s system_$S.h5 >err_output/output_period_$S.txt 2>&1 &
done


#for S in system_*.h5; do
#    nohup python3 h5_analysis.py -s $S error_output/output_period_$S.txt 2>&1 &
#done
