


for S in 10031409 10215422 10330495 10992733 11147276 11200773; do
#for S in 11232745 11704044 12004679 4285087 4352168 4678171 4773155 4839180 5652260 5802470 5871918 6029130 6927629 7987749 8356054 8543278 8984706 9353182 9532123 9971475; do
    nohup python3 h5_analysis.py -s system_$S.h5 >error_output/output_period_$S.txt 2>&1 &
done


#for S in system_*.h5; do
#    nohup python3 h5_analysis.py -s $S error_output/output_period_$S.txt 2>&1 &
#done
