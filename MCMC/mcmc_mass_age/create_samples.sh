for s in 111 128 132 18 27 33 49 5 62 68 99 
do
    nohup python3 main.py -l $s >/mnt/md0/ruskin/QstarFromTidalSynchronization/mcmc_mass_output/sampling_more_points/output/output_1_$s.txt 2>&1 &
done

