SYSTEM=$1


SLURM_DIR=/home/ruskin/projects/QstarFromTidalSynchronization/MCMC/version2_emcee/slurm_scripts

sed\
    -e 's%@@SYSTEM@@%'"$SYSTEM"'%g'\
    ${SLURM_DIR}/template.temp\
    >\
    ${SLURM_DIR}/scripts/${SYSTEM}.slurm
