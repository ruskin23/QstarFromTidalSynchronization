SYSTEM=$1
DIR=$2

SLURM_DIR=/home/$DIR/QstarFromTidalSynchronization/MCMC/version2_emcee/slurm_scripts

sed\
    -e 's%@@SYSTEM@@%'"$SYSTEM"'%g'\
    ${SLURM_DIR}/template.temp\
    >\
    ${SLURM_DIR}/scripts/${SYSTEM}.slurm
