

class directories:

    def __init__(self,home_dir):
        self.git_dir='/QstarFromTidalSynchronization/MCMC/version2_emcee'
        if home_dir=='/home/ruskin':
            self.poet_path=home_dir+'/projects/poet'
            self.mcmc_directory=home_dir+'/projects'+self.git_dir
            self.scratch_directory=home_dir+'/projects'+self.git_dir+'/h5analysis/err_output'

        if home_dir=='/home1/06850/rpatel23':
            self.poet_path=home_dir+'/poet'
            self.mcmc_directory=home_dir+self.git_dir
            self.scratch_directory='/work/06850/rpatel23/ls6/scratch'
