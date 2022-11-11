

class directories:

    def __init__(self,home_dir):
        self.git_dir='/QstarFromTidalSynchronization/MCMC/version2_emcee'
        self.binary_dir='/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code'
        if home_dir=='/home/ruskin':
            self.poet_path=home_dir+'/projects/poet'
            self.mcmc_directory=home_dir+'/projects'+self.git_dir
            self.scratch_directory=home_dir+'/projects'+self.binary_dir

        if home_dir=='/home/rxp163130':
            self.poet_path=home_dir+'/poet'
            self.mcmc_directory=home_dir+self.git_dir
            self.binary_directory=home_dir+self.binary_dir
            self.scratch_directory=home_dir+'/scratch'

        if home_dir=='/home1/06850/rpatel23':
            self.poet_path=home_dir+'/poet'
            self.mcmc_directory=home_dir+self.git_dir
            self.binary_directory=home_dir+self.binary_dir
            self.scratch_directory='/work/06850/rpatel23/ls6/scratch'
