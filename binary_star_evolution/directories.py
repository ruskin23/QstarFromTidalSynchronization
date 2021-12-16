

class directories:

    def __init__(self,home_dir):
        self.git_dir='/QstarFromTidalSynchronization/binary_star_evolution/updated_evolution_code'
        if home_dir=='/home/ruskin':
            self.poet_path=home_dir+'/projects/poet'
            self.current_directory=home_dir+'/projects'+self.git_dir
            self.results_directory=self.current_directory+'/results/kartof'
            self.scratch_directory=self.current_directory+'/output/kartof'

        if home_dir=='/home/rxp163130':
            self.poet_path=home_dir+'/poet'
            self.current_directory=home_dir+self.git_dir
            self.results_directory=self.current_directory+'/results/ganymede'
            self.scratch_directory=home_dir+'/scratch/binary_evolution'

        if home_dir=='/home1/06850/rpatel23':
            self.work_dir='/work/06850/rpatel23/stampede2'
            self.poet_path=self.work_dir+'/poet'
            self.current_directory=self.work_dir+self.git_dir
            self.results_directory=self.current_directory+'/results/stampede'
            self.scratch_directory='/scratch/06850/rpatel23/Nested_Sampling'


