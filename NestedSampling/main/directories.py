

class directories:

    def __init__(self,home_dir):
        self.git_dir='/QstarFromTidalSynchronization/NestedSampling/main'
        if home_dir=='/home/ruskin':
            self.poet_path=home_dir+'/projects/poet'
            self.current_directory=home_dir+'/projects'+self.git_dir
            self.output_directory=self.current_directory+'/results/kartof'

        if home_dir=='/home/rxp163130':
            self.poet_path=home_dir+'/poet'
            self.current_directory=home_dir+self.git_dir
            self.output_directory=self.current_directory+'/results/ganymede'

        if home_dir=='/home1/06850/rpatel23':
            self.work_dir='/work/06850/rpatel23/stampede2'
            self.poet_path=self.work_dir+'/poet'
            self.current_directory=self.work_dir+self.git_dir
            self.output_directory=self.current_directory+'/results/stampede'


