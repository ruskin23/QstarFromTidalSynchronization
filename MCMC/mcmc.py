
import sys
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')


from binary_evolution_class import evolution
from stellar_evolution.manager import StellarEvolutionManager


import scipy
from scipy.stats import norm
import numpy


class MetropolisHastings:


    def posterior_probability(self,parameter_set=None):
        """reutrns current surface spin of the star given the parameters"""

        find_spin = evolution(
                                    self.interpolator,
                                    parameter_set,
                                    self.fixed_parameters)

        Pspin_given_model = find_spin()

        prior = 1.0

        for (key_obs,value_obs),(key_parameter,value_parameter) in zip(self.observation_data.items(),parameter_set.items()):
            prior  *= scipy.stats.norm(value_obs['value'],value_obs['sigma']).pdf(value_parameter)

        likelihood = scipy.stats.norm(self.observed_Pspin['value'],self.observed_Pspin['sigma']).pdf(Pspin_given_model)

        posterior = prior*likelihood

        return posterior

    def acceptance_probability(self):


        posterior_proposed = self.posterior_probability(parameter_set=self.proposed_parameters)
        posterior_previous = self.posterior_probability(parameter_set=self.updated_parameters)
        p = posterior_proposed/posterior_previous

        return p


    def values_proposed(self):

        proposed = dict()
        for (name_obs,value_obs),(name_step,value_step) in zip(self.updated_parameters.items(),self.proposed_step.items()):
            proposed[name_obs]=scipy.stats.norm.rvs(loc=value_obs, scale=value_step)

        return proposed


    def write_on_file(self):

        # writes all the parameter names in the file

        def write_header():

            filename = ['accepted_aparemters.txt', 'rejected_parameters.txt']

            header = []
            header.append('Iteration_step')

            for key in self.updated_parameters.keys():
                header.append(key)

            for i in range(0, 1):
                with open(filename[i], 'w') as file:
                    for name in header:
                        file.write('%s ' % name)
                file.close()

        if self.iteration_step==0: write_header()

        if self.isAccepted == True: filename = 'accepted_aparemters.txt'
        else: filename = 'rejected_parameters.txt'

        with open(filename, 'w') as file:
            for key, value in self.updated_parameters.items():
                file.write('%s ' % value)
        file.close()


    def iterations(self):

        """runs the metropolis hastrings algorithm for number of iterations given"""
        if self.iteration_step==0: self.updated_parameters = self.initial_parameters
        print ('\nINITIAL PARAMETERS SET')

        while self.iteration_step<self.total_iterations:

            print ('ITERATION STEP = ',self.iteration_step)

            #draw a random value from proposal function
            self.proposed_parameters = self.values_proposed()
            print ('new values proposed')

            # calculate acceptance probablity
            print ('calculating accepatanece probability')
            p_acceptance = self.acceptance_probability()
            print ('calculated accepatanece probability')


            if p_acceptance > 1:
                self.updated_parameters = self.proposed_parameters
                self.isAccepted = True
            else:
                rand = scipy.stats.norm.rvs()
                if p_acceptance>rand : self.isAccepted=True
                else : self.isAccepted=False

            self.write_on_file()
            self.iteration_step = self.iteration_step + 1


    def __init__(
                self,
                interpolator,
                fixed_parameters,
                observation_data,
                logQ,
                proposed_step,
                total_iterations,
                observed_Pspin
                ):


        self.interpolator  = interpolator
        self.fixed_parameters = fixed_parameters
        self.observation_data = observation_data
        self.logQ = logQ
        self.proposed_step = proposed_step
        self.iteration_step = 0
        self.total_iterations= total_iterations
        self.updated_parameters = dict()
        self.proposed_parameters = dict()



        self.initial_parameters = dict()
        for key, value in self.observation_data.items():
            self.initial_parameters[key] = scipy.stats.norm.rvs(loc=value['value'], scale=value['sigma'])
        self.initial_parameters['logQ'] = numpy.random.uniform(low=self.logQ['min'],high=self.logQ['max'],size=None)

        self.isAccepted = None
        self.observed_Pspin = observed_Pspin


if __name__ == '__main__':

    serialized_dir = "/Users/ruskinpatel/Desktop/Research/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    observation_data = dict(
                        age = dict(value = 5.0, sigma = 0.1),
                        teff = dict(value = 5777.0, sigma = 50),
                        feh = dict(value = 0.0, sigma = 0.1),
                        rvk = dict(value = 615, sigma = 10),
                        Porb = dict(value=5.0, sigma=1e-5),
                        Pdisk = dict(value=2*scipy.pi / 1.4, sigma=0.1)
                    )

    observed_Pspin = dict(
                        value = 5.0,
                        sigma = 0.5
                    )

    fixed_parameters = dict(
                        disk_dissipation_age = 5e-3,
                        planet_formation_age = 5e-3,
                        wind = True,
                        wind_saturation_frequency = 2.54,
                        diff_rot_coupling_timescale = 5e-3,
                        wind_strength = 0.17,
                        inclination=scipy.pi/2

    )

    proposed_step = dict(
                        age_step = 0.1,
                        teff_step = 50,
                        feh_step = 0.1,
                        rvk_step = 10,
                        Porb_step=0.1,
                        Pdisk_step = 0.1,
                        logQ_step = 0.5
                    )


    logQ = dict(
                min = 4,
                max = 6
            )


mcmc = MetropolisHastings(interpolator,fixed_parameters,observation_data,logQ,proposed_step,1,observed_Pspin)

mcmc.iterations()
#flush()
#buffer_size =0 to write at the moment
#