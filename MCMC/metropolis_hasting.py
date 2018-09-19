
import sys
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/PythonPackage')
sys.path.append('/Users/ruskinpatel/Desktop/Research/poet/scripts')
sys.path.append('/Users/ruskinpatel/Desktop/Research/QstarFromTidalSynchronization/')


from binary_evolution.binary_evolution_class import binary_evolution
from Calculate_Mass.Deriving_Mass import DerivePrimnaryMass
from Calculate_Mass.Deriving_Mass import DeriveDeriveSecondaryMassMass

import scipy
from scipy.stats import norm
import numpy


class MetropolisHastings:


    def posterior_probability(self,parameter_set):
        """reutrns current surface spin of the star given the parameters"""

        evolution = binary_evolution(
                                    interploator = self.interpolator,
                                    observational_parameters = parameter_set,
                                    fixed_parameters = self.fixed_parameters)

        spin_given_model = evolution()

        prior = 1.0

        for (key_obs,value_obs),(key_parameter,value_parameter) in zip(self.observation_data.items(),parameter_set.items()):
            prior  *= scipy.stats.norm(value_obs['value'],value_obs['sigma']).pdf(value_parameter)

        likelihood = scipy.stats.norm(self.current_spin['value'],self.current_spin['sigma']).pdf(spin_given_model)

        posterior = prior*likelihood

        return posterior

    def acceptance_probability(self):

        posterior_proposed = self.posterior_probability(self.proposed_parameters)
        posterior_previous = self.posterior_probability(self.updated_parameters)
        p = posterior_proposed/posterior_previous

        return p


    def values_proposed(self):

        #propose values for observables
        for (name_obs,value_obs),(name_step,value_step) in zip(self.updated_parameters.items(),self.proposed_step.items()):
            proposed[name_obs]=scipy.stats.norm.rvs(loc=value_obs, scale=value_step)

        #propose value for logQ


    def write_file(self):

        data = []
        data.append(self.iteration_step)
        for value in self.updated_parameters.values():
            data.append(value)


        if self.isAccepted == True:
            f = open('accepted_aparemters.txt')

        else:
            f = open('rejected_parameters.txt')

        #f.write()
        f.close()


    def iterations(self):

        """runs the metropolis hastrings algorithm for number of iterations given"""

        file1 = open('accepted_aparemters.txt')
        file2 = open('rejected_parameters.txt'
                     )
        header = []
        header.append('Iteration_step')
        for key in self.updated_parameters.keys():
            header.append(key)


        #file1.write(header)
        #file2.write(header)

        file.close()

        if self.iteration_step==0:
            self.updated_parameters = self.initial_parameters


        while self.iteration_step < self.total_iterations:

            #draw a random value from proposal function
            self.proposed_parameters = self.values_proposed()

            # calculate acceptance probablity
            p_acceptance = acceptance_probability()

            if p_acceptance > 1:
                self.updated_parameters = self.proposed_parameters
                self.isAccepted = True
            else:
                rand = scipy.stats.norm.rvs()
                if p_acceptance>rand: self.isAccepted=True
                else:self.isAccepted=False

            self.iteration_step = self.iteration_step + 1

        return none


    def __init__(
                self,
                interpolator,
                fixed_parameters,
                observation_data,
                logQ,
                proposed_step,
                total_iterations,
                current_spin
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
            self.initial_parameters[key] = scipy.stats.norm.rvs(loc=key['value'], scale=key['sigma'])
        self.isAccepted = None
        self.current_spin = current_spin


if __name__ == '__main__':

    serialized_dir = "/Users/ruskinpatel/Desktop/Research/poet/stellar_evolution_interpolators"
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    observation_data = dict(
                        age = dict(vlaue = 1.0, sigma = 0.1),
                        teff = dict(vlaue = 1.0, sigma = 0.1),
                        Porb = dict(value = 1.0, sigma =0.1),
                        feh = dict(value = 1.0, sigma = 0.1),
                        rvk = dict(value = 1.0, sigma = 0.1),
                        inclination = dict(value = 1.0, sigma = 0.1),
                    )

    current_spin = dict(
                        value = 1.0,
                        sigma = 0.1
                        )

    fixed_parameters = dict(
                        disk_dissipation_age = 1.0,
                        disk_lock_frequency = 1.0,
                        wind = True,
                        wind_saturation_frequency = 1.0,
                        diff_rot_coupling_timescale = 1.0
                    )

    proposed_step = dict(
                        age_step = 0.1,
                        teff_step = 0.1,
                        Porb_step = 0.1,
                        feh_step = 0.1,
                        rvk_step = 0.1,
                        inclination = 0.1,
                    )


    logQ = dict(
                min = 4,
                max = 6,
                step = 0.1)




#flush()
#buffer_size =0 to write at the moment
#