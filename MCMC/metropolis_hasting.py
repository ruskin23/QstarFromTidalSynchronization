
class MetropolisHastings:


    def posterior_probability(self,parameter_set):
        """reutrns current surface spin of the star given the parameters"""

        evolution = binary_evolution(interploator = self.interpolator, parameters = parameter_set,fixed_parameters = self.fixed_parameters)
        model_spin = evolution()

        prior = 1.0

        for (key_obs,value_obs),(key_parameter,value_parameter) in zip(self.observation_data.items(),parameter_set.items()):
            prior  *= scipy.stats.norm(value_obs['value'],value_obs['sigma']).pdf(value_parameter)

        likelihood = scipy.stats.norm(self.observation_data['Pspin']['value'],self.observation_data['Pspin']['sigma']).pdf(model_spin)

        posterior = prior*likelihood

        return posterior



    def acceptance_probability(self,current_values,proposed_values):

        p = self.posterior_probability(proposed_values)/self.posterior_probability(current_values)

        return p


    def values_proposed(self,current_values):

        proposed = dict()

        for (name_obs,value_obs),(name_step,value_step) in zip(current_values.items(),self.proposed_step.items()):

            proposed[name_obs] =  scipy.stats.norm.rvs(loc=value_obs, scale=value_step)

        return proposed

    def metropolis_hasting(self):

        """runs the metropolis hastrings algorithm for number of iterations given"""

        f.open("accpted_parameters.txt")
        f.open("rejected_parameters.txt")

        # loop for number of iterations times

        while self.iteration_step < self.total_iterations:

            if self.iteration_step==0: current_values = self.start_values

            #draw a random value from proposal function
            values_proposed = self.values_proposed(current_values)

            # calculate acceptance probablity
            p_acceptance = acceptance_probability(current_values, values_proposed)

            if p_acceptance > 1:

                current_values = values_proposed
                f.open("accepted_parameters.txt")
                f.write('step = %d\t', iteration_step, 'value = %d', current_values, '\n')
                f.close()

            else:

                f.open("rejected_parameters.txt")
                f.write('step = %d\t', iteration_step, 'value = %d', proposed_values, '\n')
                f.close()

                continue

                self.iteration_step = self.iteration_step + 1

        return none


    def __init__(self,interpolator,fixed_parameters,observation_data,logQ,proposed_step,total_iterations):


        self.interpolator  = interpolator
        self.fixed_parameters = fixed_parameters
        self.observation_data = observation_data
        self.convective_phase_lag = phase_lag(logQ)
        self.proposed_step = proposed_step
        self.iteration_step = 0
        self.total_iterations= total_iterations
        self.updated_parameters = dict()
        self.start_values = dict()



if __name__ == '__main__':



    observation_data = dict(Teff = dict(vlaue = 1.0, sigma = 0.1),
                        feh = dict(value = 1.0, sigma = 0.1),
                        rvk = dict(value = 1.0, sigma = 0.1),
                        inclination = dict(value = 1.0, sigma = 0.1)
                    )

    fixed_parameters = dict(disk_dissipation_age =  dict(value = 1.0, sigma = 0.1))

    proposed_step = dict( Teff_step = 0.1,
                      feh_step = 0.1,
                      rvk_step = 0.1,
                      inclination = 0.1)



