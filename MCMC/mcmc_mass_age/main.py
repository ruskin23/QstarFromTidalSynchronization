



class metropolis_hasting:

    def posterior_probability(self,
                              parameter_set):




    def check_acceptance(self):
        self.p_acceptance = self.propose_posterior/self.current_posterior
        if self.p_acceptance>1:self.isAccepted=True
        else:
            rand=numpy.random.random_sample()
            if self.p_acceptance>rand:self.isAccepted=True
            else:self.isAccepted=False

    def propose_values(self):

        for key in self.current_parameters:
            self.proposed_parameters[key]scipy.stats.norm.rvs(local=self.current_parameters[key],scale=self.step_size[key])

        self.propose_posterior=self.posterior_probability

    def initialise_parameters(self):

        self.current_parameters['primary_mass']=numpy.random.uniform(low=self.parameters['primary_mass']['min'],high=self.parameters['primary_mass']['max'])
        self.current_parameters['feh']=scipy.stats.norm.rvs(local=)
        self.current_parameters['age']=numpy.random.uniform(low=self.parameters['age']['min'],high=self.parameters['age']['max'])

        self.current_posterior=self.posterior_probability()

        self.write_output()


    def iterations(self):
        #MCMC iterations

        #Sample the first set of parameters and calculate its posterior
        if self.iteration=1: self.initialise_parameters()

        while True:

            #propose new parameters
            self.propose_values()

            #calculate acceptance probability
            self.acceptance_probability()

            #update output files
            if self.acceptance_rate>1:self.isAccepted=True
            else self.acceptance_rate:self.isAccepted=False
            self.write_output()

            #increment iteration step:
            self.iteration_step=self.iteration_step+1




    def __init__(self,
                 interpolator,
                 parameters,
                 constraints,
                 step_size):

        self.interpolator=interpolator
        self.parameters=parameters
        self.constraints=constaints
        self.step_size=step_size

        self.updated_parameters=None






if __name__ == '__main__':

    #interpolator

    #parameters
    parameters = dict(
        primary_mass=dict(min=0.4,max=1.4),
        feh=dict(value=,sigma=),
        age=dict(min=,max=)
    )

    constaints = dict(
        teff=dict(value=,sigma=)
        logg=dict(value=,sigma=)
    )

    step_size = dict(
        feh_step_size=
    )

    #initialise class

    #metropolis_hasting


