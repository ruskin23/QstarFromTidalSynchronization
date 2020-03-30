import matplotlib.pyplot as plt
from get_cdf import GetCDF
import argparse
import scipy

def cmdline_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-l',
                        action='store',
                        dest='system',
                        help='select a system')

    parser.add_argument('-t',
                        action='store',
                        dest='test_case',
                        help='name of test model'
                        )


    return parser.parse_args()

if __name__=='__main__':


    args = cmdline_args()

    system_number=args.system
    test_case=args.test_case

    catalog_file='SpinlogQCatalog_el0.4.txt'
    solution_file='SolutionFileBreaks0.0.txt'
    samples_file='MassAgeFehSamples_'+system_number+'.txt'


    with open(catalog_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if system_number==at_system:
                Porb_value=float(x[6])
                Porb_error=float(x[7])
                eccentricity_value=float(x[8])
                eccentricity_error=float(x[9])
                primary_mass_value=float(x[15])
                age_value=float(x[16])
                feh_value=float(x[17])
                break


    with open(solution_file,'r') as f:
        next(f)
        for lines in f:
            x=lines.split()
            at_system=x[0]
            if system_number==at_system:
                logQ_value=float(x[1])
                break


    sampling_parameters=dict(Porb=dict(value=Porb_value,
                                       sigma=Porb_error,
                                       dist='Normal',
                                       step=Porb_error),

                             eccentricity=dict(value=eccentricity_value,
                                               sigma=eccentricity_error,
                                               dist='Normal',
                                               step=eccentricity_error),

                             Wdisk=dict(value=4.1,
                                        min=2*scipy.pi/14,
                                        max=2*scipy.pi/1.4,
                                        dist='Uniform',
                                        step=0.5),

                             logQ=dict(value=logQ_value,
                                       min=5.0,
                                       max=12.0,
                                       dist='Uniform',
                                       step=0.5),

                             primary_mass=dict(value=primary_mass_value,
                                               dist='Samples',
                                               step=1.0),

                             age=dict(value=age_value,
                                      dist='Samples',
                                      step=1.0),

                             feh=dict(value=feh_value,
                                      dist='Samples',
                                      step=0.3)

                             )


    model_width=dict(primary_mass=0.8,
                     age=0.8,
                     feh=0.2,
                     logQ=0.2,
                     eccentricity=0.01)


    rho=0.8

    Distrbution=GetCDF(system_number,
                       sampling_parameters,
                       model_width,
                       rho)

    fig,(ax1,ax2)=plt.subplots(1,2)

    if test_case in ['gp','gpt','gptc']:
        age_mcmc_cdf=Distrbution.get_mcmc_distribution(test_case,'age')
        age_model_cdf=Distrbution.get_model_distribution(test_case,parameter='age',phi='age')

        ax1.scatter(*zip(*age_mcmc_cdf))
        ax1.plot(*zip(*age_model_cdf))
        ax1.set_title('age')

        logQ_mcmc_cdf=Distrbution.get_mcmc_distribution(test_case,'logQ')
        logQ_model_cdf=Distrbution.get_model_distribution(test_case,parameter='logQ',phi='age')

        ax2.scatter(*zip(*logQ_mcmc_cdf))
        ax2.plot(*zip(*logQ_model_cdf))
        ax2.set_title('logQ')

        plt.show()

    else:
        eccentricity_mcmc_cdf=Distrbution.get_mcmc_distribution(test_case,'eccentricity')
        eccentricity_model_cdf=Distrbution.get_model_distribution(test_case,parameter='eccentricity',theta='eccentricity')

        ax1.scatter(*zip(*eccentricity_mcmc_cdf))
        ax1.plot(*zip(*eccentricity_model_cdf))
        ax1.set_title('eccentricity')

        logQ_mcmc_cdf=Distrbution.get_mcmc_distribution(test_case,'logQ')
        logQ_model_cdf=Distrbution.get_model_distribution(test_case,parameter='logQ',theta='eccentricity')

        ax2.scatter(*zip(*logQ_mcmc_cdf))
        ax2.plot(*zip(*logQ_model_cdf))
        ax2.set_title('logQ')

        plt.show()

