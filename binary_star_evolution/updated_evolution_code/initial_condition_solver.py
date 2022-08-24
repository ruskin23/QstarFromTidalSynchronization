#!/usr/bin/env python3

import time
import logging
import sys
from pathlib import Path
from directories import directories

home_dir=str(Path.home())
path=directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')


from orbital_evolution.transformations import phase_lag
from create_objects import BinaryObjects

import scipy
from scipy import optimize
import numpy

from lmfit import minimize, Parameters, Parameter, report_fit
from functools import lru_cache

wsun = 0.24795522138
_logger = logging.getLogger()


def check_last_nan(a):
    for i in range(len(a)):
        val=a[len(a)-i-1]
        check=numpy.isnan(val)
        if check == False:
            return(val,len(a)-i-1)




class InitialConditionSolver:
    """Find initial orbital period and eccentricity which reproduce
        current orbital period and eccentricity of a given system """

    def initial_condition_errfunc(self,initial_conditions):
        """Error function which returns the difference between final values and intial values"""

        initial_orbital_period=initial_conditions[0]
        initial_eccentricity=initial_conditions[1]
        initial_conditions=tuple([initial_orbital_period,initial_eccentricity])

        _logger.info('\nTrying Porb_initial = {!r} , e_initial = {!r}'.format(initial_orbital_period,initial_eccentricity))

        # if self.initial_guess is not None and self.function=='root':
        #     if initial_orbital_period==self.initial_guess[0] and initial_eccentricity==self.initial_guess[1]:
        #         _logger.info('final_orbital_period = {!r} , final_eccentricity = {!r}'.format(self.final_orbital_period,self.final_eccentricity))
        #         _logger.info('delta_p = {!r} , delta_e = {!r}'.format(self.delta_p,self.delta_e))
        #         _logger.info('Spin Period = %s',repr(self.spin))
        #         return numpy.array(self.err_intial_guess)

        if numpy.isnan(initial_orbital_period) or numpy.isnan(initial_eccentricity):
            _logger.warning('Solver using NaN as initial values.')
            raise ValueError('Solution cant be found as solver is trying NaN as initial values')

        if initial_eccentricity>0.80  or initial_eccentricity<0 or initial_orbital_period<0:
            _logger.warning('Encoutnered invalid initial values, returning NaN')
            if self.function=='minimize' or self.function=='nested': return scipy.nan
            else: return numpy.array([scipy.nan,scipy.nan])


        if initial_conditions not in self.solver_cache:
            binary_system=BinaryObjects(self.interpolator,self.parameters)

            binary=binary_system.create_binary_system(self.primary,
                                                    self.secondary,
                                                    initial_orbital_period=initial_orbital_period,
                                                    initial_eccentricity=initial_eccentricity,
                                                    secondary_angmom=self.secondary_angmom)


            # time_start=time.time()

            binary.evolve(
                self.age,
                self.evolution_max_time_step,
                self.evolution_precision,
                None,
                timeout=3600
            )

            # time_diff=time.time()-time_start
            
            # _logger.info('\n\nTIME TAKEN FOR EVOLUTION = {!r} mins {!r} s\n\n'.format(time_diff//60,time_diff%60))

            final_state=binary.final_state()
            if final_state.age != self.target_age:
                _logger.warning('Evolution did not reach target age, crashed at age = {!r} Gyr'.format(final_state.age))
                assert(final_state.age==self.target_age)

            
            self.final_orbital_period=binary.orbital_period(final_state.semimajor)
            self.final_eccentricity=final_state.eccentricity

            if numpy.logical_or(numpy.isnan(self.final_orbital_period),numpy.isnan(self.final_eccentricity)):
                _logger.warning('Enountered NaN in final_orbital_period={!r} or final_eccentricity={!r}'.format(self.final_orbital_period,self.final_eccentricity))
                _logger.warning('Binary was destroyed')
                if self.function=='minimize' or self.function=='nested': return scipy.nan
                else: return numpy.array([scipy.nan,scipy.nan])
            #     evolution = binary.get_evolution()
            #     self.final_eccentricity,non_nan_index=check_last_nan(evolution.eccentricity)
            #     _logger.warning('Binary system was destroyed at age = {!r} Gyr'.format(evolution.age[non_nan_index]))
            #     self.delta_p=-self.target_orbital_period-self.target_age+evolution.age[non_nan_index]

            self.delta_p=self.final_orbital_period-self.target_orbital_period
            self.delta_e=self.final_eccentricity-self.target_eccentricity

            self.spin=(2*numpy.pi*binary.primary.envelope_inertia(final_state.age)/final_state.primary_envelope_angmom)

            binary.delete()

            # if rootsum_errsq<1e-5:
            #     _logger.info('Sqrt sum of squares of delta close enough to zero. Setting dp=0,de=0')
            #     self.delta_p,self.delta_e=0.0,0.0

            self.solver_cache[initial_conditions]={'final_orbital_period':self.final_orbital_period,
                                                    'final_eccentricity':self.final_eccentricity,
                                                    'delta_p':self.delta_p,
                                                    'delta_e':self.delta_e,
                                                    'spin':self.spin}
        else:
            current_cache=self.solver_cache[initial_conditions]
            self.final_orbital_period=current_cache['final_orbital_period']
            self.final_eccentricity=current_cache['final_eccentricity']
            self.delta_p=current_cache['delta_p']
            self.delta_e=current_cache['delta_e']
            self.spin=current_cache['spin']


        _logger.info('final_orbital_period = {!r} , final_eccentricity = {!r}'.format(self.final_orbital_period,self.final_eccentricity))
        _logger.info('delta_p = {!r} , delta_e = {!r}'.format(self.delta_p,self.delta_e))
        _logger.info('Spin Period = %s',repr(self.spin))

        rootsum_errsq=numpy.sqrt(self.delta_p**2 + self.delta_e**2)

        if self.function=='minimize':
            _logger.info('Error Function = {!r}'.format(rootsum_errsq))
            return rootsum_errsq
        if self.function=='root': 
            return numpy.array([self.delta_p,self.delta_e])
        if self.function=='nested':
            return self.delta_p





    def __init__(self,
                interpolator,
                parameters,
                evolution_max_time_step=1e-3,
                evolution_precision=1e-6,
                secondary_angmom=None,
                initial_guess=None):
        """
        Args:

            - interpolator

            - parameters
                a dictionary of all the parameters

            - evolution_max_time_step:
                The maximum timestep the evolution is allowed to make.

            - evolution_precision:
                The precision to require of the evolution.

            - secondary_angmom:
                The initial angular momentum of secondary star
        """

        self.interpolator=interpolator
        self.parameters=parameters
        for item,value in parameters.items():
            setattr(self,item,value)
        if 'logQ' in parameters: self.convective_phase_lag=phase_lag(self.logQ)
        else: self.convective_phase_lag=self.phase_lag_max

        self.evolution_max_time_step=evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.secondary_angmom = secondary_angmom
        self.initial_guess=initial_guess

        self.final_orbital_period,self.final_eccentricity=scipy.nan,scipy.nan
        self.delta_p,self.delta_e=scipy.nan,scipy.nan
        self.spin=scipy.nan

        self.solver_cache=dict()
        


    def calculate_initial_simplex(self):

        Pguess=self.target_orbital_period
        e=self.target_eccentricity
        _logger.info('\nCalculating first simplex using Pguess = {} e =  {}'.format(repr(Pguess),repr(e)))
        err_fun=self.initial_condition_errfunc((Pguess,e))
        if numpy.isnan(err_fun):
            _logger.warning('Encontered nan values.')
            n=1
            while numpy.isnan(err_fun):
                n=n+1
                Pguess=Pguess*n
                err_fun=self.initial_condition_errfunc((Pguess,e))
        simplex_1=[Pguess,e]
        _logger.info('Found first simplex values at Pguess = {} e = {}'.format(repr(Pguess),repr(e)))

        Pguess=2*Pguess
        _logger.info('\nCalculating second simplex using Pguess = {} e =  {}'.format(repr(Pguess),repr(e)))
        err_fun=self.initial_condition_errfunc((Pguess,e))
        if numpy.isnan(err_fun):
            _logger.warning('Encontered nan values.')
            n=1
            while numpy.isnan(err_fun):
                n=n+1
                Pguess=Pguess*n
                err_fun=self.initial_condition_errfunc((Pguess,e))
        simplex_2=[Pguess,e]
        _logger.info('Found second simplex values at Pguess = {} e = {}'.format(repr(Pguess),repr(e)))


        e = e + 0.1 if e < 0.3 else -0.1
        _logger.info('\nCalculating third simplex using Pguess = {} e =  {}'.format(repr(Pguess),repr(e)))
        err_fun=self.initial_condition_errfunc((Pguess,e))
        if numpy.isnan(err_fun):
            _logger.warning('Encontered nan values.')
            while numpy.isnan(err_fun):
                e = e + (0.1 if e < 0.3 else -0.1)
                err_fun=self.initial_condition_errfunc((Pguess,e))
        simplex_3=[Pguess,e]
        _logger.info('Found second simplex values at Pguess = {} e = {}'.format(repr(Pguess),repr(e)))


        initial_simplex=numpy.array([simplex_1,simplex_2,simplex_3])
        _logger.info('\nInitial Simplex = {}'.format(repr(initial_simplex)))

        return initial_simplex

    def calculate_good_initial_guess(self):

        Pguess=self.target_orbital_period
        e=self.target_eccentricity
        n=0
        while True:
            n=n+1
            _=self.initial_condition_errfunc((Pguess,e))
            if numpy.isnan(self.delta_p) or numpy.isnan(self.delta_e):
                Pguess+=0.5
                continue
            else:
                self.initial_guess=[Pguess,e]
                self.err_intial_guess=[self.delta_p,self.delta_e]
                _logger.info('Found non-NaN intial guess in {!r} tries. dp={!r},de={!r}'.format(n,self.delta_p,self.delta_e))
                break


    def brent_orbital_period_func(self,Pguess,eccentricity):
        return self.initial_condition_errfunc((Pguess,eccentricity))

    def orbital_period_solver(self,
                              eccentricity,
                              P_guess=None):

        _logger.info('\nInitiating Orbial Period Solver')
        if P_guess is None:
            P_guess=self.initial_guess[0]

        _logger.info('Finding first limit for orbital Period')
        while P_guess<60:
            try:    
                dp_initial=self.brent_orbital_period_func(P_guess,eccentricity)
                if numpy.isnan(dp_initial):
                    _logger.info('Encountered NaN for initial conditions ({!r},{!r}). Increasing first orbital period limit by 5.0'.format(P_guess,eccentricity))
                    P_guess+=5.0
                else:break
            except:
                _logger.info('Solver Crashed for initial conditions ({!r},{!r}). Increasing first orbital period limit by 5.0'.format(P_guess,eccentricity))
                P_guess+=5.0

        if numpy.isnan(dp_initial):
            _logger.info('Couldnt find first limit for orbitial period. Returning Nan')
            return scipy.nan

        dp_new=dp_initial
        
        _logger.info('Finding second limit for orbital period with initial_delta_p = {!r}'.format(dp_initial))
        while dp_new*dp_initial>0 and P_guess<60 and P_guess>0:
            
            P_b=P_guess
            if dp_initial<0: P_guess+=+5.0
            else: P_guess-=2.0
            try:
                dp_new=self.brent_orbital_period_func(P_guess,eccentricity)
                if numpy.isnan(dp_new):
                    _logger.info('Encountered NaN for initial conditions ({!r},{!r}). Changing second orbitial period limit'.format(P_guess,eccentricity))
                    continue
            except:
                _logger.info('Solver Crashed for initial conditions ({!r},{!r}). Changing second orbitial period limit'.format(P_guess,eccentricity))
                continue
        
        if numpy.isnan(dp_new):
            _logger.info('Couldnt find second limit for orbitial period. Returning Nan')
            return scipy.nan
        
        if dp_new*dp_initial>0:
            _logger.info('Bounds are not opposite.No solution found for orbital_period')
            return scipy.nan

        P_a=P_guess
        
        _logger.info('Using brentq method to solve between orbital periods a={!r} and b={!r}'.format(P_a,P_b))

        return scipy.optimize.brentq(self.brent_orbital_period_func,P_a,P_b,args=(eccentricity,),xtol=1e-4,rtol=1e-5)


    def lmfit_errfunc(self,params,eccentricity):
        return self.initial_condition_errfunc((params['P_i'].value,eccentricity))


    def check_if_no_solution(self):

        P_guess=10+self.initial_guess[0]

        e_ulimit=0.70


        while True:
            try:
                _=self.initial_condition_errfunc((P_guess,e_ulimit))
                if numpy.logical_or(numpy.isnan(self.final_orbital_period),numpy.isnan(self.final_eccentricity)):
                    _logger.info('Solver returned Nan for orbital period = {!r}. Increasing orbital period by 1.0'.format(P_guess))
                    P_guess+=10.0
                    continue
                else:
                    _logger.info('Found non-Nan eccenctricity upper limit at e={!r}'.format(e_ulimit))
                    break
            except:
                _logger.info('Solver crashed for eccentricity = {!r}. Decreasing eccentricity upper limit by 0.02'.format(e_ulimit))
                P_guess+=10.0
                # e_ulimit-=0.02
                continue
            
        
        _logger.info('\nSolving for orbital period at eccentricity upper limit')


        params=Parameters()
        params.add('P_i',value=P_guess,min=0,max=60)


        # _=minimize(self.lmfit_errfunc,params,args=(e_ulimit,),method='leastsq',nan_policy='omit')
        
        _=self.orbital_period_solver(e_ulimit,P_guess=P_guess)

        _logger.info('\nSolution Check completed with Final_Eccentricity={!r} at Eccentricity_Limit={!r}'.format(self.final_eccentricity,e_ulimit))
        
        return self.final_eccentricity>self.target_eccentricity,e_ulimit,self.delta_e

    def brent_eccentricity_func(self,eccentricity,Pguess):

        _=self.orbital_period_solver(eccentricity,P_guess=Pguess)
        _logger.info('for ic=({!r},{!r}) returning delta_e={!r}'.format(eccentricity,Pguess,self.delta_e))
        return self.delta_e


    def nested_solver(self):

        self.function='nested'

        if self.final_eccentricity<1e-8:
            _logger.info('Final eccentricity = {!r} is less than 1e-8.'.format(self.final_eccentricity))
            solution_exist,e_ulimit,de_ulimit=self.check_if_no_solution()

            if solution_exist is False:
                return scipy.nan
            

            P_guess=self.initial_guess[0]
            e_llimit=0.0

            while True:
                _=self.initial_condition_errfunc((P_guess,e_llimit))
                if numpy.logical_or(numpy.isnan(self.final_orbital_period),numpy.isnan(self.final_eccentricity)):
                    _logger.info('Increasing eccentricity lower limit by 0.02')
                    e_llimit+=0.02
                    continue
                else:
                    _logger.info('Found non-Nan eccenctricity lower limit at e={!r}'.format(e_llimit))
                    break

            _=self.orbital_period_solver(e_llimit,P_guess=P_guess)
            de_llimit=self.delta_e
            
            if de_llimit*de_ulimit>0:
                _logger.info('No solution can exist between eccentricity limit. Error lower limit = {!r} Error upper limit = {!r}'.format(de_llimit,de_ulimit))
                return scipy.nan 
            else:
                e_root=scipy.optimize.brentq(self.brent_eccentricity_func,e_llimit,e_ulimit,args=(P_guess,),xtol=1e-5,rtol=1e-6)
                _logger.info('Eccentricity root found. E_ROOT={!r}'.format(e_root))
                return self.spin

        else:

            P_guess=self.initial_guess[0]
            e_llimit=self.initial_guess[1]

            _logger.info('\nFinding lower limit for eccentricity')
            p_root=self.orbital_period_solver(e_llimit,P_guess=P_guess)
            _logger.info('Found orbital period P={!r} at lower eccentricity limit e={!r}'.format(p_root,e_llimit))
            
            _logger.info('\nFinding upper limit for eccentricity')
            de_initial=self.delta_e
            de_new=de_initial
            eccentricity=e_llimit
            while de_new*de_initial>0 and eccentricity<0.75:
                eccentricity+=0.2

                while True:
                    try:
                        _=self.initial_condition_errfunc((P_guess,eccentricity))
                        break
                    except:
                        P_guess+=5.0
                        
                p_root=self.orbital_period_solver(eccentricity,P_guess=P_guess)
                _logger.info('Found orbital period P={!r} at upper eccentricity limit e={!r}'.format(p_root,eccentricity))
                de_new=self.delta_e

            e_ulimit=eccentricity
            if de_new*de_initial>0:
                _logger.info('No solution found for eccentricity')
                return scipy.nan

            e_root=scipy.optimize.brentq(self.brent_eccentricity_func,e_llimit,e_ulimit,args=(P_guess,),xtol=1e-4,rtol=1e-5)

            

    def __call__(self, primary, secondary):
        """
        Find initial conditions which reproduce the given system now.

        Args:
            - primary

            - secondary
        """

        self.primary = primary
        self.secondary = secondary

        self.target_age=self.age
        self.target_orbital_period=self.orbital_period
        self.target_eccentricity=self.eccentricity
        

        _logger.info('\nSolving for p and e using function = {} and method {}'.format(self.function,self.method))

        

        if self.function=='minimize':
            self.initial_guess=[self.target_orbital_period,self.target_eccentricity]
            initial_simplex=self.calculate_initial_simplex()
        else:
            self.calculate_good_initial_guess()
        
        # self.nested_solver()
        time_start=time.time()

        try:
            if self.function=='least_squares':
                bounds=(numpy.array([0.0,0.0]), numpy.array([numpy.inf,0.8]))
                sol = optimize.least_squares(self.initial_condition_errfunc,
                                self.initial_guess,
                                bounds=bounds,
                                method=self.method
                                )
            if self.function=='root':
                sol = optimize.root(self.initial_condition_errfunc,
                                self.initial_guess,
                                method=self.method,
                                options={'xtol': 0.0, 'ftol': 1e-05}#, 'eps':numpy.sqrt(self.evolution_precision)}
                                )
            if self.function=='minimize':
                bounds=((0.0,numpy.inf), (0.0,0.8))
                sol = optimize.minimize(self.initial_condition_errfunc,
                                self.initial_guess,
                                method=self.method,
                                bounds=bounds,
                                options={'initial_simplex' : initial_simplex}
                                )


            initial_orbital_period_sol,initial_eccentricity_sol=sol.x
        except:
            initial_orbital_period_sol,initial_eccentricity_sol=scipy.nan,scipy.nan
            self.spin=numpy.inf


        _logger.info('Solver_Results:')
        _logger.info('Intial_Orbital_Period={!r} , Initial_Eccentricity={!r}'.format(initial_orbital_period_sol,initial_eccentricity_sol))
        _logger.info('Final_Orbital_Period={!r} , Final_Eccentricity={!r}'.format(self.final_orbital_period,self.final_eccentricity))
        _logger.info('Errors: delta_p={!r} , delta_e={!r}'.format(self.delta_p,self.delta_e))
        _logger.info('Final_Spin_Period={!r}'.format(self.spin))

        time_diff=time.time()-time_start
        _logger.info('\n\nTOTAL TIME TAKEN = {!r} mins {!r} s\n\n'.format(time_diff//60,time_diff%60))

        # if numpy.isfinite(self.spin):
        #     sum_of_squares=numpy.sqrt(self.delta_p**2+self.delta_e**2)
        #     if sum_of_squares>1e-3:
        #         _logger.info('Error Margins not good enough. Finding root using Nested_Solver')
        #         self.check_if_no_solution()
        # if numpy.isnan(self.spin):
        #     _logger.warning('Spin=Nan after solver. Setting Spin=inf')
        #     self.spin=numpy.inf

        return self.spin

