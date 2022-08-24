#!/usr/bin/env python3

from cmath import isnan
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
import numpy

_logger = logging.getLogger()

class InitialConditionSolver:
    """Find initial orbital period and eccentricity which reproduce
        current orbital period and eccentricity of a given system """



    def check_if_no_solution(self):

        P_guess=10+self.initial_guess[0]

        e_ulimit=0.70

        _logger.info('\nChecking if a solutions exists at high e={!r}'.format(e_ulimit))
        while True:
            try:
                _=self.initial_condition_errfunc((P_guess,e_ulimit))
                if numpy.logical_or(numpy.isnan(self.final_orbital_period),numpy.isnan(self.final_eccentricity)):
                    _logger.info('\nSolver returned Nan for orbital period = {!r}. Increasing orbital period by 10.0'.format(P_guess))
                    P_guess+=10.0
                    continue
                else:
                    _logger.info('\nFound non-Nan eccenctricity upper limit at e={!r}'.format(e_ulimit))
                    break
            except:
                _logger.info('\nSolver crashed for eccentricity = {!r}. Decreasing eccentricity upper limit by 0.02'.format(e_ulimit))
                P_guess+=10.0
                # e_ulimit-=0.02
                continue
            
        
        _logger.info('\nSolving for orbital period at eccentricity upper limit')
        
        p_root=self.orbital_period_solver(e_ulimit,P_guess=P_guess)

        _logger.info('\nSolution Check completed with Final_Eccentricity={!r} at Eccentricity_Limit={!r}'.format(self.final_eccentricity,e_ulimit))
        
        if numpy.isnan(p_root):return scipy.nan,scipy.nan
        else: return e_ulimit,self.delta_e


    def initial_condition_errfunc(self,initial_conditions):
        """Error function which returns the difference between final values and intial values"""

        initial_orbital_period=initial_conditions[0]
        initial_eccentricity=initial_conditions[1]

        _logger.info('\nTrying Porb_initial = {!r} , e_initial = {!r}'.format(initial_orbital_period,initial_eccentricity))

        if self.initial_guess is not None and self.function=='root':
            if initial_orbital_period==self.initial_guess[0] and initial_eccentricity==self.initial_guess[1]:
                _logger.info('final_orbital_period = {!r} , final_eccentricity = {!r}'.format(self.final_orbital_period,self.final_eccentricity))
                _logger.info('delta_p = {!r} , delta_e = {!r}'.format(self.delta_p,self.delta_e))
                _logger.info('Spin Period = %s',repr(self.spin))
                return numpy.array(self.err_intial_guess[0])

        if numpy.isnan(initial_orbital_period) or numpy.isnan(initial_eccentricity):
            _logger.warning('Solver using NaN as initial values.')
            raise ValueError('Solution cant be found as solver is trying NaN as initial values')

        if initial_eccentricity>0.80  or initial_eccentricity<0 or initial_orbital_period<0:
            _logger.warning('Encoutnered invalid initial values, returning NaN')
            return scipy.nan

        if initial_conditions not in self.solver_cache:
            binary_system=BinaryObjects(self.interpolator,self.parameters)

            binary=binary_system.create_binary_system(self.primary,
                                                    self.secondary,
                                                    initial_orbital_period=initial_orbital_period,
                                                    initial_eccentricity=initial_eccentricity,
                                                    secondary_angmom=self.secondary_angmom)

            binary.evolve(
                self.age,
                self.evolution_max_time_step,
                self.evolution_precision,
                None,
                timeout=3600
            )

            final_state=binary.final_state()
            if final_state.age != self.target_age:
                _logger.warning('Evolution did not reach target age, crashed at age = {!r} Gyr'.format(final_state.age))
                assert(final_state.age==self.target_age)

            
            self.final_orbital_period=binary.orbital_period(final_state.semimajor)
            self.final_eccentricity=final_state.eccentricity

            if numpy.logical_or(numpy.isnan(self.final_orbital_period),numpy.isnan(self.final_eccentricity)):
                _logger.warning('Enountered NaN in final_orbital_period={!r} or final_eccentricity={!r}'.format(self.final_orbital_period,self.final_eccentricity))
                _logger.warning('Binary was destroyed')
                return scipy.nan

            self.delta_p=self.final_orbital_period-self.target_orbital_period
            self.delta_e=self.final_eccentricity-self.target_eccentricity

            self.spin=(2*numpy.pi*binary.primary.envelope_inertia(final_state.age)/final_state.primary_envelope_angmom)

            binary.delete()

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
                _logger.info('\nFound non-NaN intial guess in {!r} tries. dp={!r},de={!r}'.format(n,self.delta_p,self.delta_e))
                break


    def brent_orbital_period_func(self,Pguess,eccentricity):
        return self.initial_condition_errfunc((Pguess,eccentricity))

    def evaluate_dp(self,P_guess,eccentricity):

        try:
            dp=self.brent_orbital_period_func(P_guess,eccentricity)
            if numpy.isnan(dp):
                _logger.info('\nEncountered NaN for initial conditions ({!r},{!r})'.format(P_guess,eccentricity))
                dp=-P_guess
        except:
            _logger.info('\nSolver Crashed for initial conditions ({!r},{!r}).'.format(P_guess,eccentricity))
            dp=scipy.nan
        
        return dp



    def orbital_period_solver(self,
                              eccentricity,
                              P_guess=None):

        _logger.info('\nInitiating Orbial Period Solver')
        if P_guess is None:
            P_guess=self.initial_guess[0]

        _logger.info('\nFinding first limit for orbital Period')

        dp_initial=scipy.nan
        while P_guess<60:
            dp_initial=self.evaluate_dp(P_guess,eccentricity)
            if numpy.isnan(dp_initial):P_guess+=5.0
            else:break
            
        
        if numpy.isnan(dp_initial):
            _logger.info('\nCouldnt find first limit for orbitial period. Returning NaN')
            return scipy.nan

        
        if abs(dp_initial)<1e-4:
            _logger.info('\nSolution already found. Skipping orbital period solver')
            return P_guess

        dp_new=dp_initial
        P_a=P_guess
        
        last_cached_ic=list(self.solver_cache.keys())[-1]
        _logger.info('\nFinding second limit for orbital period with initial_delta_p = {!r}'.format(dp_initial))
        while dp_new*dp_initial>0 and P_guess<60 and P_guess>0:
            
            if self.solver_cache[last_cached_ic]['final_orbital_period']<last_cached_ic[0]:
                if dp_initial<0: P_guess+=+5.0
                else: P_guess-=2.0
                    # if P_guess>self.initial_guess[0]: P_guess-=2.0
                    # else:
                    #     _logger.info('P_guess went below good initial condition. Checking if solution is possible')
                    #     P_guess=P_a+0.5
                    #     dp_check=self.evaluate_dp(P_guess,eccentricity)
                    #     if numpy.isnan(dp_check):
                    #         _logger.info('No Solution for Second Limit')
                    #     else:
                    #         if dp_check*dp_initial<0:
                    #             _logger.info('Found Second Limit at P={!r}'.format(P_guess))
                    #             dp_new=dp_check
                    #         else:
                    #             dp_new=scipy.nan

                    #     break
            
            else:
                P_guess=P_guess-0.5

            
            dp_new=self.evaluate_dp(P_guess,eccentricity)
        
        if numpy.isnan(dp_new):
            _logger.info('\nCouldnt find second limit for orbitial period. Returning NaN')
            return scipy.nan
        
        if dp_new*dp_initial>0:
            _logger.info('\nBounds are not opposite.No solution found for orbital_period')
            return scipy.nan

        P_b=P_guess

        _logger.info('\nUsing brentq method to solve between orbital periods a={!r} and b={!r}'.format(P_a,P_b))
 
        return scipy.optimize.brentq(self.brent_orbital_period_func,P_a,P_b,args=(eccentricity,),xtol=1e-4,rtol=1e-5)


    def brent_eccentricity_func(self,eccentricity,Pguess):

        dp=self.orbital_period_solver(eccentricity,P_guess=Pguess)
        if numpy.isnan(dp):
            _logger('\nGot dp=Nan. returning delta_e=Nan')
            return scipy.nan
        else:
            _logger.info('\nFor initial_condictions=({!r},{!r}) Found delta_e={!r}'.format(eccentricity,Pguess,self.delta_e))
            return self.delta_e


    def nested_solver(self):

        self.function='nested'

        if self.final_eccentricity<1e-8:
            _logger.info('\nFinal eccentricity = {!r} is less than 1e-8.'.format(self.final_eccentricity))
            e_ulimit,de_ulimit=self.check_if_no_solution()
            
            if numpy.isnan(de_ulimit):
                _logger.info('Solution check crashed')
                return scipy.nan

            if self.final_eccentricity<self.target_eccentricity: 
                _logger.info('No Solution exist as Final Eccentricity={!r} < Target Eccentricity={!r} at e=0.7'.format(self.final_eccentricity,self.target_eccentricity))
                return scipy.nan
            

            P_guess=self.initial_guess[0]
            e_llimit=self.initial_guess[1]

            #should be irrelevant
            while True:
                _=self.initial_condition_errfunc((P_guess,e_llimit))
                if numpy.logical_or(numpy.isnan(self.solver_cache[(P_guess,e_llimit)]['final_orbital_period']),
                                    numpy.isnan(self.solver_cache[(P_guess,e_llimit)]['final_eccentricity'])):
                    _logger.info('\nIncreasing eccentricity lower limit by 0.02')
                    e_llimit+=0.02
                    continue
                else:
                    _logger.info('\nFound non-Nan eccenctricity lower limit at e={!r}'.format(e_llimit))
                    break
            #--------------------
            
            p_root=self.orbital_period_solver(e_llimit,P_guess=P_guess)
            if numpy.isnan(p_root):de_llimit=scipy.nan
            else: de_llimit=self.delta_e
            
            if de_llimit*de_ulimit>0 or numpy.isnan(de_llimit):
                _logger.info('\nNo solution cannot exist between eccentricity limit. Error lower limit = {!r} Error upper limit = {!r}'.format(de_llimit,de_ulimit))
                return scipy.nan 
            else:
                e_root=scipy.optimize.brentq(self.brent_eccentricity_func,e_llimit,e_ulimit,args=(P_guess,),xtol=1e-5,rtol=1e-6)
                _logger.info('\nEccentricity root found. E_ROOT={!r}'.format(e_root))
                return self.spin

        else:

            P_guess=self.initial_guess[0]
            e_llimit=self.initial_guess[1]

            _logger.info('\nFinding lower limit for eccentricity')
            p_root=self.orbital_period_solver(e_llimit,P_guess=P_guess)
            if numpy.isnan(p_root):
                _logger.info('\nSolver cannot find orbital period at lower eccentricity limit. No Solution')
                return scipy.nan
            else:
                _logger.info('\nFound orbital period P={!r} at lower eccentricity limit e={!r}'.format(p_root,e_llimit))
            
            
            last_cached_ic=list(self.solver_cache.keys())[-1]
            if numpy.logical_and(abs(self.solver_cache[last_cached_ic]['delta_p'])<1e-5,
                                 abs(self.solver_cache[last_cached_ic]['delta_e'])<1e-5):
                _logger.info('Solution already found. Skipping eccenrticity solver')
                return self.solver_cache[last_cached_ic]['spin']
            
            _logger.info('\nFinding upper limit for eccentricity')
            de_initial=self.solver_cache[last_cached_ic]['delta_e']
            de_new=de_initial
            eccentricity=e_llimit
            P_guess=p_root
            
            while de_new*de_initial>0 and eccentricity<0.75:
                
                if numpy.isnan(de_new):eccentricity+=0.1
                else: eccentricity+=min(abs(1.5*de_new),0.1)

                _logger.info('Increasing upper eccentricity limit by 0.1. e_new={!r}'.format(eccentricity))

                while True:
                    try:
                        _=self.initial_condition_errfunc((P_guess,eccentricity))
                        if numpy.logical_or(numpy.isnan(self.solver_cache[(P_guess,eccentricity)]['delta_p']),
                                            numpy.isnan(self.solver_cache[(P_guess,eccentricity)]['delta_e'])):
                            _logger.info('\nSolver returned Nan for orbital period = {!r}. Increasing orbital period by 5.0'.format(P_guess))
                            P_guess+=5.0
                        else:break
                    except:
                        _logger.info('\nSolver crashed for orbital period = {!r}. Increasing orbital period by 5.0'.format(P_guess))
                        P_guess+=5.0
                        
                p_root=self.orbital_period_solver(eccentricity,P_guess=P_guess)
                if numpy.isnan(p_root):
                    _logger.info('\nSolver cannot find orbital period at upper eccentricity limit. No Solution')
                    return scipy.nan

                last_cached_ic=list(self.solver_cache.keys())[-1]
                de_new=self.solver_cache[(last_cached_ic[0],last_cached_ic[1])]['delta_e']
                _logger.info('\nFound orbital period P={!r} at upper eccentricity limit e={!r} with de={!r}'.format(p_root,eccentricity,de_new))

            e_ulimit=eccentricity
            if de_new*de_initial>0:
                _logger.info('\nNo solution found for eccentricity')
                return scipy.nan

            _logger.info('\nFinding eccentricity root between e_a={!r} and e_b={!r}'.format(e_llimit,e_ulimit))
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

        self.calculate_good_initial_guess()
        
        self.nested_solver()

        cached_ic_tuples=list(self.solver_cache.keys())
        initial_conditions_results=cached_ic_tuples[-1]
        _logger.info('Solver_Results:')
        _logger.info('Intial_Orbital_Period={!r} , Initial_Eccentricity={!r}'.format(initial_conditions_results[0],initial_conditions_results[1]))
        _logger.info('Final_Orbital_Period={!r} , Final_Eccentricity={!r}'.format(self.final_orbital_period,self.final_eccentricity))
        _logger.info('Errors: delta_p={!r} , delta_e={!r}'.format(self.delta_p,self.delta_e))
        _logger.info('Final_Spin_Period={!r}'.format(self.spin))

        if numpy.isnan(self.spin):
            _logger.warning('Spin=Nan after solver. Setting Spin=inf')
            self.spin=numpy.inf

        return self.spin

