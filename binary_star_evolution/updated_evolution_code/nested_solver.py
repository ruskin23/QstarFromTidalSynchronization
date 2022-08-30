#!/usr/bin/env python3

import time
import logging
import sys
import scipy
import numpy

from pathlib import Path
from directories import directories

from orbital_evolution.transformations import phase_lag
from create_objects import BinaryObjects

home_dir = str(Path.home())
path = directories(home_dir)
sys.path.append(path.poet_path+'/PythonPackage')
sys.path.append(path.poet_path+'/scripts')


_logger = logging.getLogger()

class InitialConditionSolver:
    """Find initial orbital period and eccentricity which reproduce
    current orbital period and eccentricity of a given system """

    def check_if_no_solution(self):

        P_guess = 10+self.initial_guess[0]

        e_ulimit = 0.7

        _logger.info('\nChecking if a solutions exists at high e={!r}'.format(e_ulimit))

        # FIND NON-NAN INITIAL VALUE
        while P_guess < 60:
            try:
                _ = self.initial_condition_errfunc((P_guess,e_ulimit))
                if numpy.logical_or(numpy.isnan(self.final_orbital_period),
                                    numpy.isnan(self.final_eccentricity)):
                    _logger.info('\nSolver returned Nan for orbital period = {!r}.'.format(P_guess))
                    _logger.info('Increasing orbital period by 10.0')
                    P_guess+=10.0
                    continue
                else:
                    _logger.info('\nFound non-Nan eccenctricity upper limit at e={!r}'.format(e_ulimit))
                    break
            except:
                _logger.warning('\nSolver crashed for eccentricity = {!r}.'.format(e_ulimit))
                _logger.warning('Decreasing eccentricity upper limit by 0.02')
                P_guess+=10.0
                continue

        _logger.info('\nSolving for orbital period at eccentricity upper limit')

        p_root = self.orbital_period_solver(e_ulimit, P_guess=P_guess)

        if numpy.isnan(p_root):
            return scipy.nan, scipy.nan
        else: 
            _logger.info('\nSolution Check completed with Final_Eccentricity={!r} at Initial_eccentriciy={!r}'.format(self.final_eccentricity, e_ulimit))
            return e_ulimit, self.delta_e


    def initial_condition_errfunc(self,initial_conditions):
        """Error function which returns the difference between final values and intial values"""

        initial_orbital_period = initial_conditions[0]
        initial_eccentricity = initial_conditions[1]

        _logger.info('\nTrying Porb_initial = {!r} , e_initial = {!r}'.format(initial_orbital_period, initial_eccentricity))

        if self.initial_guess is not None and self.function == 'root':
            if initial_orbital_period == self.initial_guess[0] and initial_eccentricity == self.initial_guess[1]:
                _logger.info('final_orbital_period = {!r} , final_eccentricity = {!r}'.format(self.final_orbital_period, self.final_eccentricity))
                _logger.info('delta_p = {!r} , delta_e = {!r}'.format(self.delta_p, self.delta_e))
                _logger.info('Spin Period = {!r}}'.format(self.spin))
                return numpy.array(self.err_intial_guess[0])

        if numpy.isnan(initial_orbital_period) or numpy.isnan(initial_eccentricity):
            _logger.warning('Solver using NaN as initial values.')
            raise ValueError('Solution cant be found as solver is trying NaN as initial values')

        if initial_eccentricity > 0.80  or initial_eccentricity < 0 or initial_orbital_period < 0:
            _logger.warning('Encoutnered invalid initial values, returning NaN')
            return scipy.nan

        if initial_conditions not in self.solver_cache:
            binary_system = BinaryObjects(self.interpolator, self.parameters)

            binary = binary_system.create_binary_system(self.primary,
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
                assert(final_state.age == self.target_age)

            
            self.final_orbital_period=binary.orbital_period(final_state.semimajor)
            self.final_eccentricity=final_state.eccentricity
            if numpy.logical_or(numpy.isnan(self.final_orbital_period),
                                numpy.isnan(self.final_eccentricity)):
                _logger.warning('Enountered NaN in final_orbital_period={!r} or final_eccentricity={!r}'.format(self.final_orbital_period,self.final_eccentricity))
                _logger.warning('Binary was destroyed. Setting final orbital period as zero!')
                self.final_orbital_period = 0.0

            self.delta_p = self.final_orbital_period-self.target_orbital_period
            self.delta_e = self.final_eccentricity-self.target_eccentricity

            self.spin=(2
                       * numpy.pi
                       * binary.primary.envelope_inertia(final_state.age)
                       / final_state.primary_envelope_angmom
                       )

            binary.delete()

            self.solver_cache[initial_conditions] = {'final_orbital_period' : self.final_orbital_period,
                                                     'final_eccentricity' : self.final_eccentricity,
                                                     'delta_p' : self.delta_p,
                                                     'delta_e' : self.delta_e,
                                                     'spin' : self.spin}
            

        else:
            current_cache = self.solver_cache[initial_conditions]
            self.final_orbital_period = current_cache['final_orbital_period']
            self.final_eccentricity = current_cache['final_eccentricity']
            self.delta_p = current_cache['delta_p']
            self.delta_e = current_cache['delta_e']
            self.spin = current_cache['spin']


        _logger.info('final_orbital_period = {!r} , final_eccentricity = {!r}'.format(self.final_orbital_period, self.final_eccentricity))
        _logger.info('delta_p = {!r} , delta_e = {!r}'.format(self.delta_p, self.delta_e))
        _logger.info('Spin Period = %s',repr(self.spin))

        return self.delta_p





    def __init__(self,
                 interpolator,
                 parameters,
                 evolution_max_time_step=1e-2,
                 evolution_precision=1e-5,
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
        self.initial_guess = initial_guess

        self.final_orbital_period,self.final_eccentricity = scipy.nan,scipy.nan
        self.delta_p,self.delta_e = scipy.nan,scipy.nan
        self.spin = scipy.nan

        self.solver_cache = dict()


    def _calculate_good_initial_guess(self):

        Pguess=self.target_orbital_period
        e=self.target_eccentricity
        n=0
        while Pguess<60:
            n+=1
            try:
                _=self.initial_condition_errfunc((Pguess,e))
                if numpy.isnan(self.delta_p) or numpy.isnan(self.delta_e):
                    Pguess+=1.0
                    continue
                else:
                    self.initial_guess=[Pguess,e]
                    self.err_intial_guess=[self.delta_p,self.delta_e]
                    _logger.info('\nFound non-NaN intial guess in {!r} tries. dp={!r},de={!r}'.format(n,self.delta_p,self.delta_e))
                    break
            except:
                Pguess+=1.0
        
        return Pguess<60

    def brent_orbital_period_func(self,Pguess,eccentricity):
            return self.initial_condition_errfunc((Pguess,eccentricity))

    def evaluate_dp(self,P_guess,eccentricity):

        try:
            dp=self.initial_condition_errfunc((P_guess,eccentricity))
        except:
            _logger.info('\nEvolution Crashed for initial conditions ({!r},{!r}).'.format(P_guess,eccentricity))
            dp=scipy.nan
        
        return dp

    def orbital_period_solver(self,
                              eccentricity,
                              P_guess=None):

        _logger.info('\nInitiating Orbial Period Solver')
        if P_guess is None:
            P_guess=self.initial_guess[0]

        _logger.info('\nFinding first limit for orbital Period')

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
        
        
        _logger.info('\nFinding second limit for orbital period with initial_delta_p = {!r}'.format(dp_initial))
        while dp_new*dp_initial>0 or numpy.isnan(dp_new):
            
            if numpy.logical_or(P_guess>60,P_guess<0):
                _logger.info('\nReached max iteration to find orbital period limit.')
                break
            
            cached_ic_list=list(self.solver_cache.keys())
            if numpy.isnan(dp_new):
                dp_new=self.solver_cache[cached_ic_list[-1]]['delta_p']
                _logger.info('\ndp=NaN. Setting dp={!r} from previous evolution'.format(dp_new))
            else:
                final_porb=self.solver_cache[cached_ic_list[-1]]['final_orbital_period']
                init_porb=cached_ic_list[-1][0]

            if dp_new<0: 
                P_guess+=+5.0
            else:
                if final_porb<init_porb:
                    P_guess-=2.0            
                else:
                    P_guess-=0.5

            # import pdb; pdb.set_trace()
            dp_new=self.evaluate_dp(P_guess,eccentricity)
        
        if numpy.isnan(dp_new):
            _logger.info('\nCouldnt find second limit for orbitial period. Returning NaN')
            return scipy.nan
        
        if dp_new*dp_initial>0:
            _logger.info('\nBounds are not opposite.No solution found for orbital_period')
            return scipy.nan

        P_b=P_guess

        _logger.info('\nUsing brentq method to solve between orbital periods a={!r} and b={!r}'.format(P_a,P_b))

        xtol,rtol=1e-4,1e-5
        try:
            p_root=scipy.optimize.brentq(self.brent_orbital_period_func,P_a,P_b,args=(eccentricity,),xtol=xtol,rtol=rtol)
            cached_ic_list=list(self.solver_cache.keys())
            last_cached_ic=cached_ic_list[-1]
            dp=self.solver_cache[last_cached_ic]['delta_p']
            if dp>1e-1:
                _logger.warning('Bad Solution Found. dp={!r} is not small enough'.format(dp))
                p_root=scipy.nan
        except Exception as e:
            _logger.warning('\nOrbital Period Solver Crashed with Exception={!r} while using brentq method while finding solution'.format(e))
            cached_ic_list=list(self.solver_cache.keys())
            last_cached_ic=cached_ic_list[-1]
            final_e_last=self.solver_cache[last_cached_ic]['final_eccentricity']
            _logger.info('\nLast successful evolution gave final_eccentricity={!r} for initial_eccentricity={!r}'.format(final_e_last,last_cached_ic[1]))
            p_root=scipy.nan
        
        return p_root

    def brent_eccentricity_func(self,eccentricity,Pguess):

        p_root=self.orbital_period_solver(eccentricity,P_guess=Pguess)
        if numpy.isnan(p_root):
            _logger.warning('\nGot p_root=NaN for eccentricity={!r}. Exiting Solver'.format(eccentricity))
            raise ValueError()
        else:
            _logger.info('\nFor initial_condictions=({!r},{!r}) Found delta_e={!r}'.format(eccentricity,Pguess,self.delta_e))
            return self.delta_e


    def nested_solver(self):

        self.function='nested'

        if self.final_eccentricity<1e-8:
            _logger.info('\nFinal eccentricity = {!r} is less than 1e-8.'.format(self.final_eccentricity))
            e_ulimit,de_ulimit=self.check_if_no_solution()
            
            if numpy.isnan(de_ulimit):
                _logger.info('Solution check crashed as orbital period root at e=0.7 is NaN')
                last_cached_ic=list(self.solver_cache.keys())[-1]
                last_final_e=self.solver_cache[last_cached_ic]['final_eccentricity']
                _logger.info('Last evolution gave Final eccentricity={!r} and Target Eccentricity={!r}'.format(last_final_e,self.target_eccentricity))
                return scipy.nan

            if self.final_eccentricity<self.target_eccentricity: 
                _logger.info('No Solution exist as Final Eccentricity={!r} < Target Eccentricity={!r} at e=0.7'.format(self.final_eccentricity,self.target_eccentricity))
                return scipy.nan
            

            P_guess=self.initial_guess[0]
            e_llimit=self.initial_guess[1]

            #should be irrelevant
            while e_llimit<0.75:
                _=self.initial_condition_errfunc((P_guess,e_llimit))
                if numpy.logical_or(numpy.isnan(self.solver_cache[(P_guess,e_llimit)]['final_orbital_period']),
                                    numpy.isnan(self.solver_cache[(P_guess,e_llimit)]['final_eccentricity'])):
                    _logger.info('\nIncreasing eccentricity lower limit by 0.02')
                    e_llimit+=0.02
                    continue
                else:
                    _logger.info('\nFound non-Nan eccenctricity lower limit at e={!r}'.format(e_llimit))
                    break
            
            
            p_root=self.orbital_period_solver(e_llimit,P_guess=P_guess)
            if numpy.isnan(p_root):de_llimit=scipy.nan
            else: de_llimit=self.delta_e
            #--------------------

            if de_llimit*de_ulimit>0 or numpy.isnan(de_llimit):
                _logger.info('\nNo solution cannot exist between eccentricity limit. Error lower limit = {!r} Error upper limit = {!r}'.format(de_llimit,de_ulimit))
                return scipy.nan 
            try:
                e_root=scipy.optimize.brentq(self.brent_eccentricity_func,e_llimit,e_ulimit,args=(P_guess,),xtol=1e-5,rtol=1e-6)
                _logger.info('\nEccentricity root found. E_ROOT={!r}'.format(e_root))
                return self.spin
            except:
                _logger.warning('Eccenricity Solver crashed')
                return scipy.nan

        else:

            P_guess=self.initial_guess[0]
            e_llimit=self.initial_guess[1]

            _logger.info('\nFinding lower limit for eccentricity')
            p_root=self.orbital_period_solver(e_llimit,P_guess=P_guess)
            if numpy.isnan(p_root):
                _logger.info('\nSolver cannot find orbital period at lower eccentricity limit = {!r}. No Solution'.format(e_llimit))
                return scipy.nan
            else:
                _logger.info('\nFound orbital period P={!r} at lower eccentricity limit e={!r}'.format(p_root,e_llimit))
            
            
            last_cached_ic=list(self.solver_cache.keys())[-1]
            if abs(self.solver_cache[last_cached_ic]['delta_e'])<1e-4:
                _logger.info('Solution already found. Skipping eccenrticity solver')
                return self.solver_cache[last_cached_ic]['spin']
            
            _logger.info('\nFinding upper limit for eccentricity')
            de_initial=self.solver_cache[last_cached_ic]['delta_e']
            de_new=de_initial
            eccentricity=e_llimit
            P_guess=p_root
            
            while de_new*de_initial>0 or numpy.isnan(de_new): 
                
                if abs(de_new)<1e-4:
                    break

                last_cached_ic=list(self.solver_cache.keys())[-1]
                final_e=self.solver_cache[last_cached_ic]['final_eccentricity']
                if final_e<1e-4:
                    eccentricity,de_new=self.check_if_no_solution()
                    break
                
                if numpy.isnan(de_new):increment=0.1
                elif self.solver_cache[last_cached_ic]['final_eccentricity']<last_cached_ic[1]:increment=0.1#min(abs(2.5*de_new),0.1)
                else: increment=-min(abs(1.5*de_new),0.1)
                
                if eccentricity>0.65 and eccentricity<0.75:eccentricity=0.75
                else: eccentricity+=increment
                
                if numpy.logical_or(eccentricity>0.75,eccentricity<0):
                    _logger.warning('eccentricity limit reached, e={!r}, while incrementing eccentriciy to find upper limit'.format(eccentricity))
                    break

                _logger.info('Increasing upper eccentricity limit by {!r}. e_new={!r}'.format(increment,eccentricity))

                while P_guess<60:
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
                    _logger.info('\nSolver cannot find orbital period at upper eccentricity limit = {!r}. No Solution'.format(eccentricity))
                    return scipy.nan

                last_cached_ic=list(self.solver_cache.keys())[-1]
                de_new=self.solver_cache[(last_cached_ic[0],last_cached_ic[1])]['delta_e']
                _logger.info('\nFound orbital period P={!r} at upper eccentricity limit e={!r} with de={!r}'.format(p_root,eccentricity,de_new))

                
            if numpy.isnan(de_new):
                _logger.info('\ndelta_e at upper limit is NaN. No solution found for eccentricity')
                return scipy.nan
            else: e_ulimit=eccentricity

            if abs(de_new)<1e-4:
                last_cached_ic=list(self.solver_cache.keys())[-1]
                _logger.info('Solution already found. Skipping eccenrticity solver')
                return self.solver_cache[last_cached_ic]['spin']

            if de_new*de_initial>0:
                _logger.info('\nBounds are of same sign. No solution found for eccentricity')
                return scipy.nan

            _logger.info('\nFinding eccentricity root between e_a={!r} and e_b={!r}'.format(e_llimit,e_ulimit))
            try:
                e_root=scipy.optimize.brentq(self.brent_eccentricity_func,e_llimit,e_ulimit,args=(P_guess,),xtol=1e-4,rtol=1e-5)
                return self.spin
            except:
                _logger.warning('\nEccentriciy solver crashed')
                return scipy.nan

            

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
        
        age_array=numpy.linspace(5e-3,self.target_age,10000)
        for _quantity in [self.primary_mass,self.secondary_mass]:
            I=self.interpolator('iconv',_quantity,self.feh)(age_array)
            if min(I)<0:
                _logger.warning('\nConvective envelope moment of inertia goes to less than zero I_min={!r} for at age={!r}'.format(min(I),age_array[numpy.argmin(I)]))
                _logger.info('Final_Spin_Period=inf')
                return numpy.inf

        _logger.info('\nSolving for p and e using function = {} and method {}'.format(self.function,self.method))

        self._calculate_good_initial_guess()
        
        t=time.time()
        spin=self.nested_solver()
        print('time taken={}'.format(time.time()-t))

        cached_ic_tuples=list(self.solver_cache.keys())
        initial_conditions_results=cached_ic_tuples[-1]
        _logger.info('Solver_Results:')
        _logger.info('Intial_Orbital_Period={!r} , Initial_Eccentricity={!r}'.format(initial_conditions_results[0],initial_conditions_results[1]))
        _logger.info('Final_Orbital_Period={!r} , Final_Eccentricity={!r}'.format(self.final_orbital_period,self.final_eccentricity))
        _logger.info('Errors: delta_p={!r} , delta_e={!r}'.format(self.delta_p,self.delta_e))
        _logger.info('Final_Spin_Period={!r}'.format(spin))

        if self.delta_p>0.1 or self.delta_e>0.1:
            _logger.warning('Solver did not find good solutions. delta_p={!r} delta_e={!r}'.format(self.delta_p,self.delta_e))

        if numpy.isnan(spin):
            _logger.warning('Spin=Nan after solver. Setting Spin=inf')
            spin=numpy.inf

        return spin

