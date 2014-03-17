'''
Created on 23.05.2013

@author: Johannes Rieke

potentials in mV
conductances in mS
capacitances in muF
currents in muA
lengths in mum
concentrations in mM

literature:
Hodgkin / Huxley (1952): A quantitative description of membrane current and its application to conduction and excitation in nerve
Destexhe / Mainen / Sejnowski (1998): Kinetic Models of Synaptic Transmission
Jahr / Stevens (1990): Voltage  Dependence  of  NMDA-Activated  Macroscopic  Conductances Predicted  by  Single-Channel  Kinetics
'''

import numpy
import Solver
from Alphabets import greek


''' Global calculations and conversions '''
def _calculate_membrane_conductance(specific_membrane_conductance, radius, length):
    return specific_membrane_conductance * 2. * numpy.pi * radius * length

def _calculate_membrane_capacitance(specific_membrane_capacitance, radius, length):
    return specific_membrane_capacitance * 2. * numpy.pi * radius * length

def _calculate_plasma_conductance(specific_plasma_conductance, radius, length):
    return specific_plasma_conductance * numpy.pi * radius**2 / length

# 'quantities' is an iterable of tuples of the format (name, value, unit)
def _quantities_to_string(quantities, prefix='', suffix='\n'):
    s = ''
    for name, value, unit in quantities:
        s += prefix + name + ': ' + str(round(value, 2)) + ' ' + unit + suffix
    return s




'''
Base class for Rall'_s Cable Equations.
'''
class Cable:

    distance_to_soma = 0  # dummy value, is set after creation
    _resting_potential = 0.  # mV
    
    def __init__(self, radius, length, num_segments, specific_membrane_conductance, specific_membrane_capacitance, specific_plasma_conductance):
                
        self.radius = radius
        self.length = length
        self.cross_sectional_area = numpy.pi * radius**2
        self.num_segments = num_segments
        self.segment_length = length / num_segments
        
        self._segment_membrane_conductance = _calculate_membrane_conductance(specific_membrane_conductance, radius, self.segment_length)
        self._segment_membrane_capacitance = _calculate_membrane_capacitance(specific_membrane_capacitance, radius, self.segment_length)
        self._segment_plasma_conductance = _calculate_plasma_conductance(specific_plasma_conductance, radius, self.segment_length)
        
        self.potential = numpy.zeros(num_segments)
        self._outside_current = numpy.zeros(num_segments)  # Stored locally, so it doesn't have to be recreated at every simulation step.
        
    def set_resting_potential(self, resting_potential):
        self._resting_potential = resting_potential
        self.potential[:] = resting_potential
    
    # 'relative_position' between 0. (fist segment) and 1. (last segment).
    def _get_segment_index(self, relative_position):
        return round(relative_position * (self.num_segments - 1))
    
    # Takes either 'relative_position' or 'segment_index'.
    def add_current(self, current, relative_position=None, segment_index=None):
        if relative_position is not None:
            segment_index = self._get_segment_index(relative_position)
        self._outside_current[segment_index] += current
    
    def set_current(self, current, relative_position=None, segment_index=None):
        if relative_position is not None:
            segment_index = self._get_segment_index(relative_position)
        self._outside_current[segment_index] = current
        
    def remove_current(self, relative_position=None, segment_index=None):
        if relative_position is not None:
            segment_index = self._get_segment_index(relative_position)
        self._outside_current[segment_index] = 0
        
    def remove_all_currents(self):
        self._outside_current[:] = 0


    '''
    Calculates currents and potentials for one 'timeStep'.    
    'currentAtStart' and 'currentAtEnd' are the plasma currents at the respective cable ends. 
    'extraCurrent' can be caused by action potentials, e. g.   
    '''
#     def simulate(self, timeStep, currentAtStart=0, currentAtEnd=0, extraCurrent=0):
#         # TODO: maybe change to one step (performance?) <-> problem: numpy.array not easy to enhance
#         leftDifferenceE = numpy.ediff1d(self.potential, to_begin = 0)
#         rightDifferenceE = - numpy.ediff1d(self.potential, to_end = 0)
#         
#         # the current from one segment to the next one
#         # second order central difference formula
#         boundaryCurrentE = self._segment_plasma_conductance * (leftDifferenceE + rightDifferenceE)
#         boundaryCurrentE[0] += currentAtStart
#         boundaryCurrentE[-1] += currentAtEnd
#         
#         membraneCurrentE = self._segment_membrane_conductance * (self.potential - self.restingPotential)
# 
#         potentialE = self.potential + timeStep * ((self._outside_current - boundaryCurrentE - membraneCurrentE - extraCurrent) / self._segment_membrane_capacitance)
#         
#         
#         # Heun'_s
#         
#         leftDifferenceH = numpy.ediff1d(potentialE, to_begin = 0)
#         rightDifferenceH = - numpy.ediff1d(potentialE, to_end = 0)
#         
#         # the current from one segment to the next one
#         # second order central difference formula
#         boundaryCurrentH = self._segment_plasma_conductance * (leftDifferenceE + rightDifferenceE)
#         boundaryCurrentH[0] += currentAtStart
#         boundaryCurrentH[-1] += currentAtEnd
#         
#         membraneCurrentH = self._segment_membrane_conductance * (potentialE - self.restingPotential)
#         
#         self.potential = self.potential + timeStep * 0.5 * ((self._outside_current - boundaryCurrentE - membraneCurrentE - extraCurrent  +  self._outside_current - boundaryCurrentH - membraneCurrentH - extraCurrent) / self._segment_membrane_capacitance)
        

    # TODO: Use with scipy and alter names.
    def _calculate_derivatives(self, initial_potential, current_at_start=0., current_at_end=0., extra_current=0.):
        # TODO: maybe change to one step (performance?) <-> problem: numpy.array not easy to enhance
        left_difference = numpy.ediff1d(initial_potential, to_begin = 0)
        right_difference = - numpy.ediff1d(initial_potential, to_end = 0)
         
        # the current from one segment to the next one
        # second order central difference formula
        boundary_current = self._segment_plasma_conductance * (left_difference + right_difference)
        boundary_current[0] += current_at_start
        boundary_current[-1] += current_at_end
         
        membrane_current = self._segment_membrane_conductance * (initial_potential - self._resting_potential)
 
        return (self._outside_current - boundary_current - membrane_current - extra_current) / self._segment_membrane_capacitance # new potential


    def simulate(self, time_step, current_at_start=0, current_at_end=0, extra_current=0):
        self.potential = Solver.default(time_step, self._calculate_derivatives, self.potential, current_at_start, current_at_end, extra_current)

    # TODO: Use greek letters for units
    def to_string(self, leading_indentation=''):
        quantities = [('Distance to Soma', self.distance_to_soma, greek['mu'] + 'm'),
                      ('Radius', self.radius, greek['mu'] + 'm'),
                      ('Length', self.length, greek['mu'] + 'm'),
                      ('Number of Segments', self.num_segments, ''),
                      ('Segment Length', self.segment_length, greek['mu'] + 'm'),
                      ('Segment Membrane Conductance', self._segment_membrane_conductance, 'mS'),
                      ('Segment Membrane Capacitance', self._segment_membrane_capacitance, greek['mu'] + 'F'),
                      ('Segment Plasma Conductance', self._segment_plasma_conductance, 'mS')]
        return _quantities_to_string(quantities, leading_indentation)
        

class Axon(Cable):
    
    # needed because axon is added to neuron with 'Dendrite.add_successor'
    # not used in any way
    dendrite_code = ''

    # TODO: Check out use of 'leak_potential'
    def __init__(self, radius, length, num_segments, specific_membrane_conductance, specific_membrane_capacitance, specific_plasma_conductance, specific_sodium_conductance, sodium_potential, specific_potassium_conductance, potassium_potential, leak_potential):
        
        Cable.__init__(self, radius, length, num_segments, specific_membrane_conductance, specific_membrane_capacitance, specific_plasma_conductance)
        self._segment_sodium_conductance = _calculate_membrane_conductance(specific_sodium_conductance, radius, self.segment_length)
        self._sodium_potential = sodium_potential
        self._segment_potassium_conductance = _calculate_membrane_conductance(specific_potassium_conductance, radius, self.segment_length)
        self._potassium_potential = potassium_potential
        self._leak_potential = leak_potential
        self.m_gate = numpy.empty(num_segments)
        self.h_gate = numpy.empty(num_segments)
        self.n_gate = numpy.empty(num_segments)
        
        # Initial values, taken from experimental equilibrium, to prevent from firing at each start.
        self.m_gate[:] = 0.0530
        self.h_gate[:] = 0.5958
        self.n_gate[:] = 0.3178
        
        self._outgoing_synapses = []


    '''
    Calculates currents and potentials for one 'timeStep' and passes the command to all '_outgoing_synapses'.
    Axon terminal is modeled as sealed end, therefore not current at end.
    ''' 
    def _calculate_derivatives(self, initial_potential, initial_m_gate, initial_h_gate, initial_n_gate, current_at_start=0):
        u = initial_potential - self._resting_potential
            
        # c. f. Hodgkin / Huxley 1952, p. 519
        alpha_n = 0.01 * (-u + 10.) / (numpy.exp((-u + 10.) / 10.) - 1.)
        beta_n = 0.125 * numpy.exp(-u / 80.)
        alpha_m = 0.1 * (-u + 25.) / (numpy.exp((-u + 25.) / 10.) - 1.)
        beta_m = 4. * numpy.exp(-u / 18.)
        alpha_h = 0.07 * numpy.exp(-u / 20.)
        beta_h = 1. / (numpy.exp((-u + 30.) / 10.) + 1.)
        
        m_gate_derivative = alpha_m * (1. - initial_m_gate) - beta_m * initial_m_gate
        h_gate_derivative = alpha_h * (1. - initial_h_gate) - beta_h * initial_h_gate
        n_gate_derivative = alpha_n * (1. - initial_n_gate) - beta_n * initial_n_gate
            
        # TODO: Should use new gate values instead?
        sodium_current = self._segment_sodium_conductance * initial_m_gate**3 * initial_h_gate * (initial_potential - self._sodium_potential)
        potassium_current = self._segment_potassium_conductance * initial_n_gate**4 * (initial_potential - self._potassium_potential)
        leak_current = self._segment_membrane_conductance * (initial_potential - self._leak_potential)
                            
        potential_derivative = Cable._calculate_derivatives(self, initial_potential, current_at_start, 0, extra_current=(sodium_current + potassium_current + leak_current))
        
        return potential_derivative, m_gate_derivative, h_gate_derivative, n_gate_derivative
    
        
    def simulate(self, time_step, current_at_start=0):
        simulated_values = Solver.default(time_step, self._calculate_derivatives, self.potential, self.m_gate, self.h_gate, self.n_gate, current_at_start)
        for s in self._outgoing_synapses:
            s.simulate(time_step, self.potential[-1])
            
        self.potential, self.m_gate, self.h_gate, self.n_gate = simulated_values
        
            
        
    def add_outgoing_synapse(self, synapse):
        self._outgoing_synapses.append(synapse)
        
    def to_string(self, leading_indentation='', segment_index=None):
        s = leading_indentation + 'Axon\n'
        s += Cable.to_string(self, leading_indentation)

        axon_quantities = [('Segment Maximum Sodium Conductance', self._segment_sodium_conductance, 'mS'),
                           ('Segment Sodium Potential', self._sodium_potential, 'mV'),
                           ('Segment Maximum Potassium Conductance', self._segment_potassium_conductance, 'mS'),
                           ('Potassium Potential', self._potassium_potential, 'mV'),
                           ('Leak Potential', self._leak_potential, 'mV')]
        s += _quantities_to_string(axon_quantities, leading_indentation)

        # TODO: maybe unify with methods in Cable and Dendrite
        if segment_index:
            s += '\n' + leading_indentation + 'Segment ' + str(segment_index) + '\n'
            segment_quantities = [('Potential', self.potential[segment_index], 'mV'),
                                  ('m', self.m_gate[segment_index], ''),
                                  ('h', self.h_gate[segment_index], ''),
                                  ('n', self.n_gate[segment_index], '')]
            s += _quantities_to_string(segment_quantities, leading_indentation)

        return s
    

class Dendrite(Cable):
    
    _successors_cross_sectional_area = 0.  # mum^2
    dendrite_code = ''
    
    def __init__(self, radius, length, num_segments, specific_membrane_conductance, specific_membrane_capacitance, specific_plasma_conductance):
        Cable.__init__(self, radius, length, num_segments, specific_membrane_conductance, specific_membrane_capacitance, specific_plasma_conductance)
        self._successors = []
        self._incoming_synapses = []
        self._synaptic_current = numpy.zeros(num_segments)  # Stored locally, so it doesn't have to be recreated at every simulation step.
        
    def set_resting_potential(self, resting_potential):
        Cable.set_resting_potential(self, resting_potential)
        for successor in self._successors:
            successor.set_resting_potential(resting_potential)
        
    '''
    'pass_code' makes the 'dendrite' being passed to further _successors, format like '01032', count from zero.
    '''
    def add_successor(self, dendrite, pass_code=''):
        if pass_code:
            self._successors[int(pass_code[0])].add_successor(dendrite, pass_code[1:])
        else:
            # TODO
            dendrite.distance_to_soma = self.distance_to_soma + self.length
            self._successors.append(dendrite)
            self._successors_cross_sectional_area += dendrite.cross_sectional_area
            dendrite.dendrite_code += str(len(self._successors)-1)
    
    def add_incoming_synapse(self, synapse, relative_position, pass_code=''):
        if pass_code:
            self._successors[int(pass_code[0])].add_incoming_synapse(synapse, relative_position, pass_code[1:])
        else:
            synapse.segment_index = self._get_segment_index(relative_position)
            self._incoming_synapses.append(synapse)
    
    '''
    Calculates currents and potentials for one 'timeStep' and passes the command to all '_successors'.
    current at end is evaluated by number and size of _successors.
    ''' 
    def _calculate_derivatives(self, initial_potential, current_at_start=0., current_at_end=0.):
        # resetting the array to zero
        self._synaptic_current[:] = 0
        for synapse in self._incoming_synapses:
            self._synaptic_current[synapse.segment_index] += synapse.conductance * (initial_potential[synapse.segment_index] - synapse.reversal_potential)
                    
        return Cable._calculate_derivatives(self, initial_potential, current_at_start, current_at_end, extra_current=self._synaptic_current)
        
    def simulate(self, time_step, current_at_start=0.):
        
        # TODO: no good place for this, but can't be done in '_calculate_derivatives' due to several calls
        # TODO: correct for neuron? models dendrites and axon to be just on one side of the soma-cable
        current_at_end = 0.
        for s in self._successors:
            cur = 0.5 * (s._segment_plasma_conductance + s.cross_sectional_area / self._successors_cross_sectional_area * self._segment_plasma_conductance) * (- s.potential[0] + self.potential[-1])
            current_at_end += cur
            s.simulate(time_step, -cur)
        
        simulated_values = Solver.default(time_step, self._calculate_derivatives, self.potential, current_at_start, current_at_end)
        self.potential = simulated_values
            
        
                 
    # TODO: maybe abolish 'include_successors' and iterate over complete cable list in 'Neuron' <-> indentations
    def to_string(self, leading_indentation='', segment_index=None, include_successors=False):
        s = leading_indentation + 'Dendrite ' + self.dendrite_code + '\n'
        s += Cable.to_string(self, leading_indentation)
        if segment_index:
            s += '\n' + leading_indentation + 'Segment ' + str(segment_index) + '\n'
            quantities = [('Potential', self.potential[segment_index], 'mV')]
            s += _quantities_to_string(quantities, leading_indentation)
            for synapse in self._incoming_synapses:
                if synapse.segment_index == segment_index:
                    s += '\n' + synapse.to_string(leading_indentation)
        if include_successors:
            for successor in self._successors:
                s += '\n' + successor.to_string(leading_indentation+'\t')
        return s



class Neuron(Dendrite):
            
    def __init__(self, soma_radius, soma_length, specific_membrane_conductance, specific_membrane_capacitance, specific_plasma_conductance):
        self._specific_membrane_conductance = specific_membrane_conductance
        self._specific_membrane_capacitance = specific_membrane_capacitance
        self._specific_plasma_conductance = specific_plasma_conductance
        
        self.axon = None
        self.dendrites = []

        # TODO: Make easier
        self.distance_to_soma = - soma_length / 2. # needed for further calculations in dendrites and axon; "soma" is in the middle of the soma-cable
        
        Dendrite.__init__(self, soma_radius, soma_length, 1, specific_membrane_conductance, specific_membrane_capacitance, specific_plasma_conductance)
          
    def simulate(self, time_step):
        Dendrite.simulate(self, time_step, 0)


    # TODO: Add functions like add_dendrite(dendrite) etc.

    def add_dendrite(self, radius, length, num_segments, dendrite_code=''):
        d = Dendrite(radius, length, num_segments, self._specific_membrane_conductance, self._specific_membrane_capacitance, self._specific_plasma_conductance)
        d.set_resting_potential(self._resting_potential)
        if dendrite_code:
            d.dendrite_code = dendrite_code
            self.dendrites[int(dendrite_code[0])].add_successor(d, dendrite_code[1:])
        else:
            self.dendrites.append(d)
            Dendrite.add_successor(self, d)
            d.dendrite_code = str(len(self.dendrites)-1)
            #d.distance_to_soma = self.distance_to_soma - d.length

                    
    # only to be called if neuron has axon
    def add_synapse(self, synapse, post_synaptic_neuron, relative_position=0.5, dendrite_code=''):
        self.axon.add_outgoing_synapse(synapse)
        post_synaptic_neuron.add_incoming_synapse(synapse, relative_position, dendrite_code)
    
    def add_AMPA_Synapse(self, maximum_conductance, maximum_transmitter_concentration, post_synaptic_neuron, relative_position=0.5, dendrite_code=''):
        self.add_synapse(AMPASynapse(maximum_conductance, maximum_transmitter_concentration), post_synaptic_neuron, relative_position, dendrite_code)
        
    def add_NMDA_synapse(self, maximum_conductance, maximum_transmitter_concentration, post_synaptic_neuron, relative_position=0.5, dendrite_code=''):
        self.add_synapse(NMDASynapse(maximum_conductance, maximum_transmitter_concentration), post_synaptic_neuron, relative_position, dendrite_code)
        
    def add_GABA_A_synapse(self, maximum_conductance, maximum_transmitter_concentration, post_synaptic_neuron, relative_position=0.5, dendrite_code=''):
        self.add_synapse(GABA_ASynapse(maximum_conductance, maximum_transmitter_concentration), post_synaptic_neuron, relative_position, dendrite_code)
        
    def add_GABA_B_synapse(self, maximum_conductance, maximum_transmitter_concentration, post_synaptic_neuron, relative_position=0.5, dendrite_code=''):
        self.add_synapse(GABA_BSynapse(maximum_conductance, maximum_transmitter_concentration), post_synaptic_neuron, relative_position, dendrite_code)
                
    def add_incoming_synapse(self, synapse, relative_position=0.5, dendrite_code=''):
        if dendrite_code:
            self.dendrites[int(dendrite_code[0])].add_incoming_synapse(synapse, relative_position, dendrite_code[1:])
        else:
            synapse.segment_index = 0
            self._incoming_synapses.append(synapse)
            
    # Only one axon per neuron.
    def add_axon(self, radius, length, num_segments, specific_sodium_conductance, sodium_potential, specific_potassium_conductance, potassium_potential, leak_potential):
        a = Axon(radius, length, num_segments, self._specific_membrane_conductance, self._specific_membrane_capacitance, self._specific_plasma_conductance, specific_sodium_conductance, sodium_potential, specific_potassium_conductance, potassium_potential, leak_potential)
        a.set_resting_potential(self._resting_potential)
        self.axon = a
        # TODO
        Dendrite.add_successor(self, a)
        a.distance_to_soma = self.distance_to_soma + self.length


    # TODO: Change this into 'to_string' and maybe create 'soma_to_string' instead
    def neuron_to_string(self, identifier='', previous_indentation='', include_substructures=False):
        s = previous_indentation + 'Neuron'
        if identifier:
            s += ' ' + identifier
        s += '\n'
        v = [('Specific Membrane Conductane', self._specific_membrane_conductance, 'mS'), ('Specific Membrane Capacitance', self._specific_membrane_capacitance, greek['mu'] + 'F'), ('Specific Plasma Conductance', self._specific_plasma_conductance, 'mS')]
        s += _quantities_to_string(v, previous_indentation)
        
        if include_substructures:
            s += '\n'
            s += self.to_string(leading_indentation= previous_indentation + '\t') + '\n'
#             _s += Dendrite.to_string(self, previous_indentation = previous_indentation + '\t', includeSuccessors = True)
            if self.axon is not None:
                s += self.axon.to_string(previous_indentation+'\t') + '\n'
            for d in self.dendrites:
                s += d.to_string(previous_indentation+'\t', include_successors=True) + '\n'
        
        return s
    
    def to_string(self, leading_indentation='', segment_index=-1):
        s = leading_indentation + 'Soma\n'
        s += Cable.to_string(self, leading_indentation)
        if segment_index >= 0:
            s += '\n'
            v = [('Potential', self.potential[0], 'mV')]
            s += _quantities_to_string(v, leading_indentation)
        return s


'''
Base class for different synapse types.
'''
class Synapse:

    # TODO: unit
    conductance = 0.
    segment_index = 0  # dummy value, is set after creation

    # currently implemented with static values
    # c. f. Destexhe 1998
    _Vp = 2.  # mV
    _Kp = 5.  # mV
    
    def __init__(self, maximum_conductance, maximum_transmitter_concentration, reversal_potential):
        self._maximum_conductance = maximum_conductance
        self._maximum_transmitter_concentration = maximum_transmitter_concentration # mM
        self.reversal_potential = reversal_potential
        
    # Smooth step
    def get_transmitter_concentration(self, potential):
        return self._maximum_transmitter_concentration / (1. + numpy.exp(- (potential - self._Vp) / self._Kp)) # mM
    
    def to_string(self, previous_indentation=''):
        quantities = [('Reversal Potential', self.reversal_potential, 'mV'),
                      ('Conductance', self.conductance, 'mS')]
        return _quantities_to_string(quantities, previous_indentation)
            

class AMPASynapse(Synapse):
    
    _r = 0.
    
    # currently implemented with static values
    # c. f. Destexhe 1998
    _alpha_r = 1.1  # ms^-1 mM^-1
    _beta_r = 0.19  # ms^-1
    
    def __init__(self, maximum_conductance, maximum_transmitter_concentration):
        Synapse.__init__(self, maximum_conductance, maximum_transmitter_concentration, 0.)
        
    def calculate_derivatives(self, initial_r, potential):
        r_derivative = self._alpha_r * Synapse.get_transmitter_concentration(self, potential) * (1. - initial_r) - self._beta_r * initial_r
        return r_derivative
        
    def simulate(self, time_step, potential):
        simulated_values = Solver.default(time_step, self.calculate_derivatives, self._r, potential)
        self.conductance = self._r * self._maximum_conductance
        self._r = simulated_values
        
    def to_string(self, previous_indentation=''):
        s = previous_indentation + 'Synapse, AMPA\n'
        s += Synapse.to_string(self, previous_indentation)
        quantities = [('r', self._r, '')]
        s += _quantities_to_string(quantities, previous_indentation)
        return s
    
    
class NMDASynapse(Synapse):
    
    _r = 0.
    _magnesium_block = 0.  # Instantiated here to include in 'to_string' method
    
    _alpha_r = 0.072  # ms^-1 mM^-1
    _beta_r = 0.0066  # ms^-1
    _extracellular_magnesium_concentration = 0.8 # mM, in physiological conditions between 1 and 2 mM, c. f. Destexhe 1998
    
    def __init__(self, maximum_conductance, maximum_transmitter_concentration):
        Synapse.__init__(self, maximum_conductance, maximum_transmitter_concentration, 0.)
        
    def calculate_derivatives(self, initial_r, potential):
        r_derivative = self._alpha_r * Synapse.get_transmitter_concentration(self, potential) * (1. - initial_r) - self._beta_r * initial_r
        return r_derivative
 
    def simulate(self, time_step, potential):
        simulated_values = Solver.default(time_step, self.calculate_derivatives, self._r, potential)
        self.conductance = self._r * self._magnesium_block * self._maximum_conductance
        # c. f. Jahr 1990
        self._magnesium_block = 1. / (1. + numpy.exp(-0.062 * potential) * (self._extracellular_magnesium_concentration / 3.57))
        self._r = simulated_values
        
    def to_string(self, previous_indentation=''):
        s = previous_indentation + 'Synapse, NMDA\n'
        s += Synapse.to_string(self, previous_indentation)
        quantities = [('r', self._r, ''),
                      ('Magnesium Block', self._magnesium_block, '')]
        s += _quantities_to_string(quantities, previous_indentation)
        return s
        

class GABA_ASynapse(Synapse):
    
    _r = 0.
    
    _alpha_r = 5. # ms^-1 mM^-1
    _beta_r = 0.18 # ms^-1
    
    def __init__(self, maximum_conductance, maximum_transmitter_concentration):
        Synapse.__init__(self, maximum_conductance, maximum_transmitter_concentration, -80.)
        
    # equivalent to AMPA type
    def calculate_derivatives(self, initial_r, potential):
        r_derivative = self._alpha_r * Synapse.get_transmitter_concentration(self, potential) * (1. - initial_r) - self._beta_r * initial_r
        return r_derivative
        
    def simulate(self, timeStep, potential):
        simulated_values = Solver.default(timeStep, self.calculate_derivatives, self._r, potential)
        self.conductance = self._r * self._maximum_conductance
        self._r = simulated_values
        
    def to_string(self, previous_indentation=''):
        s = previous_indentation + 'Synapse, GABA_A\n'
        s += Synapse.to_string(self, previous_indentation)
        quantities = [('r', self._r, '')]
        s += _quantities_to_string(quantities, previous_indentation)
        return s
        
        
class GABA_BSynapse(Synapse):
    
    _r = 0.
    _s = 0.  # muM
    
    _Kd = 100.  # muM^4
    _K1 = 0.09  # ms^-1 mM^-1
    _K2 = 0.0012  # ms^-1
    _K3 = 0.18  # muM ms^-1
    _K4 = 0.034  # ms^-1
    _n = 4  # binding sites
    
    def __init__(self, maximum_conductance, maximum_transmitter_concentration):
        Synapse.__init__(self, maximum_conductance, maximum_transmitter_concentration, -95.)
    
    def calculate_derivatives(self, initial_r, initial_s, potential):
        r_derivative = self._K1 * Synapse.get_transmitter_concentration(self, potential) * (1. - initial_r) - self._K2 * initial_r
        s_derivative = self._K3 * initial_r - self._K4 * initial_s  # uses muM for calculation!!
        return r_derivative, s_derivative
        
    def simulate(self, time_step, potential):
        simulated_values = Solver.default(time_step, self.calculate_derivatives, self._r, self._s, potential)
        self.conductance = self._s**self._n / (self._s**self._n + self._Kd) * self._maximum_conductance
        self._r, self._s = simulated_values
        
    def to_string(self, previous_indentation=''):
        s = previous_indentation + 'Synapse, GABA_B\n'
        s += Synapse.to_string(self, previous_indentation)
        quantities = [('r', self._r, ''),
                      ('s', self._s, greek['mu'] + 'M')]
        s += _quantities_to_string(quantities, previous_indentation)
        return s
            
