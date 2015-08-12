NeuroSim
========

Simple neuron simulator using NumPy and matplotlib.<br>
Can simulate multi-compartment neurons with soma, axon, dendrites and synapses. 
Implements Hodgkin-Huxley equations (<a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/pdf/jphysiol01442-0106.pdf">Hodgkin / Huxley 1952</a>), Rall's cable theory (<a href="http://stg.rutgers.edu/courses/old/CompNeuro07/Handouts/Rall%20-%20Core%20conductor%20theory.pdf">Rall 1977</a>) and several synapse models (<a href="http://cnl.salk.edu/~alain/abstracts/KSchap96.html">Destexhe / Mainen / Sejnowsky 1998</a>). 
Supports Euler's and Heun's algorithms as differential equation solvers.<br>
Visualizes the potentials in an animated matplotlib graph and allows direct user interaction (e. g. insert current, show segment parameters).
<br><br>
Run src/simulation.py to start. Requries Python 3.3, NumPy and matplotlib. Published under the MIT License.
<br><hr>
User interface with voltage graph and data browser:<br><br>
<img src="/gui.png">
