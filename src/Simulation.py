"""
Created on 23.05.2013

@author: Johannes Rieke
"""

import copy
import numpy
import matplotlib.pyplot as plot
import matplotlib.animation as animation
from src.Neuron import *
import tkinter as tk
import time

neurons = []

################################################################
# neuron setup
################################################################
n1 = Neuron(25., 25., 0.3, 1., 20.)
n1.set_resting_potential(-70.)
n1.add_axon(10., 300., 50, 120., 35., 36., -82., -58.997)
n1.add_dendrite(20., 20., 5)
n1.add_dendrite(20., 20., 5, '0')
n1.add_dendrite(10., 30., 5)
neurons.append(n1)

n2 = Neuron(25., 25., 0.3, 1., 20.)
n2.set_resting_potential(-70.)
n2.add_axon(10., 200., 40, 120., 35., 36., -82., -58.997)
n2.add_dendrite(10., 50., 10)
neurons.append(n2)

n1.add_AMPA_Synapse(10000., 10., n2, 0, '0')
n1.add_NMDA_synapse(10000., 10., n2, 0.2, '0')
n1.add_GABA_A_synapse(10000., 10., n2, 0.4, '0')
n1.add_GABA_B_synapse(10000., 10., n2, 0., '0')
################################################################


elapsedTime = 0.
startTime = 0.
timeStep = 0.04
timeStepsBetweenPlots = 5

paused = False

current = 100e3  # muA, positive charges
# only for 'plotTemporalProgress' and 'simulate'
currentStartTime = 20. # ms
currentEndTime = 24. # ms


# create window to display additional information for selected data points
master = tk.Tk()
master.title('Data Browser')


def create_log():
    s = ''
    for i, n in enumerate(neurons):
        s += n.neuron_to_string(previous_indentation='', identifier=str(i), include_substructures=True)
        s += '\n'
    return s

'''
Print log with the current neuron setup to a text file or to the console
@param filename - empty for console
'''
def print_log(filename=''):
    if filename:
        txt = open(filename, 'w')
        for i, n in enumerate(neurons):
            txt.write(n.neuron_to_string(previous_indentation='', identifier=str(i), include_substructures=True))
        txt.close()
    else:
        for i, n in enumerate(neurons):
            print(n.neuron_to_string(previous_indentation='', identifier=str(i), include_substructures=True))

print_log()  # 'C:\\Users\\Johannes\\Desktop\\Log.txt')


axes = []
f1 = plot.figure('Membrane Potential')

# create one subplot for each neuron and append it to 'axes'
num_plots = len(neurons)
for i in range(num_plots):
    ax = f1.add_subplot(num_plots, 1, i+1)
    ax.set_ylabel('Voltage / mV')
    ax.set_xlabel('Distance from Soma / mum')
    ax.set_ylim([-100, 40])
    ax.set_xlim([-100, 300])
    axes.append(ax)


# TODO: Print on figure, not axes
time_text = axes[0].text(-80, 20, '', verticalalignment='top')
current_text = axes[0].text(-80, 0, '', color = 'r', verticalalignment='top')
speed_text = axes[0].text(-80, 10, 'Speed (time between frames): ', verticalalignment='top')

def update_speed_text():
    speed_text.set_text('Speed (time between frames): ' + str(round(timeStepsBetweenPlots * timeStep, 2)) + ' ms')

update_speed_text()


picked_element = None
picked_segment = -1
pick_tolerance = 5

picked_neuron_description = tk.StringVar()
picked_neuron_description_default = 'Pick data point\nto display further information\nRight click to unpick'
picked_neuron_description.set(picked_neuron_description_default)
tk.Label(master, textvariable=picked_neuron_description).pack()

picked_element_description = tk.StringVar()
picked_element_description_default = 'Space: Pause\nUp/Down: Change speed\nT: Single Step\nC: Current injection\n(into picked segment)\nLeft/Right: Pick next segment'
picked_element_description.set(picked_element_description_default)
tk.Label(master, textvariable=picked_element_description).pack()

tk.Label(master, text='Current injection (in mA):').pack()
current_description = tk.StringVar()
current_description.set('100')
tk.Entry(master, textvariable=current_description).pack()


# store one list for each neuron that contains all its cables (without hierarchy)
cable_matrix = []

# store one list for each neuron that contains the matplotlib lines for all its cables; indices as in 'cable_matrix'
line_matrix = []

# reuse colors periodically
colors = ['b', 'g', 'r', 'c', 'm', 'y']
while len(colors) < len(neurons):
    colors.extend(colors)

# plot initial setup
for n, ax, c in zip(neurons, axes, colors):
    lines = []
    cables = []
    
    somaLine, = ax.plot(0, n.potential[0], c+'o', picker = pick_tolerance)
    lines.append(somaLine)
    cables.append(n)
    
    if n.axon:
        axonDistances = numpy.arange(n.axon.num_segments) * n.axon.segment_length #)[:] + 0.5 * _n.axon.segmentLength
        axonDistances[:] += 0.5 * n.axon.segment_length    # take middle of segment
        axonDistances[:] += n.axon.distance_to_soma
        axonLine, = ax.plot(axonDistances, n.axon.potential, c+'.', picker = pick_tolerance)
        lines.append(axonLine)
        cables.append(n.axon)
    
    # TODO: best way?
    maxCrossSectionalArea = 0.
    for d in n.dendrites:
        maxCrossSectionalArea = max(maxCrossSectionalArea, d.cross_sectional_area)
    
    # depth-first
    dendritesToVisit = copy.copy(n.dendrites)
    while dendritesToVisit: # returns True if list isn't empty
        currentDendrite = dendritesToVisit.pop(0)
        
        dendriteDistances = numpy.arange(0, -currentDendrite.num_segments, -1) * currentDendrite.segment_length
        dendriteDistances[:] -= 0.5 * currentDendrite.segment_length    # take middle of segment
        dendriteDistances[:] -= currentDendrite.distance_to_soma
        dendriteLine, = ax.plot(dendriteDistances, currentDendrite.potential, c+'.', picker = pick_tolerance, alpha = currentDendrite.cross_sectional_area / maxCrossSectionalArea)
        lines.append(dendriteLine)
        cables.append(currentDendrite)
        
        dendritesToVisit[:0] = currentDendrite._successors
        
    line_matrix.append(lines)
    cable_matrix.append(cables)
    
while elapsedTime < startTime:
    elapsedTime += timeStep
    for n in neurons:
        n.simulate(timeStep)

def singleStep(i):
    if not paused:
        global elapsedTime
        global insertCurrent
        global removeCurrent
        
        for step in range(timeStepsBetweenPlots):
            elapsedTime += timeStep
            for n in neurons:
                n.simulate(timeStep)
                              
        time_text.set_text('Time: ' + str(round(elapsedTime, 1)) + ' ms')
        
        if picked_element:
            updatePickedElementDescription()
            
        # plot new values
        for lines, cables in zip(line_matrix, cable_matrix):
            for line, cable in zip(lines, cables):
                line.set_ydata(cable.potential)
                
def insertCurrent():
    try:
        c = 1000. * float(current_description.get())   # current is entered as mA but needed for computation as muA
        if picked_element:
            picked_element.set_current(c, segment_index=picked_segment)
        else:   # default
            neurons[0].set_current(c, 0)
        current_text.set_text('Current Injection: ' + str(round(c / 1000., 1)) + ' mA')
    except ValueError:  # non-number digit was entered
        tk.messagebox.showerror('Current Injection', 'Please enter a valid number!')
        
def removeCurrent():
    if picked_element:
        picked_element.remove_current(segment_index=picked_segment)
    else:   # default
        neurons[0].remove_current(0)
    current_text.set_text('')

def onKeyPressed(event):
    global paused
    global timeStepsBetweenPlots
    global picked_segment

    print('key pressed', event.key)
    
    if event.key == ' ':
        paused = not paused
    elif event.key == 't' and paused:
        paused = False
        singleStep(0)
        paused = True
    elif event.key == 'right' and picked_element:
        removeCurrent()
        picked_segment = max(0, min(picked_segment + 1, picked_element.num_segments-1))
        updatePickedElementDescription()
    elif event.key == 'left' and picked_element:
        removeCurrent()
        picked_segment = max(0, min(picked_segment - 1, picked_element.num_segments-1))
        updatePickedElementDescription()
    elif event.key == 'up':
        timeStepsBetweenPlots += 1
        update_speed_text()
    elif event.key == 'down':
        timeStepsBetweenPlots = max(1, timeStepsBetweenPlots-1)
        update_speed_text()
    elif event.key == 'c':
        insertCurrent()
        
def onKeyRealeased(event):
    global removeCurrent
    if event.key == 'alt+c':
        removeCurrent()
        
        
def updatePickedElementDescription():
    picked_element_description.set(picked_element.to_string(segment_index=picked_segment))

def onPick(event):
    pickedLine = event.artist
    pickedPointIndex = event.ind[0]
    
    global picked_element
    global picked_segment
    global pickedIdentifier
    
    for neuronIndex, lines in enumerate(line_matrix):
        try:
            lineIndex = lines.index(pickedLine)
            
            # only reached when index function finds something            
            # removes the inserted current to prevent trouble with changing 'picked_element'
            removeCurrent()            
            picked_element = cable_matrix[neuronIndex][lineIndex]
            picked_segment = pickedPointIndex
            picked_neuron_description.set(neurons[neuronIndex].neuron_to_string(identifier=str(neuronIndex)))
            updatePickedElementDescription()            
            
            # leaves the loop as soon as picked element is found
            break
        except ValueError:
            # could not find picked element
            pass
    
def onButtonPressed(event):
    global picked_element
    global picked_segment
    if event.button == 3:   # right mouse button
        removeCurrent()
        picked_element = None
        picked_segment = -1
        picked_neuron_description.set(picked_neuron_description_default)
        picked_element_description.set(picked_element_description_default)
        
 
f1.canvas.mpl_connect('key_press_event', onKeyPressed)
f1.canvas.mpl_connect('key_release_event', onKeyRealeased)
f1.canvas.mpl_connect('button_press_event', onButtonPressed)
f1.canvas.mpl_connect('pick_event', onPick)               
ani = animation.FuncAnimation(f1, singleStep, interval = 1, repeat = True)
        

plot.show()
tk.mainloop()


'''
Simulate the neuron setup for a duration of 'runtime'; used with 'plotTemporalProgress'.
@return - time series arrays storing several quantities for each neuron
'''
def simulate(runtime, h):
    numSteps = int(round(runtime / h))
    numNeurons = neurons.__len__()
        
    elapsedTime = numpy.zeros(numSteps)
    somaVoltages = numpy.zeros((numSteps, numNeurons))
    
    axonVoltages = numpy.zeros((numSteps, numNeurons, 3))    
    axonNGates = numpy.zeros((numSteps, numNeurons, 3))    
    axonMGates = numpy.zeros((numSteps, numNeurons, 3))    
    axonHGates = numpy.zeros((numSteps, numNeurons, 3))      
    dendriteVoltages = numpy.zeros((numSteps, numNeurons, 3))   
    
    timeNew = 0.
    for step in range(numSteps):
        elapsedTime[step] = timeNew
                                
        for i in range(numNeurons):
            somaVoltages[step, i] = neurons[i].potential[0]
            
            for x, segment in enumerate([0, round(neurons[i].axon.numSegments / 2.), -1]):
                axonVoltages[step, i, x] = neurons[i].axon.potential[segment]
                axonMGates[step, i, x] = neurons[i].axon.mGate[segment]
                axonHGates[step, i, x] = neurons[i].axon.hGate[segment]
                axonNGates[step, i, x] = neurons[i].axon.nGate[segment]
                
            if len(neurons[i].dendrites) > 0:
                for y, segment in enumerate([0, round(neurons[i].dendrites[0].numSegments / 2.), -1]):
                    dendriteVoltages[step, i, y] = neurons[i].dendrites[0].potential[segment]
            
            
        if timeNew >= currentStartTime and timeNew < currentEndTime:
            neurons[0].dendrites[0].set_current(current, 0.1)
        elif timeNew >= currentEndTime:
            neurons[0].dendrites[0].remove_all_currents()
        
        for n in neurons:
#             p = Process(target=_n.simulate, args=(h,))
#             p.start()
#             p.join()
            n.simulate(h)
                        
        timeNew = elapsedTime[step] + h
            
    return elapsedTime, somaVoltages, axonVoltages, axonMGates, axonHGates, axonNGates, dendriteVoltages


'''
Simulate current neuron setup with short current injection for 100 ms with a 0.04 ms time step and print the time needed for computation.
'''
def performance_test():
    time_to_simulate = 100  # ms
    simulation_time_step = 0.04  # ms
    print('Simulating neuron setup and current injection for', time_to_simulate, 'ms with', simulation_time_step, 'ms time step...')
    time.clock()
    simulate(time_to_simulate, simulation_time_step)
    end_time = time.clock()
    print('Result:', end_time, '_s')

    
# performance_test()


def plotTemporalProgress(t, sv, av, am, ah, an, dv):
    f1 = plot.figure('Membrane Potential')
    ax1 = f1.add_subplot(111)
    ax1.set_ylabel('Voltage / mV')
    ax1.set_xlabel('Time / ms')
    ax1.axvspan(currentStartTime, currentEndTime, facecolor = 'b', linewidth = 0., alpha = 0.05)
        
    f2 = plot.figure('Gates')
    axM = f2.add_subplot(221)
    axN = f2.add_subplot(223)
    axH = f2.add_subplot(222)
    for ax, title in zip([axM, axH, axN], ['m', 'h', '_n']):
        ax.set_title(title)
        ax.set_ylim([0, 1])
        ax.axvspan(currentStartTime, currentEndTime, facecolor = 'b', linewidth = 0., alpha = 0.05)
        
    neuronColors = ['b', 'g']
    for neuron, color in zip(range(len(neurons)), neuronColors):
        ax1.plot(t, dv[:, neuron, 0], color+':', alpha = 0.6)
        ax1.plot(t, dv[:, neuron, 1], color+':', alpha = 0.4)
        ax1.plot(t, dv[:, neuron, 2], color+':', alpha = 0.2)
        
        ax1.plot(t, sv[:, neuron], color, alpha = 1.0)
        
        ax1.plot(t, av[:, neuron, 0], color, alpha = 0.45)
        ax1.plot(t, av[:, neuron, 1], color, alpha = 0.3)
        ax1.plot(t, av[:, neuron, 2], color, alpha = 0.15)
        
        for ax, a in zip([axM, axH, axN], [am, ah, an]):    
            ax.plot(t, a[:, neuron, 0], color, alpha = 0.45)
            ax.plot(t, a[:, neuron, 1], color, alpha = 0.3)
            ax.plot(t, a[:, neuron, 2], color, alpha = 0.15)
    
    plot.show()

# t, sv, av, am, ah, an, dv = simulate(100., 0.02)
# plotTemporalProgress(t, sv, av, am, ah, an, dv)

