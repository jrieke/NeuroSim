'''
Methods that implement different numerical solvers for ordinary differential equations.

Created on 30.10.2013
@author: Johannes Rieke
'''

def euler(time_step, derivatives_func, *args):
    euler_derivatives = derivatives_func(*args)
    # Check whether multiple values are returned (as a tuple).
    if isinstance(euler_derivatives, tuple):
        # 'zip' trims 'args' to the length of 'euler_derivatives' and thus to the number of variable values.
        euler = tuple([initial + time_step * derivative for initial, derivative in zip(args, euler_derivatives)])
    else:
        euler = args[0] + time_step * euler_derivatives
        
    return euler


def heun(time_step, derivatives_func, *args):
    euler_derivatives = derivatives_func(*args)
    # Check whether multiple values are returned (as a tuple).
    if isinstance(euler_derivatives, tuple):
        # 'zip' trims 'args' to the length of 'euler_derivatives' and thus to the number of variable values.
        euler = tuple([initial + time_step * derivative for initial, derivative in zip(args, euler_derivatives)])
        heun_derivatives = derivatives_func(*(euler + args[len(euler_derivatives):]))
        heun = tuple([initial + time_step * 0.5 * (eulerDerivative + heunDerivative) for initial, eulerDerivative, heunDerivative in zip(args, euler_derivatives, heun_derivatives)])
    else:
        euler = args[0] + time_step * euler_derivatives
        heun_derivatives = derivatives_func(euler, *args[1:])
        heun = args[0] + time_step * 0.5 * (euler_derivatives + heun_derivatives)
        
    return heun


'''
Convenient way to switch solvers in the whole code at once.
'''
def default(timeStep, derivativesFunc, *args):
    return euler(timeStep, derivativesFunc, *args)
    
