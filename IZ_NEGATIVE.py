#%%
from tkinter import Y
from brian2 import *
from scipy import stats
import matplotlib.pyplot as plt

plt.style.use("seaborn")
#----------------------------------SECOND_BIFURCATION_PLOT-------------------------
defaultclock.dt = 0.03*ms
N = 200

# define izhikevich model
model = '''dvm/dt = (0.04/ms/mV)*vm**2+(5/ms)*vm+140*mV/ms-w + I :  volt
          dw/dt = a*(b*vm-w) :                                      volt/second
          d :                                                       volt/second
          I :                                                       volt/second '''

threshold = "vm >= 30*mV"

reset = "vm = c; w += d"

#model parameters associated with the system for analysis 
I = -99*volt/second
a = 0.2 / ms
b = 2 / ms
c = -55 * mV

#Starting the simulation 
start_scope()
neuron_bifurcation_negative = NeuronGroup(N, model=model, threshold=threshold,
                    reset=reset, method='euler')


neuron_bifurcation_negative.d = np.linspace(-18*volt/second, -1.0*volt/second, N)
neuron_bifurcation_negative.I = I

Initial_run_time_frame = 1*second
run(Initial_run_time_frame, report='text')

states_of_system = StateMonitor(neuron_bifurcation_negative, ["w", "vm"], record=True, when='start')
spikes_in_system = SpikeMonitor(neuron_bifurcation_negative)
run(1*second, report='text')

# Get the values of V and u for each spike
D = neuron_bifurcation_negative.d[spikes_in_system.i]
w = states_of_system.w[spikes_in_system.i, int_((spikes_in_system.t-Initial_run_time_frame)/defaultclock.dt)]
vm = states_of_system.w[spikes_in_system.i, int_((spikes_in_system.t-Initial_run_time_frame)/defaultclock.dt)]

plt.figure(figsize=(16,9))
plt.scatter(D / volt/second, w / volt/second, marker=".", color="k", linewidths=0.1)
plt.xlabel('d (volt/second)', fontsize = 14)
plt.ylabel('w [volt/second]', fontsize = 14)
plt.title("Bifurcation diagram of d vs u on system dynamics", fontsize = 18)
plt.show()

#%%
#---------------------------------STABILITY_ANALYSIS_-------------------------------
#%%
#imports
from brian2 import mV, ms, volt, second, umetre, ufarad, siemens, cm, msiemens, amp, uA, nA
from brian2 import start_scope, NeuronGroup, StateMonitor, run

import numpy as np
from sympy import symbols, solve, nsolve, lambdify, sympify, dsolve, Eq, solveset, linear_eq_to_matrix, nonlinsolve, Matrix, diff, sqrt, exp
import sympy as smp

#%%

def examine_stable_points(expression_for_v, expression_for_u , x_var, y_var, solutions):
  # First, calculate the jacobian matrix for our equations
  equation_matrix = Matrix([expression_for_v, expression_for_u])
  var_mat = Matrix([x_var, y_var])
  jacobian = equation_matrix.jacobian(var_mat)
  # Set up list that contains all stable points
  stable_points = []
  # Calculate eigenvalues for each of the stablepoints
  for stable_point in solutions:
    # Eigenvalue calculation
    eqmat = jacobian.subs([(x_var, stable_point[0]), (y_var, stable_point[1])])
    eigenvalues = list(eqmat.eigenvals().keys())
    # Check the eigenvalues to determine type of stable point
    if eigenvalues[0].is_real:
        if eigenvalues[0] > 0 and eigenvalues[1] > 0:
            stable_point_type = 'Unstable Node'
        elif eigenvalues[0] < 0 and eigenvalues[1] < 0:
            stable_point_type = 'Stable Node'
        elif (eigenvalues[0] < 0 and eigenvalues[1] > 0) or (eigenvalues[0] > 0 and eigenvalues[1] < 0):
            stable_point_type = 'Saddle Point'
    else:
        if eigenvalues[0].args[0] > 0:
            stable_point_type = 'Unstable Focus'
        if eigenvalues[0].args[0] < 0:
            stable_point_type = 'Stable Focus'
    # Add tuple for each stable point to list
    stable_points.append((stable_point, stable_point_type))
  # Return list of stable points
  return stable_points

def solve_dynamical_system(expression_for_v, expression_for_u
, x_range, y_range, x_var, y_var): 
  # Convert our equations to functions
  f1 = lambdify((x_var, y_var), expression_for_v)
  f2 = lambdify((x_var, y_var), expression_for_u)
  # Define range that we will examine the system in
  start_x = x_range[0]
  end_x = x_range[1]
  start_y = y_range[0]
  end_y = y_range[1]
  # Set up list of x and y values inside this range 
  x_range = np.linspace(start_x, end_x,1000)
  y_range = np.linspace(start_y, end_y,1000)
  x1_range = np.linspace(start_x, end_x)
  y1_range = np.linspace(start_y, end_y)
  total_range = np.linspace(min(start_x, start_y), max(end_x, end_y),1000)
  # Compute quivers by calculating the expression for a combination of x and y values
  f1_val = [[f1(x_cur, y_cur) for x_cur in x1_range] for y_cur in y1_range];
  f2_val = [[f2(x_cur, y_cur) for x_cur in x1_range] for y_cur in y1_range];
  # Solve analytically using sympy
  solutions = smp.solve([smp.Eq(expression_for_v, 0), smp.Eq(expression_for_u
  , 0)], (x_var, y_var)) 
  # Calculate nullclines
  x_nullclines = calculate_nullclines(expression_for_v, y_var, x_var, total_range)
  y_nullclines = calculate_nullclines(expression_for_u
  , y_var, x_var, total_range)

  # Get stable points and check their type
  stable_points = examine_stable_points(expression_for_v, expression_for_u
  , x_var, y_var, solutions)
  # Return relevant results
  return x_range, y_range, x1_range, y1_range, f1_val, f2_val, x_nullclines, y_nullclines, stable_points

def calculate_nullclines(eq, solvar, plotvar, inputrange):
  # Set our equations to zero and solve them
  eq = Eq(eq, 0)
  sol = solve(eq, solvar)
  nullclines = []
  # Calculate y-values for each of the solutions
  for s in sol:
      f = lambdify((plotvar), sol)
      nullclines.append([f(input) for input in inputrange])
  return nullclines

def plot_dynamical_system(xlimit, ylimit,x_range, y_range, x1_range, y1_range, f1_val, f2_val, x_nullclines, y_nullclines, stable_points, x_behavior, y_behavior, x_var_name, y_var_name):
  # Create plotting area and set axis limits
  fig, ax = plt.subplots(figsize=(10,10))
  ax.set_xlim([x_range[0], x_range[-1]])
  ax.set_ylim([y_range[0], y_range[-1]])  
  # Set axis lables
  ax.set_xlabel(f"{x_var_name} [V]", fontsize = 12)
  ax.set_ylabel(f"{y_var_name} [V/s]", fontsize = 12)
  ax.set_title("Stability analysis of Neuron Network", fontsize = 14)
  # Plot quivers
  ax.quiver(x1_range, y1_range, f1_val, f2_val, alpha=.5, headwidth = 5, headlength = 3, headaxislength = 3)
  # Plot nullclines
  [ax.plot(x_range, nullcline, c = 'b', alpha=.7, label= x_var_name + ' Nullcline') for nullcline in x_nullclines]
  [ax.plot(x_range, nullcline, c = 'r', alpha=.7, label= y_var_name + ' Nullcline') for nullcline in y_nullclines]
  # Plot stable points
  for stable_point in stable_points:
    # Extract information
    stable_point_x = stable_point[0][0]
    stable_point_y = stable_point[0][1]
    stable_point_type = stable_point[1]
    # Plot stable points with different color for each type 
    try: # We have to use try to ignore errors that occur when trying to plot complex value
      if stable_point_type == 'Unstable Node':
        ax.scatter(stable_point_x, stable_point_y, marker = '.', label = stable_point_type, s = 150, c = 'coral')
      elif stable_point_type == 'Stable Node':
        ax.scatter(stable_point_x, stable_point_y, marker = '.', label = stable_point_type, s = 150, c = 'cyan')
      elif stable_point_type == 'Saddle Point':
        ax.scatter(stable_point_x, stable_point_y, marker = '.', label = stable_point_type, s = 150, c = 'violet')
      elif stable_point_type == 'Unstable Focus':
        ax.scatter(stable_point_x, stable_point_y, marker = '.', label = stable_point_type, s = 150, c = 'orange')
      elif stable_point_type == 'Stable Focus':
        ax.scatter(stable_point_x, stable_point_y, marker = '.', label = stable_point_type, s = 150, c = 'springgreen')  
    except:
      # Do nothing for complex values
      pass
  # Plot the actual behavior of the variables over time


  desired_range = []
  for i in range(10000, len(y_behavior)):
      if  -130*volt/second < y_behavior[i] < -80*volt/second:
        desired_range.append(i)

  u_range = []
  vm_range = []

  for i in desired_range:
    u_range.append(y_behavior[i])
    vm_range.append(x_behavior[i])


  ax.plot(vm_range, u_range, c = 'mediumseagreen', label = 'Model trajectory')
  # Add legend
  ax.set_xlim(xlimit[0],xlimit[1])
  ax.set_ylim(ylimit[0],ylimit[1])
  ax.legend(fontsize = 14)
  # Show plot
  return ax

# %%
# Brain 2 implementation of Izhikevich neuron
def create_izhikevich_neuron(v_max):
    IZeq = '''   
        dv/dt = I + (0.04/ms/mV)*v**2 + (5/ms)*v + 140*mV/ms - u : volt

        du/dt = a*(b*v-u) : volt/second
        
        I : volt/second 
        '''
    reset = ''' 
        v = c
        u += d
        '''
    threshold = 'v >= {}*mV'.format(v_max)

    neuron_stability_analysis = NeuronGroup(1, IZeq, threshold = threshold, reset = reset, method = 'euler')

    # Return NeuronGroup object
    return neuron_stability_analysis

# Create function that creates a neuron and plots its behavior based on the given parameters
def plot_izhikevich_dynamical_sytem(boolean, xlimit, ylimit,I_ext, a_input, b_input, c_input, d_input, v_max):

  # Regular simulation of Izhikevich model using brian2
  # Start the scope to register all activity
  defaultclock.dt = 0.01*ms
  start_scope()
  # Define the neuron
  neuron = create_izhikevich_neuron(v_max)      
  # Set neuron parameters
  a = a_input/ms
  b = b_input/ms
  c = c_input * mV 
  d = d_input * volt/second  
  # Start monitoring the neurons state
  statemon = StateMonitor(source = neuron, variables = ['v', 'u'], record = True)
  # Run neuron simulation for 100ms without input

  # Set input current to neuron
  neuron.I = I_ext * volt / second
  # Run 500ms with input
  run(500*ms)
  # Remove input current to neuron

  if boolean == "True":

    fig, ax1 = plt.subplots(figsize=(16,9))
    fig, ax2 = plt.subplots(figsize=(16,9))

    ax1.plot(statemon.t/ms, statemon.v[0])
    ax1.set_title("Volatge behaviour of neuron network", fontsize = 16)
    ax1.set_xlabel('Time [ms]', fontsize = 14)
    ax1.set_xlim(0, 200)
    ax1.set_ylabel("V [mV]", fontsize = 14)

  
    ax2.plot(statemon.t/ms, statemon.u[0])
    ax2.set_title("Volatge behaviour of neuron network", fontsize = 16)
    ax2.set_xlabel('Time [ms]', fontsize = 14)
    ax2.set_ylim(-130, -90)
    ax2.set_xlim(0, 200)
    ax2.set_ylabel("V [mV]", fontsize = 14)

  # Define model for sympy and calculate nullclines  
  # Define the symbols
  v, u = symbols('v u')
  # First expression
  expression_for_v = I_ext + 0.04*v**2 + 5*v + 140 - u
  # Second expression
  expression_for_u = a_input*((b_input*v) - u)  

  # Solve dynamical system

  x_range, y_range,x1_range, y1_range, f1_val, f2_val, x_nullclines, y_nullclines, stable_points = solve_dynamical_system(expression_for_v, expression_for_u
  , x_range = [-180, 80] , y_range = [-140, -10], x_var = v, y_var = u)
  # Plot results
  plot_dynamical_system(xlimit, ylimit, x_range, y_range, x1_range, y1_range, f1_val, f2_val, x_nullclines, y_nullclines, stable_points, x_behavior = statemon.v[0]/mV, y_behavior = statemon.u[0], x_var_name = 'v', y_var_name = 'u')

# %%
#-----------------------------FIRST_VALUE_ANALYSIS------------------------------
#a zoomed out view 
I_ext_def =-99
a_def = 0.2
b_def = 2
c_def = -56.
d_def = -10
vmax_def = 30.

x_limit = [-100, 80]
y_limit = [-130, -10]

plotone = plot_izhikevich_dynamical_sytem("True",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)

#%%
#a zoomed in view 
#-----------------------------SECOND_VALUE_ANALYSIS------------------------------
I_ext_def =-99
a_def = 0.2
b_def = 2
c_def = -56.
d_def = -12
vmax_def = 30.

x_limit = [-100, 80]
y_limit = [-130, -80]

plotone = plot_izhikevich_dynamical_sytem("True",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)

# %%
#a zoomed  view 
#-----------------------------THIRD_VALUE_ANALYSIS------------------------------
I_ext_def =-99
a_def = 0.2
b_def = 2
c_def = -56.
d_def = -13
vmax_def = 30.

x_limit = [-100, 80]
y_limit = [-130, -80]

plotone = plot_izhikevich_dynamical_sytem("True",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)

#%%
#a zoomed out view 
#-----------------------------FINAL_VALUE_ANALYSIS------------------------------
I_ext_def =-99
a_def = 0.2
b_def = 2
c_def = -56.
d_def = -16
vmax_def = 30.

x_limit = [-100, 80]
y_limit = [-130, -80]

plotone = plot_izhikevich_dynamical_sytem("True",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)

#%%
