#%%
from tkinter import Y
from brian2 import *
from scipy import stats
import matplotlib.pyplot as plt

plt.style.use("seaborn")

#%%
#bifurcation plot of the system 

#----------------------------creating the bifurcation diagram ---------------------
defaultclock.dt = 0.05*ms

#creating a synapse with 500 nuerons 
N = 500

# define izhikevich model equations, threshold and parameters 
model = '''dv/dt = (0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u + I : volt
du/dt = a*(b*v-u) : volt/second
d : volt/second
I : volt/second '''
threshold = "v >= 30*mV"
reset = "v = c; u += d"

#the constant variables according to research
I = 10 * volt/second
a = 0.02 / ms
b = 0.2 / ms
c = -55 * mV

# do the simulation
start_scope()
neuron = NeuronGroup(N, model=model, threshold=threshold,
                    reset=reset, method='euler')


neuron.d = linspace(0.50*volt/second, 1.5*volt/second, N)
neuron.I = I

init_time = 1*second
run(init_time, report='text')
#we discard the first spikes, if we run it without doing this it makes it messy


states = StateMonitor(neuron, "u", record=True, when='start')
spikes = SpikeMonitor(neuron)
run(1*second, report='text')

# Get the values of d and w(or u)for each spike
D = neuron.d[spikes.i]
u = states.u[spikes.i, int_((spikes.t-init_time)/defaultclock.dt)]

plt.figure(figsize=(16,9))
plt.scatter(D / volt/second, u / volt/second, marker=".", color="k", linewidths=0.1)
plt.xlabel('d (volt/second)', fontsize = 14)
plt.ylabel('u[volt/second]', fontsize = 14)
# plt.xlim(0.83, 0.94)
plt.title("Bifurcation diagram of d vs u on system dynamics", fontsize = 18)
plt.show()

#---------------------------END OF BIFURCATION ------------------------------------------
# %%
#testing Regular spiking behaviour 

#--------------------------START_OF_REGULAR_SPIKING-------------------------------
I_ext= 10 * volt/second
a = 0.02 / ms
b = 0.2 / ms
c = -65 * mV

start_scope()
neuronRS = NeuronGroup(1, model=model, threshold=threshold,
                    reset=reset, method='euler')



statemon = StateMonitor(source = neuronRS, variables = ['vm', 'w', 'I'], record = True)

run(100*ms)
neuronRS.I = I_ext
neuronRS.d = 8*volt/second
run(500*ms)

fig, ax1 = plt.subplots(figsize=(16,9))
fig, ax2 = plt.subplots(figsize=(16,9))

#just found out that command d selects multiple values of the cell
#eg command d and then plot will keepm selecting the plot labels in the cell
#gamechanger 

ax1.plot(statemon.t/ms, statemon.vm[0])
ax1.set_title("Regular spiking behaviour of neuron network", fontsize = 15)
ax1.set_xlabel('Time [ms]', fontsize = 14)
ax1.set_xlim(2, 600)
ax1.set_ylabel('V [mV]', fontsize = 14)


ax2.plot(statemon.t/ms, statemon.w[0])
ax2.set_title("Regular spiking behaviour of neuron network", fontsize = 15)
ax2.set_xlabel('Time [ms]', fontsize = 14)
ax2.set_xlim(2, 600)
ax2.set_ylabel('u [mV]', fontsize = 14)

#%%
#displaying single period oscilations
desired_range = []
for i in range(0, len(statemon.w[0])):
  if  -8*volt/second < statemon.w[0][i] < 0.5*volt/second:
    desired_range.append(i)

u_range = []
vm_range = []

for i in desired_range:
  u_range.append(statemon.w[0][i])
  vm_range.append(statemon.vm[0][i])

plt.figure(figsize=(16,9))
plt.plot( vm_range, u_range)
plt.title("Phase plane diagram of single periodic Voltage vs Induced current", fontsize =16)
plt.xlabel('V [mV]', fontsize = 14)
plt.ylabel('u [v/s]', fontsize = 14)

#-------------------------END_OF_REGULAR_SPIKING---------------------------------------------
# %%
#intrinsicly bursting 
#------------------------START_OF_INTRINSICALLY_BURSTING-----------------------
I = 10* volt/second
a = 0.02 / ms
b = 0.2 / ms
c = -55 * mV

start_scope()
neuronIS = NeuronGroup(1, model=model, threshold=threshold,
                    reset=reset, method='euler')



statemon = StateMonitor(source = neuronIS, variables = ['vm', 'w', 'I'], record = True)

run(100*ms)
neuronIS.I = I
neuronIS.d = 4*volt/second
run(500*ms)

fig, ax1 = plt.subplots(figsize=(16,9))
fig, ax2 = plt.subplots(figsize=(16,9))

ax1.plot(statemon.t/ms, statemon.vm[0])
ax1.set_title("IB spiking behaviour of neuron network", fontsize = 15)
ax1.set_xlabel('Time [ms]', fontsize = 14)
ax1.set_xlim(2, 600)
ax1.set_ylabel('V [mV]', fontsize = 14)


ax2.plot(statemon.t/ms, statemon.w[0])
ax2.set_title("IB spiking behaviour of neuron network", fontsize = 15)
ax2.set_xlabel('Time [ms]', fontsize = 14)
ax2.set_xlim(2, 600)
ax2.set_ylabel('u [mV]', fontsize = 14)


#%%
#phase plane diagram
#can see initial trajectory and then it settles into a single period
desired_range = []
for i in range(0, len(statemon.w[0])):
  if  -8*volt/second < statemon.w[0][i] < 0.5*volt/second:
    desired_range.append(i)

u_range = []
vm_range = []

for i in desired_range:
  u_range.append(statemon.w[0][i])
  vm_range.append(statemon.vm[0][i])

plt.figure(figsize=(16,9))
plt.plot( vm_range, u_range)
plt.title("Phase plane diagram of Intrinsicly bursting Voltage vs Induced current", fontsize =16)
plt.xlabel('V[mV]', fontsize = 14)
plt.ylabel('u[v/s]', fontsize = 14)
#------------------END_OF_INTRINSICLY_BURSTING -----------------------------
# %%
#Chatering bursting 
#----------------------START_OF_CHATTERING_BURSTING---------------------------
I = 10* volt/second
a = 0.02 / ms
b = 0.2 / ms
c = -50 * mV

start_scope()
neuronCH = NeuronGroup(1, model=model, threshold=threshold,
                    reset=reset, method='euler')



statemon = StateMonitor(source = neuronCH, variables = ['vm', 'w', 'I'], record = True)

run(100*ms)
neuronCH.I = I
neuronCH.d = 2*volt/second
run(500*ms)

fig, ax1 = plt.subplots(figsize=(16,9))
fig, ax2 = plt.subplots(figsize=(16,9))

ax1.plot(statemon.t/ms, statemon.vm[0])
ax1.set_title("CB spiking behaviour of neuron network", fontsize = 15)
ax1.set_xlabel('Time [ms]', fontsize = 14)
ax1.set_xlim(2, 600)
ax1.set_ylabel('V [mV]', fontsize = 14)

ax2.plot(statemon.t/ms, statemon.w[0])
ax2.set_title("CB spiking behaviour of neuron network", fontsize = 15)
ax2.set_xlabel('Time [ms]', fontsize = 14)
ax2.set_xlim(2, 600)
ax2.set_ylabel('u [mV]', fontsize = 14)

#%%
#Phase plane diagram
#can see multiple periods associated with the system
desired_range = []
for i in range(0, len(statemon.w[0])):
  if  -8*volt/second < statemon.w[0][i] < 0.5*volt/second:
    desired_range.append(i)

u_range = []
vm_range = []

for i in desired_range:
  u_range.append(statemon.w[0][i])
  vm_range.append(statemon.vm[0][i])

plt.figure(figsize=(16,9))
plt.plot( vm_range, u_range)
plt.title("Phase plane diagram of Chattering Voltage vs Induced current", fontsize =16)
plt.xlabel('V[mV]', fontsize = 14)
plt.ylabel('u [v/s]', fontsize = 14)

#----------------------END_CHATERING_BURSTING -----------------------------

# %%
#Tonic spiking 
# The most common type of excitatory neurons in the mammalian neocortex are 
# pyramidal cells that fire spikes with decreasing frequency. When presented with 
# a prolonged stimulus, the neurons fire a few spikes with a short inter-spike
# period and then the period increases. This is called spike frequency adapation 
# in biology. It means that the frequency is relatively high at the onset of the stimulation and then it adapts.

I = 10* volt/second
a = 0.02 / ms
b = 0.2 / ms
c = -65* mV

start_scope()
neuronCH = NeuronGroup(1, model=model, threshold=threshold,
                    reset=reset, method='euler')



statemon = StateMonitor(source = neuronCH, variables = ['vm', 'w', 'I'], record = True)

run(10*ms)
neuronCH.I = I
neuronCH.d = 6*volt/second
run(100*ms)

plt.figure(figsize=(16,9))
plt.plot(statemon.t/ms, statemon.vm[0])
plt.title("tonic bursting behaviour of neuron network", fontsize = 16)
plt.xlabel('Time [ms]', fontsize = 14)
plt.ylabel("V[mV]", fontsize = 14)


#%%
#tonic Spiking bursting 

# Some neurons, such as the chattering neurons in cat neocortex
# [7], fire periodic bursts of spikes when stimulated, as in Fig. 1(c).
# The interburst (i.e., between bursts) frequency may be as high
# as 50 Hz, and it is believed that such neurons contribute to the
# gamma-frequency oscillations in the brain.

I = 10* volt/second
I2 = 8* volt/second
a = 0.02 / ms
b = 0.25 / ms
c = -65* mV

defaultclock.dt = 0.05*ms

start_scope()
neuronCH = NeuronGroup(500, model=model, threshold=threshold,
                    reset=reset, method='euler')



statemon = StateMonitor(source = neuronCH, variables = ['vm', 'w', 'I'], record = True)

run(50*ms)
neuronCH.I = I
neuronCH.d = 2*volt/second
run(400*ms)

plt.figure(figsize=(16,9))
plt.plot(statemon.t/ms, statemon.vm[0])
plt.title("tonic bursting behaviour of neuron network", fontsize = 16)
plt.xlabel('Time [ms]', fontsize = 14)
plt.xlim(0,150)
plt.ylabel("V[mV]", fontsize = 14)

#%%
I = 15* volt/second
a = 0.02 / ms
b = 0.25 / ms
c = -55* mV

start_scope()
neuronCH = NeuronGroup(1, model=model, threshold=threshold,
                    reset=reset, method='euler')



statemon = StateMonitor(source = neuronCH, variables = ['vm', 'w', 'I'], record = True)

run(50*ms)
neuronCH.I = I
neuronCH.d = 2*volt/second
run(300*ms)

plt.figure(figsize=(16,9))
plt.plot(statemon.t/ms, statemon.vm[0])
plt.title("tonic bursting behaviour of neuron network", fontsize = 16)
plt.xlabel('Time [ms]', fontsize = 14)
plt.ylabel("V[mV]", fontsize = 14)
#%%
#----------------------START_OF_STABILITY---------------------------------------
# Main packages
import numpy as np
from sympy import symbols, solve, nsolve, lambdify, sympify, dsolve, Eq, solveset, linear_eq_to_matrix, nonlinsolve, Matrix, diff, sqrt, exp
import sympy as smp
# Brian2 package
# Unit definitions
from brian2 import mV, ms, volt, second, umetre, ufarad, siemens, cm, msiemens, amp, uA, nA
# Other stuff
from brian2 import start_scope, NeuronGroup, StateMonitor, run

# Plotting stuff
import matplotlib.pyplot as plt
import seaborn as sns

# Interactive widgets
import ipywidgets as widgets
from ipywidgets import interact, interactive, fixed, HBox, VBox, Layout

#%%
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

def examine_stable_points(expr1, expr2, x_var, y_var, solutions):
  # First, calculate the jacobian matrix for our equations
  equation_matrix = Matrix([expr1, expr2])
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

def solve_dynamical_system(expr1, expr2, x_range, y_range, x_var, y_var): 
  # Convert our equations to functions
  f1 = lambdify((x_var, y_var), expr1)
  f2 = lambdify((x_var, y_var), expr2)
  # Define range that we will examine the system in
  start_x = x_range[0]
  end_x = x_range[1]
  start_y = y_range[0]
  end_y = y_range[1]
  # Set up list of x and y values inside this range 
  x_range = np.linspace(start_x, end_x,500)
  y_range = np.linspace(start_y, end_y,500)
  x1_range = np.linspace(start_x, end_x)
  y1_range = np.linspace(start_y, end_y)
  total_range = np.linspace(min(start_x, start_y), max(end_x, end_y),500)
  # Compute quivers by calculating the expression for a combination of x and y values
  f1_val = [[f1(x_cur, y_cur) for x_cur in x1_range] for y_cur in y1_range];
  f2_val = [[f2(x_cur, y_cur) for x_cur in x1_range] for y_cur in y1_range];
  # Solve analytically using sympy
  solutions = smp.solve([smp.Eq(expr1, 0), smp.Eq(expr2, 0)], (x_var, y_var)) 
  # Calculate nullclines
  x_nullclines = calculate_nullclines(expr1, y_var, x_var, total_range)
  y_nullclines = calculate_nullclines(expr2, y_var, x_var, total_range)

  # Get stable points and check their type
  stable_points = examine_stable_points(expr1, expr2, x_var, y_var, solutions)
  # Return relevant results
  return x_range, y_range, x1_range, y1_range, f1_val, f2_val, x_nullclines, y_nullclines, stable_points

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
  ax.quiver(x1_range, y1_range, f1_val, f2_val, alpha=.5, headwidth = 3, headlength = 2, headaxislength = 2)
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
  for i in range(7500, len(y_behavior)):
      if  -5*volt/second < y_behavior[i] < -3.5*volt/second:
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
  ax.legend()
  # Show plot
  return ax

# %%
# Brain 2 implementation of Izhikevich neuron
def create_izhikevich_neuron(v_max):
    """Creates a brian2 NeuronGroup that contains a single izhikevich neuron"""
    # Define differential equation for izhikevich neuron
    eqs = '''   
        dv/dt = I + (0.04/ms/mV)*v**2 + (5/ms)*v + 140*mV/ms - u : volt

        du/dt = a*(b*v-u) : volt/second
        
        I : volt/second 
        '''
    # Define reset function
    reset = ''' 
        v = c
        u += d
        '''
    # Define threshold
    threshold = 'v > {}*mV'.format(v_max)

    neuron = NeuronGroup(1, eqs, threshold = threshold, reset = reset, method = 'euler')

    # Return NeuronGroup object
    return neuron

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

  final = 0

  if boolean == "True":
    plt.figure(figsize=(16,9))
    plt.plot(statemon.t/ms, statemon.v[0])
    plt.title("Volatge behaviour of neuron network", fontsize = 16)
    plt.xlabel('Time [ms]', fontsize = 14)
    plt.xlim(2,180)
    plt.ylabel("V [mV]", fontsize = 14)

  # Define model for sympy and calculate nullclines  
  # Define the symbols
  v, u = symbols('v u')
  # First expression
  expr1 = I_ext + 0.04*v**2 + 5*v + 140 - u
  # Second expression
  expr2 = a_input*((b_input*v) - u)  

  # Solve dynamical system
  x_range, y_range,x1_range, y1_range, f1_val, f2_val, x_nullclines, y_nullclines, stable_points = solve_dynamical_system(expr1, expr2, x_range = [-160, 60] , y_range = [-65, 65], x_var = v, y_var = u)
  xlim = xlimit
  ylim = ylimit
  # Plot results
  plot_dynamical_system(xlim, ylim, x_range, y_range,x1_range, y1_range, f1_val, f2_val, x_nullclines, y_nullclines, stable_points, x_behavior = statemon.v[0]/mV, y_behavior = statemon.u[0], x_var_name = 'v', y_var_name = 'u')

#%%
#---------------------------THE_ROUTE_TO_CHAOS------------------------------
#---------------------------SINGLE_PERIOD_STABILITY-------------------------
#a zoomed out view 
I_ext_def =10
a_def = 0.02
b_def = 0.2
c_def = -55.
d_def = 0.80
vmax_def = 30.

x_limit = [-160, 60]
y_limit = [-20, 10]

plotone = plot_izhikevich_dynamical_sytem("True",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)
#%%
#a zoomed in more detailed view
x_limit = [-80, 35]
y_limit = [-7, -2]

plotone = plot_izhikevich_dynamical_sytem("false",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)

#%%
#--------------------DOUBLE_PERIOD_STABILTY-----------------------------------
I_ext_def =10
a_def = 0.02
b_def = 0.2
c_def = -55.
d_def = 0.85
vmax_def = 30.

x_limit = [-160, 60]
y_limit = [-20, 10]

plotone = plot_izhikevich_dynamical_sytem("True",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)
#%%
#a zoomed in more detailed view
x_limit = [-80, 35]
y_limit = [-5.2, -3.5]

plotone = plot_izhikevich_dynamical_sytem("false",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)

#%%
#------------------------------QUAD_PERIOD_STABILITY-----------------------------------
I_ext_def =10
a_def = 0.02
b_def = 0.2
c_def = -55.
d_def = 0.89
vmax_def = 30.

x_limit = [-160, 60]
y_limit = [-20, 10]

plotone = plot_izhikevich_dynamical_sytem("True",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)
#%%
#a zoomed in more detailed view
x_limit = [-80, 30]
y_limit = [-5.5, -3.5]

plotone = plot_izhikevich_dynamical_sytem("false",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)

#%%
#------------------------------CHAOS_PERIOD_STABILITY-----------------------------------
I_ext_def =10
a_def = 0.02
b_def = 0.2
c_def = -55.
d_def = 0.93
vmax_def = 30.

x_limit = [-160, 60]
y_limit = [-20, 10]

plotone = plot_izhikevich_dynamical_sytem("True",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)
#%%
#a zoomed in more detailed view
x_limit = [-80, 30]
y_limit = [-5.3, -3.5]

plotone = plot_izhikevich_dynamical_sytem("false",x_limit,y_limit, I_ext_def, a_def, b_def, c_def, d_def, vmax_def)

#%%
#---------------APPLYING_SCIPY_INTEGRATION_TO_DOUBLE_CHECK_RESULTS---------------
# %%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from scipy.integrate import solve_ivp
plt.close('all')

a = 0.02
b = 0.2
c = -55
d = 0.8
i = 10

p = [a,b,c,d,i]

def fun(t, u):
    du = [0,0]
    if u[0] < 30: #Checking if the threshold has been reached
        du[0] = (0.04*u[0] + 5)*u[0] + 140 - u[1] + p[4]
        du[1] = p[0]*(p[1]*u[0]-u[1])
    else:
        u[0] = p[2] #reset to -65    
        u[1] = u[1] + p[3] 

    return du

y0 = [0,0]

tspan = (0,100)
sol = solve_ivp(fun, tspan, y0)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)         
plt.plot(sol.t,sol.y[0, :],'k',linewidth = 5)
plt.plot(sol.t,sol.y[1, :],'r',linewidth = 5)
myleg = plt.legend(['v','u'],loc='upper right',prop = {'size':28,'weight':'bold'}, bbox_to_anchor=(1,0.9))
# %%

x = np.linspace(-150,-30, 1000)
x2 = np.linspace(-105,100, 1000)
a = 0.04
b = 5
c = 150
y = ((a*x*x) + (b*x) + c)
y1 = 0.2*x

plt.plot(sol.y[0, :], sol.y[1, :])
plt.plot(x,y)
plt.plot(x2,y1)

# %%
