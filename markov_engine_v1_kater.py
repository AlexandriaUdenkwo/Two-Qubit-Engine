
# Define initial conditions
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 13:36:52 2024

@author: katermurch
"""

import numpy as np
import matplotlib.pyplot as plt

# Define initial conditions
state = np.array([0, 1,0, 0])  # Initial state is 00 10 01 11
work = np.array([0])
probabilities = np.array([0., 1., 0., 0.])  # P_{00}, P_{01}, P_{10}, P_{11}
variance =2
readout  = np.array([np.random.normal(1, variance),np.random.normal(-1, variance)])
fidelity_matrix = np.array([[0.9, 0.05, 0.03, 0.02],  # Example fidelity matrix
                            [0.05, 0.9, 0.02, 0.03],
                            [0.03, 0.02, 0.9, 0.05],
                            [0.02, 0.03, 0.05, 0.9]])
P_theta = .2 # Probability of bit exchange
num_cycles = 10
P_T1_QA = 0.1
P_T1_QB = 0.1 #Probaility that the qubit decays during readout.

# To store the history
state_history = [state.copy()]
prob_history = [probabilities.copy()]
work_history = [work]
readout_history = [readout]

def apply_T1_decay(whichqubit,state):
    if whichqubit == "A":
        if np.array_equal(state, np.array([0, 1, 0, 0])):
             return np.array([1, 0, 0, 0])
        elif np.array_equal(state, np.array([0, 0, 0, 1])):
             return np.array([0, 0 , 1, 0])
        else: return state
    if whichqubit == "B":
        if np.array_equal(state, np.array([0, 0, 1, 0])):
             return np.array([1, 0, 0, 0])
        elif np.array_equal(state, np.array([0, 0, 0, 1])):
             return np.array([0, 1 , 0, 0])
        else: return state
        


def update_probabilities(fidelity, state):
    return np.dot(fidelity, state)

def apply_not_gate(state):  # pi pulses on QA and QB
    if np.array_equal(state, np.array([0, 1, 0, 0])):
        return np.array([0, 0, 1, 0])
    elif np.array_equal(state, np.array([0, 0, 1, 0])):
        return np.array([0, 1, 0, 0])
    elif np.array_equal(state, np.array([0, 0, 0, 1])):
        return np.array([1, 0, 0, 0])
    elif np.array_equal(state, np.array([1, 0, 0, 0])):
        return np.array([0, 0, 0, 1])
    else:
        print("error")
        return state  # Ensure we return the current state if none of the conditions are met

    
def get_readouts(state):
    readout_temp = np.zeros(2)
    if state[1] == 1 or state[3] == 1:
        readout_temp[0] = np.random.normal(1, variance)
    else:
        readout_temp[0] = np.random.normal(-1, variance)
    
    if state[2] == 1 or state[3] == 1:
         readout_temp[1] = np.random.normal(1, variance)
    else:# state[2] == 0 or state[3] == 0:
        readout_temp[1] = np.random.normal(-1, variance)
    return readout_temp
    

for cycle in range(num_cycles):
    # Step 1: Bit exchange with probability P_theta
    if np.random.rand() < P_theta:
        state = apply_not_gate(state)
    
    # Step 2: apply readout
    
    readout = get_readouts(state)
    
    #probabilities = update_probabilities(fidelity_matrix, state)
    
    #Step 3: allow a chance of T1 decay
    if np.random.rand() < P_T1_QA:
        state = apply_T1_decay("A",state)
    if np.random.rand() < P_T1_QB:
        state = apply_T1_decay("B",state)
    
    # Step 3: Apply conditional NOT gate if P_{01} > 0.5
    if readout[1] > 0: #if we think QB is in excited state
        state = apply_not_gate(state) #apply the pi pulses
        #modify, because if we apply the pi pulse when state is 11, state[3]=1, then we also get work
        if state[2] == 0 or state[0]==0:  # This is where work is extracted state[2] = 01
            work[0] += 1
        if state[2] == 1 or state[3] == 1:  # If QB is in the ground state we pay work to excite it
            work[0] -= 1

        
            

    
    # Store the history
    state_history.append(state.copy())
    prob_history.append(probabilities.copy())
    work_history.append(work.copy())
    readout_history.append(readout)

# Convert histories to numpy arrays for plotting
state_history = np.array(state_history)
prob_history = np.array(prob_history)
work_history = np.array(work_history)

# # Plot the results
# fig, axs = plt.subplots(3, 1, figsize=(10, 15))
plt.figure(figsize=(10, 4))
state_str_history = [''.join(map(str, s)) for s in state_history]
plt.plot(state_str_history, marker='o')
plt.title('State of the Bits at Each Step')
plt.xlabel('Cycle')
plt.ylabel('State')
plt.show()

state_mapping = {
    (1, 0, 0, 0): 0,
    (0, 1, 0, 0): 1,
    (0, 0, 1, 0): 2,
    (0, 0, 0, 1): 3
}

# Convert state_history to y-values using the mapping
y_values = [state_mapping[tuple(state)] for state in state_history]

# Plot the results
plt.figure(figsize=(10, 4))
plt.plot(y_values, marker='o')

# Set custom y-axis labels
plt.yticks([0, 1, 2, 3], ['[1, 0, 0, 0]', '[0, 1, 0, 0]', '[0, 0, 1, 0]', '[0, 0, 0, 1]'])
plt.title('State of the Bits at Each Step')
plt.xlabel('Cycle')
plt.ylabel('State')
plt.grid(True)
plt.show()

# # Plot the probabilities of the states at each step
# plt.figure(figsize=(10, 4))
# for i, states in enumerate(['00', '01', '10', '11']):
#     plt.plot(prob_history[:, i], label=f'P_{{{states}}}')
# plt.title('Probabilities of the States at Each Step')
# plt.xlabel('Cycle')
# plt.ylabel('Probability')
# plt.legend()
# plt.show()

# Plot the work history
plt.figure(figsize=(10, 4))
plt.plot(work_history, marker='o', color='r')
plt.title('Work at Each Step')
plt.xlabel('Cycle')
plt.ylabel('Work')
plt.show()

# Plot the readouts history
plt.figure(figsize=(10, 4))

plt.plot(readout_history)
plt.title('readout at Each Step')
plt.xlabel('Cycle')
plt.ylabel('readout')
plt.show()




