# Modeling and Simulation

## Sonar System Simulation

This project simulates a distributed sonar system, which consists of various independent nodes interconnected to function collaboratively. Maintaining a synchronized clock signal across these nodes is vital for the system's effective operation. One node assumes the role of a primary clock signal distributor, while the others function as receivers.

### System Overview

- **Node States**: Nodes can transition between primary and secondary roles based on system requirements and node conditions. The system can adaptively select a new primary node if the existing one encounters issues.
- **Redundancy**: Even if only `k` nodes work synchronously, the system can function correctly. Having more nodes than `k` provides redundancy, enhancing the system's survival rate and longevity.
- **Switch Failures**: Both switches A and B in nodes can exhibit different types of failures. These failures influence node states and the overall system state.

### System Structure

1. **Constants & Parameters**: Both the Julia and MATLAB codes initiate the system with predefined settings. These constants, like `NUM_SYSTEM`, `NUM_NODE`, and `LIFE_LIMIT`, set the scope of the simulation. Switch-specific constants also dictate their failure probabilities.
2. **Simulation Core**: The core logic assesses the system's state by evaluating the state of switches and nodes. If required, the primary node can be reselected to maintain system functionality.
3. **Output**: Post-simulation, the results detailing system lifetimes are recorded. Visualizations, such as histograms, provide insights into system reliability and MTTF (Mean Time To Failure).

### Usage

**Julia**:

- Launch the Julia script `simulation.jl` in an editor or the Julia REPL.
- Run the simulation, and the results will be saved as CSV files.
- Inspect the MTTF histogram.

**MATLAB**:

- Execute the MATLAB script `simulation_varia_timestep.m`.
- On completing the simulation for a series of nodes, MTTF and reliability data will be presented in the MATLAB console.

### Notice

- For tailored simulations, modify parameters like the total systems, node counts, or life limits.
- As random number generation is integral to these simulations, different runs might give slightly varying results.
- Differences between the Julia and MATLAB codes primarily stem from language-specific nuances.
