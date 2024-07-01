# Floquet Dynamics
[This repository](https://github.com/TanerTure/Lambda_RWA_Code) contains code that uses [recently derived Magnus Expansion propagators](https://pubs.acs.org/doi/full/10.1021/acs.jpca.3c07866) for the time evolution of a lambda system. Because the time evolution is periodic, the dynamics can be sped up exponentially similar to the average Hamiltonian theory or Floquet-Magnus expansion. However, our method has the additional benefit of including the intraperiod motion. For theory and general description of the system, [see the corresponding submitted paper](). 
## File Structure

### Producing Data
The three .py files beginning with **data** (`data_dynamics.py`, `data_wp.py`, `data_parallel_3D.py`) are scripts to run the dynamics and store either steady state values or the dynamical data. Files can be run by using `python data_dynamics.py` or `python3 data_dynamics.py`.
The data is stored as binary files in the `Data` folder. These binary files are required to produce the plots, as described in the next section.

The most time consuming calculations are in `data_parallel_3D.py`; to reduce the computational time, the code is written to perform calculations in parallel.


### Making Plots
Plots are made using the corresponding .py file that begin with **plotting**: (plotting_dynamics.py, plotting_wp.py, data_parallel_3D.py). Any of these can be run by using
    `python plotting_dynamics.py` or `python3 data_dynamics.py`
    
### Helper Files
The code concerning the Magnus Expansion propagators can be found in `Magnus_Expansion.py`. The system parameters are defined in `ME_params.py`. Other commonly used functions for the dynamics (such as taking advantage of the Floquet time evolution, or transformation to the rotating frame) are defined in `dynamics_functions.py`. Helper functions for making the plots can be found in `plotting_functions.py`. 

## Code Convention
The density operator in the code has the convention $|1\rangle$ as the excited state, $|2\rangle$ as the state coupled by the control pulse $\omega_c$, and $|3\rangle$ as the state coupled by the probe pulse $\omega_p$. This is a less used convention, and different from the [main paper](). These differences are summarized in the following table. 

|Code|Paper|
|-----|-----|
|$\|1\rangle$|$\|2\rangle$|
|$\|2\rangle$|$\|3\rangle$|
|$\|3\rangle$|$\|1\rangle$|
