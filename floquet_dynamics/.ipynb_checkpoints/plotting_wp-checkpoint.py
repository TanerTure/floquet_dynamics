import dynamics_functions as dfuncs
import pickle
import plotting_functions as plot

if __name__ == "__main__":
    fontsize=18
    with open("Data/A_wp","rb") as file:
        times,steady_state_data = pickle.load(file)
    with open("Data/A_wp_RWA","rb") as file:
        times,steady_state_data_RWA = pickle.load(file)
    plot.steady_state_wp_plots(steady_state_data_RWA,steady_state_data,Hamiltonian="lambda",save_name="TPR_A",params="TPR_A",legend_fontsize=fontsize)
    
    with open("Data/B_wp","rb") as file:
        times,steady_state_data = pickle.load(file)
    with open("Data/B_wp_RWA","rb") as file:
        times,steady_state_data_RWA = pickle.load(file)
    plot.steady_state_wp_plots(steady_state_data_RWA,steady_state_data,Hamiltonian="lambda",save_name="TPR_B",params="TPR_B",legend_fontsize=fontsize)
    
    with open("Data/C_wp","rb") as file:
        times,steady_state_data = pickle.load(file)
    with open("Data/C_wp_RWA","rb") as file:
        times,steady_state_data_RWA = pickle.load(file)
    plot.steady_state_wp_plots(steady_state_data_RWA,steady_state_data,Hamiltonian="lambda",save_name="TPR_C",params="TPR_C",legend_fontsize=fontsize)
    
    
    