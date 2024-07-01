import dynamics_functions as dfuncs
import pickle

if __name__ == "__main__":        
    times, steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],params="A")
    with open("Data/A_wp","wb") as file:
        pickle.dump((times,steady_state_data), file)
    
    times, RWA_steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda_RWA",methods=['M_14+M_23+M_32+M_41'],params="A")
    with open("Data/A_wp_RWA","wb") as file:
        pickle.dump((times,RWA_steady_state_data), file)
    
    times, steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],params="B")
    with open("Data/B_wp","wb") as file:
        pickle.dump((times,steady_state_data), file)
    
    times, RWA_steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda_RWA",methods=['M_14+M_23+M_32+M_41'],params="B")
    with open("Data/B_wp_RWA","wb") as file:
        pickle.dump((times,RWA_steady_state_data), file)
        
    times, steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],params="C")
    with open("Data/C_wp","wb") as file:
        pickle.dump((times,steady_state_data), file)
    
    times, RWA_steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda_RWA",methods=['M_14+M_23+M_32+M_41'],params="C")
    with open("Data/C_wp_RWA","wb") as file:
        pickle.dump((times,RWA_steady_state_data), file)
    
    

    
    
    