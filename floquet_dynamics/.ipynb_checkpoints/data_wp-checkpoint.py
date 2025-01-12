import dynamics_functions as dfuncs
import pickle

if __name__ == "__main__":
    #method = 'M_14+M_23+M_32+M_41'
    method = 'M_12+M_21'
    times, steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda",methods=[method],params="A")
    with open("Data/A_wp","wb") as file:
        pickle.dump((times,steady_state_data), file)
    
    times, RWA_steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda_RWA",methods=[method],params="A")
    with open("Data/A_wp_RWA","wb") as file:
        pickle.dump((times,RWA_steady_state_data), file)
    
    times, steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda",methods=[method],params="B")
    with open("Data/B_wp","wb") as file:
        pickle.dump((times,steady_state_data), file)
    
    times, RWA_steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda_RWA",methods=[method],params="B")
    with open("Data/B_wp_RWA","wb") as file:
        pickle.dump((times,RWA_steady_state_data), file)
        
    times, steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda",methods=[method],params="C")
    with open("Data/C_wp","wb") as file:
        pickle.dump((times,steady_state_data), file)
    
    times, RWA_steady_state_data = dfuncs.get_steadystate_data(Hamiltonian="lambda_RWA",methods=[method],params="C")
    with open("Data/C_wp_RWA","wb") as file:
        pickle.dump((times,RWA_steady_state_data), file)
    
    

    
    
    