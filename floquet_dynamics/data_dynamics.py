import Magnus_Expansion as ME
import ME_params
import dynamics_functions as funcs
import pickle
import numpy as np

if __name__ == "__main__":
    params=["A","B","C"]
    for param in params:
        times,data = funcs.run_dynamics("TPR_{}".format(param),t_final=2*np.pi,Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],num_point=2**13)
        times,RWA_data = funcs.run_dynamics("TPR_{}".format(param),t_final=2*np.pi,Hamiltonian="lambda_RWA",methods=['M_14+M_23+M_32+M_41'],num_point=2**13)
        with open("Data/{}_dynamics".format(param),"wb") as file:
            pickle.dump((times,data),file)
        with open("Data/{}_dynamics_RWA".format(param),"wb") as file:
            pickle.dump((times,RWA_data),file)
    
        times,data = funcs.run_dynamics("TPR_{}_closed".format(param),t_final=2*np.pi,Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],num_point=2**13)
        times,RWA_data = funcs.run_dynamics("TPR_{}_closed".format(param),t_final=2*np.pi,Hamiltonian="lambda_RWA",methods=['M_14+M_23+M_32+M_41'],num_point=2**13)
        with open("Data/{}_dynamics_closed".format(param),"wb") as file:
            pickle.dump((times,data),file)
        with open("Data/{}_dynamics_RWA_closed".format(param),"wb") as file:
            pickle.dump((times,RWA_data),file)
    
    
    