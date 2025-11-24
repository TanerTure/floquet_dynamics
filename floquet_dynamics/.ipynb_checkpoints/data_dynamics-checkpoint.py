'''
    Runs the dynamics using either the sixth or fourth-order ME to the desired t_final. Produces data files 
    using pickle, which stores the list times with dimension [num_points], and the list of lists data with dimension 
    [1,1, (num_points,9,9)]. The parentheses represent the ndarray within. 

'''
import numpy as np
import pickle

import Magnus_Expansion as ME
import ME_params
import dynamics_functions as funcs

#num_point = 2**13
#t_final = 16 *np.pi
num_point = 2**12 #2**1
t_final = 2*np.pi

if __name__ == "__main__":
    params=["A","B","C"]
    method,method_name = 'M_14+M_23+M_32+M_41', "sixth_order"
    #method, method_name = "M_12+M_21","fourth_order"
    for param in params:
        times,data = funcs.run_dynamics("TPR_{}".format(param),t_final=t_final,Hamiltonian="lambda",methods=[method],num_point=num_point)
        times,RWA_data = funcs.run_dynamics("TPR_{}".format(param),t_final=t_final,Hamiltonian="lambda_RWA",methods=[method],num_point=num_point)
        with open("Data/{}_dynamics_{}".format(param,method_name),"wb") as file:
            pickle.dump((times,data),file)
        with open("Data/{}_dynamics_RWA_{}".format(param,method_name),"wb") as file:
            pickle.dump((times,RWA_data),file)
    
        times,data = funcs.run_dynamics("TPR_{}_closed".format(param),t_final=t_final,Hamiltonian="lambda",methods=[method],num_point=num_point)
        times,RWA_data = funcs.run_dynamics("TPR_{}_closed".format(param),t_final=t_final,Hamiltonian="lambda_RWA",methods=[method],num_point=num_point)
        with open("Data/{}_dynamics_closed_{}".format(param,method_name),"wb") as file:
            pickle.dump((times,data),file)
        with open("Data/{}_dynamics_RWA_closed_{}".format(param,method_name),"wb") as file:
            pickle.dump((times,RWA_data),file)
    
    
    