import Magnus_Expansion as ME
import ME_params
import dynamics_functions as funcs
import matplotlib.pyplot as plt
import sympy as sp
import numpy as np
import timeit
import multiprocessing as mp
import pickle
plt.rcParams["mathtext.fontset"]="cm"


def get_steadystate_data_3d(Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],params="A",initial_condition=8):
    '''Uses 8 cores to run the 825 calculations in batches of 8'''
    total_time_1=timeit.time.time()
    params_name="TPR_{}".format(params)
    steady_state_data = np.zeros((33,25,9,9,1),dtype=np.complex128)
    convergence_time_data = np.zeros((33,25),dtype=np.complex128)
    
    full_params_name = params_name+"_{}_{}".format(32,24)
    times,data = funcs.run_dynamics(full_params_name,Hamiltonian=Hamiltonian,methods=methods,t_final=8*np.pi,num_point=2**15)
    n = funcs.find_convergence_time(times,data[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=initial_condition)
    new_times,new_data = funcs.make_periodic_data(n,times,data[0][0][-1])
    convergence_time_data[32][24] = new_times[0]
    steady_state_data[32][24] = funcs.get_averages(new_times[:-1],data[0][0][:-1]@new_data,ME_params.params[full_params_name])
    
    for counter in range(103):
        time_1=timeit.time.time()
        indices = [counter*8+k for k in range(8)]
        j_indices = [index//33 for index in indices]
        i_indices = [index % 33 for index in indices]
        params_names = [params_name+"_{}_{}".format(i_indices[k],j_indices[k]) for k in range(8)]
        #full_params_names_list = [params_name+"_{}_{}".format(i,k) for k in range(j,j+8)] 
            #full_params_name = params_name+"_{}_{}".format(i,j)
        function_args = [params_name,Hamiltonian,methods,False,0,8*np.pi,2**15]
        function_args_list = [tuple(adjust_list(function_args,params_names[k])) for k in range(8)]
        with mp.Pool(processes=8) as pool:
            a = pool.starmap(funcs.run_dynamics, function_args_list)
      
        for k in range(8):
            times=a[k][0]
            data=a[k][1]
            i = i_indices[k]
            j = j_indices[k]
            n = funcs.find_convergence_time(times,data[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=initial_condition)
            new_times,new_data = funcs.make_periodic_data(n,times,data[0][0][-1])
            convergence_time_data[i][j]=new_times[0]
            steady_state_data[i][j] = funcs.get_averages(new_times[:-1],data[0][0][:-1]@new_data,ME_params.params[params_names[k]])
        time_2=timeit.time.time()

        print("Counter is {}/102 ".format(counter) + "Time took :",time_2-time_1)


    total_time_2=timeit.time.time()
    
    print("Calculations complete, time is :", total_time_2-total_time_1)
    
    return convergence_time_data,steady_state_data




def adjust_list(function_args,params_name):
    ''' params_names,function_args is a list'''
    new_function_args = function_args.copy()
    new_function_args[0]=params_name
    return new_function_args



if __name__ == "__main__":
    times,data =  get_steadystate_data_3d(Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],params="A",initial_condition=8)    
    with open("Data/steady_state_data_A","wb")as file:
        pickle.dump((times,data),file)    
    times,data =  get_steadystate_data_3d(Hamiltonian="lambda_RWA",methods=['M_14+M_23+M_32+M_41'],params="A",initial_condition=8)    
    with open("Data/steady_state_data_A_RWA","wb")as file:
        pickle.dump((times,data),file)    
    times,data =  get_steadystate_data_3d(Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],params="B",initial_condition=8)    
    with open("Data/steady_state_data_B","wb")as file:
        pickle.dump((times,data),file)    
    times,data =  get_steadystate_data_3d(Hamiltonian="lambda_RWA",methods=['M_14+M_23+M_32+M_41'],params="B",initial_condition=8)    
    with open("Data/steady_state_data_B_RWA","wb")as file:
        pickle.dump((times,data),file)  
    # times,data =  get_steadystate_data_3d(Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],params="C",initial_condition=8)    
    # with open("Data/steady_state_data_C","wb")as file:
    #     pickle.dump((times,data),file)
    # times,data =  get_steadystate_data_3d(Hamiltonian="lambda_RWA",methods=['M_14+M_23+M_32+M_41'],params="C",initial_condition=8)    
    # with open("Data/steady_state_data_C_RWA","wb")as file:
    #     pickle.dump((times,data),file)
    #times,data =  get_steadystate_data_3d(Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],params="D",initial_condition=8)    
    #with open("Data/steady_state_data_D","wb")as file:
    #    pickle.dump((times,data),file)    
    #times,data =  get_steadystate_data_3d(Hamiltonian="lambda_RWA",methods=['M_14+M_23+M_32+M_41'],params="D",initial_condition=8)    
    #with open("Data/steady_state_data_D_RWA","wb")as file:
    #    pickle.dump((times,data),file)    
