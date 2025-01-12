import plotting_functions as plot
import dynamics_functions as funcs
import pickle
import scipy
import numpy as np


if __name__ == "__main__":
    steady_state = False
    #params=['A','B','C']
    params = ['A','B','C']
    #method_names = ['fourth_order', 'CFME_equal_opt','CFME_equal','CFME_gauss_opt','CFME_gauss','RK4','AM4','AM4_alt','CN2']
    method_names = ['fourth_order','CFME_equal','RK4','AM4_alt']
   # method_names = ['CFME_equal','CFME']
    file_strings = [["_open",""],["_open","_RWA"],["_closed",""],["_closed","_RWA"]]
    #file_strings_sixth_order = [["",""],["","_RWA",["_closed",""],["_closed","_RWA"]
    for file_string in file_strings:
        #file_strings related to RWA vs. non-RWA, open vs. closed
        
        for param in params:
            if(file_string[0] == "_open"):
                with open("Data/{}_dynamics{}{}_sixth_order".format(param,file_string[1],""),"rb") as file:
                    times, comparison = pickle.load(file)
            else:
                with open("Data/{}_dynamics{}{}_sixth_order".format(param,file_string[1],file_string[0]),"rb") as file:
                    times, comparison = pickle.load(file)
           # print(type(comparison),len(comparison),"for comparison")
           # print(type(comparison[0]),len(comparison[0]), "for comparison [0]")
           # print(type(comparison[0][0]),len(comparison[0][0]), "for comparison [0][0]")
            #comparison = [comparison[-1]]
           # print(type(comparison))
           # print(comparison[0][0][-1].shape,"shape")
            comparison = [comparison[0][0][-1]]
            # #begin steady state code
            if(steady_state == True):
                print("before find convergence time")
                if(file_string[0] == "_open"):
                    n = funcs.find_convergence_time(times,comparison[0],max_cycle=10000000,print_out=False,initial_condition=8)
                else:
                    n = 40
                print("after find convergence time" , n, "is n")
                comparison = [np.linalg.matrix_power(comparison[0],n+1)]
                
            # #end steady state code
            data = []
        
            for method_name in method_names:
                with open("Data/{}_dynamics{}{}_stepsizes_{}".format(param,file_string[1],file_string[0],method_name),"rb") as file:
                    stepsizes,data_method = pickle.load(file)
                    #print(len(data_method), " is len data_method")
                    #print(len(data_method[0]),"is len dat_method[0]")
                    #print(len(stepsizes), "is length stepsizes")
                    #print(stepsizes.shape)
                    #begin steady state code
                    #times = np.linspace(0,2*np.pi, 2**13)
                    #n = funcs.find_convergence_time(times,data_method[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=8)
                    if(steady_state == True):
                        for i in range(4):
                            data_method[0][i][-1] = np.linalg.matrix_power(data_method[0][i][-1],n+1) #only change the last matrix, because that is what is used for the error
                    #end steady state code

                    data.append(data_method[0])
                    
                    #code for steady state
                    
            print(len(data_method),"data_method")
            print(len(data_method[0]))
            print(len(data_method[0][0]))
            print(len(data_method[0][0][0]),"last")
            print("before plot func")
           
            #plot.make_error_plot(data, comparison, stepsizes, "Lambda", method_names,save_name=param+file_string[1]+file_string[0],legend_fontsize=11)
            file_name = param+file_string[1]+file_string[0]
            if(steady_state == True):
                file_name = file_name + "ss" 
            
            plot.make_error_plot(data, comparison, stepsizes, "Lambda", method_names,save_name=file_name,legend_fontsize=11)
            print("after plot func")
#         plot.make_3d_plot_times(times_3D,save_name="{}_time".format(param))
#         plot.make_3d_plot_matrix_element(data_3D,i=0,j=8,save_name="{}_pop_1".format(param))
        
#         with open("Data/steady_state_data_{}_RWA".format(param),"rb") as file:
#             times_3D_RWA, data_3D_RWA = pickle.load(file)
        
#         plot.make_3d_plot_times(times_3D_RWA,save_name="{}_time_RWA".format(param))
#         plot.make_3d_plot_matrix_element(data_3D_RWA,i=0,j=8,save_name="{}_pop_1_RWA".format(param))
        
       
        
        
#          with open("Data/{}_dynamics_stepsizes_{}".format(param,method_name),"wb") as file:
#             pickle.dump((stepsizes,data),file)
#         with open("Data/{}_dynamics_stepsizes_RWA_{}".format(param,method_name),"wb") as file:
#             pickle.dump((stepsizes,RWA_data),file)
    
#         stepsizes,data = funcs.run_dynamics_stepsizes("TPR_{}_closed".format(param),t_final=2*np.pi,Hamiltonian="lambda",methods=[method],num_point=2**13,num_stepsizes = num_stepsizes,return_stepsize=True)
#         stepsizes,RWA_data = funcs.run_dynamics_stepsizes("TPR_{}_closed".format(param),t_final=2*np.pi,Hamiltonian="lambda_RWA",methods=[method],num_point=2**13,num_stepsizes = num_stepsizes,return_stepsize=True)
#         with open("Data/{}_dynamics_closed_stepsizes_{}".format(param,method_name),"wb") as file:
#             pickle.dump((stepsizes,data),file)
#         with open("Data/{}_dynamics_RWA_closed_stepsizes_{}".format(param,method_name),"wb") as file:
#             pickle.dump((stepsizes,RWA_data),file)

#def make_error_plot(data, comparison,stepsizes,Hamiltonian, methods,
#                   indices=None,save_name=None,title=None):