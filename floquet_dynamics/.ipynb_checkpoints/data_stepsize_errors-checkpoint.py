import Magnus_Expansion as ME
import ME_params
import dynamics_functions as funcs
import pickle
import numpy as np

def main():
    num_point=2**10
    num_stepsizes=4
    t_final = 2*np.pi
    
    params=["A","B","C"]
    method_data = [['M_14+M_23+M_32+M_41', "sixth_order"],
                    ["M_12+M_21","fourth_order"]
                   ]
    for method_datum in method_data:
        method = method_datum[0]
        method_name = method_datum[1]
#     #method,method_name = 'M_14+M_23+M_32+M_41', "sixth_order"
#     method, method_name = "M_12+M_21","fourth_order"

        #def run_dynamics(params_name,Hamiltonian="lambda_RWA",methods = ['M_14+M_23+M_32+M_41'],mapping=False,
        #             t_initial=0,t_final=4*np.pi,num_point=2**14,save_func='save_all',num_stepsizes=1,return_stepsize=False):
        for param in params:
            stepsizes,data = funcs.run_dynamics_stepsizes("TPR_{}".format(param),t_final=t_final,Hamiltonian="lambda",methods=[method],num_point=num_point,num_stepsizes = num_stepsizes,return_stepsize=True)
            stepsizes,RWA_data = funcs.run_dynamics_stepsizes("TPR_{}".format(param),t_final=t_final,Hamiltonian="lambda_RWA",methods=[method],num_point=num_point,num_stepsizes = num_stepsizes,return_stepsize=True)
            with open("Data/{}_dynamics_open_stepsizes_{}".format(param,method_name),"wb") as file:
                pickle.dump((stepsizes,data),file)
            with open("Data/{}_dynamics_RWA_open_stepsizes_{}".format(param,method_name),"wb") as file:
                pickle.dump((stepsizes,RWA_data),file)

            stepsizes,data = funcs.run_dynamics_stepsizes("TPR_{}_closed".format(param),t_final=t_final,Hamiltonian="lambda",methods=[method],num_point=num_point,num_stepsizes = num_stepsizes,return_stepsize=True)
            stepsizes,RWA_data = funcs.run_dynamics_stepsizes("TPR_{}_closed".format(param),t_final=t_final,Hamiltonian="lambda_RWA",methods=[method],num_point=num_point,num_stepsizes = num_stepsizes,return_stepsize=True)
            with open("Data/{}_dynamics_closed_stepsizes_{}".format(param,method_name),"wb") as file:
                pickle.dump((stepsizes,data),file)
            with open("Data/{}_dynamics_RWA_closed_stepsizes_{}".format(param,method_name),"wb") as file:
                pickle.dump((stepsizes,RWA_data),file)
    # method_data = [['CFME_equal', 'CFME_equal'], ['CFME_equal_opt', 'CFME_equal_opt'],['CFME_gauss', 'CFME_gauss'], ['CFME_gauss_opt', 'CFME_gauss_opt']]
    # method_data = [['CN2','CN2'],
    #                ['AM4_alt','AM4_alt'],
    #                ['AM4', "AM4"],
    #                ["RK4","RK4"],]
    method_data = [["CFME_equal","CFME_equal"],
                   ['AM4_alt','AM4_alt'],
                   ['AM4', "AM4"],
                   ["RK4","RK4"],
                   ]
    
    for method_datum in method_data:
        method = method_datum[0]
        method_name = method_datum[1]
    #method, method_name = "RK4", "RK4"
        for param in params:
            stepsizes, data = funcs.run_dynamics_stepsizes_class("TPR_{}".format(param),Hamiltonian="lambda", methods = [method],
                                                t_initial = 0, t_final = t_final, num_point=num_point, num_stepsizes=num_stepsizes, return_stepsize=True)
            stepsizes, RWA_data = funcs.run_dynamics_stepsizes_class("TPR_{}".format(param),Hamiltonian="lambda_RWA", methods = [method],
                                                t_initial = 0, t_final = t_final, num_point=num_point, num_stepsizes=num_stepsizes, return_stepsize=True)
            with open("Data/{}_dynamics_open_stepsizes_{}".format(param,method_name),"wb") as file:
                pickle.dump((stepsizes,data),file)
            with open("Data/{}_dynamics_RWA_open_stepsizes_{}".format(param,method_name),"wb") as file:
                pickle.dump((stepsizes,RWA_data),file)

            stepsizes,data = funcs.run_dynamics_stepsizes_class("TPR_{}_closed".format(param),t_final=t_final,Hamiltonian="lambda",methods=[method],num_point=num_point,num_stepsizes = num_stepsizes,return_stepsize=True)
            stepsizes,RWA_data = funcs.run_dynamics_stepsizes_class("TPR_{}_closed".format(param),t_final=t_final,Hamiltonian="lambda_RWA",methods=[method],num_point=num_point,num_stepsizes = num_stepsizes,return_stepsize=True)
            with open("Data/{}_dynamics_closed_stepsizes_{}".format(param,method_name),"wb") as file:
                pickle.dump((stepsizes,data),file)
            with open("Data/{}_dynamics_RWA_closed_stepsizes_{}".format(param,method_name),"wb") as file:
                pickle.dump((stepsizes,RWA_data),file)
        
if __name__ == "__main__":
    main()