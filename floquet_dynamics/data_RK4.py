from functools import partial
import ME_params
import Magnus_Expansion as ME
import numpy as np
import matplotlib.pyplot as plt
import timeit
import pickle
import dynamics_functions as funcs


def get_func(param_name = "TPR_A",Liouvillian=ME.L_lambda):
    parameters_func = partial(Liouvillian, **ME_params.params[param_name]["np"])
    def new_func (t,y):
        return (parameters_func(t)@(y.reshape((9,9)))).reshape(-1)
    return new_func

def get_data(times,param_name = "TPR_A",Liouvillian=ME.L_lambda):
    time_1 = timeit.time.time()
    func = get_func(param_name, Liouvillian)
    y_initial = np.eye(9,dtype=np.complex128).reshape(-1)
    
    runge_kutta_object = funcs.RungeKutta(times, func, y_initial)
    runge_kutta_object.propagate()
    time_2 = timeit.time.time()
    print("Time took for RK4 is ", time_2-time_1)
    return runge_kutta_object.y_vals.reshape((len(times),9,9))



def main():
    params = ['A','C','B']
    times = np.linspace(0,2*np.pi,2**13)
    for param in params:
        param_name = "TPR_{}".format(param)
        Liouvillian = ME.L_lambda
        y_vals = get_data(times,param_name,Liouvillian)
        with open("Data/"+param+"_dynamics_RK4","wb") as file:
            pickle.dump((times,y_vals),file)
            
        param_name = "TPR_{}_closed".format(param)
        Liouvillian = ME.L_lambda
        y_vals = get_data(times,param_name,Liouvillian)
        with open("Data/"+param+"_dynamics_closed_RK4","wb") as file:
            pickle.dump((times,y_vals),file)
          
        param_name = "TPR_{}".format(param)
        Liouvillian = ME.L_lambda_RWA
        y_vals = get_data(times,param_name,Liouvillian)
        with open("Data/"+param+"_dynamics_RWA_RK4","wb") as file:
            pickle.dump((times,y_vals),file)
            
        param_name = "TPR_{}_closed".format(param)
        Liouvillian = ME.L_lambda_RWA
        y_vals = get_data(times,param_name,Liouvillian)
        with open("Data/"+param+"_dynamics_RWA_closed_RK4","wb") as file:
            pickle.dump((times,y_vals),file)

if __name__ == "__main__":
    main()
    
    