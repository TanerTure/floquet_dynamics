import numpy as np
import scipy
import pickle
import ME_params
from functools import partial
import Magnus_Expansion as ME
import timeit

def scipy_propagate(integrator):
    times = []
    y_vals = []
    for _ in range(81):
        y_vals.append([])
       # interp_y_vals.append([])
    
    for i in range(10000):
        times.append(integrator.t)
        if(i!=0):
            interpolant = integrator.dense_output()
        for j in range(81):
            y_vals[j].append(integrator.y[j])
            #if(i!=0):
               # interp_y_vals[j].append(interpolant.__call__(np_times[i])[j])
        try:
            integrator.step()
        except Exception as error:
            print(error)
            break
    return times,y_vals

def get_integrator(param_name = "TPR_A",Liouvillian=ME.L_lambda):
    new_L_func = partial(Liouvillian, **ME_params.params[param_name]["np"])
    def rhs_func(t,y):
            return  (new_L_func(t)@y.reshape(9,9)).reshape(-1)
    y_initial = np.eye(9,dtype=np.complex128).reshape(-1)
    integrator = scipy.integrate.RK45(rhs_func,0, y_initial, 2*np.pi,first_step=(2*np.pi)/(2**13-1),
                                          max_step=(2*np.pi)/(2**13-1),rtol=1e9,atol=1e9)
    return integrator

def get_data(param_name = "TPR_A",Liouvillian=ME.L_lambda):
    time_1 = timeit.time.time()
    integrator = get_integrator(param_name, Liouvillian)    
    times,y_vals = scipy_propagate(integrator)
    time_2 = timeit.time.time()
    print("Time took for RK45 is ", time_2-time_1)
    return times,y_vals

def main():
    #params = ['A']
    params = ["A","B","C"]
    for param in params:
        param_name = "TPR_{}".format(param)
        Liouvillian = ME.L_lambda
        times,y_vals = get_data(param_name,Liouvillian)
        with open("Data/"+param+"_dynamics_RK5","wb") as file:
            pickle.dump((times,y_vals),file)
       
      #for closed system:      
        param_name = "TPR_{}_closed".format(param)
        Liouvillian = ME.L_lambda
        times,y_vals = get_data(param_name,Liouvillian)
        with open("Data/"+param+"_dynamics_closed_RK5","wb") as file:
            pickle.dump((times,y_vals),file)
      
        param_name = "TPR_{}".format(param)
        Liouvillian = ME.L_lambda_RWA
        times,y_vals = get_data(param_name,Liouvillian)
        with open("Data/"+param+"_dynamics_RWA_RK5","wb") as file:
            pickle.dump((times,y_vals),file)
            
        param_name = "TPR_{}_closed".format(param)
        Liouvillian = ME.L_lambda_RWA
        times,y_vals = get_data(param_name,Liouvillian)
        with open("Data/"+param+"_dynamics_RWA_closed_RK5","wb") as file:
            pickle.dump((times,y_vals),file)
        
            
if __name__ == "__main__":
    main()


        