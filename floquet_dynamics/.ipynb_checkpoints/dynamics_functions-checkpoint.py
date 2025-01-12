import Magnus_Expansion as ME
import ME_params
import numpy as np
import timeit
from functools import partial
import scipy

class RungeKutta:
    #order:int = NotImplemented
    def __init__(self,times,func,y_initial):
        self.times = times
        self.func = func
       # self.y_vals = np.zeros(len(times),len(y_initial),dtype=type(y_initial[0]))
        self.y_vals = self.initiate_y_vals(y_initial)
        
    def initiate_y_vals(self,y_initial):
        y_vals = np.zeros((len(self.times),len(y_initial)),dtype=np.complex128)
        y_vals[0] = y_initial[:]
        return y_vals
    
    def propagate(self):
        for i in range(len(self.times)-1):
            dt = self.times[i+1]-self.times[i]
            k_1 = self.func(self.times[i],self.y_vals[i])
            k_2 = self.func(self.times[i]+dt/2,self.y_vals[i]+k_1*dt/2)
            k_3 = self.func(self.times[i]+dt/2,self.y_vals[i]+k_2*dt/2)
            k_4 = self.func(self.times[i]+dt,self.y_vals[i]+k_3*dt)
    
            self.y_vals[i+1] = self.y_vals[i]+dt/6*(k_1+2*k_2+2*k_3+k_4)

class AdamsMoulton:
    #order:int = NotImplemented
    def __init__(self,times,func,y_initial):
        self.times = times
        self.dt = times[1]-times[0] # for constant stepsize
        self.func = func
       # self.y_vals = np.zeros(len(times),len(y_initial),dtype=type(y_initial[0]))
        self.y_vals = self.initiate_y_vals(y_initial)
        
    def initiate_y_vals(self,y_initial):
        y_vals = np.zeros((len(self.times),len(y_initial)),dtype=np.complex128)
        y_vals[0] = y_initial[:]
        return y_vals
    
    def RK4_initialize(self):
        #gets the first four points
        for i in range(3):
            k_1 = self.func(self.times[i],self.y_vals[i])
            k_2 = self.func(self.times[i]+self.dt/2,self.y_vals[i]+k_1*self.dt/2)
            k_3 = self.func(self.times[i]+self.dt/2,self.y_vals[i]+k_2*self.dt/2)
            k_4 = self.func(self.times[i]+self.dt,self.y_vals[i]+k_3*self.dt)
            
            self.y_vals[i+1] = self.y_vals[i]+self.dt/6*(k_1+2*k_2+2*k_3+k_4)
        return
    
    def Adams_Bashforth(self,i):
        return self.y_vals[i+3] + self.dt/24*(55*self.func(self.times[i+3],self.y_vals[i+3])
                                         -59*self.func(self.times[i+2],self.y_vals[i+2])
                                         +37*self.func(self.times[i+1],self.y_vals[i+1])
                                         -9*self.func(self.times[i],self.y_vals[i])
                                        )
    def Adams_Moulton(self,i,prediction):
        self.y_vals[i+3] = self.y_vals[i+2] + self.dt/24*(9*self.func(self.times[i+3],prediction)
                                                          +19*self.func(self.times[i+2],self.y_vals[i+2])
                                                          -5*self.func(self.times[i+1],self.y_vals[i+1])
                                                          +self.func(self.times[i],self.y_vals[i])
                                                         )
        
    def propagate(self):
        self.RK4_initialize()
        for i in range(len(self.times)-4):
            prediction = self.Adams_Bashforth(i)
            self.Adams_Moulton(i+1,prediction)
        return
    
#AM4_alt is also different from the rest because you shouldn't use f(t,y);
#Therefore, func here stands for function that returns the Liouvillian matrix (my definition) 
class AdamsMoulton2:
    def __init__(self,times,func,y_initial):
        self.times = times
        self.dt = times[1]-times[0] # for constant stepsize
        self.func = func
       # self.y_vals = np.zeros(len(times),len(y_initial),dtype=type(y_initial[0]))
        self.y_vals = self.initiate_y_vals(y_initial)
        
    def initiate_y_vals(self,y_initial):
        y_vals = np.zeros((len(self.times),len(y_initial)),dtype=np.complex128)
        y_vals[0] = y_initial[:]
        return y_vals
    
    def RK4_initialize(self):
        n = self.get_matrix_length()
        #gets the first four points
        for i in range(2):
            k_1 = self.func(self.times[i])@self.y_vals[i].reshape(n,n)
            k_2 = self.func(self.times[i]+self.dt/2)@(self.y_vals[i].reshape(n,n)+k_1*self.dt/2)
            k_3 = self.func(self.times[i]+self.dt/2)@(self.y_vals[i].reshape(n,n)+k_2*self.dt/2)
            k_4 = self.func(self.times[i]+self.dt)@(self.y_vals[i].reshape(n,n)+k_3*self.dt)
            
            self.y_vals[i+1] = self.y_vals[i]+self.dt/6*(k_1+2*k_2+2*k_3+k_4).reshape(-1)
        return
    
    def get_matrix_length(self):
        return round(np.sqrt(len(self.y_vals[0])))
    # def Adams_Bashforth(self,i):
    #     return self.y_vals[i+3] + self.dt/24*(55*self.func(self.times[i+3],self.y_vals[i+3])
    #                                      -59*self.func(self.times[i+2],self.y_vals[i+2])
    #                                      +37*self.func(self.times[i+1],self.y_vals[i+1])
    #                                      -9*self.func(self.times[i],self.y_vals[i])
    #                                     )
    def adams_moulton(self,i):
        n = self.get_matrix_length()
        A = np.eye(n,dtype=np.complex128) - self.dt*9/24*self.func(self.times[i+3])
        B = self.y_vals[i+2].reshape((n,n))+self.dt/24*(19*self.func(self.times[i+2])@self.y_vals[i+2].reshape((n,n))
                                                        - 5 * self.func(self.times[i+1])@self.y_vals[i+1].reshape((n,n))
                                                        + 1 * self.func(self.times[i])@self.y_vals[i].reshape((n,n))
                                                       )
        return np.linalg.solve(A,B).reshape(-1)                                                 
        
        # self.y_vals[i+3] = self.y_vals[i+2] + self.dt/24*(9*self.func(self.times[i+3],prediction)
        #                                                   +19*self.func(self.times[i+2],self.y_vals[i+2])
        #                                                   -5*self.func(self.times[i+1],self.y_vals[i+1])
        #                                                   +self.func(self.times[i],self.y_vals[i])
        #                                                  )
        
    def propagate(self):
        self.RK4_initialize()
        for i in range(len(self.times)-3):
            self.y_vals[i+3] = self.adams_moulton(i)
        return

#CranckNicolson is different from the rest because you shouldn't use f(t,y);
#Therefore, func here stands for function that returns the Liouvillian matrix (my definition) 
class CrankNicolson:
    def __init__(self,times,func,y_initial):
        self.times = times
        self.dt = times[1]-times[0] # for constant stepsize
        self.func = func
       # self.y_vals = np.zeros(len(times),len(y_initial),dtype=type(y_initial[0]))
        self.y_vals = self.initiate_y_vals(y_initial)
        
    def initiate_y_vals(self,y_initial):
        y_vals = np.zeros((len(self.times),len(y_initial)),dtype=np.complex128)
        y_vals[0] = y_initial[:]
        return y_vals
    
    def get_matrix_length(self):
        return round(np.sqrt(len(self.y_vals[0])))
        
    def crank_nicolson(self,i):
        n = self.get_matrix_length()
        A = np.eye(n,dtype=np.complex128) - self.dt/2*self.func(self.times[i+1])
        B = (np.eye(n,dtype=np.complex128)+self.dt/2*self.func(self.times[i]))@self.y_vals[i].reshape(n,n)
        return np.linalg.solve(A,B).reshape(-1)
    
    def propagate(self):
        for i in range(len(self.times)-1):
            self.y_vals[i+1] = self.crank_nicolson(i)
        return
    
class CFME_gauss():
    def __init__(self,times,func,y_initial):
        self.order = 4
        self.legend_name = "Thalhammer2006"
        self.times = times
        self.dt = times[1]-times[0] #for constant stepsize
        self.func = func
        self.y_vals = self.initiate_y_vals(y_initial)
        self.num_points,self.points = self.initiate_points()
        self.num_exponentials,self.coefficients = self.initiate_coefficients()

        
    def initiate_y_vals(self,y_initial):
        y_vals = np.zeros((len(self.times),len(y_initial)),dtype=np.complex128)
        y_vals[0] = y_initial[:]
        return y_vals
    
    def initiate_points(self):
        num_points = 2 
        points = np.zeros(num_points,dtype=np.complex128)
        points[0] = .5-np.sqrt(3)/6
        points[1] = .5 + np.sqrt(3)/6
        return num_points, points
    
    def initiate_coefficients(self):
        num_exponentials = 2
        coefficients = np.zeros((num_exponentials,self.num_points),dtype=np.complex128)
        coefficients[0][0] = .25 + np.sqrt(3)/6
        coefficients[0][1] = .25 - np.sqrt(3)/6
        coefficients[1][0] = coefficients[0][1]
        coefficients[1][1] = coefficients[0][0]
        return num_exponentials,coefficients
        
    def propagate(self):
        length = round(np.sqrt(self.y_vals.shape[1]))
        for i in range(len(self.times)-1):
            matrices = np.zeros((self.num_points, length, length),dtype=np.complex128)
            y_val = self.y_vals[i] 
            for j in range(self.num_points):
                t = self.times[i]+self.dt*self.points[j]
                matrices[j] = self.func(t).reshape(length,length)
            for j in range(self.num_exponentials):
                exponent_sum = np.zeros((length,length),dtype=np.complex128) 
                for k in range(self.num_points):
                    exponent_sum += self.dt*self.coefficients[j][k]*matrices[k]
                y_val = (scipy.linalg.expm(exponent_sum)@y_val.reshape(length,length)).reshape(-1)
            self.y_vals[i+1] = y_val
        return

class CFME_gauss_opt():
    def __init__(self,times,func,y_initial):
        self.order = 4
        self.legend_name = "Thalhammer2006"
        self.times = times
        self.dt = times[1]-times[0] #for constant stepsize
        self.func = func
        self.y_vals = self.initiate_y_vals(y_initial)
        self.num_points,self.points = self.initiate_points()
        self.num_exponentials,self.coefficients = self.initiate_coefficients()

        
    def initiate_y_vals(self,y_initial):
        y_vals = np.zeros((len(self.times),len(y_initial)),dtype=np.complex128)
        y_vals[0] = y_initial[:]
        return y_vals
    
    def initiate_points(self):
        num_points = 3 
        points = np.zeros(num_points,dtype=np.complex128)
        points[0] = .5 - np.sqrt(15)/10
        points[1] = .5
        points[2] = .5 + np.sqrt(15)/10
        return num_points, points
    
    def initiate_coefficients(self):
        num_exponentials = 3
        coefficients = np.zeros((num_exponentials,self.num_points),dtype=np.complex128)
        coefficients[0][0] = 10*np.sqrt(15)/261+37/240
        coefficients[0][1] = -1/30
        coefficients[0][2] = -10*np.sqrt(15)/261+37/240
        coefficients[1][0] = -11/360
        coefficients[1][1] = 23/45
        coefficients[1][2] = -11/360
        coefficients[2][0] = -10 * np.sqrt(15)/261 + 37/240
        coefficients[2][1] = -1/30
        coefficients[2][2] = 10*np.sqrt(15)/261+37/240
        return num_exponentials,coefficients
        
    def propagate(self):
        length = round(np.sqrt(self.y_vals.shape[1]))
        for i in range(len(self.times)-1):
            matrices = np.zeros((self.num_points, length, length),dtype=np.complex128)
            y_val = self.y_vals[i] 
            for j in range(self.num_points):
                t = self.times[i]+self.dt*self.points[j]
                matrices[j] = self.func(t).reshape(length,length)
            for j in range(self.num_exponentials):
                exponent_sum = np.zeros((length,length),dtype=np.complex128) 
                for k in range(self.num_points):
                    exponent_sum += self.dt*self.coefficients[j][k]*matrices[k]
                y_val = (scipy.linalg.expm(exponent_sum)@y_val.reshape(length,length)).reshape(-1)
            self.y_vals[i+1] = y_val
        return 

class CFME_equal():
    def __init__(self,times,func,y_initial):
        self.order = 4
        self.legend_name = "CFME"
        self.times = times
        self.dt = times[1]-times[0] #for constant stepsize
        self.func = func
        self.y_vals = self.initiate_y_vals(y_initial)
        self.num_points,self.points = self.initiate_points()
        self.num_exponentials,self.coefficients = self.initiate_coefficients()

        
    def initiate_y_vals(self,y_initial):
        y_vals = np.zeros((len(self.times),len(y_initial)),dtype=np.complex128)
        y_vals[0] = y_initial[:]
        return y_vals
    
    def initiate_points(self):
        num_points = 3 
        points = np.zeros(num_points,dtype=np.complex128)
        points[0] = 0
        points[1] = .5
        points[2] = 1
        return num_points, points
    
    def initiate_coefficients(self):
        num_exponentials = 2
        coefficients = np.zeros((num_exponentials,self.num_points),dtype=np.complex128)
        coefficients[0][0] = 3/12
        coefficients[0][1] = 4/12
        coefficients[0][2] = -1/12
        coefficients[1][0] = -1/12
        coefficients[1][1] = 4/12
        coefficients [1][2] = 3/12
        # coefficients[0][0] = -1/12
        # coefficients[0][1] = 4/12
        # coefficients[0][2] = 3/12
        # coefficients[1][0] = 3/12
        # coefficients[1][1] = 4/12
        # coefficients [1][2] = -1/12
        return num_exponentials,coefficients
    
    def propagate(self):
        length = round(np.sqrt(self.y_vals.shape[1]))
        for i in range(len(self.times)-1):
            matrices = np.zeros((self.num_points, length, length),dtype=np.complex128)
            y_val = self.y_vals[i] 
            for j in range(self.num_points):
                t = self.times[i]+self.dt*self.points[j]
                matrices[j] = self.func(t).reshape(length,length)
            for j in range(self.num_exponentials):
                exponent_sum = np.zeros((length,length),dtype=np.complex128) 
                for k in range(self.num_points):
                    exponent_sum += self.dt*self.coefficients[j][k]*matrices[k]
                y_val = (scipy.linalg.expm(exponent_sum)@(y_val.reshape(length,length))).reshape(-1)
            self.y_vals[i+1] = y_val
        return

class CFME_equal_opt():
    def __init__(self,times,func,y_initial):
        self.order = 4
        self.legend_name = "CFME_equal_opt"
        self.times = times
        self.dt = times[1]-times[0] #for constant stepsize
        self.func = func
        self.y_vals = self.initiate_y_vals(y_initial)
        self.num_points,self.points = self.initiate_points()
        self.num_exponentials,self.coefficients = self.initiate_coefficients()

        
    def initiate_y_vals(self,y_initial):
        y_vals = np.zeros((len(self.times),len(y_initial)),dtype=np.complex128)
        y_vals[0] = y_initial[:]
        return y_vals
    
    def initiate_points(self):
        num_points = 3 
        points = np.zeros(num_points,dtype=np.complex128)
        points[0] = 0
        points[1] = .5
        points[2] = 1
        return num_points, points
    
    def initiate_coefficients(self):
        num_exponentials = 3
        coefficients = np.zeros((num_exponentials,self.num_points),dtype=np.complex128)
        coefficients[0][0] = 1931/6960
        coefficients[0][1] = -1/20
        coefficients[0][2] = 331/6960
        coefficients[1][0] = -19/120
        coefficients[1][1] = 23/30
        coefficients[1][2] = -19/120
        coefficients[2][0] = 331/6960
        coefficients[2][1] = -1/20
        coefficients[2][2] = 1931/6960
        # coefficients[0][0] = -1/12
        # coefficients[0][1] = 4/12
        # coefficients[0][2] = 3/12
        # coefficients[1][0] = 3/12
        # coefficients[1][1] = 4/12
        # coefficients [1][2] = -1/12
        #print(coefficients)
        return num_exponentials,coefficients
        
    def propagate(self):
        length = round(np.sqrt(self.y_vals.shape[1]))
        for i in range(len(self.times)-1):
            matrices = np.zeros((self.num_points, length, length),dtype=np.complex128)
            y_val = self.y_vals[i] 
            for j in range(self.num_points):
                t = self.times[i]+self.dt*self.points[j]
                matrices[j] = self.func(t).reshape(length,length)
            for j in range(self.num_exponentials):
                exponent_sum = np.zeros((length,length),dtype=np.complex128) 
                for k in range(self.num_points):
                    exponent_sum += self.dt*self.coefficients[j][k]*matrices[k]
                y_val = (scipy.linalg.expm(exponent_sum)@(y_val.reshape(length,length))).reshape(-1)
            self.y_vals[i+1] = y_val
        return
                                    

    
def get_func(param_name = "TPR_A",Liouvillian=ME.L_lambda):
    parameters_func = partial(Liouvillian, **ME_params.params[param_name]["np"])
    def new_func (t,y):
        return (parameters_func(t)@(y.reshape((9,9)))).reshape(-1)
    return new_func

def get_func_implicit(param_name = "TPR_A",Liouvillian=ME.L_lambda):
    parameters_func = partial(Liouvillian, **ME_params.params[param_name]["np"])
    return parameters_func

#for implicit functions, it is better to have the function be the Liouvillian

# def get_data(times,param_name = "TPR_A",Liouvillian=ME.L_lambda, method_func=RK4):
#     time_1 = timeit.time.time()
#     func = get_func(param_name, Liouvillian)
#     y_initial = np.eye(9,dtype=np.complex128).reshape(-1)
    
#     runge_kutta_object = RungeKutta(times, func, y_initial)
#     runge_kutta_object.propagate()
#     time_2 = timeit.time.time()
#     print("Time took for RK4 is ", time_2-time_1)
#     return runge_kutta_object.y_vals.reshape((len(times),9,9))

def get_data(times,param_name = "TPR_A",Liouvillian=ME.L_lambda, method_class=RungeKutta):
    time_1 = timeit.time.time()
    if( method_class == RungeKutta or method_class == AdamsMoulton):
        func = get_func(param_name, Liouvillian)
    else:
        func = get_func_implicit(param_name,Liouvillian)
    y_initial = np.eye(9,dtype=np.complex128).reshape(-1)
    
    #runge_kutta_object = RuneKutta(times, func, y_initial)
    ODE_solver_object = method_class(times,func,y_initial)
    ODE_solver_object.propagate()
    time_2 = timeit.time.time()
    print("Time took for {} is ".format(method_class), time_2-time_1)
    return ODE_solver_object.y_vals.reshape((len(times),9,9))
        
def get_stepsizes(t_initial,t_final,num_point,num_stepsizes):
    stepsizes = np.zeros(num_stepsizes)
    num_points = np.zeros(num_stepsizes)
    num_points[0] = num_point
    t = np.linspace(t_initial,t_final,int(num_points[0]))
    stepsizes[0] = t[1]-t[0]
    for i in range(num_stepsizes-1):
        num_points[i+1] = num_points[i]/2
        t = np.linspace(0,t_final,int(num_points[i+1]))
        stepsizes[i+1]=t[1]-t[0]
    return stepsizes, num_points

def run_dynamics_stepsizes(params_name,Hamiltonian="lambda_RWA",methods = ['M_14+M_23+M_32+M_41'],mapping=False,
                 t_initial=0,t_final=4*np.pi,num_point=2**14,save_func='save_all',num_stepsizes=1,return_stepsize=False):
    ''' methods: should be a list containing strings'''
    #t_initial = 0
    #num_stepsizes = 1
    stepsizes = np.zeros(num_stepsizes)
    num_points = np.zeros(num_stepsizes)
    #num_points[0] = 2**14-31
    num_points[0] = num_point
    #um_points[0]=2
    t = np.linspace(0,t_final,int(num_points[0]))
    stepsizes[0] = t[1]-t[0]
    for i in range(num_stepsizes-1):
        num_points[i+1] = num_points[i]/2
        print(num_points)
        t = np.linspace(0,t_final,int(num_points[i+1]))
        stepsizes[i+1]=t[1]-t[0]
    #print(stepsizes)

    data = [ [] for method in methods]
    for i in range(len(stepsizes)):
        for j in range(len(methods)):
            times = np.linspace(t_initial,t_final,int(num_points[i]))
            time_1=timeit.time.time()

            data[j].append(ME.propagate(Hamiltonian,methods[j],save_func,
                                        ME_params.params[params_name],times,mapping=mapping)
                            )
            #save_all
            time_2=timeit.time.time()                    
            print("Time took for {}: ".format(methods[j]),time_2-time_1)
    if(return_stepsize == False):
        return t,data
    else:
        return stepsizes, data
class2methods = {"RK4": RungeKutta,
                 "AM4": AdamsMoulton,
                 "AM4_alt":AdamsMoulton2,
                 "CN2":CrankNicolson,
                 "CFME_gauss":CFME_gauss,
                 "CFME_equal":CFME_equal,
                 "CFME_gauss_opt":CFME_gauss_opt,
                 "CFME_equal_opt":CFME_equal_opt,
                }
#Hamiltonian2func = {"lambda_RWA":
def run_dynamics_stepsizes_class(params_name,Hamiltonian="lambda_RWA", methods = ["RK4"],
                                            t_initial = 0, t_final = 4*np.pi, num_point=2**14, num_stepsizes=1, return_stepsize=False):
    stepsizes,num_points = get_stepsizes(t_initial, t_final, num_point, num_stepsizes)
    data = [ [] for method in methods]
    for i in range(len(stepsizes)):
        for j in range(len(methods)):
            times = np.linspace(t_initial,t_final,int(num_points[i]))
            #propagator_object = class2methods[methods[j]](times,func,y_initial)
            Liouvillian = ME.H_name_to_func[Hamiltonian][0]
            #print(type(Liouvillian))
           # print(Liouvillian)
            method_class = class2methods[methods[j]]
            data[j].append(get_data(times,params_name, Liouvillian,method_class))
    
    if(return_stepsize == False):
        return t,data
    else:
        return stepsizes, data

def run_dynamics(params_name,Hamiltonian="lambda_RWA",methods = ['M_14+M_23+M_32+M_41'],mapping=False,
                 t_initial=0,t_final=4*np.pi,num_point=2**14):
    ''' methods: should be a list containing strings'''
    #t_initial = 0
    num_stepsizes = 1
    stepsizes = np.zeros(num_stepsizes)
    num_points = np.zeros(num_stepsizes)
    #num_points[0] = 2**14-31
    num_points[0] = num_point
    #um_points[0]=2
    t = np.linspace(0,t_final,int(num_points[0]))
    stepsizes[0] = t[1]-t[0]
    for i in range(num_stepsizes-1):
        num_points[i+1] = num_points[i]/2
        print(num_points)
        t = np.linspace(0,t_final,int(num_points[i+1]))
        stepsizes[i+1]=t[1]-t[0]
    #print(stepsizes)

    data = [ [] for method in methods]
    for i in range(len(stepsizes)):
        for j in range(len(methods)):
            times = np.linspace(t_initial,t_final,int(num_points[i]))
            #time_1=timeit.time.time()

            data[j].append(ME.propagate(Hamiltonian,methods[j],'save_all',
                                        ME_params.params[params_name],times,mapping=mapping)
                            )
            #save_all
            #time_2=timeit.time.time()                    
            #print("Time took for {}: ".format(methods[j]),time_2-time_1)
    
    return t,data

def make_periodic_data(num_periods, times, matrix):
    new_matrix = np.linalg.matrix_power(matrix,num_periods)
    new_times =times + num_periods*times[-1]
    #print(new_matrix
    return new_times, new_matrix


def find_convergence_time(times,matrix,max_cycle=1000,print_out=False,initial_condition=None):
    '''initial_condition:the number of the column for the initial condition
        e.g. 8 means initial condition is in the ground state
        if None, then uses the Frobenius norm of the unitary operator
    '''
    for i in range(max_cycle):
        if(check_periodic(times,matrix,n=i,print_out=print_out,initial_condition=initial_condition) == True):
            #return make_periodic_data(i,times,matrix)[0][0]
            return i
    raise ValueError('exceeded max_cycles for convergence: {}'.format(max_cycle))
            #returns new_times, new_matrix
        
def check_periodic(times,matrix,n=10,print_out=False,initial_condition=None):
    '''initial_condition:the number of the column for the initial condition
        e.g. 8 means initial condition is in the ground state
        if None, then uses the Frobenius norm of the unitary operator
    '''
    new_times_1,new_matrix_1 = make_periodic_data(n,times,matrix)
    new_times_2,new_matrix_2 = make_periodic_data(n+1,times,matrix)
    if(initial_condition is None):
        error = ME.error_matrix(new_matrix_1,new_matrix_2)
    else:
        size = round(np.sqrt(new_matrix_1.shape[0]))
        error= ME.error_matrix(new_matrix_1[:,initial_condition].reshape(size,size),new_matrix_2[:,initial_condition].reshape(size,size))
    if(error > 1e-9):
        if(print_out==True):
            print("Not periodic, error is {}".format(error))
        return False
    return True

def get_averages(times,data,params):
    averages = np.zeros((9,1),dtype=np.complex128)
    transform = transformation_vector(times,params)
    #print(np.shape(transform))
    #print(np.shape(averages))
    transformed_data = np.einsum('ijl,ijk->ijkl',transform,data)
    #print(np.shape(transformed_data))
    averages = np.average(transformed_data,axis=0)
    #returns full matrix
    return averages

def transformation_vector(t,params):
    #t must be one-dimensional 
    vector = np.ones((len(t),9,1),dtype=np.complex128)
    w_p = params["np"]["w_p"]
    w_c = params["np"]["w_c"]
    vector[:,1,0] = np.exp(1j*w_c*t)
    vector[:,3,0] = np.conj(vector[:,1,0])
    vector[:,2,0] = np.exp(1j*w_p*t)
    vector[:,6,0] = np.conj(vector[:,2,0])
    vector[:,5,0] = np.exp(1j*(w_p-w_c)*t)
    vector[:,7,0] = np.conj(vector[:,5,0])
    return vector

def get_steadystate_data(Hamiltonian="lambda",methods=['M_14+M_23+M_32+M_41'],params="A",n=20,print_out=False,initial_condition=8):
    params_name = "TPR_{}".format(params)
    steady_state_data = np.zeros((33,9,9,1),dtype=np.complex128)
    for i in range(33): 
        full_params_name = params_name+"_{}".format(i)
        times,data = run_dynamics(full_params_name,Hamiltonian=Hamiltonian,methods=methods,t_final=8*np.pi,num_point=2**15)
        n = find_convergence_time(times,data[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=initial_condition)
        #find_convergence_time(times,matrix,max_cycle=1000,print_out=False,initial_condition=None)
        #check_periodic(times,data[0][0][-1],n-1,print_out=print_out)
        new_times,new_data = make_periodic_data(n,times,data[0][0][-1])
        steady_state_data[i] = get_averages(new_times[:-1],data[0][0][:-1]@new_data,ME_params.params[full_params_name])
    return new_times,steady_state_data

