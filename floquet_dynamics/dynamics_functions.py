import Magnus_Expansion as ME
import ME_params
import numpy as np
import timeit

def run_dynamics(params_name,Hamiltonian="lambda_RWA",methods = ['M_14+M_23+M_32+M_41'],mapping=False,
                 t_initial=0,t_final=4*np.pi,num_point=2**14):
    ''' methods: should be a list containing strings'''
    t_initial = 0
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
            time_1=timeit.time.time()

            data[j].append(ME.propagate(Hamiltonian,methods[j],'save_all',
                                        ME_params.params[params_name],times,mapping=mapping)
                            )
            #save_all
            time_2=timeit.time.time()                    
            print("Time took for {}: ".format(methods[j]),time_2-time_1)
    
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

