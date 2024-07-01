import numpy as np
import sympy as sp
import scipy.linalg #for scipi.linalg.expm

#delta_t, t_k, t_k1 = sp.symbols('\delta\t t_k,t_k+1',real=True)
delta_t = sp.symbols('\delta\t',real=True)


def comm(a,b):
    return a @ b - b@a
def comm_sp(a,b):
    return a*b - b*a

def H_three_state(t, w_c=1,w_1=2,w_2=1,
                  w_p =2,w_3 =0,
                  Omega_c=0.56, Omega_p=0.50,
                 ):
    
    #assume hbar=1
    H = np.zeros((3,3),dtype=np.complex128) 
    H[0][0]=w_1
    H[1][1]=w_2
    H[2][2]=w_3
    H[0][1]=-Omega_c*2*np.cos(w_c*t)
    H[1][0]=H[0][1]
    H[0][2]=-Omega_p*2*np.cos(w_p*t)
    H[2][0]=H[0][2]
    H*=-1j
    return H



def L_lambda (t, w_c=1,w_1=2,w_2=1,
                  w_p =2,w_3 =0,
                  Omega_c=0.56, Omega_p=0.50,
                  gamma_12=0.9, gamma_13=1,
                 n_12 = 0,n_13 = 0
                ):
    L_sys_int = (np.kron(np.eye(3),H_three_state(t,w_c=w_c,w_1=w_1,
                                                w_2=w_2,w_p =w_p,w_3 =w_3,
                                                Omega_c=Omega_c, Omega_p=Omega_p,
                                               ).T)
    - np.kron(H_three_state(t,w_c=w_c,w_1=w_1,
            w_2=w_2,w_p =w_p,w_3 =w_3,
            Omega_c=Omega_c, Omega_p=Omega_p,
            ), np.eye(3))
                )*-1
    L_lambda_dissipative(L_sys_int,gamma_12=gamma_12, gamma_13=gamma_13,
                 n_12 = n_12,n_13 = n_13
                        )
    
    return L_sys_int

def H_three_state_RWA(t, w_c=1,w_1=2,w_2=1,
                  w_p =2,w_3 =0,
                  Omega_c=0.56, Omega_p=0.50,
                 ):
    
    #assume hbar=1
    H = np.zeros((3,3),dtype=np.complex128)
    H[0][0]=w_1
    H[1][1]=w_2
    H[2][2]=w_3
    H[0][1]=-Omega_c*np.exp(-1j*w_c*t)
    H[1][0]=np.conj(H[0][1])
    H[0][2]=-Omega_p*np.exp(-1j*w_p*t)
    H[2][0]=np.conj(H[0][2])
    H*=-1j
    return H



def U_RWA_closed(t,w_c=1,w_1=2,w_2=1,
       w_p =2,w_3 =0,
       Omega_c=0.56, Omega_p=0.50
      ): #Only good at TPR!
    U = np.zeros((3,3),dtype=np.complex128)
    val = np.sqrt(Omega_c**2+Omega_p**2)
    U[2][2] = (Omega_c**2+Omega_p**2*np.cos(val*t))/val**2
    U[2][0] = 1j*Omega_p*np.sin(val*t)/(2*val)
    U[2][1] = Omega_c*Omega_p*(np.cos(val*t)-1)/val**2
    U[0][0] = np.cos(val*t)
    U[0][1] = 1j*Omega_c*np.sin(val*t)/(2*val)
    U[1][1] = (Omega_c**2*np.cos(val*t)+Omega_p**2)/(val**2)
    U[0][2] = U[2][0]
    U[1][2] = U[2][1]
    U[1][0] = U[0][1]
   
    return U

def W_RWA(t,w_c=1,w_1=2,w_2=1,
      w_p =2,w_3 =0,
      Omega_c=0.56, Omega_p=0.50
         ):#only good at TPR, no environment
    W = np.zeros((3,3),dtype=np.complex128)
    W[0][0] = np.exp(-1j*w_p*t)
    W[1][1] = np.exp(-1j*(w_p-w_c)*t)
    W[2][2] = 1
    return W

def L_lambda_RWA_closed (t, w_c=1,w_1=2,w_2=1,
                  w_p =2,w_3 =0,
                  Omega_c=0.56, Omega_p=0.50,
                  gamma_12=0.9, gamma_13=1,
                 n_12 = 0,n_13 = 0
                ):
    full_U = W_RWA(t,w_c=w_c,w_1=w_1,
                   w_2=w_2,w_p =w_p,w_3 =w_3,
                  Omega_c=Omega_c, Omega_p=Omega_p,
                  )@U_RWA_closed(t,w_c=w_c,w_1=w_1,
                   w_2=w_2,w_p =w_p,w_3 =w_3,
                  Omega_c=Omega_c, Omega_p=Omega_p,
                                )
    return np.kron(full_U,np.conj(full_U.T).T)


def L_lambda_RWA (t,w_c=1,w_1=2,w_2=1,
                  w_p =2,w_3 =0,
                  Omega_c=0.56, Omega_p=0.50,
                  gamma_12=0.9, gamma_13=1,
                 n_12 = 0,n_13 = 0
                 ):
    L_sys_int = (np.kron(np.eye(3),H_three_state_RWA(t,w_c=w_c,w_1=w_1,
                                                w_2=w_2,w_p =w_p,w_3 =w_3,
                                                Omega_c=Omega_c, Omega_p=Omega_p,
                                               ).T)
    - np.kron(H_three_state_RWA(t,w_c=w_c,w_1=w_1,
            w_2=w_2,w_p =w_p,w_3 =w_3,
            Omega_c=Omega_c, Omega_p=Omega_p,
            ), np.eye(3))
                )*-1 #since the factor of -i is included in the Hamiltonian
    
    L_lambda_dissipative(L_sys_int,gamma_12=gamma_12, gamma_13=gamma_13,
                 n_12 = n_12,n_13 = n_13
                        )
    return L_sys_int

    
def L_lambda_dissipative( L_sys_int,gamma_12=0.9, gamma_13=1,
                 n_12 = 0,n_13 = 0):
    
    a = gamma_12*(1+n_12)
    b = gamma_12*n_12
    c = gamma_13*(1+n_13)
    d = gamma_13*n_13
    L_sys_int[0][0] += -1*(a+c)
    L_sys_int[1][1] += -1/2*(a+b+c)
    L_sys_int[2][2] += -1/2*(a+c+d)
    L_sys_int[3][3] += -1/2*(a+b+c)
    L_sys_int[4][4] += -b
    L_sys_int[5][5] += -1/2*(b+d)
    L_sys_int[6][6] += -1/2*(a+c+d)
    L_sys_int[7][7] += -1/2*(b+d)
    L_sys_int[8][8] += -d
    L_sys_int[4][0] += a
    L_sys_int[8][0] += c
    L_sys_int[0][4] += b
    L_sys_int[0][8] += d
    
    return

def H_three_state(t, w_c=1,w_1=2,w_2=1,
                  w_p =2,w_3 =0,
                  Omega_c=0.56, Omega_p=0.50,
                 ):
    
    #assume hbar=1
    H = np.zeros((3,3),dtype=np.complex128)
    H[0][0]=w_1
    H[1][1]=w_2
    H[2][2]=w_3
    H[0][1]=-Omega_c*np.cos(w_c*t)
    H[1][0]=H[0][1]
    H[0][2]=-Omega_p*np.cos(w_p*t)
    H[2][0]=H[0][2]
    H*=-1j
    return H 
 


def M_14_M_23_M_32_M_41(points,dt):
    M_14 = dt/90*(7*points[0]+32*points[1]+12*points[3]+32*points[5]+7*points[6])
    M_23 = (dt**2/6720*(comm(232/39*points[0]+1152/13*points[2]+72*points[4],
                            2*points[0]+81/8*points[2]-13/8*points[6]))
           +dt**2/180*comm(points[6],points[0])
           )
    M_32 = dt**3/15120*(64*(comm(points[3]+points[6],comm(points[3],points[0]))
                            +comm(points[3]+points[0],comm(points[3],points[6])))
                        +44*(comm(points[0],comm(points[0],points[3]))
                             +comm(points[6],comm(points[6],points[3])))
                        +9*comm(points[6]-points[0],comm(points[6],points[0]))
                       )
    #M_32=dt**3/240*comm(points[6]-points[0],comm(points[6],points[0]))
    #M_32=0
    c = -(5+np.sqrt(21))/2
    M_41 =dt**4/5040*comm(1/c*points[0]-points[6],comm(points[6]-c*points[0],comm(points[6],points[0])))
    #M_41=0
    return scipy.linalg.expm(M_14+M_23+M_32+M_41)

def M_12_M_21(points,dt):
    M_12 = dt/6*(points[0]+4*points[1]+points[2])
    M_21 = .5*dt**2/6*comm(points[2],points[0])
   # M_2 = 0
    return scipy.linalg.expm(M_12+M_21)

def M_12_M_22(points,dt):
    M_12 = dt/6*(points[0]+4*points[1]+points[2])
    M_22 = dt**2/15*comm(points[0]/4+points[1],points[0]-points[2])
    return scipy.linalg.expm(M_12+M_22)


def save_final(U,i,saved_data,total=0):
    if(i == 0):
        length = len(U[0])
        data = np.zeros((1,length,length),dtype=np.complex128)
        data[0] = U
        return data
    else:
        saved_data[0] = U
        return

        
def save_all(U,i,saved_data,total=0):
    length = len(U[0])
    if(i ==0 ):
        data = np.zeros((total,length,length),dtype=np.complex128)
        data[0] = U
        return data
    else:
        saved_data[i] = U
        return

def get_funcs_from_names(H_name,method_name,save_name):
    method_func = method_name_to_func[method_name]
    save_func = save_name_to_func[save_name]
    H_func = H_name_to_func[H_name]
    return H_func,method_func, save_func

def propagate(H_name,method_name,save_name,H_parameters,times,mapping=False):
    H_func,method_func,save_func = get_funcs_from_names(H_name,method_name,save_name)
    dt = times[1]-times[0]
   # print(dt)
    if(mapping==False):
        if (method_func[0] == "save_point"):
            return propagate_save_point(H_func,method_func,save_func, H_parameters["np"],dt,times)
    if(mapping==True):
        if (method_func[0] == "save_point"):
            return propagate_save_point_mapping(H_func,method_func,save_func, H_parameters["np"],dt,times)


        
def propagate_save_point_mapping(H_func,method_func,save_func, H_parameters,dt,times):
    points = [np.complex128(point.subs(delta_t,dt)) for point in method_func[1]]
    points = np.asarray(points,dtype=np.complex128)
    H_function = H_func[0] #H_func also stores additional info, size of matrix
    method_function = method_func[2] #method_function stores additional info
    
    H_points = [H_function(point+times[0],**H_parameters) for point in points]

    
    U = np.eye(H_func[1],dtype=np.complex128)
    saved_data = save_func(U,0,None,len(times))

    for i in range(len(times)-1):
        
      
        U = method_function(H_points,dt)
        save_func(U,i+1,saved_data)
        H_points[0] = H_points[-1]
        H_points[1:] =[H_function(point+times[i+1],**H_parameters) for point in points[1:]]
        
      
    return saved_data
    
                                                                    
                                    
    
def propagate_save_point(H_func,method_func,save_func, H_parameters,dt,times):
    points = [np.complex128(point.subs(delta_t,dt)) for point in method_func[1]]
    points = np.asarray(points,dtype=np.complex128)
    H_function = H_func[0] #H_func also stores additional info, size of matrix
    method_function = method_func[2] #method_function stores additional info
    #H_points = np.zeros((len(points),H_func[1],H_func[1]),dtype=np.complex128)
    H_points = [H_function(point+times[0],**H_parameters) for point in points]

    #print(H_points.shape)
    #for i in range(len(points)):
        #H_points[i,:,:] = H_function(points[i],**H_parameters)
   # np.array((H_function(point,**H_parameters) for point in points),dtype=np.complex128)
    U = np.eye(H_func[1],dtype=np.complex128)
    saved_data = save_func(U,0,None,len(times))

    for i in range(len(times)-1):
        
       # print(U_dt)
        U = method_function(H_points,dt)@U
        save_func(U,i+1,saved_data)
        H_points[0] = H_points[-1]
        #for i in range(1,len(points)):
            #H_points[i,:,:] = H_function(points[i],**H_parameters)
        H_points[1:] =[H_function(point+times[i+1],**H_parameters) for point in points[1:]]
        
      
    return saved_data

def Froebinus_norm(matrix):
    sum=0
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            sum += np.abs(matrix[i,j])**2
    return np.sqrt(sum)

H_name_to_func = {
                  "lambda":[L_lambda,9],
                  "lambda_RWA":[L_lambda_RWA,9],
                 }
#H_evaluation function, the size of row of H matrix, the integral function

method_name_to_func = {
                       "M_12+M_21":["save_point",[sp.Integer(0),delta_t/2,delta_t],M_12_M_21],
                       
                       "M_12+M_22":["save_point",[sp.Integer(0),delta_t/2,delta_t],M_12_M_22],
                       
                       "M_14+M_23+M_32+M_41":["save_point",[sp.Integer(0),delta_t/4,delta_t/3,delta_t/2,
                                                            2*delta_t/3,3*delta_t/4,delta_t],M_14_M_23_M_32_M_41],

}

save_name_to_func = {"U_final":save_final,
                     "save_all":save_all #saves all unitary in interval
}

def Froebinus_norm(matrix):
    sum=0
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            sum += np.abs(matrix[i,j])**2
    return np.sqrt(sum)

def vector_norm(vec):
    return np.sum(np.abs(vec)**2)
    
                   
def error_matrix(exp,obs):
    return Froebinus_norm(exp-obs)/Froebinus_norm(exp)


folder_name = {
               #"lambda":"C:/Users/Taner/Desktop/GIT/Research/Notes/Remarks/lambda_RWA_paper/Figures/lambda/",
               "lambda":"",
               #"lambda":"C:\\Users\\Taner\\Desktop\\GIT\\Research\\Notes\\Remarks\\lambda_RWA_paper\\Figures\\lambda\\",

               "lambda_RWA":"C:/Users/Taner/Desktop/GIT/Research/Notes/Remarks/lambda_RWA_paper/Figures/lambda_RWA/",
               }


import pickle
def save_data(filename,data,Hamiltonian,comparison=None):
    filename = folder_name[Hamiltonian]+filename
    print(filename)
    with open(filename,'wb') as file:
        pickle.dump(data,file)
        if(comparison is not None):
            pickle.dump(comparison,file)

def get_data(filename,Hamiltonian,comparison=False):
    #input filename should be just the filename; 
    #must be placed in the folder corresponding to the Hamiltonian
    filename=folder_name[Hamiltonian]+filename
    with open(filename,'rb') as file:
        data = pickle.load(file)
        if(comparison==True):
            comparison_data = pickle.load(file)
    
    if(comparison==True):
        return data,comparison_data
    else:
        return data