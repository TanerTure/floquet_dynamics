#parameter lists for Hamiltonian used
import sympy as sp
import numpy as np
import copy 

TPR_A_params = {"w_c":4, 
       "w_1" : 6, 
       "w_2":2,
       "w_p":6,
       "w_3" : 0,
       "Omega_c" : 0.56, 
       "Omega_p": 0.50,
       "gamma_12":0.90,
       "gamma_13": 1.00,
       "n_12": 0,
       "n_13": 0
    }

TPR_A_closed_params = {"w_c":4, 
       "w_1" : 6, 
       "w_2":2,
       "w_p":6,
       "w_3" : 0,
       "Omega_c" : 0.56, 
       "Omega_p": 0.50,
       "gamma_12":0,
       "gamma_13": 0,
       "n_12": 0,
       "n_13": 0
    }

TPR_B_params = {"w_c":4, 
       "w_1" : 6, 
       "w_2":2,
       "w_p":6,
       "w_3" : 0,
       "Omega_c" : 0.1, 
       "Omega_p": 0.1,
       "gamma_12":0.90,
       "gamma_13": 1.00,
       "n_12": 0,
       "n_13": 0
    }

TPR_B_closed_params = {"w_c":4, 
       "w_1" : 6, 
       "w_2":2,
       "w_p":6,
       "w_3" : 0,
       "Omega_c" : 0.1, 
       "Omega_p": 0.1,
       "gamma_12":0,
       "gamma_13": 0,
       "n_12": 0,
       "n_13": 0
    }

TPR_C_params = {"w_c":4, 
       "w_1" : 6, 
       "w_2":2,
       "w_p":6,
       "w_3" : 0,
       "Omega_c" : 5, 
       "Omega_p": 0.5,
       "gamma_12":0.90,
       "gamma_13": 1.00,
       "n_12": 0,
       "n_13": 0
    }

TPR_C_closed_params = {"w_c":4, 
       "w_1" : 6, 
       "w_2":2,
       "w_p":6,
       "w_3" : 0,
       "Omega_c" : 5, 
       "Omega_p": 0.5,
       "gamma_12":0,
       "gamma_13": 0,
       "n_12": 0,
       "n_13": 0
    }

TPR_D_params = {"w_c":4, 
       "w_1" : 6, 
       "w_2":2,
       "w_p":6,
       "w_3" : 0,
       "Omega_c" : 3, 
       "Omega_p": 3,
       "gamma_12":0.90,
       "gamma_13": 1.00,
       "n_12": 0,
       "n_13": 0
    }

TPR_D_closed_params = {"w_c":4, 
       "w_1" : 6, 
       "w_2":2,
       "w_p":6,
       "w_3" : 0,
       "Omega_c" : 3, 
       "Omega_p": 3,
       "gamma_12":0,
       "gamma_13": 0,
       "n_12": 0,
       "n_13": 0
    }

# STIRAP_closed_params = {"w_c":(2-4/6)*np.pi*2.5*10**6, 
#        "w_1" : 4*np.pi*2.5*10**6, 
#        "w_2":4/6*np.pi*2.5*10**6,
#        "w_p":2*np.pi*2.5*10**6,
#        "w_3" : 0,
#        "Omega_c" : 2*np.pi*5*10**3/2, 
#        "Omega_p": 2*np.pi*5*10**3/2,
#        "gamma_12":0,
#        "gamma_13": 0,
#        "n_12": 0,
#        "n_13": 0
#     }

# STIRAP_params = {"w_c":(2-4/6)*np.pi*2.5*10**6, 
#        "w_1" : 4*np.pi*2.5*10**6, 
#        "w_2":4/6*np.pi*2.5*10**6,
#        "w_p":2*np.pi*2.5*10**6,
#        "w_3" : 0,
#        "Omega_c" : 2*np.pi*5*10**3/2, 
#        "Omega_p": 2*np.pi*5*10**3/2,
#        "gamma_12":0.90,
#        "gamma_13": 1.00,
#        "n_12": 0,
#        "n_13": 0
#     }

params = {
    "TPR_A":{"sp":{},"np":TPR_A_params},
    "TPR_A_closed":{"sp":{},"np":TPR_A_closed_params},
    "TPR_B":{"sp":{},"np":TPR_B_params},
    "TPR_B_closed":{"sp":{},"np":TPR_B_closed_params},
    "TPR_C":{"sp":{},"np":TPR_C_params},
    "TPR_C_closed":{"sp":{},"np":TPR_C_closed_params},
    "TPR_D":{"sp":{},"np":TPR_D_params},
    "TPR_D_closed":{"sp":{},"np":TPR_D_closed_params},
    "empty":{"sp":{},"np":{}},
    }

def make_w_p_lambda_params(string="A"):
    '''Function for creating params for varying w_p figures '''
    params_name="TPR_{}".format(string)
    w_ps = np.linspace(2,10,33,dtype=np.float64)
    for i in range(len(w_ps)):
        parameters = copy.deepcopy(params[params_name])
        parameters["np"]["w_p"] = w_ps[i]
        params[params_name+"_{}".format(i)] = copy.deepcopy(parameters)
    return

make_w_p_lambda_params("A")
make_w_p_lambda_params("B")
make_w_p_lambda_params("C")

def make_w_p_w_c_lambda_params(string="A"):
    '''Function for creating params for 3d plots (w_p,w_c vary)'''
    params_name="TPR_{}".format(string)
    w_ps = np.linspace(2,10,33,dtype=np.float64)
    w_cs = np.linspace(1,7,25,dtype=np.float64)
    for i in range(len(w_ps)):
        for j in range(len(w_cs)):
            parameters = copy.deepcopy(params[params_name])
            parameters["np"]["w_p"] = w_ps[i]
            parameters["np"]["w_c"] = w_cs[j]
            params[params_name+"_{}_{}".format(i,j)] = copy.deepcopy(parameters)
    return

make_w_p_w_c_lambda_params("A")
make_w_p_w_c_lambda_params("B")
make_w_p_w_c_lambda_params("C")
make_w_p_w_c_lambda_params("D")


# def create_lambda_params(w_1,w_2,w_3,w_c,w_p,
#                          d_12,
#                          d_13,
#                          E_c,
#                          E_p,
#                         ):
    
#     sp_parameters = { 
#             "w_1":w_1,
#             "w_2":w_2,
#             "w_3":w_3,
#             "w_c":w_c,
#             "w_p":w_p,
#             "d_12":d_12,
#             "d_13":d_13,
#             "E_c":E_c,
#             "E_p":E_p,
#         }
    
#     np_parameters = {
#             "w_1":np.complex128(sp_parameters["w_1"]),
#             "w_2":np.complex128(sp_parameters["w_2"]),
#             "w_3":np.complex128(sp_parameters["w_3"]),
#             "w_c":np.complex128(sp_parameters["w_c"]),
#             "w_p":np.complex128(sp_parameters["w_p"]),
#             "d_12":np.array(sp_parameters["d_12"]).astype(np.complex128),
#             "d_13":np.array(sp_parameters["d_13"]).astype(np.complex128),
#             "E_c":np.array(sp_parameters["E_c"]).astype(np.complex128),
#             "E_p":np.array(sp_parameters["E_p"]).astype(np.complex128), 
#         }
#     parameters = {
#             "sp":sp_parameters,
#             "np":np_parameters
#     }
#     return parameters


# empty_params = {}