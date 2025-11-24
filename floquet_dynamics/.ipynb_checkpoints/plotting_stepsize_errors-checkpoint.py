import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pickle
import scipy

import dynamics_functions as funcs
import Magnus_Expansion as ME
import plotting_functions as plot

params = ['A','B','C']
method_names = ['fourth_order','CFME_equal','RK4','AM4_alt']
name_to_legend = {'fourth_order': 'Magnus Expansion',
                  'CFME_equal': 'Commutator-free ME',
                  'RK4': 'Runge-Kutta',
                  'AM4_alt':'Adams-Moulton'
                 }
file_strings = [["_open",""],["_open","_RWA"],["_closed",""],["_closed","_RWA"]]
for file_string in file_strings:
    for param in params:
        if(file_string[0] == "_open"):
            with open(f"Data/{param}_dynamics{file_string[1]}{''}_sixth_order","rb") as file:
                times, comparison = pickle.load(file)
        else:
            with open("Data/{}_dynamics{}{}_sixth_order".format(param,file_string[1],file_string[0]),"rb") as file:
                times, comparison = pickle.load(file)
                    
        comparison = comparison[0][0][-1]
        n_ss = 8
        comparison_ss = np.linalg.matrix_power(comparison,n_ss+1)
        x_vals_plots = [[],[]]
        y_vals_plots = [[],[]]
        y_vals_plots_ss = [[],[]]
        for i,method_name in enumerate(method_names):
            x_vals_plots[0].append([])
            x_vals_plots[1].append([])
            y_vals_plots[0].append([])
            y_vals_plots[1].append([])
            with open(f"Data/{param}_dynamics{file_string[1]}{file_string[0]}_stepsizes_{method_name}","rb") as file:
                stepsizes, data_method = pickle.load(file)
                for j in range(4):
                    x_vals_plots[0][-1].append(stepsizes[j])
                    y_vals_plots[0][-1].append(ME.error_matrix(comparison,data_method[0][j][-1]))
                    x_vals_plots[1][-1].append(stepsizes[j])
                    y_vals_plots[1][-1].append(ME.error_matrix(comparison_ss,np.linalg.matrix_power(data_method[0][j][-1], n_ss + 1)))

        
        file_name = param+file_string[1]+file_string[0]+"_both"
        
        legend_names = [name_to_legend[name] for name in method_names]
        plot.loglog_lob_2(x_vals_plots, y_vals_plots, legend_names, file_name=file_name)
        plot.make_text_files_stepsizes(x_vals_plots, y_vals_plots, method_names, file_name)