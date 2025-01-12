import dynamics_functions as funcs
import pickle
import ME_params
import plotting_functions as plot

if __name__ == "__main__":
    params = ['A','C','B']
    #method = "sixth_order"
    method = "fourth_order"
    #begin plots for closed system dynamics
    for param in params:
        
        with open("Data/{}_dynamics_closed_{}".format(param,method),"rb") as file:
            times,data = pickle.load(file)
        with open("Data/{}_dynamics_RWA_closed_{}".format(param,method),"rb") as file:
            times,RWA_data = pickle.load(file)
        
        n=0
        new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
        new_times,new_RWA_matrix = funcs.make_periodic_data(n,times,RWA_data[0][0][-1])
        
        
        save_name="TPR_both_{}_closed_{}".format(param,n+1)
        plot.make_pop_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18)
        #plot.make_off_diag_plots_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18)
        plot.make_off_diag_plots_tilde_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18)
        
        
        with open("Data/{}_dynamics_{}".format(param, method),"rb") as file:
            times,data = pickle.load(file)
        with open("Data/{}_dynamics_RWA_{}".format(param, method),"rb") as file:
            times,RWA_data = pickle.load(file)
        
        n=0
        new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
        new_times,new_RWA_matrix = funcs.make_periodic_data(n,times,RWA_data[0][0][-1])
        
        
        save_name="TPR_both_{}_open_{}".format(param,n+1)
        plot.make_pop_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18)
        plot.make_off_diag_plots_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18)
        plot.make_off_diag_plots_tilde_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18)
                  
        
        n = funcs.find_convergence_time(times,data[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=8)
        print("n is ",n)
        n_RWA = funcs.find_convergence_time(times,RWA_data[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=8)
        print("n_RWA is ", n_RWA)
        new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
        new_times,new_RWA_matrix = funcs.make_periodic_data(n_RWA,times,RWA_data[0][0][-1])
        
        save_name="TPR_both_{}_open_{}".format(param,n+1)
        if(param == "B"):
            x_fontsize=20
        else:
            x_fontsize=25
        plot.make_pop_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18,x_fontsize=x_fontsize)
        plot.make_off_diag_plots_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18,x_fontsize=x_fontsize)
        plot.make_off_diag_plots_tilde_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18,x_fontsize=x_fontsize)
        #new_times,new_data = make_periodic_data(n,times,data[0][0][-1])
    
    n=144
    x_fontsize=25
    new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
    new_times,new_RWA_matrix = funcs.make_periodic_data(n,times,RWA_data[0][0][-1])
    
    save_name="TPR_both_{}_open_{}".format(param,n+1)
    plot.make_pop_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18,x_fontsize=x_fontsize)
    plot.make_off_diag_plots_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18,x_fontsize=x_fontsize)
    plot.make_off_diag_plots_tilde_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18,x_fontsize=x_fontsize)
# n=0
# initial_condition=8
# fontsize=25
# #fontsize=20

# check_periodic(times,data[0][0][-1],n=n)
# check_periodic(times,RWA_data[0][0][-1],n=n)
# name = "TPR_both_C_closed_{}".format(n+1)
# #name = "Test"
# new_times,new_matrix = make_periodic_data(n,times,data[0][0][-1])
# new_times_RWA,new_matrix_RWA=make_periodic_data(n,times,RWA_data[0][0][-1])





# #new_times_RWA,new_matrix_RWA=make_periodic_data(n,analytical_times,RWA_data[0][0][-1])
# #for i in range(len(new_times_RWA)):
# #    analytical_data[i] = ME.L_lambda_RWA_closed(new_times_RWA[i],**ME_params.params["TPR_B"]["np"])
# #new_times,new_matrix = make_periodic_data(n,analytical_times,analytical_data[-1])



# make_pop_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_matrix_RWA,Hamiltonian,save_name=name,title=None,
#                    print_avg=True,x_fontsize=fontsize,initial_condition=initial_condition)
# make_off_diag_plots_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_matrix_RWA,Hamiltonian,save_name=name,title=None,
#                          print_avg=True,x_fontsize=fontsize,initial_condition=initial_condition)
# make_off_diag_plots_tilde_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_matrix_RWA,Hamiltonian,save_name=name,title=None,
#                                print_avg=True,x_fontsize=fontsize,initial_condition=initial_condition)
                  
                  
                  
                  
    #begin plots for open system dynamics
    
    
    
    #begin plots for steady state open system dynamics