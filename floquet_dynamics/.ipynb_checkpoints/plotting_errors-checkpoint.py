import dynamics_functions as funcs
import pickle
import ME_params
import plotting_functions as plot


def main():
    #params = ['A','B','C']
    params = ['A','B','C']
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
        plot.make_error_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,
                                  save_name=save_name,title=None, print_avg=True,x_fontsize=25, initial_condition = 8)
        
        with open("Data/{}_dynamics_{}".format(param, method),"rb") as file:
            times,data = pickle.load(file)
        with open("Data/{}_dynamics_RWA_{}".format(param, method),"rb") as file:
            times,RWA_data = pickle.load(file)
        
        n=0
        new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
        new_times,new_RWA_matrix = funcs.make_periodic_data(n,times,RWA_data[0][0][-1])
        
        save_name="TPR_both_{}_open_{}".format(param,n+1)
        plot.make_error_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,
                                  save_name=save_name,title=None, print_avg=True,x_fontsize=25, initial_condition = 8)
        
        n = funcs.find_convergence_time(times,data[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=8)
        print("n is ",n)
        n_RWA = funcs.find_convergence_time(times,RWA_data[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=8)
        print("n_RWA is ", n_RWA)
        new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
        new_times,new_RWA_matrix = funcs.make_periodic_data(n_RWA,times,RWA_data[0][0][-1])
        save_name="TPR_both_{}_open_{}".format(param,n+1)
        plot.make_error_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,
                                  save_name=save_name,title=None, print_avg=True,x_fontsize=25, initial_condition = 8)
            
if __name__ == "__main__":
    main()
    
    
# make_error_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,save_name=name,title=None,
#                    print_avg=True,x_fontsize=fontsize,initial_condition=initial_condition)


# n=9
# initial_condition=8
# fontsize=25
# name = "TPR_both_C_open_{}".format(n+1)
# #name = "Test"
# new_times,new_matrix = make_periodic_data(n,times,data[0][0][-1])
# new_times_RWA,new_matrix_RWA=make_periodic_data(n,times,RWA_data[0][0][-1])

# make_error_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_matrix_RWA,Hamiltonian,save_name=name,title=None,
#                    print_avg=True,x_fontsize=fontsize,initial_condition=initial_condition)