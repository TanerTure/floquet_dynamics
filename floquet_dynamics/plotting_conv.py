import dynamics_functions as funcs
import pickle
import ME_params
import plotting_functions as plot

def main():
    #params = ['A','B','C']
    params = ['A']
    #method = "sixth_order"
    #method = "fourth_order"
    #begin plots for closed system dynamics
    for param in params:
        
        with open("Data/{}_dynamics_closed_fourth_order".format(param),"rb") as file:
            times,data_fourth_order = pickle.load(file)
        with open("Data/{}_dynamics_closed_sixth_order".format(param),"rb") as file:
            times, data_sixth_order = pickle.load(file)
        # with open("Data/{}_dynamics_RWA_closed_{}".format(param,method),"rb") as file:
        #     times,RWA_data = pickle.load(file)
        
        n=0
        new_times,new_matrix_fourth = funcs.make_periodic_data(n,times,data_fourth_order[0][0][-1])
        new_times,new_matrix_sixth = funcs.make_periodic_data(n,times,data_sixth_order[0][0][-1])
        
        
        save_name="Error_{}_closed_{}".format(param,n+1)
        plot.make_pop_plot_error(new_times[:-1],data_fourth_order[0][0][:-1]@new_matrix_fourth, data_sixth_order[0][0][:-1]@new_matrix_sixth, Hamiltonian = None,
                                  save_name=save_name,title=None, print_avg=True,x_fontsize=25, initial_condition = 8)
        
        
        with open("Data/{}_dynamics_fourth_order".format(param),"rb") as file:
            times,data_fourth_order = pickle.load(file)
        with open("Data/{}_dynamics_sixth_order".format(param),"rb") as file:
            times,data_sixth_order = pickle.load(file)
        
        n=0
        new_times,new_matrix_fourth = funcs.make_periodic_data(n,times,data_fourth_order[0][0][-1])
        new_times,new_matrix_sixth = funcs.make_periodic_data(n,times, data_sixth_order [0][0][-1]) 
        
        save_name="Error_{}_open_{}".format(param,n+1)
        plot.make_pop_plot_error(new_times[:-1],data_fourth_order[0][0][:-1]@new_matrix_fourth, data_sixth_order[0][0][:-1]@new_matrix_sixth, Hamiltonian = None,
                                  save_name=save_name,title=None, print_avg=True,x_fontsize=25, initial_condition = 8)
        
        n = funcs.find_convergence_time(times,data_fourth_order[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=8)
        print("n is ",n)
        #n_RWA = funcs.find_convergence_time(times,RWA_data[0][0][-1],max_cycle=10000000,print_out=False,initial_condition=8)
        #print("n_RWA is ", n_RWA)
        new_times,new_matrix_fourth = funcs.make_periodic_data(n,times,data_fourth_order[0][0][-1])
        new_times,new_matrix_sixth = funcs.make_periodic_data(n,times,data_sixth_order[0][0][-1])
        save_name="Error_{}_open_{}".format(param,n+1)
        plot.make_pop_plot_error(new_times[:-1],data_fourth_order[0][0][:-1]@new_matrix_fourth,data_sixth_order[0][0][:-1]@new_matrix_sixth, Hamiltonian = None,
                                  save_name=save_name,title=None, print_avg=True,x_fontsize=25, initial_condition = 8)
            
if __name__ == "__main__":
    main()