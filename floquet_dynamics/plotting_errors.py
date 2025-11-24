import dynamics_functions as funcs
import pickle
import ME_params
import plotting_functions as plot
#import Magnus_Expansion as ME

def main():
    #params = ['A','B','C']
    params = ['B','A','C']
    #method = "sixth_order"
    method = "fourth_order"
    #begin plots for closed system dynamics
    all_times = []
    all_data = []
    all_RWA_data = []
    all_times_closed = []
    all_data_closed = []
    all_RWA_data_closed = []
    for param in params:
        
        with open("Data/{}_dynamics_closed_{}".format(param,method),"rb") as file:
            times,data = pickle.load(file)
        with open("Data/{}_dynamics_RWA_closed_{}".format(param,method),"rb") as file:
            times,RWA_data = pickle.load(file)
        
        n=0
        new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
        new_times,new_RWA_matrix = funcs.make_periodic_data(n,times,RWA_data[0][0][-1])
        all_times_closed.append(new_times[:-1])
        all_data_closed.append(data[0][0][:-1]@new_matrix)
        all_RWA_data_closed.append(RWA_data[0][0][:-1]@new_RWA_matrix)
        
        #save_name="TPR_both_{}_closed_{}".format(param,n+1)
        #make_error_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,
        #                          save_name=save_name,title=None, print_avg=True,x_fontsize=25, initial_condition = 8)
        
        with open("Data/{}_dynamics_{}".format(param, method),"rb") as file:
            times,data = pickle.load(file)
        with open("Data/{}_dynamics_RWA_{}".format(param, method),"rb") as file:
            times,RWA_data = pickle.load(file)
        
        n=0
        new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
        new_times,new_RWA_matrix = funcs.make_periodic_data(n,times,RWA_data[0][0][-1])
        
        all_times.append(new_times[:-1])
        all_data.append(data[0][0][:-1]@new_matrix)
        all_RWA_data.append(RWA_data[0][0][:-1]@new_RWA_matrix)

    plot.make_error_plot_both(all_times, all_data, all_RWA_data, save_name="open",text="II")
    plot.make_error_plot_both(all_times_closed, all_data_closed, all_RWA_data_closed, save_name="closed",text="I")
if __name__ == "__main__":
    main()
