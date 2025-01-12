import plotting_functions as plot
import pickle


if __name__ == "__main__":
    params=['A','B','C','D']
    #params=['A']
    for param in params:
        with open("Data/steady_state_data_{}".format(param),"rb") as file:
            times_3D, data_3D = pickle.load(file)
        
        plot.make_3d_plot_times(times_3D,save_name="{}_time".format(param))
        plot.make_3d_plot_matrix_element(data_3D,i=0,j=8,save_name="{}_pop_1".format(param))
        
        with open("Data/steady_state_data_{}_RWA".format(param),"rb") as file:
            times_3D_RWA, data_3D_RWA = pickle.load(file)
        
        plot.make_3d_plot_times(times_3D_RWA,save_name="{}_time_RWA".format(param))
        plot.make_3d_plot_matrix_element(data_3D_RWA,i=0,j=8,save_name="{}_pop_1_RWA".format(param))
        
        #("Data/steady_state_data_B","wb")as