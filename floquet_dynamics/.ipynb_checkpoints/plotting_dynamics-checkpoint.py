import dynamics_functions as funcs
import pickle
import ME_params
import plotting_functions as plot
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    #params = ['A','C','B']
    params = ['A','C','B']
    #method = "sixth_order"
    method = "fourth_order"
    #begin plots for closed system dynamics
    fig_to_legend_params = {"default":{"legend_fontsize":20, "loc":"best", "bbox_to_anchor":None},
                            ("A","closed", 0): {"legend_fontsize":19, "loc":"upper left", "bbox_to_anchor":(.21,1.02)},
                           }
    fig_to_off_diag_legend_params = { "default":{"legend_fontsize":[20,20,20],"loc":["best","best","best"], "bbox_to_anchor":[None, None, None],
                                                "extend":[0,0,0],"ncol":[1,1,1]},
                                     ("A","closed", 0):{"legend_fontsize":[20,18,20],"loc":["best","best","best"], "bbox_to_anchor":[None, None, None],
                                                "extend":[0,0,0]},
                                     ("C","closed",0):{"legend_fontsize":[18,18,18],"loc":["best","upper left","best"], "bbox_to_anchor":[None, (.3, .285), None],
                                                "extend":[5,0,0]},
                                     ("B","open",144):{"legend_fontsize":[17,20,17],"loc":["upper left","best","upper left"], "bbox_to_anchor":[None, None, None],
                                                "extend":[15,0,15],"ncol":[2,1,2]},
                                     ("A","open",10): {"legend_fontsize":[18,20,18], "loc":["upper left", "best", "upper left"], "bbox_to_anchor":[None,None,None],
                                                    "extend":[15,0,15],"ncol":[2,1,2]},
                                     ("C","open",7):{"legend_fontsize":[18,18,18],"loc":["upper left","upper left","upper left"], "bbox_to_anchor":[None, None, None],
                                                "extend":[15,15,15],"ncol":[2,2,2]},
                                     ("C","open",0):{"legend_fontsize":[19,20,20],"loc":["best","best","best"], "bbox_to_anchor":[None, None, None],
                                                "extend":[0,0,0],"ncol":[2,2,1]},
                                     ("B","open",143):{"legend_fontsize":[18,20,18],"loc":["upper left","best","upper left"], "bbox_to_anchor":[None, None, None],
                                                "extend":[15,0,15],"ncol":[2,1,2]}
                                     
                                    }
    for param in params:
        
        with open("Data/{}_dynamics_closed_{}".format(param,method),"rb") as file:
            times,data = pickle.load(file)
        with open("Data/{}_dynamics_RWA_closed_{}".format(param,method),"rb") as file:
            times,RWA_data = pickle.load(file)
        
        n=0
        new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
        new_times,new_RWA_matrix = funcs.make_periodic_data(n,times,RWA_data[0][0][-1])
        
        
        save_name="TPR_both_{}_closed_{}".format(param,n+1)
        fig_id = (param, "closed", n)
        if fig_id not in fig_to_legend_params:
            fig_id = "default"
        plot.make_pop_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,**fig_to_legend_params[fig_id])
        #make_off_diag_plots_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=18)
        fig_id = (param, "closed", n)
        if fig_id not in fig_to_off_diag_legend_params:
            fig_id = "default"
        plot.make_off_diag_plots_tilde_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,
                                       **fig_to_off_diag_legend_params[fig_id])
        
        
        with open("Data/{}_dynamics_{}".format(param, method),"rb") as file:
            times,data = pickle.load(file)
        with open("Data/{}_dynamics_RWA_{}".format(param, method),"rb") as file:
            times,RWA_data = pickle.load(file)
        
        n=0
        new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
        new_times,new_RWA_matrix = funcs.make_periodic_data(n,times,RWA_data[0][0][-1])
        
        
        save_name="TPR_both_{}_open_{}".format(param,n+1)
        fig_id = (param, "open", n)
        if fig_id not in fig_to_legend_params:
            fig_id = "default"
        plot.make_pop_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,**fig_to_legend_params[fig_id])
        #make_off_diag_plots_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=20)
        fig_id = (param, "open", n)
        if fig_id not in fig_to_off_diag_legend_params:
            fig_id = "default"
        plot.make_off_diag_plots_tilde_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,**fig_to_off_diag_legend_params[fig_id])
                  
        
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
        fig_id = (param, "open", n)
        if fig_id not in fig_to_legend_params:
            fig_id = "default"
        plot.make_pop_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,x_fontsize=x_fontsize, **fig_to_legend_params[fig_id])
        #make_off_diag_plots_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=20,x_fontsize=x_fontsize)
        fig_id = (param, "open", n)
        if fig_id not in fig_to_off_diag_legend_params:
            fig_id = "default"
        plot.make_off_diag_plots_tilde_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,x_fontsize=x_fontsize, **fig_to_off_diag_legend_params[fig_id])
        #new_times,new_data = make_periodic_data(n,times,data[0][0][-1])
    
    n=144
    x_fontsize=25
    new_times,new_matrix = funcs.make_periodic_data(n,times,data[0][0][-1])
    new_times,new_RWA_matrix = funcs.make_periodic_data(n,times,RWA_data[0][0][-1])
    
    save_name="TPR_both_{}_open_{}".format(param,n+1)
    fig_id = ('B', "open", n)
    if fig_id not in fig_to_legend_params:
        fig_id = "default"
    plot.make_pop_plot_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,x_fontsize=x_fontsize, **fig_to_legend_params[fig_id])
    
    #make_off_diag_plots_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,legend_fontsize=20,x_fontsize=x_fontsize)
    fig_id = ('B', "open", n)
    if fig_id not in fig_to_off_diag_legend_params:
        fig_id = "default"
    plot.make_off_diag_plots_tilde_both(new_times[:-1],data[0][0][:-1]@new_matrix,RWA_data[0][0][:-1]@new_RWA_matrix,"lambda",save_name=save_name,x_fontsize=x_fontsize, **fig_to_off_diag_legend_params[fig_id])