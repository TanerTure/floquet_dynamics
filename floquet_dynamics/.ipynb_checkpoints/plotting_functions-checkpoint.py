import steady_state as ss
import Magnus_Expansion as ME
import ME_params
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
plt.rcParams["mathtext.fontset"] = "cm"

def make_3d_plot_times(convergence_times_3d,save_name="Test_conv_times"):
    matplotlib.rc('xtick', labelsize=20) 
    matplotlib.rc('ytick', labelsize=20)
    dw_ps = np.linspace(2,10,33)-6
    dw_cs = np.linspace(1,7,25)-4
    dW_Ps,dW_Cs = np.meshgrid(dw_ps,dw_cs)
    fig,ax = plt.subplots(1,1)
    mappable = ax.contourf(dW_Ps,dW_Cs,np.einsum('ij->ji',convergence_times_3d),levels=100)
    fig.colorbar(mappable)
    ax.set_xlabel(r'$\delta\omega_p$',fontsize=25)
    ax.set_ylabel(r'$\delta\omega_c$',fontsize=25)
    fig.savefig("Figures/"+"3D/"+save_name,dpi=300,bbox_inches='tight')

def make_3d_plot_matrix_element(steady_state_3d,i=0,j=0,save_name="Test_matrix"):
    fig,ax = plt.subplots(1,1)
    dw_ps = np.linspace(2,10,33)-6
    dw_cs = np.linspace(1,7,25)-4
    dW_Ps,dW_Cs = np.meshgrid(dw_ps,dw_cs)
    mappable = plt.contourf(dW_Ps,dW_Cs,np.einsum('ij->ji',steady_state_3d[:,:,i,j,0]),levels=100)
    fig.colorbar(mappable)
    ax.set_xlabel(r'$\delta\omega_p$',fontsize=25)
    ax.set_ylabel(r'$\delta\omega_c$',fontsize=25)
    fig.savefig("Figures/"+"3D/"+save_name,dpi=300,bbox_inches='tight')


#Using new conventions |1> in ground state |2> in excited state, etc.
def make_pop_plot_both(times,data,data_RWA,Hamiltonian,save_name=None,title=None,print_avg=True,x_fontsize=25,
                       legend_fontsize=18,initial_condition=8):
    #data should only have data from one method
    plt.plot(times,data[:,8,initial_condition],label=r"$\rho_{11}$")
    print(np.max(np.abs(np.imag(data[:,8,initial_condition]))))
    print(np.max(np.abs(np.imag(data[:,0,initial_condition]))))
    print(np.max(np.abs(np.imag(data[:,4,initial_condition]))))


    plt.plot(times,data[:,0,initial_condition],label=r"$\rho_{22}$")
    plt.plot(times,data[:,4,initial_condition],label=r"$\rho_{33}$")
    plt.plot(times,data_RWA[:,8,initial_condition],"C0--",linewidth=2,label=r"$\rho_{11}$RWA")
    plt.plot(times,data_RWA[:,0,initial_condition],"C1--",linewidth=2,label=r"$\rho_{22}$RWA")
    plt.plot(times,data_RWA[:,4,initial_condition],"C2--",linewidth=2,label=r"$\rho_{33}$RWA")
    print(np.max(np.abs(np.imag(data_RWA[:,8,initial_condition]))))
    print(np.max(np.abs(np.imag(data_RWA[:,0,initial_condition]))))
    print(np.max(np.abs(np.imag(data_RWA[:,4,initial_condition]))))
    
    plt.xlabel(r"$t$",fontsize=20)
    if(print_avg==True):
        print(np.average(data[:,8,initial_condition]))
        print(np.average(data[:,0,initial_condition]))
        print(np.average(data[:,4,initial_condition]))
        print("RWA_pop_1",np.average(np.real(data_RWA[:,8,initial_condition])))
        print("RWA_pop_2",np.average(np.real(data_RWA[:,0,initial_condition])))
        print("RWA_pop_3",np.average(np.real(data_RWA[:,4,initial_condition])))


              
    plt.legend(fontsize=legend_fontsize)
    plt.xticks(fontsize=x_fontsize)
    plt.yticks(fontsize=25)
    if save_name is None:
        pass
    else:
        plt.savefig("Figures/"+"dynamics/"+save_name+"_pops_both",dpi=300,bbox_inches='tight')
        #ME.save_data(save_name,data,Hamiltonian,comparison=comparison)
    plt.show()

def make_off_diag_plots_both(times,data,data_RWA,Hamiltonian,save_name=None,title=None,print_avg=True,x_fontsize=25,
                            legend_fontsize=18,initial_condition=8):
    plt.plot(times,np.real(data[:,6,initial_condition]),"C3",label=r"$\rho_{12}^R$")
    plt.plot(times,np.imag(data[:,6,initial_condition]),"C4",label=r"$\rho_{12}^I$")
    plt.plot(times,np.real(data_RWA[:,6,initial_condition]),"C3--",linewidth=2,label=r"$\rho_{12}^RRWA$")
    plt.plot(times,np.imag(data_RWA[:,6,initial_condition]),"C4--",linewidth=2,label=r"$\rho_{12}^IRWA$")
    if(print_avg==True):
        print(np.average(np.real(data[:,6,initial_condition])))
        print(np.average(np.imag(data[:,6,initial_condition])))
        print("RWA_real_avg",np.average(np.real(data_RWA[:,6,initial_condition])))
        print("RWA_imag_avg",np.average(np.imag(data_RWA[:,6,initial_condition])))

    plt.xlabel(r"$t$",fontsize=20)
    #plt.axes
    plt.legend(fontsize=legend_fontsize)
    plt.xticks(fontsize=x_fontsize)
    plt.yticks(fontsize=25)
    if save_name is None:
        pass
    else:
        plt.savefig("Figures/"+"dynamics/"+save_name+"_p12_both",dpi=300,bbox_inches='tight')
    plt.show()
    plt.plot(times,np.real(data[:,7,initial_condition]),"C3",label=r"$\rho_{13}^R$")
    plt.plot(times,np.imag(data[:,7,initial_condition]),"C4",label=r"$\rho_{13}^I$")
    plt.plot(times,np.real(data_RWA[:,7,initial_condition]),"C3--",linewidth=2,label=r"$\rho_{13}^RRWA$")
    plt.plot(times,np.imag(data_RWA[:,7,initial_condition]),"C4--",linewidth=2,label=r"$\rho_{13}^IRWA$")
    if(print_avg==True):
        print(np.average(np.real(data[:,7,initial_condition])))
        print(np.average(np.imag(data[:,7,initial_condition])))
        print("RWA_real_avg",np.average(np.real(data_RWA[:,7,initial_condition])))
        print("RWA_imag_avg",np.average(np.imag(data_RWA[:,7,initial_condition])))

    plt.legend(fontsize=legend_fontsize)
    plt.xticks(fontsize=x_fontsize)
    plt.yticks(fontsize=25)
    plt.xlabel(r"$t$",fontsize=20)

    if save_name is None:
        pass
    else:
        plt.savefig("Figures/"+"dynamics/"+save_name+"_p13_both",dpi=300,bbox_inches='tight')
    plt.show()
    
    plt.plot(times,np.real(data[:,1,initial_condition]),"C3",label=r"$\rho_{23}^R$")
    plt.plot(times,np.imag(data[:,1,initial_condition]),"C4",label=r"$\rho_{23}^I$")
    plt.plot(times,np.real(data_RWA[:,1,initial_condition]),"C3--",linewidth=2,label=r"$\rho_{23}^R$RWA")
    plt.plot(times,np.imag(data_RWA[:,1,initial_condition]),"C4--",linewidth=2,label=r"$\rho_{23}^I$RWA")
    if(print_avg==True):
        print(np.average(np.real(data[:,1,initial_condition])))
        print(np.average(np.imag(data[:,1,initial_condition])))
        print("RWA_real_avg",np.average(np.real(data_RWA[:,1,initial_condition])))
        print("RWA_imag_avg",np.average(np.imag(data_RWA[:,1,initial_condition])))
    plt.legend(fontsize=legend_fontsize)
    plt.xticks(fontsize=x_fontsize)
    plt.yticks(fontsize=25)
    plt.xlabel(r"$t$",fontsize=20)
    if save_name is None:
        pass
    else:
        plt.savefig("Figures/"+"dynamics/"+save_name+"_p23_both",dpi=300,bbox_inches='tight')
    plt.show()
    
def make_off_diag_plots_tilde_both(times,data,data_RWA,Hamiltonian,save_name=None,w_p=6,w_c=4,title=None,print_avg=True,x_fontsize=25,
                                  legend_fontsize=18,initial_condition=8):
    plt.plot(times,np.real(data[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C3",label=r"$\tilde{\rho}_{12}^R$")
    plt.plot(times,np.imag(data[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C4",label=r"$\tilde{\rho}_{12}^I$")
    plt.plot(times,np.real(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C3--",linewidth=2,label=r"$\tilde{\rho}_{12}^R$RWA")
    plt.plot(times,np.imag(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C4--",linewidth=2,label=r"$\tilde{\rho}_{12}^I$RWA")
    if(print_avg==True):
        print(np.average(np.real(data[:,6,initial_condition]*np.exp(-1j*w_p*times))))
        print(np.average(np.imag(data[:,6,initial_condition]*np.exp(-1j*w_p*times))))
        print("RWA_real_avg",np.average(np.real(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times))))
        print("RWA_imag_avg",np.average(np.imag(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times))))
    plt.legend(fontsize=legend_fontsize)
    plt.xticks(fontsize=x_fontsize)
    plt.yticks(fontsize=25)
    plt.xlabel(r"$t$",fontsize=20)

    if save_name is None:
        pass
    else:
        plt.savefig("Figures/"+"dynamics/"+save_name+"_p12tilde_both",dpi=300,bbox_inches='tight')
    plt.show()
    plt.plot(times,np.real(data[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times)),"C3",label=r"$\tilde{\rho}_{13}^R$")
    plt.plot(times,np.imag(data[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times)),"C4",label=r"$\tilde{\rho}_{13}^I$")
    plt.plot(times,np.real(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times)),"C3--",linewidth=2,label=r"$\tilde{\rho}_{13}^R$RWA")
    plt.plot(times,np.imag(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times)),"C4--",linewidth=2,label=r"$\tilde{\rho}_{13}^I$RWA")
    print(np.average(np.real(data[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print(np.average(np.imag(data[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_real_avg",np.average(np.real(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_imag_avg",np.average(np.imag(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_real_max",np.max(np.real(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_imag_max",np.max(np.imag(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_real_min",np.min(np.real(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_imag_min",np.min(np.imag(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))


    
    plt.legend(fontsize=legend_fontsize)
    plt.xticks(fontsize=x_fontsize)
    plt.yticks(fontsize=25)
    plt.xlabel(r"$t$",fontsize=20)

    if save_name is None:
        pass 
    else:
        plt.savefig("Figures/"+"dynamics/"+save_name+"_p13tilde_both",dpi=300,bbox_inches='tight')
    plt.show()
    
    plt.plot(times,np.real(data[:,1,initial_condition]*np.exp(1j*w_c*times)),"C3",label=r"$\tilde{\rho}_{23}^R$")
    plt.plot(times,np.imag(data[:,1,initial_condition]*np.exp(1j*w_c*times)),"C4",label=r"$\tilde{\rho}_{23}^I$")
    plt.plot(times,np.real(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times)),"C3--",linewidth=2,label=r"$\tilde{\rho}_{23}^R$RWA")
    plt.plot(times,np.imag(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times)),"C4--",linewidth=2,label=r"$\tilde{\rho}_{23}^I$RWA")
    print(np.average(np.real(data[:,1,initial_condition]*np.exp(1j*w_c*times))))
    print(np.average(np.imag(data[:,1,initial_condition]*np.exp(1j*w_c*times))))
    print("RWA_real_avg",np.average(np.real(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times))))
    print("RWA_real_max",np.max(np.real(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times))))
    print("RWA_real_min",np.min(np.real(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times))))

    
    # plt.plot(times,np.real(data[:,1,initial_condition]*np.exp(-1j*w_c*times)),label=r"$\tilde{\rho}_{23}^R$")
    # plt.plot(times,np.imag(data[:,1,initial_condition]*np.exp(-1j*w_c*times)),label=r"$\tilde{\rho}_{23}^I$")
    # plt.plot(times,np.real(data_RWA[:,1,initial_condition]*np.exp(-1j*w_c*times)),"r--",linewidth=2,label=r"$\tilde{\rho}_{23}^R$RWA")
    # plt.plot(times,np.imag(data_RWA[:,1,initial_condition]*np.exp(-1j*w_c*times)),"m--",linewidth=2,label=r"$\tilde{\rho}_{23}^I$RWA")
    # print(np.average(np.real(data[:,1,initial_condition]*np.exp(-1j*w_c*times))))
    # print(np.average(np.imag(data[:,1,initial_condition]*np.exp(-1j*w_c*times))))
    # print("RWA_real_avg",np.average(np.real(data_RWA[:,1,initial_condition]*np.exp(-1j*w_c*times))))
    # print("RWA_real_max",np.max(np.real(data_RWA[:,1,initial_condition]*np.exp(-1j*w_c*times))))
    # print("RWA_real_min",np.min(np.real(data_RWA[:,1,initial_condition]*np.exp(-1j*w_c*times))))


    print("RWA_imag_avg",np.average(np.imag(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times))))
    print("RWA_imag_max",np.max(np.imag(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times))))
    print("RWA_imag_min",np.min(np.imag(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times))))
    plt.legend(fontsize=legend_fontsize)
    plt.xticks(fontsize=x_fontsize)
    plt.yticks(fontsize=25)
    plt.xlabel(r"$t$",fontsize=20)
    if save_name is None:
        pass
    else:
        plt.savefig("Figures/"+"dynamics/"+save_name+"_p23tilde_both",dpi=300,bbox_inches='tight')
    plt.show()
    



def steady_state_wp_plots(RWA_steady_state_data,steady_state_data,Hamiltonian="lambda",save_name=None,
                         params="TPR_A",legend_fontsize=14):
    plt.rcParams["mathtext.fontset"] = "cm"
    Omega_p = ME_params.params[params]["np"]["Omega_p"]
    Omega_c = ME_params.params[params]["np"]["Omega_c"]
    data = (steady_state_data,RWA_steady_state_data)
    ss_x_axis = np.linspace(2,10,999)
    ss_RWA = ss.get_rhos(w_p = ss_x_axis,Omega_c = Omega_c,Omega_p=Omega_p)
    x_axis = np.linspace(2,10,33)-6
    ss_x_axis -= 6

    plt.plot(x_axis,np.real(steady_state_data[:,8,0,0]),"x",label=r"$\rho_{11}$")
    plt.plot(x_axis,np.real(steady_state_data[:,0,0,0]),"x",label=r"$\rho_{22}$")
    plt.plot(x_axis,np.real(steady_state_data[:,4,0,0]),"x",label=r"$\rho_{33}$" )
    plt.plot(ss_x_axis,ss_RWA[0],"C1--")
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,8,0,0]),"C0.",label=r"$\rho_{11}$RWA")
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,0,0,0]),"C1.",label=r"$\rho_{22}$RWA")
    plt.plot(ss_x_axis,ss_RWA[1],"C2--")
    plt.plot(ss_x_axis,ss_RWA[2],"C0--")


    plt.plot(x_axis,np.real(RWA_steady_state_data[:,4,0,0]),"C2.",label=r"$\rho_{33}$RWA" )
    plt.xlabel(r"$\delta \omega_p$",fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(fontsize=legend_fontsize)
    

    if save_name is not None:
        #plt.savefig(ME.folder_name[Hamiltonian]+save_name+"steady_state_w_p_pops",dpi=300,bbox_inches='tight')
        plt.savefig("Figures/"+"w_p/"+save_name+"steady_state_w_p_pops",dpi=300,bbox_inches='tight')

        #ME.save_data("steady_state_w_p",data,Hamiltonian)
        
    plt.show()
    
    
    plt.plot(x_axis,np.real(steady_state_data[:,6,0,0]),"C3x",label=r"$\tilde{\rho}_{12}^R$")
    plt.plot(ss_x_axis,ss_RWA[5],"C3--")
    plt.plot(ss_x_axis,ss_RWA[6]*-1,"C4--")
    plt.plot(x_axis,np.imag(steady_state_data[:,6,0,0]),"C4x",label=r"$\tilde{\rho}_{12}^I$")
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,6,0,0]),"C3.",label=r"$\tilde{\rho}_{12}^R$RWA")
    plt.plot(x_axis,np.imag(RWA_steady_state_data[:,6,0,0]),"C4.",label=r"$\tilde{\rho}_{12}^I$RWA")
    plt.xlabel(r"$\delta \omega_p$",fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(fontsize=legend_fontsize)
    if save_name is not None:
        plt.savefig("Figures/"+"w_p/"+save_name+"steady_state_w_p_rho_12tilde",dpi=300,bbox_inches='tight')
    plt.show()




    plt.plot(x_axis,np.real(steady_state_data[:,7,0,0]),"C3x",label=r"$\tilde{\rho}_{13}^R$")
    plt.plot(x_axis,np.imag(steady_state_data[:,7,0,0]),"C4x",label=r"$\tilde{\rho}_{13}^I$")
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,7,0,0]),"C3.",label=r"$\tilde{\rho}_{13}^R$RWA")
    plt.plot(x_axis,np.imag(RWA_steady_state_data[:,7,0,0]),"C4.",label=r"$\tilde{\rho}_{13}^I$RWA")
    
    plt.plot(ss_x_axis,ss_RWA[7],"C3--")
    plt.plot(ss_x_axis,ss_RWA[8]*-1,"C4--")
    plt.xlabel(r"$\delta \omega_p$",fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    #plt.legend(fontsize=13,loc="lower center",bbox_to_anchor=(.55, 0))
    plt.legend(fontsize=legend_fontsize)

    if save_name is not None:
        #plt.savefig("steady_state_w_p_rho_13tilde",dpi=300,bbox_inches='tight')
        plt.savefig("Figures/"+"w_p/"+save_name+"steady_state_w_p_rho_13tilde",dpi=300,bbox_inches='tight')

    plt.show()

    plt.plot(x_axis,np.real(steady_state_data[:,1,0,0]),"C3x",label=r"$\tilde{\rho}_{23}^R$")
    plt.plot(x_axis,np.imag(steady_state_data[:,1,0,0]),"C4x",label=r"$\tilde{\rho}_{23}^I$")
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,1,0,0]),"C3.",label=r"$\tilde{\rho}_{23}^R$RWA")
    plt.plot(x_axis,np.imag(RWA_steady_state_data[:,1,0,0]),"C4.",label=r"$\tilde{\rho}_{23}^I$RWA")
    plt.plot(ss_x_axis,ss_RWA[3],"C3--")
    plt.plot(ss_x_axis,ss_RWA[4],"C4--")
    plt.xlabel(r"$\delta \omega_p$",fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(fontsize=legend_fontsize)
    if save_name is not None:
        plt.savefig("Figures/"+"w_p/"+save_name+"steady_state_w_p_rho_23tilde",dpi=300,bbox_inches='tight')
    plt.show()
    
    #---------
# size_ratios = [
#      1.3,
#      1.1,
#     .9,
#     .5, # bullshit
#     1.75,
#     .75,
#     .75,
#      .75,
#     1.75,
# ]
size_ratios = [
     1,
     1,
    1,
    1, # bullshit
    1,
    1,
    1,
     1,
    1,
]
legend_handles = []
legend_names = []
base_size = 45
graph_dicts = [{'marker':'s',
             'color': 'blue',
             'markerfacecolor':'none',
             'markersize':np.sqrt(base_size*size_ratios[0]),
             's': base_size*size_ratios[0],
             'linestyle':'-',
             'linewidth':1,
                },
               {'marker':'^',
             'color': 'red',
             'markerfacecolor':'none',
             'markersize':np.sqrt(base_size*size_ratios[1]),
             's': base_size*size_ratios[1],
             'linestyle':'-',
             'linewidth':1,
                },
               {'marker':'s',
             'color': 'C1', #orange
             'markerfacecolor':'none',
             'markersize':np.sqrt(base_size*size_ratios[2]),
             's': base_size*size_ratios[2],
             'linestyle':'-',
             'linewidth':1,
                },
                {'marker':'o', #bullshit one don't edit
             'color': 'black',
             'markerfacecolor':'none',
             'markersize':np.sqrt(base_size*size_ratios[3]),
             's': base_size*size_ratios[3],
             'linestyle':'-',
             'linewidth':1,
                },
               {'marker':'o',
             'color': 'C9', #cyan
             'markerfacecolor':'none',
             'markersize':np.sqrt(base_size*size_ratios[4]),
             's': base_size*size_ratios[4],
             'linestyle':'-',
             'linewidth':1,
                },
               {'marker':'o',
             'color': 'magenta',
             'markerfacecolor':'magenta',
             'markersize':np.sqrt(base_size*size_ratios[5]),
             's': base_size*size_ratios[5],
             'linestyle':'--',
             'linewidth':1,
                },
               {'marker':'o',
             'color': 'green',
             'markerfacecolor':'green',
             'markersize':np.sqrt(base_size*size_ratios[6]),
             's': base_size*size_ratios[6],
             'linestyle':'-.',
             'linewidth':1,
                },
               {'marker':'o',
             'color': 'yellow',
             'markerfacecolor':'none',
             'markersize':np.sqrt(base_size*size_ratios[7]),
             's': base_size*size_ratios[7],
             'linestyle':'--',
             'linewidth':1,
                },
               {'marker':'^',
             'color': 'C7',
             'markerfacecolor':'C7',
             'markersize':np.sqrt(base_size*size_ratios[8]),
             's': base_size*size_ratios[8],
             'linestyle':'-.',
             'linewidth':1,
                }
            ]
method_to_legend_string = {
             'CFME_equal_opt':'Fourth-order CFME(equal,opt)',
             'CFME_gauss_opt':'Fourth-order CFME (gauss,opt)',
             'CFME_equal':"Fourth-order CFME(equal)",
             'CFME_gauss':"Fourth-order CFME (gauss)",
             'fourth_order':"Fourth-order Magnus Expansion",
             'sixth_order':"Sixth-order Magnus Expansion",
             'RK4':"Fourth-order Runge-Kutta",
             'AM4':"Fourth-order Adams-Moulton PC",
             'AM4_alt':"Fourth-order Adams-Moulton",
             'CN2':"Second-order Crank-Nicolson",
             'M_11':r"$M_{1}^{(1)}$",
             'M_12':r"$M_{1}^{(2)}$",              
             'M_12+M_21':r"Eq. (18)",
             'M_12+M_22':r"Eq. (25)",
             'M_1':r"$M_{1}$",
             'M_1+M_21(One point)':r"$M_{1}$+one point(third order)",
             'M_1+M_21':r"$M_{1}+M_2^{(1)}$",
             'M_1+M_22':r"$M_{1}+M_2^{(2)}$",
             'M_1+M_2':r"$M_{1}+M_2$",
             'M_1_sp':r"$M_{1}$sp",
             'M_14+M_23+M_32+M_41':r"Eq. (27)",
             'Blanes_6th':'Blanes 6th-order (gauss)',
             "Iserles_4th":"Iserles 4th-order (gauss)",
             "Blanes_4th_equal":"Blanes 4th-order",
             "Blanes_4th_gauss":"Blanes 4th-order (gauss)",
             "M_12+M_22+M_31":r"Eq. (22)",
             "M_12+M_22+M_31+M_41":r"$M_1^{(2)}+M_2^{(2)}+M_3^{(1)}+M_4^{(1)}$" 

}
#--------------------------------------------------------------------------------------------
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def make_error_plot(data, comparison,stepsizes,Hamiltonian, methods,
                    indices=None,save_name=None,title=None,legend_fontsize=14.2):
    #plots only the given indices
    if(indices == None):
        indices = [i for i in range(len(data))]
    num_stepsizes = len(stepsizes)
    #for i in range(len(graph_dicts)):
    for i in range(len(methods)):
        graph_dicts[i]["legend_string"]=method_to_legend_string[methods[i]]
    legend_handles = []
    legend_names = []
    if(save_name==None):
        save_name=title
    fig,ax = plt.subplots(1,1)
    #begin inset stuff
    #x1, x2, y1, y2 = -1.93, -1.60, -8.9, -7.3 # works for case I
    #x1, x2, y1, y2 = -1.92, -1.60, -8.4, -6.7 # works for case II
    #x1, x2, y1, y2 = -1.92, -1.60, -6.8, -4.95 # works for case III
    #x1, x2, y1, y2 = -1.92, -1.60, -6.73, -4.88 # works for case IV
    #sf = 1.2
    #extent = [x1,x2,y1,y2]
   # axins = ax.inset_axes(
    #     [0, 0.71, 0.47, 0.29], #can extend a little for case III
    #     xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
    #axins.set_xticks([])
    #axins.set_yticks([])

    #axins.imshow(Z2, extent=extent, origin="lower")
    #end inset stuff
    
  
    for i in indices:
        gd = graph_dicts[i]
        error = np.zeros(num_stepsizes,dtype=np.float64)
        for j in range(num_stepsizes):
            error[j] = ME.error_matrix(comparison[0], data[i][j][-1])
            #print(comparison[0],"is comparison[0]")
            #print(data[i][j][0]," is data[i][j][0]")
            
        log_x_data = np.log10(stepsizes[:])
        log_y_data = np.log10(error)
        if(gd["odd"]==True):
            scatter = ax.scatter(log_x_data[1::2], log_y_data[1::2], color=gd['color'], marker=gd['marker'], facecolor=gd['markerfacecolor'], 
                              s=gd['s'])
            #x_scatter_idx = np.where((log_x_data[1::2]< x2) & (log_x_data[1::2] > x1))
           # y_scatter_idx = np.where((log_y_data[1::2]<y2) & (log_y_data[1::2] > y1 ))
           # scatter_idx =[index for index in x_scatter_idx if index in y_scatter_idx]
            #axins.scatter(log_x_data[x_scatter_idx],log_y_data[x_scatter_idx],color=gd['color'], marker=gd['marker'], facecolor=gd['markerfacecolor'], 
            #                  s=gd['s'])
            #axins.scatter(log_x_data[1::2], log_y_data[1::2], color=gd['color'], marker=gd['marker'], facecolor=gd['markerfacecolor'], 
            #                  s=gd['s']*sf)
#origin="lower"??
        if(gd["even"]==True):
            scatter = ax.scatter(log_x_data[::2], log_y_data[::2], color=gd['color'], marker=gd['marker'], facecolor=gd['markerfacecolor'], 
                              s=gd['s'])
            #x_scatter_idx = np.where((log_x_data[::2]< x2) & (log_x_data[::2] > x1))
           # y_scatter_idx = np.where((log_y_data[::2]<y2) & (log_y_data[::2] > y1 ))
           # scatter_idx =[index for index in x_scatter_idx if index in y_scatter_idx]
           # axins.scatter(log_x_data[::2],log_y_data[::2],color=gd['color'], marker=gd['marker'], facecolor=gd['markerfacecolor'], 
            #                  s=gd['s']*sf)
        print(log_x_data,log_y_data)
        m,b = np.polyfit(log_x_data, log_y_data, 1)
        print("slope,y-intercept is: ",m,b,methods[i])
        plt.plot(log_x_data, b + m * log_x_data, color=gd['color'], linestyle=gd['linestyle'], linewidth=gd['linewidth'])
        #axins.plot(log_x_data, b + m * log_x_data, color=gd['color'], linestyle=gd['linestyle'], linewidth=gd['linewidth']*sf)
    
        #legend stuff
        legend_handle = Line2D([0], [0], marker=gd['marker'], color=gd['color'], markerfacecolor=gd['markerfacecolor'], markersize=gd['markersize'],
                               linestyle=gd['linestyle'], linewidth=gd['linewidth'])
        legend_handles.append(legend_handle)
        legend_names.append(gd["legend_string"])
        
    plt.xlabel(r"$\log_{10}$(stepsize/$t_c$)",fontsize=25)
    plt.ylabel(r"$\log_{10}$(error)",fontsize=25)
    #plt.xlim(right=-1.0) for II
    #plt.ylim(bottom=-30) for II

    if(title!=None):
        plt.title(title,fontsize=30)
    plt.legend(legend_handles,legend_names,fontsize=legend_fontsize,ncols=2)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    #inset_rect = Rectangle((x1, y1), x2 - x1, y2 - y1, linewidth=2, edgecolor='black', facecolor='none')
    #ax.add_patch(inset_rect)
    ax.spines['top'].set_linewidth(2)     # Set the top spine linewidth
    ax.spines['bottom'].set_linewidth(2)  # Set the bottom spine linewidth
    ax.spines['left'].set_linewidth(2)    # Set the left spine linewidth
    ax.spines['right'].set_linewidth(2)
    #axins.spines['top'].set_linewidth(2)     # Set the top spine linewidth
    #axins.spines['bottom'].set_linewidth(2)  # Set the bottom spine linewidth
    #axins.spines['left'].set_linewidth(2)    # Set the left spine linewidth
    #axins.spines['right'].set_linewidth(2)
    fig.set_size_inches(12.9, 4.8)

    #ax.indicate_inset_zoom(axins, edgecolor="black") #for the original box, and the lines denoting smaller box
    if(save_name!=None):
        print("Figures/"+"stepsize_errors/"+save_name+"_errors")
        #plt.savefig("Figures/"+"stepsize_errors/"+"testing")
        plt.savefig("Figures/"+"stepsize_errors/"+save_name+"_errors.png",dpi=300,bbox_inches='tight')
        #plt.savefig(ME.folder_name[Hamiltonian]+save_name,dpi=300,bbox_inches='tight')
        #ME.save_data(save_name,data,Hamiltonian,comparison=comparison)
        #save_data(filename,data,Hamiltonian,comparison=None)
    # inset axes....
   


    #plt.show() remove comment to see figure
#def make_plot_error(stepsizes, y_vals, comparison_y_vals):
    #create all legend handles here
#THE MOST RECENT CODE HERE ASSUMES DATA IS IN THE ORDER OF THE PAPER.
#See the example code for plotting: make_plot(data,comparison,[0,1,2,4,7,5,6,8],Hamiltonian,"test")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
# size_ratios = [
#     .079,
#     .080,
#     .087,
#     .161, # bullshit
#     .161,
#     .080,
#     .080,
#     .086,
#     .101
# ]
size_ratios = [
     1.3,
     1.1,
    .9,
    .5, # bullshit
    1.75,
    .75,
    .75,
     .75,
    1.75,
]


# plot_type = [
#         {'marker': 'o', 'linestyle': 'solid', 'markerfacecolor': 'none','color':'green'},
#         {'marker': 'o', 'linestyle': 'dotted', 'markerfacecolor': 'none','color':'red'},
#         {'marker': 'o', 'linestyle': 'dashed', 'markerfacecolor': 'none','color':'black'},
#         {'marker': 'o', 'linestyle': 'dashdot', 'markerfacecolor': 'none','color':'red'},
#         {'marker': '^', 'linestyle': 'solid','color':'magenta'},
#         {'marker': '^', 'linestyle': 'dotted','color':'red'},
#         {'marker': '^', 'linestyle': 'dashed','color':'blue'},
#         {'marker': '.', 'linestyle': 'dashdot','color':'yellow'},
#         {'marker': '*', 'linestyle': 'dashdot','color':'blue'},
#     ]


#2,8
    #odd = [0,7,2,4,6,8]     #works for case I
    #even = [0,1,2,3,4,5,8]  #works for case I
    #odd = [0,5,4,6,8,7]   #works for case II
    #even = [0,1,2,3,4,8,7] #works for case II
    #odd = [0,2,5,4,6,8,7] #works for case III
    #even = [0,1,2,3,4,5,6,8] #works for case III
#odd = [0,2,5,4,6,8,7] #works for case IV
#even = [0,1,2,3,4,5,6,8] #works for case IV
odd = [0,1,2,3,4,5,6,7,8]
even = [0,1,2,3,4,5,6,7,8]

for i in range(len(graph_dicts)):
    if (i in odd):
        graph_dicts[i]["odd"]=True
      #  print(i)
    else:
        graph_dicts[i]["odd"]=False
    if (i in even):
        graph_dicts[i]["even"]=True
      #  print(i)
    else:
        graph_dicts[i]["even"]=False
       
def make_error_plot_both(times,data,data_RWA,save_name=None,title=None,print_avg=True,x_fontsize=25,
                      initial_condition=8):
    #data should only have data from one method
    errors=np.zeros(len(times))
    for j in range(len(times)):
        errors[j] = ME.error_matrix(data[j,:,initial_condition].reshape(3,3),data_RWA[j,:,initial_condition].reshape(3,3))
    plt.plot(times,errors,label=r"$\rho_{11}$")
    print(np.max(np.abs(np.imag(data[:,8,initial_condition]))))
    print(np.max(np.abs(np.imag(data[:,0,initial_condition]))))
    print(np.max(np.abs(np.imag(data[:,4,initial_condition]))))


#     plt.plot(times,data[:,0,initial_condition],label=r"$\rho_{22}$")
#     plt.plot(times,data[:,4,initial_condition],label=r"$\rho_{33}$")
#     plt.plot(times,data_RWA[:,8,initial_condition],"--",linewidth=2,label=r"$\rho_{11}$RWA")
#     plt.plot(times,data_RWA[:,0,initial_condition],"--",linewidth=2,label=r"$\rho_{22}$RWA")
#     plt.plot(times,data_RWA[:,4,initial_condition],"--",linewidth=2,label=r"$\rho_{33}$RWA")
#     print(np.max(np.abs(np.imag(data_RWA[:,8,initial_condition]))))
#     print(np.max(np.abs(np.imag(data_RWA[:,0,initial_condition]))))
#     print(np.max(np.abs(np.imag(data_RWA[:,4,initial_condition]))))
    
    plt.xlabel(r"$t$",fontsize=20)
    plt.ylabel(r"Error",fontsize=20)
    # if(print_avg==True):
    #     print(np.average(data[:,8,initial_condition]))
    #     print(np.average(data[:,0,initial_condition]))
    #     print(np.average(data[:,4,initial_condition]))
    #     print("RWA_pop_1",np.average(np.real(data_RWA[:,8,initial_condition])))
    #     print("RWA_pop_2",np.average(np.real(data_RWA[:,0,initial_condition])))
    #     print("RWA_pop_3",np.average(np.real(data_RWA[:,4,initial_condition])))


              
   # plt.legend(fontsize=14)
    plt.xticks(fontsize=x_fontsize)
    plt.yticks(fontsize=25)
    if save_name is None:
        pass
    else:
        #plt.savefig(ME.folder_name[Hamiltonian]+save_name+"_error_both",dpi=300,bbox_inches='tight')
        plt.savefig("Figures/errors/"+save_name+"_error_both",dpi=300,bbox_inches='tight')
        #ME.save_data(save_name,data,Hamiltonian,comparison=comparison)
    plt.show()    
    
        
    


