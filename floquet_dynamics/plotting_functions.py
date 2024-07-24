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
