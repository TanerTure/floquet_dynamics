import steady_state as ss
import Magnus_Expansion as ME
import ME_params
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
plt.rcParams["mathtext.fontset"] = "cm"
#plt.locator_params(nbins = 6)

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
                       legend_fontsize=16,initial_condition=8, loc = "best", bbox_to_anchor = None, extend=0):
    fig, ax = plt.subplots(1,1, constrained_layout=True)
    #data should only have data from one method
    #plt.locator_params(nbins = 6)
    ax.plot(times,data[:,8,initial_condition],label=r"$\rho_{11}$")
    print(np.max(np.abs(np.imag(data[:,8,initial_condition]))))
    print(np.max(np.abs(np.imag(data[:,0,initial_condition]))))
    print(np.max(np.abs(np.imag(data[:,4,initial_condition]))))


    ax.plot(times,data[:,0,initial_condition],label=r"$\rho_{22}$")
    ax.plot(times,data[:,4,initial_condition],label=r"$\rho_{33}$")
    ax.plot(times,data_RWA[:,8,initial_condition],"C0--",linewidth=2)
    ax.plot(times,data_RWA[:,0,initial_condition],"C1--",linewidth=2)
    ax.plot(times,data_RWA[:,4,initial_condition],"C2--",linewidth=2)
    print(np.max(np.abs(np.imag(data_RWA[:,8,initial_condition]))))
    print(np.max(np.abs(np.imag(data_RWA[:,0,initial_condition]))))
    print(np.max(np.abs(np.imag(data_RWA[:,4,initial_condition]))))
    
    ax.set_xlabel(r"$t$",fontsize=25)
    if(print_avg==True):
        print(np.average(data[:,8,initial_condition]))
        print(np.average(data[:,0,initial_condition]))
        print(np.average(data[:,4,initial_condition]))
        print("RWA_pop_1",np.average(np.real(data_RWA[:,8,initial_condition])))
        print("RWA_pop_2",np.average(np.real(data_RWA[:,0,initial_condition])))
        print("RWA_pop_3",np.average(np.real(data_RWA[:,4,initial_condition])))
    if extend != 0:
        bottom_limit, upper_limit = ax.get_ylim()
        ax.set_ylim(bottom_limit, upper_limit + (upper_limit - bottom_limit)*extend/100)
    ax.legend(fontsize=legend_fontsize, loc=loc, bbox_to_anchor= bbox_to_anchor)
    # if legend_flag == 1:
    #     ax.legend(fontsize=legend_fontsize, loc="upper left", bbox_to_anchor = (.21, 1.02))
    # else:
    #     ax.legend(fontsize=legend_fontsize)
    ax.tick_params(axis='x', which="major",labelsize=x_fontsize)
    ax.tick_params(axis='y', which="major",labelsize=25)
    #plt.yticks(fontsize=25)
    if save_name is None:
        pass
    else:
        #fig.savefig("test")
        fig.savefig("Figures/"+"dynamics/"+save_name+"_pops_both",dpi=300,bbox_inches='tight')
        #ME.save_data(save_name,data,Hamiltonian,comparison=comparison)
    fig.show()

def make_pop_plot_error(times,data_fourth,data_sixth,Hamiltonian,save_name=None,title=None,print_avg=True,x_fontsize=25,
                       legend_fontsize=18,initial_condition=8):
    #data should only have data from one method
    #plt.locator_params(nbins = 6)
    data = data_sixth-data_fourth
    plt.plot(times,data[:,8,initial_condition],label=r"$\rho_{11}$")
    print(np.max(np.abs(np.imag(data[:,8,initial_condition]))))
    print(np.max(np.abs(np.imag(data[:,0,initial_condition]))))
    print(np.max(np.abs(np.imag(data[:,4,initial_condition]))))


    plt.plot(times,data[:,0,initial_condition],label=r"$\rho_{22}$")
    plt.plot(times,data[:,4,initial_condition],label=r"$\rho_{33}$")
    # plt.plot(times,data_RWA[:,8,initial_condition],"C0--",linewidth=2,label=r"$\rho_{11}$RWA")
    # plt.plot(times,data_RWA[:,0,initial_condition],"C1--",linewidth=2,label=r"$\rho_{22}$RWA")
    # plt.plot(times,data_RWA[:,4,initial_condition],"C2--",linewidth=2,label=r"$\rho_{33}$RWA")
    # print(np.max(np.abs(np.imag(data_RWA[:,8,initial_condition]))))
    # print(np.max(np.abs(np.imag(data_RWA[:,0,initial_condition]))))
    # print(np.max(np.abs(np.imag(data_RWA[:,4,initial_condition]))))
    
    plt.xlabel(r"$t$",fontsize=20)
    if(print_avg==True):
        print(np.average(data[:,8,initial_condition]))
        print(np.average(data[:,0,initial_condition]))
        print(np.average(data[:,4,initial_condition]))
        # print("RWA_pop_1",np.average(np.real(data_RWA[:,8,initial_condition])))
        # print("RWA_pop_2",np.average(np.real(data_RWA[:,0,initial_condition])))
        # print("RWA_pop_3",np.average(np.real(data_RWA[:,4,initial_condition])))


              
    plt.legend(fontsize=legend_fontsize)
    plt.xticks(fontsize=x_fontsize)
    plt.yticks(fontsize=25)
    if save_name is None:
        pass
    else:
        plt.savefig("Figures/"+"conv/"+save_name+"_pops",dpi=300,bbox_inches='tight')
        #ME.save_data(save_name,data,Hamiltonian,comparison=comparison)
    plt.show()

def make_off_diag_plots_both(times,data,data_RWA,Hamiltonian,save_name=None,title=None,print_avg=True,x_fontsize=25,
                            legend_fontsize=18,initial_condition=8):
    #plt.locator_params(nbins = 6)
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
                                  legend_fontsize=18,initial_condition=8, loc="best", bbox_to_anchor=None, extend = [0,0,0], ncol=[1,1,1]):
    fig, ax = plt.subplots(1,1, constrained_layout=True)
    plt.locator_params(nbins = 6)
    ax.plot(times,np.real(data[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C3",label=r"$\tilde{\rho}_{12}^R$")
    ax.plot(times,np.imag(data[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C4",label=r"$\tilde{\rho}_{12}^I$")
    #plt.plot(times,np.real(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C3--",linewidth=2,label=r"$\tilde{\rho}_{12}^R$RWA")
    ax.plot(times,np.real(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C3--",linewidth=2)
    #plt.plot(times,np.imag(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C4--",linewidth=2,label=r"$\tilde{\rho}_{12}^I$RWA")
    ax.plot(times,np.imag(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times)),"C4--",linewidth=2)

    if(print_avg==True):
        print(np.average(np.real(data[:,6,initial_condition]*np.exp(-1j*w_p*times))))
        print(np.average(np.imag(data[:,6,initial_condition]*np.exp(-1j*w_p*times))))
        print("RWA_real_avg",np.average(np.real(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times))))
        print("RWA_imag_avg",np.average(np.imag(data_RWA[:,6,initial_condition]*np.exp(-1j*w_p*times))))
    ax.legend(fontsize=legend_fontsize[0], loc=loc[0], bbox_to_anchor=bbox_to_anchor[0],ncol=ncol[0])

    #ax.xticks(fontsize=x_fontsize)
    #ax.yticks(fontsize=25)
    ax.set_xlabel(r"$t$",fontsize=25)
    ax.tick_params(axis='x', which="major", labelsize=x_fontsize)
    ax.tick_params(axis='y', which="major", labelsize=25)
    if extend[0] != 0:
        bottom_limit, upper_limit = ax.get_ylim()
        ax.set_ylim(bottom_limit, upper_limit + (upper_limit - bottom_limit)*extend[0]/100)
    if save_name is None:
        pass
    else:
        fig.savefig("Figures/"+"dynamics/"+save_name+"_p12tilde_both",dpi=300,bbox_inches='tight')
    #plt.show()
    fig, ax = plt.subplots(1,1, constrained_layout=True)
    ax.plot(times,np.real(data[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times)),"C3",label=r"$\tilde{\rho}_{13}^R$")
    ax.plot(times,np.imag(data[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times)),"C4",label=r"$\tilde{\rho}_{13}^I$")
    ax.plot(times,np.real(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times)),"C3--",linewidth=2)#label=r"$\tilde{\rho}_{13}^R$RWA"
    ax.plot(times,np.imag(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times)),"C4--",linewidth=2) #,label=r"$\tilde{\rho}_{13}^I$RWA")
    print(np.average(np.real(data[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print(np.average(np.imag(data[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_real_avg",np.average(np.real(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_imag_avg",np.average(np.imag(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_real_max",np.max(np.real(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_imag_max",np.max(np.imag(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_real_min",np.min(np.real(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    print("RWA_imag_min",np.min(np.imag(data_RWA[:,7,initial_condition]*np.exp(-1j*(w_p-w_c)*times))))
    

    
    ax.legend(fontsize=legend_fontsize[1], loc=loc[1], bbox_to_anchor=bbox_to_anchor[1],ncol=ncol[1])
    ax.tick_params(axis="x", which="major", labelsize=x_fontsize)
    ax.tick_params(axis="y", which="major", labelsize=25)
    ax.set_xlabel(r"$t$",fontsize=20)
    if extend[1] != 0:
        bottom_limit, upper_limit = ax.get_ylim()
        ax.set_ylim(bottom_limit, upper_limit + (upper_limit - bottom_limit)*extend[1]/100)
        
    if save_name is None:
        pass 
    else:
        fig.savefig("Figures/"+"dynamics/"+save_name+"_p13tilde_both",dpi=300,bbox_inches='tight')
    #plt.show()
    fig, ax = plt.subplots(1,1)
    ax.plot(times,np.real(data[:,1,initial_condition]*np.exp(1j*w_c*times)),"C3",label=r"$\tilde{\rho}_{23}^R$")
    ax.plot(times,np.imag(data[:,1,initial_condition]*np.exp(1j*w_c*times)),"C4",label=r"$\tilde{\rho}_{23}^I$")
    ax.plot(times,np.real(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times)),"C3--",linewidth=2)#,label=r"$\tilde{\rho}_{23}^R$RWA")
    ax.plot(times,np.imag(data_RWA[:,1,initial_condition]*np.exp(1j*w_c*times)),"C4--",linewidth=2)#,label=r"$\tilde{\rho}_{23}^I$RWA")
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
    ax.legend(fontsize=legend_fontsize[2], loc=loc[2], bbox_to_anchor=bbox_to_anchor[2],ncol=ncol[2])
    ax.tick_params(axis="x", which="major", labelsize=x_fontsize)
    ax.tick_params(axis="y", which="major", labelsize=25)
    #plt.yticks(fontsize=25)
    ax.set_xlabel(r"$t$",fontsize=20)
    if extend[2] != 0:
        bottom_limit, upper_limit = ax.get_ylim()
        ax.set_ylim(bottom_limit, upper_limit + (upper_limit - bottom_limit)*extend[2]/100)
    if save_name is None:
        pass
    else:
        fig.savefig("Figures/"+"dynamics/"+save_name+"_p23tilde_both",dpi=300,bbox_inches='tight')

    #plt.show()    


def steady_state_wp_plots(RWA_steady_state_data,steady_state_data,Hamiltonian="lambda",save_name=None,
                         params="TPR_A",legend_fontsize=14, ncol=[1,1,1,1]):
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
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,8,0,0]),"C0.")
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,0,0,0]),"C1.")
    plt.plot(ss_x_axis,ss_RWA[1],"C2--")
    plt.plot(ss_x_axis,ss_RWA[2],"C0--")


    plt.plot(x_axis,np.real(RWA_steady_state_data[:,4,0,0]),"C2.")
    plt.xlabel(r"$\delta \omega_p$",fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(fontsize=legend_fontsize,ncol=ncol[0])
    

    if save_name is not None:
        #plt.savefig(ME.folder_name[Hamiltonian]+save_name+"steady_state_w_p_pops",dpi=300,bbox_inches='tight')
        plt.savefig("Figures/"+"w_p/"+save_name+"steady_state_w_p_pops",dpi=300,bbox_inches='tight')

        #ME.save_data("steady_state_w_p",data,Hamiltonian)
        
    plt.show()
    
    
    plt.plot(x_axis,np.real(steady_state_data[:,6,0,0]),"C3x",label=r"$\tilde{\rho}_{12}^R$")
    plt.plot(ss_x_axis,ss_RWA[5],"C3--")
    plt.plot(ss_x_axis,ss_RWA[6]*-1,"C4--")
    plt.plot(x_axis,np.imag(steady_state_data[:,6,0,0]),"C4x",label=r"$\tilde{\rho}_{12}^I$")
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,6,0,0]),"C3.")
    plt.plot(x_axis,np.imag(RWA_steady_state_data[:,6,0,0]),"C4.")
    plt.xlabel(r"$\delta \omega_p$",fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(fontsize=legend_fontsize, ncol=ncol[1])
    if save_name is not None:
        plt.savefig("Figures/"+"w_p/"+save_name+"steady_state_w_p_rho_12tilde",dpi=300,bbox_inches='tight')
    plt.show()




    plt.plot(x_axis,np.real(steady_state_data[:,7,0,0]),"C3x",label=r"$\tilde{\rho}_{13}^R$")
    plt.plot(x_axis,np.imag(steady_state_data[:,7,0,0]),"C4x",label=r"$\tilde{\rho}_{13}^I$")
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,7,0,0]),"C3.")
    plt.plot(x_axis,np.imag(RWA_steady_state_data[:,7,0,0]),"C4.")
    
    plt.plot(ss_x_axis,ss_RWA[7],"C3--")
    plt.plot(ss_x_axis,ss_RWA[8]*-1,"C4--")
    plt.xlabel(r"$\delta \omega_p$",fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    #plt.legend(fontsize=13,loc="lower center",bbox_to_anchor=(.55, 0))
    plt.legend(fontsize=legend_fontsize, ncol=ncol[2])

    if save_name is not None:
        #plt.savefig("steady_state_w_p_rho_13tilde",dpi=300,bbox_inches='tight')
        plt.savefig("Figures/"+"w_p/"+save_name+"steady_state_w_p_rho_13tilde",dpi=300,bbox_inches='tight')

    plt.show()

    plt.plot(x_axis,np.real(steady_state_data[:,1,0,0]),"C3x",label=r"$\tilde{\rho}_{23}^R$")
    plt.plot(x_axis,np.imag(steady_state_data[:,1,0,0]),"C4x",label=r"$\tilde{\rho}_{23}^I$")
    plt.plot(x_axis,np.real(RWA_steady_state_data[:,1,0,0]),"C3.")
    plt.plot(x_axis,np.imag(RWA_steady_state_data[:,1,0,0]),"C4.")
    plt.plot(ss_x_axis,ss_RWA[3],"C3--")
    plt.plot(ss_x_axis,ss_RWA[4],"C4--")
    plt.xlabel(r"$\delta \omega_p$",fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(fontsize=legend_fontsize,ncol=ncol[3])
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

       
def make_error_plot_both(times,data,data_RWA,save_name=None,title=None,print_avg=True,x_fontsize=25,
                      initial_condition=8,text="I"):
    fig, ax = plt.subplots(3,1, constrained_layout=True, sharex=True, figsize=(5,9))
    ax[2].set_xlabel(r"$t$",fontsize=25)
    #data should only have data from one method
    letters = ["A","B","C"]
    for i in range(3):
        errors=np.zeros(len(times[i]))
        for j in range(len(times[i])):
            errors[j] = ME.error_matrix(data[i][j,:,initial_condition].reshape(3,3),data_RWA[i][j,:,initial_condition].reshape(3,3))
        ax[i].plot(times[i],errors,"b")
        ax[i].text(-.25, 1.08, "("+ letters[i] +"-" + text + ")", transform = ax[i].transAxes, fontsize=25)
        print(np.max(np.abs(np.imag(data[i][:,8,initial_condition]))))
        print(np.max(np.abs(np.imag(data[i][:,0,initial_condition]))))
        print(np.max(np.abs(np.imag(data[i][:,4,initial_condition]))))

        ax[i].set_ylabel(r"Error",fontsize=20)

        ax[i].tick_params(axis="x", which="major", labelsize=27)
        ax[i].tick_params(axis="y", which="major", labelsize=25)
        #plt.xticks(fontsize=x_fontsize)

    if save_name is None:
        pass
    else:
        fig.savefig("Figures/errors/"+save_name+"_error_both",dpi=300,bbox_inches='tight')
    plt.show()
    
    
def loglog_lob_2(x_vals_plots, y_vals_plots, names, colors = ['r','g', 'b', 'm'], markers = ['s','s','^','o'],
                 filled = ['none','g','none', 'none'], file_name="test"):
    '''
        Makes log log graph of data for two graphs which share the same x-axes. Draws the line of best fit for both graphs.
        Saves the figure using the keyword argument file_name, and saves the data points into separate .txt files  
        
        params:
        x_vals_plots:list of lists containing x_values, with dimensions [2, num_lines, num_points]
        names: list containing names to be used in legend of dimension [num_lines]
        colors = list containing colors to be used for points and lines of graph, of dimension [num_lines]
        markers = list containing shapes to be used for points on the graph, of dimension [num_lines]
        
        returns:
        fig, ax of graph
    '''
    from matplotlib.ticker import LogLocator
    plt.rcParams['mathtext.fontset'] = 'cm'
    fig, ax = plt.subplots(2, 1, figsize = (8, 10), sharex= True, constrained_layout=True)
    for i in range(2):
        if i == 0:
            letter = '(a)'
        else:
            letter = '(b)'
        ax[i].text(-.1, 1, letter ,transform=ax[i].transAxes, fontsize=20)
        ax[i].set_xscale('log')
        ax[i].set_yscale('log')
        ax[i].set_ylabel("Error", fontsize=25)
        ax[i].tick_params(axis="both", labelsize=25, size=11)
        ax[i].tick_params(axis="both",which="minor", labelsize=17, size = 4)
        ax[i].minorticks_on()
        #ax[i].xaxis.set_major_locator(LogLocator(base=10.0, subs=[1,2,4,8]))
        ax[i].xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2,10)*0.1))
        y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = numpy.arange(1.0, 10.0) * 0.1, numticks = 10)
        
        for j in range(len(names)):
            x_data, y_data  = np.log10(x_vals_plots[i][j]), np.log10(y_vals_plots[i][j])
            m,b = np.polyfit(x_data, y_data, 1)
            ax[i].plot(10**x_data, 10**(b + m * x_data), color = colors[j], linewidth = 1)
            ax[i].plot(x_vals_plots[i][j], y_vals_plots[i][j], color=colors[j], marker=markers[j], markersize=7, markerfacecolor=filled[j], linestyle='none' )
        for key in ax[i].spines:
            ax[i].spines[key].set_linewidth(2)
        legend_names = []
        legend_handles = []
        if i == 1:
            ax[1].set_xlabel(r"$\delta t           $", fontsize=25)
        if i == 0:
            for j in range(len(names)):
                legend_names.append(names[j])
                legend_handles.append(matplotlib.lines.Line2D([0],[0],color=colors[j], marker = markers[j], markersize=7, markerfacecolor=filled[j]))
            ax[0].legend(legend_handles, legend_names, fontsize=15)
    #fig.tight_layout()
    fig.savefig(f"Figures/stepsize_errors_both/{file_name}", dpi=300, bbox_inches = "tight")
    return fig, ax

def make_text_files_stepsizes(x_vals_plots, y_vals_plots, names, file_name):
    with open(f"Figures/Figures_txt/{file_name}" + "_short_time", "w") as file:
        for name in names:
            file.write(f"stepsize {name}_error ")
        file.write("\n")
        for i in range(len(x_vals_plots[0][0])):
            for j in range(len(x_vals_plots[0])):
                file.write(str(x_vals_plots[0][j][i]) + " " + str(y_vals_plots[0][j][i]) + " ")
            file.write("\n")
    with open(f"Figures/Figures_txt/{file_name}" + "_long_time","w") as file:
        for name in names:
            file.write(f"stepsize {name}_error")
        file.write("\n")
        for i in range(len(x_vals_plots[1][0])):
            for j in range(len(x_vals_plots[1])):
                file.write(str(x_vals_plots[1][j][i]) + " " + str(y_vals_plots[1][j][i]) + " ")
            file.write("\n")
    return
    
        
    


