def get_rhos(w_c=4,w_1=6,w_2=2,
                  w_p =6,w_3 =0,
                  Omega_c=0.56, Omega_p=0.50,
                  gamma_12=0.9, gamma_13=1,
                 n_12 = 0,n_13 = 0
                ):
    if(n_12 != 0 or n_13 != 0):
        print("Warning; the results here are only good in the limit n_12 = 0 and n_13 = 0")
    gamma = gamma_12/gamma_13
    #print(gamma)
    dw_p = w_p-(w_1-w_3)
    
    D = (4*Omega_c**2*dw_p**4+
    (gamma*(gamma+1)**2*Omega_p**2+(gamma+1)*(gamma+1+8*Omega_p**2)*Omega_c**2-8*Omega_c**4)*dw_p**2+
    4*(Omega_c**2+Omega_p**2)**2 * (Omega_c**2+gamma*Omega_p**2)
        )
    p11_tilde = 4*(gamma+1)*Omega_c**2*Omega_p**2*dw_p**2/D
    p22_tilde = Omega_p**2*(gamma*((gamma+1)**2+4*Omega_c**2)*dw_p**2+4*(Omega_c**2+Omega_p**2)*(Omega_c**2+gamma*Omega_p**2))/D
    p33_tilde = Omega_c**2*(4*dw_p**4+((gamma+1)**2-8*Omega_c**2+4*Omega_p**2)*dw_p**2+4*(Omega_c**2+Omega_p**2)*(Omega_c**2+gamma*Omega_p**2))/D
    p12_R = -4*Omega_c*Omega_p**2*(Omega_c**2+gamma*Omega_p**2)*dw_p/D
    p12_I = 2*gamma*(gamma+1)*Omega_c*Omega_p**2*dw_p**2/D
    p13_R = 4*Omega_c**2*Omega_p*(Omega_c**2+gamma*Omega_p**2-dw_p**2)*dw_p/D
    p13_I = 2*(gamma+1)*Omega_c**2*Omega_p*dw_p**2/D
    p23_R = 4*Omega_c*Omega_p*(Omega_c**2*dw_p**2 - (Omega_c**2+Omega_p**2)*(Omega_c**2+gamma*Omega_p**2))/D
    p23_I = -2*(gamma+1)*(Omega_c**2+gamma*Omega_p**2)*Omega_c*Omega_p*dw_p/D
    return [p11_tilde,p22_tilde,p33_tilde,p12_R,p12_I,
            p13_R,p13_I,p23_R,p23_I
           ]