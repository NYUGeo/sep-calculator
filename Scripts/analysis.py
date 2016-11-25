import numpy as np

#        omega = retaining wall inclination
#         beta = backfill of surface slope
#        theta = angle to the vertical
#        alpha = angle normal to wall
#  alpha_alpha = obliquity of the earth thrust
#        gamma = unit weight (kN/m3)
#          phi = internal friction angle
#            c = cohesion (kPa)
#            g = gravitational acceleration (9.80665 m/s2)
#           kh = horizontal seismic coefficient
#           kv = vertical seismic coefficient
#      K_alpha = active lateral earth pressure coefficient
#      alpha_h = horizontal pseudo-static acceleration
#      alpha_v = vertical pseudo-static acceleration
#      g_theta = modified gravitational acceleration
#  gamma_theta = modified unit weight
#            z = depth below sloping ground surface
#      z_theta = modified axis along theta
#  sigma_theta = modified effective stress
#  sigma_alpha = stress acting on the wall
#    sigma_AEH =
#      J_alpha = active condition


def figure5(kh, kv, omega_deg, beta_deg, phi_deg, alpha_deg, c, gamma, H, z_w):
    ### Gravitational & pseudo-static accelerations
    g = 9.80665
    alpha_h = kh * g
    alpha_v = kv * g

    ### Convert degrees to radians
    omega = np.radians(omega_deg)
    beta = np.radians(beta_deg)
    phi = np.radians(phi_deg)
    alpha = np.radians(alpha_deg)

    ### Heights
    # z_w = H ???
    z_l = z_w/np.cos(omega)
    z = z_w * (np.cos(beta-omega))/(np.cos(beta)*np.cos(omega))

    ### Equation 1
    g_theta = np.sqrt((g + alpha_v)**2 + alpha_h**2)

    ### Equation 2
    # theta = np.arctan(alpha_h / (g + alpha_v))
    # or
    theta = np.deg2rad(np.arctan(kh/(1 + kv)))

    ### Equation 3
    gamma_theta = (gamma * (1 + kv))/(np.cos(theta))

    ### Equation 4
    sigma_theta = gamma * z * (1 + kv) * (np.cos(beta)/np.cos(theta))

    ### Equation 13
    J_alpha = (1/(np.cos(phi)**2)) * (
                (
                 ((gamma*z*np.cos(beta)*np.cos(beta+theta)*(1+kv))/np.cos(theta))
                 + c * np.cos(phi) * np.sin(phi)
                )
                - np.sqrt(
                 (gamma**2 * z**2 * np.cos(beta)**2 * (1+kv)**2 *
                    ((np.cos(beta+theta)**2 - np.cos(phi)**2)/np.cos(theta)**2))
                 + c**2 + np.cos(phi)**2
                 + ((2*c*gamma*z*np.cos(phi)*np.sin(phi)*np.cos(beta)*np.cos(beta+theta)*(1+kv))/np.cos(theta))
                )
    )

    ### Equation 16
    K_alpha =   (
                (np.cos(beta)*(1+kv)*(np.sin(theta+omega)**2 - np.cos(beta-omega)**2))/(np.cos(alpha)*np.cos(beta+theta)*np.cos(theta))
                ) + (
                ((2*(J_alpha/(gamma*z))*np.cos(beta-omega)**2)/np.cos(alpha))
    )

    ### Equation 15
    sigma_alpha = gamma * z * K_alpha

    ### Equation 18
    alpha_alpha = np.arctan(
                (
                 ((2*np.cos(theta)*np.cos(beta+theta)*J_alpha)/(np.cos(beta)*(1+kv)*gamma*z)-1)
                 * np.sin(beta-omega)**2 + np.sin(theta+omega)**2
                )/2*(
                 ((2*np.cos(theta)*np.cos(beta+theta)*J_alpha)/(np.cos(beta)*(1+kv)*gamma*z)-1)
                 * np.cos(beta-omega)**2 + np.sin(theta+omega)**2
                )
    )

    alpha_alpha_deg = np.rad2deg(alpha_alpha)

    ### Equation 19
    sigma_AEH = sigma_alpha * np.cos(alpha_alpha + omega)

    ### Results
    #return z_w, z_l, z, J_alpha, alpha_alpha, K_alpha, sigma_alpha, sigma_AEH
    print (round(z_w,4),round(z_l,4),round(z,4),round(J_alpha,2),round(alpha_alpha_deg,2),\
        round(K_alpha,2),round(sigma_alpha,2),round(sigma_AEH,2))

### Test with:
### figure5(0.2, 0.1, 20, 15, 30, 5, 20, 23, 15, 0.0001)
