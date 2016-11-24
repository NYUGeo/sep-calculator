import numpy as np

#        omega = retaining wall inclination
#         beta = backfill of surface slope
#        theta = angle to the vertical
#        alpha = angle normal to wall
#  alpha_alpha = obliquity of the earth thrust
#        gamma = unit weight
#          phi = internal friction angle
#            c = cohesion
#            g = gravitational acceleration (9.80665 m/s2)
#           kh = horizontal seismic coefficient
#           kv = vertical seismic coefficient
#      K_alpha = active lateral earth pressure coefficient
#      alpha_h = horizontal pseudo-static acceleration
#      alpha_v = vertical pseudo-static acceleration
#      g_theta = modified gravitational acceleration
#  gamma_theta = modified unit weight
#         zeta = depth below sloping ground surface
#   zeta_theta = modified axis along theta
#  sigma_theta = modified effective stress
#  sigma_alpha = stress acting on the wall
#    sigma_AEH =
#      J_alpha = active condition


def calculate(kh, kv, ):
    ### Gravitational & pseudo-static accelerations
    g = 9.80665
    alpha_h = kh * g
    alpha_v = kv * g
    ### Equation 1
    g_theta = np.sqrt((g + alpha_v)**2 + alpha_h**2)

    ### Equation 2
    # theta = np.arctan(alpha_h / (g + alpha_v))
    # or
    theta = np.arctan(kh/(1 + kv))

    ### Equation 3
    gamma_theta = (gamma * (1 + kv))/(np.cos(theta))

    ### Equation 4
    sigma_theta = gamma * zeta * (1 + kv) * (np.cos(beta)/np.cos(theta))

    ### Equation 13
    J_alpha = (1/np.cos(phi)**2) * (
                (
                 ((gamma*zeta*np.cos(beta)*np.cos(beta+theta)(1+kv))/np.cos(theta))
                 + c * np.cos(phi) * np.sin(phi)
                )
                - np.sqrt(
                 (gamma**2 * zeta**2 * np.cos(beta)**2 * (1+kv)**2 *
                    ((np.cos(beta+theta)**2 - np.cos(phi)**2)/np.cos(theta)**2))
                 + c**2 + np.cos(phi)**2
                 + ((2*c*gamma*zeta*np.cos(phi)*np.sin(phi)*np.cos(beta)*np.cos(beta+theta)*(1+kv))/np.cos(theta))
                )
    )

    ### Equation 16
    K_alpha =   (
                (np.cos(beta)*(1+kv)*(np.sin(theta+omega)**2 - np.cos(beta-omega)**2))/(np.cos(alpha)*np.cos(beta+theta)*np.cos(theta))
                ) + (
                ((2*(J_alpha/(gamma*zeta))*np.cos(beta-omega)**2)/np.cos(alpha))
    )

    ### Equation 15
    sigma_alpha = gamma * zeta * K_alpha

    ### Equation 18
    alpha_alpha = np.arctan(
                (
                 ((2*np.cos(theta)*np.cos(beta+theta)*J_alpha)/(np.cos(beta)*(1+kv)*gamma*zeta)-1)
                 * np.sin(beta-omega)**2 + np.sin(theta+omega)**2
                )/2*(
                 ((2*np.cos(theta)*np.cos(beta+theta)*J_alpha)/(np.cos(beta)*(1+kv)*gamma*zeta)-1)
                 * np.cos(beta-omega)**2 + np.sin(theta+omega)**2
                )
    )

    ### Equation 19
    sigma_AEH = sigma_alpha * np.cos(alpha_alpha + omega)
