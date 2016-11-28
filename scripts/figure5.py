import numpy as np

kh = 0.2
kv = -0.1
omega_deg = 20
beta_deg = 15
phi_deg = 30
c = 20
gamma = 23
H = 15
zw = 6

### Convert degrees to radians
omega = np.radians(omega_deg)
beta = np.radians(beta_deg)
phi = np.radians(phi_deg)

### Gravitational & pseudo-static accelerations
g = 9.80665
alpha_h = kh * g
alpha_v = kv * g

### Heights
zl = zw/np.cos(omega)
z = zw * (np.cos(beta-omega))/(np.cos(beta)*np.cos(omega))

### Equation 1
g_theta = np.sqrt((g + alpha_v)**2 + alpha_h**2)

### Equation 2
theta = np.arctan(kh/(1 + kv))
theta_deg = np.rad2deg(theta)

### Equation 3
gamma_theta = (gamma * (1 + kv))/(np.cos(theta))

### Equation 13
J_part1 = ((gamma*z*np.cos(beta)*np.cos(beta+theta)*(1+kv))/np.cos(theta)) \
            + c*np.cos(phi)*np.sin(phi)
J_part2 = (gamma**2)*(z**2)*((np.cos(beta))**2)*((1+kv)**2) \
            * ((((np.cos(beta+theta))**2)-((np.cos(phi))**2))/((np.cos(theta))**2))
J_part3 = (c**2)*((np.cos(phi))**2)
J_part4 = (2*c*gamma*z*np.cos(phi)*np.sin(phi)*np.cos(beta)*np.cos(beta+theta)*(1+kv))/(np.cos(theta))

Ja = (1/((np.cos(phi))**2))*(J_part1 - np.sqrt(J_part2+J_part3+J_part4))

### Equation 18
alpha_a_p1 = (((2*np.cos(theta)*np.cos(beta+theta))/(np.cos(beta)*(1+kv)))*(Ja/(gamma*z))-1)
alpha_a_p2 = 2*(np.sin(beta-omega))*(np.cos(beta-omega))
alpha_a_p3 = 2*(np.sin(theta+omega))*(np.cos(theta+omega))
alpha_a = np.arctan(((alpha_a_p1*alpha_a_p2)+alpha_a_p3) \
            /(2*(alpha_a_p1*((np.cos(beta-omega))**2)+((np.sin(theta+omega))**2))))
alpha_a_deg = np.rad2deg(alpha_a)

### Equation 16
Ka = (
      (np.cos(beta)*(1+kv)*(np.sin(theta+omega)**2 - np.cos(beta-omega)**2)) \
       /(np.cos(alpha_a)*np.cos(beta+theta)*np.cos(theta))
      ) + (
      ((2*(Ja/(gamma*z))*np.cos(beta-omega)**2)/np.cos(alpha_a))
)

### Equation 15
sigma_a = gamma * z * Ka

### Equation 19
sigma_AEH = sigma_a * np.cos(alpha_a + omega)

### Equation 21
Hz = H * (np.cos(beta-omega))/(np.cos(beta)*np.cos(omega))

print ('zw =',zw,'\n'
       'zl =',zl,'\n'
       ' z =',z,'\n'
       'g_theta =',g_theta,'\n'
       'theta_deg =',theta_deg,'\n'
       #'theta =',theta,'\n'
       #'gamma_theta =',gamma_theta,'\n'
       #'J_part1 =',J_part1,'\n'
       #'J_part2 =',J_part2,'\n'
       #'J_part3 =',J_part3,'\n'
       #'J_part4 =',J_part4,'\n'
       'Ja =',Ja,'\n'
       #'alpha_a =',alpha_a,'\n'
       'alpha_a_deg =',alpha_a_deg,'\n'
       'Ka =',Ka,'\n'
       'sigma_a =',sigma_a,'\n'
       'sigma_AEH =',sigma_AEH)