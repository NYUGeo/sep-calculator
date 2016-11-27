import numpy as np

class sep (object):
    """
    Add a description here...
    """

    def __init__(self, kh, kv, omega, beta, phi, gamma, c, H, zw, g=9.80665):
        self.kh = kh
        self.kv = kv
        self.g = g     # Gravitational acceleration (m/s2)
        self.omega = np.radians(omega)
        self.beta = np.radians(beta)
        self.phi = np.radians(phi)
        self.gamma = gamma    # unit weight (kN/m3)
        self.c = c     # cohesion (kPa)
        self.H = H     # wall height (m)
        self.zw = zw

    def alpha_h(self):
        """
        Horizontal component of pseudo-static acceleration
        """
        return self.kh * self.g

    def alpha_v(self):
        """
        Vertical component of pseudo-static acceleration
        """
        return self.kv * self.g

    def zl(self):
        """
        Depth along length of retaining wall
        """
        return self.zw/np.cos(self.omega)

    def z(self):
        """
        Depth below the sloping ground surface
        """
        return self.zw * (np.cos(self.beta-self.omega))/(np.cos(self.beta) \
                * np.cos(self.omega))

    def g_theta(self):
        """
        Equation 1: Modified gravitational acceleration
        """
        return np.sqrt((self.g + self.alpha_v())**2 + self.alpha_h()**2)

    def theta(self):
        """
        Equation 2: Angle to vertical
        """
        return np.arctan(self.kh/(1 + self.kv))

    def gamma_theta(self):
        """
        Equation 3: Modified unit weight
        """
        return (self.gamma * (1 + self.kv))/(np.cos(self.theta()))

    def Ja(self):
        """
        Equation 13: Active condition
        """
        J_p1 = ((self.gamma*self.z()*np.cos(self.beta)*np.cos(self.beta \
                + self.theta())*(1+self.kv))/np.cos(self.theta())) \
                + self.c*np.cos(self.phi)*np.sin(self.phi)
        J_p2 = (self.gamma**2)*(self.z()**2)*((np.cos(self.beta))**2) \
                * ((1+self.kv)**2) * ((((np.cos(self.beta+self.theta()))**2) \
                - ((np.cos(self.phi))**2))/((np.cos(self.theta()))**2))
        J_p3 = (self.c**2)*((np.cos(self.phi))**2)
        J_p4 = (2*self.c*self.gamma*self.z()*np.cos(self.phi)*np.sin(self.phi) \
                * np.cos(self.beta) * np.cos(self.beta+self.theta())*(1+self.kv)) \
                / (np.cos(self.theta()))

        return (1/((np.cos(self.phi))**2))*(J_p1 - np.sqrt(J_p2+J_p3+J_p4))

    def alpha_a(self, degrees=False):
        """
        Equation 18: Obliquity
        """
        alpha_a_p1 = (((2*np.cos(self.theta())*np.cos(self.beta+self.theta())) \
                        / (np.cos(self.beta)*(1+self.kv)))*(self.Ja() \
                        / (self.gamma*self.z()))-1)
        alpha_a_p2 = 2*(np.sin(self.beta-self.omega))*(np.cos(self.beta-self.omega))
        alpha_a_p3 = 2*(np.sin(self.theta()+self.omega))*(np.cos(self.theta()+self.omega))

        alpha_a = np.arctan(((alpha_a_p1*alpha_a_p2)+alpha_a_p3) \
                    / (2*(alpha_a_p1*((np.cos(self.beta-self.omega))**2) \
                    + ((np.sin(self.theta()+self.omega))**2))))

        if degrees:
            return np.rad2deg(alpha_a)
        else:
            return alpha_a

    def Ka(self):
        """
        Equation 16: Active lateral earth pressure coefficient
        """
        return ((np.cos(self.beta)*(1+self.kv)*(np.sin(self.theta()+self.omega)**2 \
                - np.cos(self.beta-self.omega)**2)) \
               /(np.cos(self.alpha_a())*np.cos(self.beta+self.theta()) \
               * np.cos(self.theta()))
               ) + (
               ((2*(self.Ja()/(self.gamma*self.z()))*np.cos(self.beta-self.omega)**2) \
               / np.cos(self.alpha_a())))

    def sigma_a(self):
        """
        Equation 15: stress acting on the wall
        """
        return self.gamma * self.z() * self.Ka()

    def sigma_AEH(self):
        """
        Equation 19:
        """
        return self.sigma_a() * np.cos(self.alpha_a() + self.omega)

    def Hz(self):
        """
        Equation 21: Vertical distance between the heel of the retaining wall
                     and its backfill slope surface
        """
        return self.H * (np.cos(self.beta-self.omega)) \
                / (np.cos(self.beta)*np.cos(self.omega))

    def Hl(self):
        """
        Equation 23: Length of retaining wall
        """
        return self.H/np.cos(self.omega)


class sep2 (object):
    """
    Add a description here...
    """

    def __init__(self, kh, kv, omega, beta, phi, gamma, c, H, zw, g=9.80665):
        self.kh = kh
        self.kv = kv
        self.g = g     # Gravitational acceleration (m/s2)
        self.omega = np.radians(omega)
        self.beta = np.radians(beta)
        self.phi = np.radians(phi)
        self.gamma = gamma    # unit weight (kN/m3)
        self.c = c     # cohesion (kPa)
        self.H = H     # wall height (m)
        self.zw = zw

    def alpha_h(self):
        """
        Horizontal component of pseudo-static acceleration
        """
        return self.kh * self.g

    def alpha_v(self):
        """
        Vertical component of pseudo-static acceleration
        """
        return self.kv * self.g

    def zl(self):
        """
        Depth along length of retaining wall
        """
        return self.zw/np.cos(self.omega)

    def z(self):
        """
        Depth below the sloping ground surface
        """
        return self.zw * (np.cos(self.beta-self.omega))/(np.cos(self.beta) \
                * np.cos(self.omega))

    def g_theta(self):
        """
        Equation 1: Modified gravitational acceleration
        """
        return np.sqrt((self.g + self.alpha_v())**2 + self.alpha_h()**2)

    def theta(self):
        """
        Equation 2: Angle to vertical
        """
        return np.arctan(self.kh/(1 + self.kv))

    def gamma_theta(self):
        """
        Equation 3: Modified unit weight
        """
        return (self.gamma * (1 + self.kv))/(np.cos(self.theta()))

    def Ja(self, z=None):
        """
        Equation 13: Active condition
        """
        if z:
            z = z
        else:
            z = self.z()

        J_p1 = ((self.gamma*z*np.cos(self.beta)*np.cos(self.beta \
                + self.theta())*(1+self.kv))/np.cos(self.theta())) \
                + self.c*np.cos(self.phi)*np.sin(self.phi)
        J_p2 = (self.gamma**2)*(z**2)*((np.cos(self.beta))**2) \
                * ((1+self.kv)**2) * ((((np.cos(self.beta+self.theta()))**2) \
                - ((np.cos(self.phi))**2))/((np.cos(self.theta()))**2))
        J_p3 = (self.c**2)*((np.cos(self.phi))**2)
        J_p4 = (2*self.c*self.gamma*z*np.cos(self.phi)*np.sin(self.phi) \
                * np.cos(self.beta) * np.cos(self.beta+self.theta())*(1+self.kv)) \
                / (np.cos(self.theta()))

        return (1/((np.cos(self.phi))**2))*(J_p1 - np.sqrt(J_p2+J_p3+J_p4))

    def alpha_a(self, z=None, degrees=False):
        """
        Equation 18: Obliquity
        """
        if z:
            z = z
        else:
            z = self.z()

        alpha_a_p1 = (((2*np.cos(self.theta())*np.cos(self.beta+self.theta())) \
                        / (np.cos(self.beta)*(1+self.kv)))*(self.Ja() \
                        / (self.gamma*z))-1)
        alpha_a_p2 = 2*(np.sin(self.beta-self.omega))*(np.cos(self.beta-self.omega))
        alpha_a_p3 = 2*(np.sin(self.theta()+self.omega))*(np.cos(self.theta()+self.omega))

        alpha_a = np.arctan(((alpha_a_p1*alpha_a_p2)+alpha_a_p3) \
                    / (2*(alpha_a_p1*((np.cos(self.beta-self.omega))**2) \
                    + ((np.sin(self.theta()+self.omega))**2))))

        if degrees:
            return np.rad2deg(alpha_a)
        else:
            return alpha_a

    def Ka(self, z=None):
        """
        Equation 16: Active lateral earth pressure coefficient
        """
        if z:
            z = z
        else:
            z = self.z()

        return ((np.cos(self.beta)*(1+self.kv)*(np.sin(self.theta()+self.omega)**2 \
                - np.cos(self.beta-self.omega)**2)) \
               /(np.cos(self.alpha_a())*np.cos(self.beta+self.theta()) \
               * np.cos(self.theta()))
               ) + (
               ((2*(self.Ja()/(self.gamma*z))*np.cos(self.beta-self.omega)**2) \
               / np.cos(self.alpha_a())))

    def sigma_a(self, z=None):
        """
        Equation 15: stress acting on the wall
        """
        if z:
            z = z
        else:
            z = self.z()

        return self.gamma * z * self.Ka()

    def sigma_AEH(self):
        """
        Equation 19:
        """
        return self.sigma_a() * np.cos(self.alpha_a() + self.omega)

    def Hz(self):
        """
        Equation 21: Vertical distance between the heel of the retaining wall
                     and its backfill slope surface
        """
        return self.H * (np.cos(self.beta-self.omega)) \
                / (np.cos(self.beta)*np.cos(self.omega))

    def Hl(self):
        """
        Equation 23: Length of retaining wall
        """
        return self.H/np.cos(self.omega)

    def zc(self):
        """
        Equation 20: Depth of tension crack
        """
        #var1 = self.Hl()
        #var2 = self.Ja(z=var1)

        #return self.Hl() * (1 - ()/())



### Test with:
### test = sep(0.2, -0.1, 20, 15, 30, 23, 20, 15, 6)
