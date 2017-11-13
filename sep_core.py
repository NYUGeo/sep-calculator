import numpy as np

class sep (object):
    """
    Add a description here...
    """

    def __init__(self, kh, kv, omega, beta, phi, gamma, c, H, g=9.80665):
        self.kh = kh
        self.kv = kv
        self.g = g     # Gravitational acceleration (m/s2)
        self.omega = np.radians(omega)
        self.beta = np.radians(beta)
        self.phi = np.radians(phi)
        self.gamma = gamma    # unit weight (kN/m3)
        self.c = c     # cohesion (kPa)
        self.H = H     # wall height (m)

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

    def zl(self, zw):
        """
        Depth along length of retaining wall
        """
        return zw/np.cos(self.omega)

    def z(self, zw):
        """
        Depth below the sloping ground surface
        """
        return zw * (np.cos(self.beta-self.omega))/(np.cos(self.beta) \
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

    def Ja(self, zw):
        """
        Equation 13: Active condition
        """
        J_p1 = ((self.gamma*self.z(zw)*np.cos(self.beta)*np.cos(self.beta \
                + self.theta())*(1+self.kv))/np.cos(self.theta())) \
                + self.c*np.cos(self.phi)*np.sin(self.phi)
        J_p2 = (self.gamma**2)*(self.z(zw)**2)*((np.cos(self.beta))**2) \
                * ((1+self.kv)**2) * ((((np.cos(self.beta+self.theta()))**2) \
                - ((np.cos(self.phi))**2))/((np.cos(self.theta()))**2))
        J_p3 = (self.c**2)*((np.cos(self.phi))**2)
        J_p4 = (2*self.c*self.gamma*self.z(zw)*np.cos(self.phi)*np.sin(self.phi) \
                * np.cos(self.beta) * np.cos(self.beta+self.theta())*(1+self.kv)) \
                / (np.cos(self.theta()))

        return (1/((np.cos(self.phi))**2))*(J_p1 - np.sqrt(np.absolute(J_p2+J_p3+J_p4)))

    def alpha_a(self, zw, degrees=False):
        """
        Equation 18: Obliquity
        """
        alpha_a_p1 = (((2*np.cos(self.theta())*np.cos(self.beta+self.theta())) \
                        / (np.cos(self.beta)*(1+self.kv)))*(self.Ja(zw) \
                        / (self.gamma*self.z(zw)))-1)
        alpha_a_p2 = 2*(np.sin(self.beta-self.omega))*(np.cos(self.beta-self.omega))
        alpha_a_p3 = 2*(np.sin(self.theta()+self.omega))*(np.cos(self.theta()+self.omega))

        alpha_a = np.arctan(((alpha_a_p1*alpha_a_p2)+alpha_a_p3) \
                    / (2*(alpha_a_p1*((np.cos(self.beta-self.omega))**2) \
                    + ((np.sin(self.theta()+self.omega))**2))))

        if degrees:
            return np.rad2deg(alpha_a)
        else:
            return alpha_a

    def Ka(self, zw):
        """
        Equation 16: Active lateral earth pressure coefficient
        """
        return ((np.cos(self.beta)*(1+self.kv)*(np.sin(self.theta()+self.omega)**2 \
                - np.cos(self.beta-self.omega)**2)) \
               /(np.cos(self.alpha_a(zw))*np.cos(self.beta+self.theta()) \
               * np.cos(self.theta()))
               ) + (
               ((2*(self.Ja(zw)/(self.gamma*self.z(zw)))*np.cos(self.beta-self.omega)**2) \
               / np.cos(self.alpha_a(zw))))

    def sigma_a(self, zw, user_Ka=999):
        """
        Equation 15: stress acting on the wall

        Allow manual user input for Ka to facilitate multiple layers.
        """
        if user_Ka == 999:
            return self.gamma * self.z(zw) * self.Ka(zw)
        else:
            return self.gamma * self.z(zw) * user_Ka

    def sigma_AEH(self, zw, user_Ka=999, user_alpha_a=999):
        """
        Equation 19:

        Allow manual user input for Ka and alpha_a to facilitate multiple layers.
        """
        if user_Ka == 999:
            return self.sigma_a(zw) * np.cos(self.alpha_a(zw) + self.omega)
        else:
            return self.sigma_a(zw,user_Ka) * np.cos(user_alpha_a + self.omega)

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
        var1 = self.Hz()
        var2 = self.sigma_AEH(self.H)
        var3 = self.sigma_AEH(0.1*self.H)

        return var1*(1-(0.9*var2)/(var2-var3))

    def Hzc(self):
        """
        add...
        """
        return self.zc()*(np.cos(self.beta)/np.cos(self.beta-self.omega))

    def Hp1(self):
        """
        add...
        """
        return (self.Hl() - self.Hzc())/3

    def P_AEH1(self):
        """
        Equation 22: Horizontal lateral earth thrust
        """
        return 0.5*self.sigma_AEH(self.H)*(self.Hl() - self.zc() \
                * (np.cos(self.beta)/np.cos(self.beta-self.omega)))

    def Hp2(self):
        """
        add...
        """
        return self.Hl()/3

    def P_AEH2(self):
        """
        Equation 24:
        """
        return 0.5*self.sigma_AEH(self.H)*self.Hl()


### Test with:
### figure5 = sep(0.2, -0.1, 20, 15, 30, 23, 20, 15)
### figure13 = sep(0.2, 0.1, 10, 15, 0, 20, 100, 10)


def line_circle_intersect(h,k,r,angle,c):

    # Slope
    a = np.tan(np.radians(angle))

    # Solving the system for x
    x_inter_pos = ((np.sqrt((a**2 + 1) * r**2 - c**2 - 2*c*(a*h - k)
                            - a**2 * h**2 + 2*a*h*k - k**2)
                    - c * a + a * k + h
                    )/(a**2 + 1))

    x_inter_neg = (-(np.sqrt((a**2 + 1) * r**2 - c**2 - 2*c*(a*h - k)
                            - a**2 * h**2 + 2*a*h*k - k**2)
                    + c * a - a * k - h
                    )/(a**2 + 1))

    x_intersect = np.array([x_inter_pos, x_inter_neg])

    y_intersect = a * x_intersect + c

    return x_intersect, y_intersect
