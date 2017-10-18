import numpy as np

h = 279.49854871962066
k = 0
r = 157.06978243549909
deg = 0
c = 20

# Slope
a = np.tan(np.radians(deg))

x = np.arange(0,500,0.001)

# Equation of a circle
p_sqr = np.sqrt(-h**2 + 2*h*x + r**2 - x**2)
n_sqr = - p_sqr

# Equation of line
tan_line = a*x + c

# Solving the system for x
x_inter_pos = ((np.sqrt((a**2 + 1) * r**2 - c**2 - 2*c*(a*h - k) - a**2 * h**2
                        + 2*a*h*k - k**2)
                - c * a + a * k + h
                )
               /(a**2 + 1))

x_inter_neg = (-(np.sqrt((a**2 + 1) * r**2 - c**2 - 2*c*(a*h - k) - a**2 * h**2
                        + 2*a*h*k - k**2)
                + c * a - a * k - h
                )
               /(a**2 + 1))

x_intersect = np.array([x_inter_pos, x_inter_neg])

y_intersect = a * x_intersect + c

print(x_intersect, y_intersect)


from bokeh.plotting import figure, output_file, show

output_file("line.html")

p = figure(plot_width=400, plot_height=400)

# add a line renderer
p.line(x, p_sqr + k, line_width=2)
p.line(x, n_sqr + k, line_width=2)

p.line(x, tan_line, line_width=2, color='orange')

p.circle(x_intersect,y_intersect, color='red', size=7)

show(p)
