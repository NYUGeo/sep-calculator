import numpy as np
from analysis import sep

from bokeh.plotting import figure
from bokeh.io import output_file, show, curdoc
from bokeh.charts import HeatMap, bins, output_file, show
from bokeh.models import ColumnDataSource, Slider
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models.glyphs import Patch
from bokeh.layouts import widgetbox, row, column

### Default input values
kh = 0.2
kv = -0.1
omega = 20
beta = 15
phi = 30
gamma = 23
c = 20
H = 15

### Retaining Wall Coordinates
xA = 0
yA = 0
angleA = np.radians(88)

xB = H / np.tan(angleA)
yB = H

xC = xB + 3
yC = H

xD = xC + H * np.tan(np.radians(omega))
yD = 0


### Backfill Coordinates
xE = xC + 15
yE = yC + 15 * np.tan(np.radians(beta))

xF = xE
yF = 0

### List of coordinates
x_wall = [xA, xB, xC, xD]
y_wall = [yA, yB, yC, yD]
x_earth = [xC, xE, xF, xD]
y_earth = [yC, yE, yF, yD]

### Calculating example values
example = sep(kh, kv, omega, beta, phi, gamma, c, H)

sigma_y = np.arange(0.0001, example.Hl(), 0.1)
sigma_x = example.sigma_AEH(sigma_y)

x_sigma = sigma_x.tolist()
y_sigma = sigma_y.tolist()
x_sigma.extend([0,0])
y_sigma.extend([example.Hl(),0])
source_example = ColumnDataSource(data=dict(x=x_sigma, y=y_sigma))

plot_sigma = figure(x_axis_label='sigma_AEH (kPa)', y_axis_label= \
                    "Depth Along Wall Length 'zl' (m)", y_range=(example.Hl(),0), \
                    plot_width=350, plot_height=350)
sigma_plot = Patch(x='x', y='y', fill_color = '#EEEEEE', line_color = 'black')
plot_sigma.add_glyph(source_example, sigma_plot)


### Define Data Sources
source_wall = ColumnDataSource(data=dict(x=x_wall,y=y_wall))
source_earth = ColumnDataSource(data=dict(x=x_earth,y=y_earth))

# Add patches to figure
plot = figure(x_axis_label='x-axis', y_axis_label='Height (m)', y_range=(0,23), \
                plot_width=350, plot_height=350)
plot.xaxis.visible = False
wall = Patch(x='x', y='y', fill_color = '#DCDDDE', line_color = 'black')
plot.add_glyph(source_wall, wall)
earth = Patch(x='x', y='y', fill_color = '#FFFF99', line_color = 'black')
plot.add_glyph(source_earth, earth)

zwi=[0.0001,(example.H/5),2*(example.H/5),3*(example.H/5),4*(example.H/5),(example.H)]
### Tabulating data
source_table = ColumnDataSource(data=dict(
        zw=[0.0001,(example.H/5),2*(example.H/5),3*(example.H/5),4*(example.H/5),(example.H)],
        zl=[example.zl(i) for i in zwi],
        z=[example.z(i) for i in zwi],
        Ja=[example.Ja(i) for i in zwi],
        a_a=[example.alpha_a(i, degrees=True) for i in zwi],
        Ka=[example.Ka(i) for i in zwi],
        s_a=[example.sigma_a(i) for i in zwi],
        s_AEH=[example.sigma_AEH(i) for i in zwi]
))

columns = [
        TableColumn(field='zw',title='Zw (m)'),
        TableColumn(field='zl',title='Zl (m)'),
        TableColumn(field='z',title='Z (m)'),
        TableColumn(field='Ja',title='Ja'),
        TableColumn(field='a_a',title='a (deg)'),
        TableColumn(field='Ka',title='Ka'),
        TableColumn(field='s_a',title='Sa (kPa)'),
        TableColumn(field='s_AEH',title='S_AEH (kPa)')
]

data_table = DataTable(source=source_table, columns=columns, width=700)

# Callback function that updates the plot
def update_plot(attr, old, new):
    omega = omega_slider.value
    beta = beta_slider.value
    phi = phi_slider.value
    H = H_value.value
    c = c_slider.value
    gamma = gamma_slider.value

    ### New Retaining Wall Coordinates
    new_xA = 0
    new_yA = 0
    new_angleA = np.radians(88)

    new_xB = H / np.tan(new_angleA)
    new_yB = H

    new_xC = new_xB + 3
    new_yC = H

    new_xD = new_xC + H * np.tan(np.radians(omega))
    new_yD = 0

    ### New Backfill Coordinates
    new_xE = new_xC + 15
    new_yE = new_yC + 15 * np.tan(np.radians(beta))

    new_xF = new_xE
    new_yF = 0

    ### New list of coordinates
    new_x_wall = [new_xA, new_xB, new_xC, new_xD]
    new_y_wall = [new_yA, new_yB, new_yC, new_yD]
    new_x_earth = [new_xC, new_xE, new_xF, new_xD]
    new_y_earth = [new_yC, new_yE, new_yF, new_yD]

    ### New example values
    new_example = sep(kh, kv, omega, beta, phi, gamma, c, H)
    new_sigma_y = np.arange(0.0001, new_example.Hl(), 0.1)
    new_sigma_x = new_example.sigma_AEH(new_sigma_y)
    new_x_sigma = new_sigma_x.tolist()
    new_y_sigma = new_sigma_y.tolist()
    new_x_sigma.extend([0,0])
    new_y_sigma.extend([new_example.Hl(),0])

    ### Update the data
    source_wall.data = dict(x=new_x_wall, y=new_y_wall)
    source_earth.data = dict(x=new_x_earth, y=new_y_earth)
    source_example.data = dict(x=new_x_sigma, y=new_y_sigma)


### Sliders & Inputs
H_value = Slider(start=10,end=25,step=1,value=15,title='Retaining Wall Height (m)')
omega_slider = Slider(start=0,end=25,step=1,value=20,title='Omega Angle (degrees)')
beta_slider = Slider(start=0,end=20,step=1,value=15,title='Beta Angle (degrees)')
phi_slider = Slider(start=0,end=45,step=1,value=30,title='Phi Angle (degrees)')
c_slider = Slider(start=0,end=100,step=5,value=20,title='Cohesion (kPa)')
gamma_slider = Slider(start=20,end=30,step=1,value=23,title='Unit Weight (kN/m3)')

# Attach the callback to the 'value' property of slider
H_value.on_change('value', update_plot)
omega_slider.on_change('value', update_plot)
beta_slider.on_change('value', update_plot)
phi_slider.on_change('value', update_plot)
c_slider.on_change('value', update_plot)
gamma_slider.on_change('value', update_plot)

# Make a row layout of widgetbox(slider) and plot and add it to the current document
layout = column(row(widgetbox(H_value,omega_slider,beta_slider,phi_slider,c_slider,gamma_slider),\
        plot, plot_sigma),data_table)

curdoc().add_root(layout)


# run with:
# bokeh serve --show web-app.py




# Specify the name of the output file and show the result
#output_file('web-plot.html')
#show(plot)
