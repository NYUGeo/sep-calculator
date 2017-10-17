from os.path import dirname, join

import numpy as np
from analysis import sep

from bokeh.client import push_session
from bokeh import models
from bokeh.plotting import figure
from bokeh.io import output_file, show, curdoc
from bokeh.models import ColumnDataSource, Slider, Arrow, OpenHead, Label, Div
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models.glyphs import Patch
from bokeh.layouts import widgetbox, row, column, layout


# Default input values
kh = 0.2
kv = -0.1
omega = 20
beta = 15
phi = 30
gamma = 23
c = 20
H = 15

# Retaining Wall Coordinates
xA = 0
yA = 0
angleA = np.radians(88)

xB = H / np.tan(angleA)
yB = H

xC = xB + 3
yC = H

xD = xC + H * np.tan(np.radians(omega))
yD = 0


# Backfill Coordinates
xE = xC + 15
yE = yC + 15 * np.tan(np.radians(beta))

xF = xE
yF = 0


# List of coordinates
x_wall = [xA, xB, xC, xD]
y_wall = [yA, yB, yC, yD]
x_earth = [xC, xE, xF, xD]
y_earth = [yC, yE, yF, yD]

# Calculating example values
example = sep(kh, kv, omega, beta, phi, gamma, c, H)

sigma_y = np.arange(0.0001, example.Hl(), 0.1)
sigma_x = example.sigma_AEH(sigma_y)

x_sigma = sigma_x.tolist()
y_sigma = sigma_y.tolist()
x_sigma.extend([0,0])
y_sigma.extend([example.Hl(),0])
source_example = ColumnDataSource(data=dict(x=x_sigma, y=y_sigma))

# Arrow data
source_arrow1 = ColumnDataSource(data={
        'x0': [example.sigma_AEH(example.H)],
        'y0': [example.Hl()-example.Hp1()],
        'y1': [example.Hl()-example.Hp1()]
})


##############################################################################
###                              STRESS PLOT                               ###
##############################################################################

# Tools in toolbox
TOOLS = 'pan,box_zoom,reset,save'

plot_sigma = figure(#title="Horizontal Pseudo-Static Lateral Earth Pressure",
                    #title_location="above",
                    x_axis_label="\u03C3'\u1D00\u1D07\u029C (kPa)", # sigma_AEH
                    y_axis_label="Depth Along Wall Length 'Zl' (m)",
                    y_range=(0.99*example.Hl(),example.Hl()-30),
                    plot_width=250,
                    plot_height=400,
                    #toolbar_location="above",
                    toolbar_location=None,
                    toolbar_sticky=False,
                    tools=TOOLS)

sigma_plot = Patch(x='x', y='y', fill_color = '#EEEEEE', line_color = 'black')
plot_sigma.add_glyph(source_example, sigma_plot)

arrow1 = Arrow(end=OpenHead(line_color="black",
                            line_width=3,
                            line_join='bevel'),
               line_width=3,
               x_start='x0',
               y_start='y0',
               x_end=0,
               y_end='y1',
               source=source_arrow1)
plot_sigma.add_layout(arrow1)

load1 = Label(x=0.5*example.sigma_AEH(example.H),
              y=example.Hl()-example.Hp1(),
              text="{:.0f} kN".format(example.P_AEH1()),
              text_font_style='bold',
              border_line_width=2,
              text_font_size='16pt',
              text_color='red')
plot_sigma.add_layout(load1)

h_load1 = Label(x=0.5*example.sigma_AEH(example.H),
                y=example.Hl()-example.Hp1(),
                text="@ {:.2f} m".format(example.Hp1()),
                text_font_style='bold',
                border_line_width=2,
                text_font_size='12pt',
                text_color='red',
                y_offset=-20)
plot_sigma.add_layout(h_load1)
plot_sigma.min_border_left = 50



##############################################################################
###                               WALL PLOT                                ###
##############################################################################

# Define Data Sources
source_wall = ColumnDataSource(data=dict(x=x_wall,y=y_wall))
source_earth = ColumnDataSource(data=dict(x=x_earth,y=y_earth))


# Add patches to figure
plot = figure(#title="Retaining Wall and Backfill Properties",
              #title_location="above",
              x_axis_label='x-axis',
              y_axis_label='Height (m)',
              plot_width=350,
              plot_height=400,
              y_range=(0,30),
              #toolbar_location="above",
              toolbar_location=None,
              toolbar_sticky=False,
              tools=TOOLS)
plot.xaxis.visible = False

wall = Patch(x='x', y='y', fill_color = '#DCDDDE', line_color = 'black')
plot.add_glyph(source_wall, wall)

earth = Patch(x='x', y='y', fill_color = '#FFFF99', line_color = 'black')
plot.add_glyph(source_earth, earth)
plot.min_border_left = 50


##############################################################################
###                              MOHR CIRCLE                               ###
##############################################################################

# Slope of failure line
fail_slope = np.tan(np.radians(phi))

# Slope of conjugate line
conj_slope = np.tan(np.radians(beta) + example.theta())

mohr_line_data = ColumnDataSource(data=dict(
                x_fail = [0, 1.5*example.Ja(H)],
                y_fail = [c, c + (np.tan(np.radians(phi))*1.5*example.Ja(H))],
                x_conj = [0, 1.7*example.Ja(H)],
                y_conj = [0, (np.tan(np.radians(beta) + example.theta())
                             ) * 1.7 * example.Ja(H)],
                )
            )


mohr_circle_data = ColumnDataSource(data=dict(
                y = [0],
                center = [example.Ja(H)],
                # Circle radius from equation 11
                radius = [(c * (1/np.tan(np.radians(phi)))
                           + example.Ja(H)) * np.sin(np.radians(phi))],
                )
            )



mohr_plot = figure(x_axis_label="\u03C3' (kPa???)",
                   y_axis_label='\u03C4 (kPa)',
                   plot_width=400,
                   plot_height=400,
                   toolbar_location=None,
                   x_range=(0, 2.0 * example.Ja(H)),
                   y_range=(0, 2.03 * example.Ja(H)),
                   #match_aspect=True,
                   #aspect_scale = 1.0
                   )
# mohr_plot.y_range.start = 0
# mohr_plot.x_range.start = 0
# mohr_plot.y_range.end = 2.0 * example.Ja(H)
# mohr_plot.x_range.end = 2.0 * example.Ja(H)
#mohr_plot.xaxis.visible = False
#mohr_plot.yaxis.visible = False


fail_line = mohr_plot.line(
                x='x_fail',
                y='y_fail',
                source=mohr_line_data,
                line_width=2,
                legend='Effective stress failure envelope')

conj_line = mohr_plot.line(
                x='x_conj',
                y='y_conj',
                source=mohr_line_data,
                line_width=3,
                color='orange',
                legend='Conjugate stress failure envelope')

mohr_plot.circle(x='center',
                 y='y',
                 radius='radius',
                 fill_color=None,
                 line_width=2,
                 color='black',
                 #radius_dimension='y',
                 source=mohr_circle_data)



##############################################################################
###                                 TABLE                                  ###
##############################################################################


### Tabulating data
zwi=[0.0001,(example.H/5),2*(example.H/5),3*(example.H/5),4*(example.H/5),(example.H)]
source_table = ColumnDataSource(data=dict(
                    zw  = zwi,
                    zl  = [example.zl(i) for i in zwi],
                    z   = [example.z(i) for i in zwi],
                    Ja  = [example.Ja(i) for i in zwi],
                    a_a = [example.alpha_a(i, degrees=True) for i in zwi],
                    Ka  = [example.Ka(i) for i in zwi],
                    s_a = [example.sigma_a(i) for i in zwi],
                    s_AEH=[example.sigma_AEH(i) for i in zwi]
))

columns = [
        TableColumn(field='zw',
                    title='Zw (m)',
                    formatter=models.NumberFormatter(format='0.000')),
        TableColumn(field='zl',
                    title='Zl (m)',
                    formatter=models.NumberFormatter(format='0.000')),
        TableColumn(field='z',
                    title='Z (m)',
                    formatter=models.NumberFormatter(format='0.000')),
        TableColumn(field='Ja',
                    title='J\u03B1',
                    formatter=models.NumberFormatter(format='0.000')),
        TableColumn(field='a_a',
                    title='\u03B1\u2090 (deg)',
                    formatter=models.NumberFormatter(format='0.000')),
        TableColumn(field='Ka',
                    title='K\u03B1',
                    formatter=models.NumberFormatter(format='0.000')),
        TableColumn(field='s_a',
                    title="\u03C3'\u03B1 (kPa)",
                    formatter=models.NumberFormatter(format='0.000')),
        TableColumn(field='s_AEH',
                    title="\u03C3'\u1D00\u1D07\u029C (kPa)",
                    formatter=models.NumberFormatter(format='0.000'))
]

data_table = DataTable(source=source_table,
                       columns=columns,
                       width=1050,
                       height=200)



##############################################################################
###                               CALLBACK                                 ###
##############################################################################


# Callback function that updates the plot
def update_plot(attr, old, new):
    omega = omega_slider.value
    beta = beta_slider.value
    phi = phi_slider.value
    H = H_value.value
    c = c_slider.value
    gamma = gamma_slider.value
    kh = kh_slider.value
    kv = kv_slider.value


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

    mohr_plot.x_range.end = 2.0 * new_example.Ja(H)
    mohr_plot.y_range.end = 2.03 * new_example.Ja(H)


    # error = Label(x=100,
    #               y=200,
    #               x_units='screen',
    #               y_units='screen',
    #               text="Inadmissible Condition!",
    #               text_font_style='bold',
    #               text_align='center',
    #               border_line_width=2,
    #               text_font_size='12pt',
    #               text_color='red',
    #               text_alpha=0)
    # plot_sigma.add_layout(error)
    #
    #
    ### IF clause for beta + theta < phi
    if beta + np.degrees(new_example.theta()) > phi:
        #conj_line.visible = False
        sigma_plot.fill_alpha = 0
        sigma_plot.line_alpha = 0
        arrow1.visible = False
        load1.visible = False
        h_load1.visible = False

    else:
        sigma_plot.fill_alpha = 1
        sigma_plot.line_alpha = 1
        arrow1.visible = True
        load1.visible = True
        h_load1.visible = True



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
    plot_sigma.y_range.start=0.99*new_example.Hl()
    plot_sigma.y_range.end=(new_example.Hl()-30)

    # Update Mohr circle
    mohr_line_data.data = dict(
        x_fail = [0, 1.5*new_example.Ja(H)],
        y_fail = [c, c + (np.tan(np.radians(phi))*1.5*new_example.Ja(H))],
        x_conj = [0, 1.7*new_example.Ja(H)],
        y_conj = [0, (np.tan(np.radians(beta) + new_example.theta())
                     ) * 1.7 * new_example.Ja(H)],
        )

    mohr_circle_data.data = dict(
        y = [0],
        center = [new_example.Ja(H)],
        # Circle radius from equation 11
        radius = [(c * (1/np.tan(np.radians(phi)))
                   + new_example.Ja(H)) * np.sin(np.radians(phi))],
        )





    ### Update tabulated data
    new_zwi=[0.0001,(new_example.H/5),2*(new_example.H/5),3*(new_example.H/5), \
             4*(new_example.H/5),(new_example.H)]
    source_table.data = dict(
            zw=new_zwi,
            zl=[new_example.zl(i) for i in new_zwi],
            z=[new_example.z(i) for i in new_zwi],
            Ja=[new_example.Ja(i) for i in new_zwi],
            a_a=[new_example.alpha_a(i, degrees=True) for i in new_zwi],
            Ka=[new_example.Ka(i) for i in new_zwi],
            s_a=[new_example.sigma_a(i) for i in new_zwi],
            s_AEH=[new_example.sigma_AEH(i) for i in new_zwi]
    )

    ### Update arrow data
    source_arrow1.data = {
            'x0': [new_example.sigma_AEH(new_example.H)],
            'y0': [new_example.Hl()-new_example.Hp1()],
            'y1': [new_example.Hl()-new_example.Hp1()]
    }
    load1.text="{:.0f} kN".format(new_example.P_AEH1())
    load1.x=0.5*new_example.sigma_AEH(new_example.H)
    load1.y=new_example.Hl()-new_example.Hp1()
    h_load1.text="@ {:.2f} m".format(new_example.Hp1())
    h_load1.x=0.5*new_example.sigma_AEH(new_example.H)
    h_load1.y=new_example.Hl()-new_example.Hp1()


# Sliders
H_value = Slider(
                start=10,
                end=25,
                step=0.5,
                value=15,
                title='Retaining wall height, H, (m)')
omega_slider = Slider(
                start=0,
                end=30,
                step=1,
                value=20,
                title='Wall inclination, \u03C9, (degrees)')
beta_slider = Slider(
                start=-30,
                end=30,
                step=1,
                value=15,
                title='Surface slope, \u03B2, (degrees)')
phi_slider = Slider(
                start=0,
                end=45,
                step=1,
                value=30,
                title='Ιnternal friction angle, \u03C6, (degrees)')
c_slider = Slider(
                start=0,
                end=100,
                step=5,
                value=20,
                title='Cohesion, c, (kPa)')
gamma_slider = Slider(
                start=16,
                end=22,
                step=1,
                value=20,
                title='Unit weight, \u03B3, (kN/m\u00B3)')
kh_slider = Slider(
                start=0,
                end=0.4,
                step=0.1,
                value=0.2,
                title='Horizontal Seismic Coefficient, kh')
kv_slider = Slider(
                start=-0.4,
                end=0.4,
                step=0.1,
                value=-0.1,
                title='Vertical Seismic Coefficient, kv')
v_slider = Slider(
                start=-0.4,
                end=0.4,
                step=0.1,
                value=-0.1,
                title='Vertical Slider',
                orientation="vertical",
                height=300)


# Attach the callback to the 'value' property of slider
H_value.on_change('value', update_plot)
omega_slider.on_change('value', update_plot)
beta_slider.on_change('value', update_plot)
phi_slider.on_change('value', update_plot)
c_slider.on_change('value', update_plot)
gamma_slider.on_change('value', update_plot)
kh_slider.on_change('value', update_plot)
kv_slider.on_change('value', update_plot)


page_header = Div(text=open(join(dirname(__file__), "page_header.html")).read(), width=1050)
page_footer = Div(text=open(join(dirname(__file__), "page_footer.html")).read(), width=1050)
wall_controls = [H_value,omega_slider,beta_slider]
wall_inputs = widgetbox(*wall_controls, width=260)

seismic_controls = [kh_slider,kv_slider]
seismic_inputs = widgetbox(*seismic_controls, width=260)

soil_controls = [phi_slider,c_slider,gamma_slider]
soil_inputs = widgetbox(*soil_controls, width=250)

# The layout function replaces the row and column functions
page_layout = layout([
                [page_header],
                [Div(text="<h3>Wall Properties:</h3>", width=260),
                 Div(text="<h3>Soil Properties:</h3>", width=250)],
                [wall_inputs,soil_inputs],
                [Div(text="<hr>", width=1050)],
                [Div(text="<h3>Seismic Coefficients:</h3>", width=260)],
                [seismic_inputs],
                [Div(text="<h4></h4>", width=40),
                 Div(text="<h4>Retaining Wall and Backfill Properties</h4>",
                     width=350),
                 Div(text="<h4>Horizontal Pseudo-Static<br>"
                          "Lateral Earth Pressure</h4>",
                     width=240),
                 Div(text="<h4>Mohr circle corresponding to conjugate stress "
                          "state at depth of the horizontal component of the "
                          "lateral earth thrust</h4>",
                     width=390)],
                [plot,plot_sigma,mohr_plot],
                [v_slider],
                [Div(text="<h3>Calculated values at several vertical depths "
                          "from the top of wall surface, Zw</h3>",
                     width=600)],
                [data_table],
                [page_footer]
])


curdoc().add_root(page_layout)
curdoc().title = "SEP Calculator"

### test with:
### bokeh serve --show sep-calculator.py

### run forever on server with:
### nohup bokeh serve sep-calculator.py --allow-websocket-origin cue3.engineering.nyu.edu:5006 --host cue3.engineering.nyu.edu:5006
