from os.path import dirname, join

import numpy as np
from analysis import sep, line_circle_intersect

from bokeh.client import push_session
from bokeh import models
from bokeh.plotting import figure
from bokeh.io import output_file, show, curdoc
from bokeh.models import ColumnDataSource, Slider, Arrow, OpenHead, Label, LabelSet, Div
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models.glyphs import Patch
from bokeh.layouts import widgetbox, row, column, layout, Spacer


# Default input values
kh = 0.2
kv = -0.1
omega = 20
beta = 15
phi = 30
gamma = 23
c = 20
H = 15
ewt = -H
gamma_w = 9.8


##############################################################################
###                     INITIAL WALL AND SOIL PATCHES                      ###
##############################################################################

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


# EWT Coordinates
xK = xC + abs(ewt) * np.tan(np.radians(omega))
yK = H + ewt    # Note: ewt is negative

xL = xE
yL = H + ewt    # Note: ewt is negative


# List of coordinates
x_wall = [xA, xB, xC, xD]
y_wall = [yA, yB, yC, yD]
x_earth = [xC, xE, xF, xD]
y_earth = [yC, yE, yF, yD]
x_ewt = [xK, xL, xF, xD]
y_ewt = [yK, yL, yF, yD]


##############################################################################
###                                LAYERS                                  ###
##############################################################################

# Calculating values
example = sep(kh, kv, omega, beta, phi, gamma, c, H)

layer_dry = sep(kh, kv, omega, beta, phi, gamma, c, abs(ewt))
layer_wet = sep(kh, kv, omega, beta, phi, gamma-gamma_w, c, H+ewt)

all_layer_Hl = layer_dry.Hl() + layer_wet.Hl()

print(layer_dry.theta(), layer_wet.theta())
# Shoelace formula
# https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
def PolyArea(x,y):
    # for this problem, keep only the positive stresses
    x_pos = []
    y_pos = []

    for k, l in zip(x, y):
        if k >= 0:
            x_pos.append(k)
            y_pos.append(l)
        else:
            pass

    return 0.5 * np.abs(np.dot(x_pos, np.roll(y_pos, 1))
                        - np.dot(y_pos, np.roll(x_pos, 1)))



##############################################################################
###                              STRESS PLOT                               ###
##############################################################################

# Old
sigma_y = np.arange(0.0001, example.Hl(), 0.1)
sigma_x = example.sigma_AEH(sigma_y)

x_sigma = sigma_x.tolist()
y_sigma = sigma_y.tolist()
x_sigma.extend([0,0])
y_sigma.extend([example.Hl(),0])
source_example = ColumnDataSource(data=dict(x=x_sigma, y=y_sigma))


# New
sigma_y_dry = np.arange(0.0001, layer_dry.Hl(), 0.1)
sigma_x_dry = layer_dry.sigma_AEH(sigma_y_dry)
sigma_y_wet = np.arange(0.0001, layer_wet.Hl(), 0.1)
sigma_x_wet = (layer_wet.sigma_AEH(sigma_y_wet)
               + layer_dry.sigma_AEH(layer_dry.Hl()))

x_sigma_all = sigma_x_dry.tolist()
x_sigma_all.extend(sigma_x_wet.tolist())
x_sigma_all.extend([0, 0])
y_sigma_all = sigma_y_dry.tolist()
y_sigma_all.extend((sigma_y_wet + layer_dry.Hl()).tolist())
y_sigma_all.extend([all_layer_Hl, 0])
sigma_data = ColumnDataSource(data=dict(x=x_sigma_all, y=y_sigma_all))

load_height_top = all_layer_Hl - (all_layer_Hl - layer_dry.Hzc()) / 3
load_height_bot = (all_layer_Hl - layer_dry.Hzc()) / 3



plot_sigma = figure(#title="Horizontal Pseudo-Static Lateral Earth Pressure",
                    x_axis_label="\u03C3'\u1D00\u1D07\u029C (kPa)", # sigma_AEH
                    y_axis_label="Depth Along Wall Length 'Zl' (m)",
                    y_range=(0.99*example.Hl(),example.Hl()-30),
                    plot_width=200,
                    plot_height=400,
                    toolbar_location=None)

sigma_plot = Patch(x='x', y='y', fill_color = '#EEEEEE', line_color = 'black')
plot_sigma.add_glyph(source_example, sigma_plot)


# Arrow data
source_arrow1 = ColumnDataSource(data={
    'x0': [example.sigma_AEH(example.H)],
    'y0': [example.Hl() - example.Hp1()],
    'y1': [example.Hl() - example.Hp1()]
})


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



load1 = Label(x=0.3 * example.sigma_AEH(example.H),
              y=example.Hl() - example.Hp1(),
              text="{:.0f} kN".format(example.P_AEH1()),
              text_font_style='bold',
              border_line_width=2,
              text_font_size='16pt',
              text_color='red')
plot_sigma.add_layout(load1)

h_load1 = Label(x=0.3 * example.sigma_AEH(example.H),
                y=example.Hl() - example.Hp1(),
                text="@ {:.2f} m".format(example.Hp1()),
                text_font_style='bold',
                border_line_width=2,
                text_font_size='12pt',
                text_color='red',
                y_offset=-20)
plot_sigma.add_layout(h_load1)
plot_sigma.min_border_left = 50


error = Label(x=115,
              y=300,
              x_units='screen',
              y_units='screen',
              render_mode='css',
              text="Inadmissible\nCondition!",
              text_font_style='bold',
              text_align='center',
              border_line_width=2,
              text_font_size='12pt',
              text_color='red')
plot_sigma.add_layout(error)
error.visible = False







sigma_figure = figure(x_axis_label="\u03C3'\u1D00\u1D07\u029C (kPa)",  # sigma_AEH
                      y_axis_label="Depth Along Wall Length 'Zl' (m)",
                      y_range=(0.99 * all_layer_Hl, all_layer_Hl - 30),
                      plot_width=200,
                      plot_height=400,
                      toolbar_location=None,
                      toolbar_sticky=False)


sigma_patch = Patch(x='x', y='y', fill_color='#EEEEEE', line_color='black')
sigma_figure.add_glyph(sigma_data, sigma_patch)


force = PolyArea(x_sigma_all, y_sigma_all) * np.cos(np.radians(omega))

arrow_data = ColumnDataSource(data=dict(
    x0 = [max(x_sigma_all)],
    y0 = [load_height_top],
    y1 = [load_height_top]
))


force_arrow = Arrow(end=OpenHead(line_color="black",
                                 line_width=3,
                                 line_join='bevel'),
                    line_width=3,
                    x_start='x0',
                    y_start='y0',
                    x_end=0,
                    y_end='y1',
                    source=arrow_data)
sigma_figure.add_layout(force_arrow)


arrow_load = Label(x=0.27 * max(x_sigma_all),
                   y=load_height_top,
                   text='{:.0f} kN'.format(force),
                   text_font_style='bold',
                   border_line_width=2,
                   text_font_size='16pt',
                   text_color='red')
sigma_figure.add_layout(arrow_load)

arrow_height = Label(x=0.27 * max(x_sigma_all),
                     y=load_height_top,
                     text='@ {:.2f} m'.format(load_height_bot),
                     text_font_style='bold',
                     border_line_width=2,
                     text_font_size='12pt',
                     text_color='red',
                     y_offset=-20)
sigma_figure.add_layout(arrow_height)
sigma_figure.min_border_left = 50


sigma_error = Label(x=115,
              y=300,
              x_units='screen',
              y_units='screen',
              render_mode='css',
              text="Inadmissible\nCondition!",
              text_font_style='bold',
              text_align='center',
              border_line_width=2,
              text_font_size='12pt',
              text_color='red')
sigma_figure.add_layout(sigma_error)
sigma_error.visible = False








##############################################################################
###                               WALL PLOT                                ###
##############################################################################

# Define Data Sources
source_wall = ColumnDataSource(data=dict(x=x_wall,y=y_wall))
source_earth = ColumnDataSource(data=dict(x=x_earth,y=y_earth))
source_ewt = ColumnDataSource(data=dict(x=x_ewt,y=y_ewt))

# Add patches to figure
plot = figure(#title="Retaining Wall and Backfill Properties",
              #title_location="above",
              x_axis_label='x-axis',
              y_axis_label='Height (m)',
              plot_width=350,
              plot_height=400,
              y_range=(0,30),
              toolbar_location=None,
              background_fill_alpha=0.1)
plot.xaxis.visible = False

wall = Patch(x='x',
             y='y',
             fill_color = '#DCDDDE',
             line_color = 'black',
             line_width=1.5)
plot.add_glyph(source_wall, wall)

earth = Patch(x='x',
              y='y',
              fill_color = '#FFFF99',
              line_color = 'black',
              line_width=1.5)
plot.add_glyph(source_earth, earth)
plot.min_border_left = 50

water = Patch(x='x',
              y='y',
              fill_color = 'LightSeaGreen',
              fill_alpha = 0.075,
              line_color = 'LightSeaGreen',
              line_width=1.5)
plot.add_glyph(source_ewt, water)


wall_label_data = ColumnDataSource(data=dict(
                        x=[100,100,100,100],
                        y=[100-i*15 for i in range(4)],
                        names=['\u03B1h: {:.1f}g'.format(kh),
                               '\u03B1v: {:.1f}g'.format(kv),
                               '\u03C9: {:.0f}\u1d52'.format(omega),
                               '\u03B2: {:.0f}\u1d52'.format(beta)]))

wall_labels = LabelSet(x='x',
                       y='y',
                       x_units='screen',
                       y_units='screen',
                       text='names',
                       text_font_size='9pt',
                       text_color='black',
                       text_font_style='bold',
                       text_align='center',
                       #background_fill_color='white',
                       source=wall_label_data)
plot.add_layout(wall_labels)


soil_label_data = ColumnDataSource(data=dict(
                        x=[227,210,227,227,227],
                        y=[180-i*15 for i in range(5)],
                        names=['H: {:.1f} m'.format(H),
                               'EWT: {:.1f} m'.format(ewt),
                               '\u03C6: {:.0f}\u1d52'.format(phi),
                               'c: {:.0f} kPa'.format(c),
                               '\u03B3: {:.0f}\u1d52'.format(gamma)]))

soil_labels = LabelSet(x='x',
                       y='y',
                       x_units='screen',
                       y_units='screen',
                       text='names',
                       text_font_size='9pt',
                       text_color='black',
                       text_font_style='bold',
                       text_align='left',
                       source=soil_label_data)
plot.add_layout(soil_labels)

error_c = Label(x=150,
              y=300,
              x_units='screen',
              y_units='screen',
              render_mode='css',
              text="Set c = 0 or drop EWT below wall",
              text_font_style='bold',
              text_align='center',
              border_line_width=2,
              text_font_size='12pt',
              text_color='red')
plot.add_layout(error_c)
error_c.visible = False


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
                x_conj = [0, 1.9*example.Ja(H)],
                y_conj = [0, (np.tan(np.radians(beta) + example.theta())
                             ) * 1.9 * example.Ja(H)],
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


mohr_plot = figure(x_axis_label="\u03C3' (kPa)",
                   y_axis_label='\u03C4 (kPa)',
                   plot_width=400,
                   plot_height=400,
                   toolbar_location=None,
                   x_range=(0, 2.0 * example.Ja(H)),
                   y_range=(0, 2.03 * example.Ja(H)),
                   background_fill_alpha=0.1)


fail_line = mohr_plot.line(
                x='x_fail',
                y='y_fail',
                source=mohr_line_data,
                line_width=2,
                legend='Effective stress M-C failure envelope')

conj_line = mohr_plot.line(
                x='x_conj',
                y='y_conj',
                source=mohr_line_data,
                line_width=3,
                color='orange',
                legend='Conjugate stress line')

mohr_plot.circle(x='center',
                 y='y',
                 radius='radius',
                 fill_color=None,
                 line_width=2,
                 color='black',
                 source=mohr_circle_data)


# Calculate intersection points for effective stress envelope
x_line_inter, y_line_inter = line_circle_intersect(
                                h=mohr_circle_data.data['center'][0],
                                k=0,
                                r=mohr_circle_data.data['radius'][0],
                                angle=phi,
                                c=c)

# Calculate intersection points for conjugate stress line
x_conj_inter, y_conj_inter = line_circle_intersect(
                                h=mohr_circle_data.data['center'][0],
                                k=0,
                                r=mohr_circle_data.data['radius'][0],
                                angle=beta+np.degrees(example.theta()),
                                c=0)

# Calculate intersection points at zero (sigma1, sigma3)
x_sigma_inter, y_sigma_inter = line_circle_intersect(
                                h=mohr_circle_data.data['center'][0],
                                k=0,
                                r=mohr_circle_data.data['radius'][0],
                                angle=0,
                                c=0)


# Store data in a ColumnDataSource
intersect_data = ColumnDataSource(data=dict(
                x_line_inter = x_line_inter,
                y_line_inter = y_line_inter,
                x_conj_inter = x_conj_inter,
                y_conj_inter = y_conj_inter,
                x_sigma_inter = x_sigma_inter,
                y_sigma_inter = y_sigma_inter
                )
            )

# Plot intersection points for effective stress envelope
# mohr_plot.circle(x='x_line_inter',
#                  y='y_line_inter',
#                  line_color='#1F77B4',
#                  line_width=2,
#                  fill_color='white',
#                  size=7,
#                  source=intersect_data)

# Plot intersection points for conjugate stress line
mohr_plot.circle(x='x_conj_inter',
                 y='y_conj_inter',
                 line_color='orange',
                 line_width=2,
                 fill_color='white',
                 size=7,
                 source=intersect_data)

# Plot intersection points at zero (sigma1, sigma3)
# mohr_plot.circle(x='x_sigma_inter',
#                  y='y_sigma_inter',
#                  line_color='black',
#                  line_width=2,
#                  fill_color='white',
#                  size=5,
#                  source=intersect_data)


# arc_data = ColumnDataSource(data=dict(
#                 x_phi = [-c/np.radians(phi)],
#                 y_phi = [0],
#                 phi_radius = [0.15 * 2.0 * example.Ja(H)],
#                 phi_end_angle = [np.radians(phi)],
#                 x_beta = [0],
#                 y_beta = [0],
#                 beta_radius = [0.16 * 2.0 * example.Ja(H)],
#                 beta_end_angle = [np.radians(beta)+example.theta()]
#                 )
#             )
#
# mohr_plot.arc(x='x_phi',
#               y='y_phi',
#               radius='phi_radius',
#               start_angle=0,
#               end_angle='phi_end_angle',
#               line_width=2,
#               color="#1F77B4",
#               source=arc_data)
#
# mohr_plot.arc(x='x_beta',
#               y='y_beta',
#               radius='beta_radius',
#               start_angle=0,
#               end_angle='beta_end_angle',
#               line_width=2,
#               color="orange",
#               source=arc_data)


mohr_bold_label_data = ColumnDataSource(data=dict(
                        x=[70,77.8,66,70,70],
                        y=[300-i*17 for i in range(5)],
                        names=['Zw: {:.1f} m'.format(H),
                               '\u03D5: {:.0f}\u1d52'.format(phi),
                               '\u03B2+\u03B8: {:.0f}\u1d52'.format(beta +
                                            np.degrees(example.theta())),
                               "\u03C3'\u03B2: {:.0f} kPa".format(
                                    intersect_data.data['x_conj_inter'][1]),
                               "\u03C3'\u03B8: {:.0f} kPa".format(
                                    intersect_data.data['x_conj_inter'][0]),
                               ]))

mohr_bold_labels = LabelSet(x='x',
                       y='y',
                       x_units='screen',
                       y_units='screen',
                       text='names',
                       text_font_size='9pt',
                       text_color='black',
                       text_font_style='bold',
                       text_align='left',
                       background_fill_color='white',
                       source=mohr_bold_label_data)
mohr_plot.add_layout(mohr_bold_labels)


mohr_sigma_label_data = ColumnDataSource(data=dict(
                        x=intersect_data.data['x_conj_inter'],
                        y=intersect_data.data['y_conj_inter'],
                        names=["\u03C3'\u03B8",
                               "\u03C3'\u03B2",
                               ]))

mohr_sigma_labels = LabelSet(x='x',
                       y='y',
                       text='names',
                       text_font_size='9pt',
                       text_color='black',
                       text_font_style='bold',
                       text_align='left',
                       background_fill_color='white',
                       background_fill_alpha=0.6,
                       x_offset=8,
                       #y_offset=5,
                       source=mohr_sigma_label_data)
mohr_plot.add_layout(mohr_sigma_labels)

mohr_plot.legend.location = "top_left"




circle_center = layer_dry.Ja(min(H,abs(ewt))) + layer_wet.Ja(max(0,H+ewt))

mohr_plot2 = figure(x_axis_label="\u03C3' (kPa)",
                   y_axis_label='\u03C4 (kPa)',
                   plot_width=400,
                   plot_height=400,
                   toolbar_location=None,
                   x_range=(0, 2.0 * circle_center),
                   y_range=(0, 2.03 * circle_center),
                   background_fill_alpha=0.1)


mohr2_circle_data = ColumnDataSource(data=dict(
                y = [0],
                center = [circle_center],
                # Circle radius from equation 11
                radius = [(c * (1/np.tan(np.radians(phi)))
                           + circle_center) * np.sin(np.radians(phi))],
                )
            )


mohr_plot2.circle(x='center',
                 y='y',
                 radius='radius',
                 fill_color=None,
                 line_width=2,
                 color='black',
                 source=mohr2_circle_data)










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


# Callback function that updates all plots
def update_plot(attr, old, new):
    omega = omega_slider.value
    beta = beta_slider.value
    phi = phi_slider.value
    H = H_value.value
    c = c_slider.value
    gamma = gamma_slider.value
    kh = kh_slider.value
    kv = kv_slider.value
    ewt = ewt_slider.value
    if ewt == 0:
        ewt = 0.0000001

    ewt_slider.start = -H
    zw_slider.end = H
    zw_slider.value = H

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

    # New EWT Coordinates
    new_xK = new_xC + abs(ewt) * np.tan(np.radians(omega))
    new_yK = H + ewt    # Note: ewt is negative

    new_xL = new_xE
    new_yL = H + ewt    # Note: ewt is negative

    ### New list of coordinates
    new_x_wall = [new_xA, new_xB, new_xC, new_xD]
    new_y_wall = [new_yA, new_yB, new_yC, new_yD]
    new_x_earth = [new_xC, new_xE, new_xF, new_xD]
    new_y_earth = [new_yC, new_yE, new_yF, new_yD]
    new_x_ewt = [new_xK, new_xL, new_xF, new_xD]
    new_y_ewt = [new_yK, new_yL, new_yF, new_yD]

    ### New calculated values
    new_example = sep(kh, kv, omega, beta, phi, gamma, c, H)
    new_layer_dry = sep(kh, kv, omega, beta, phi, gamma, c, min(H,abs(ewt)))
    new_layer_wet = sep(kh, kv, omega, beta, phi, gamma-gamma_w, c, max(0,H+ewt))
    new_all_layer_Hl = new_layer_dry.Hl() + new_layer_wet.Hl()

    mohr_plot.x_range.end = 2.0 * new_example.Ja(H)
    mohr_plot.y_range.end = 2.03 * new_example.Ja(H)

    # Old
    new_sigma_y = np.arange(0.0001, new_example.Hl(), 0.1)
    new_sigma_x = new_example.sigma_AEH(new_sigma_y)
    new_x_sigma = new_sigma_x.tolist()
    new_y_sigma = new_sigma_y.tolist()
    new_x_sigma.extend([0,0])
    new_y_sigma.extend([new_example.Hl(),0])

    # New
    new_sigma_y_dry = np.arange(0.0001, new_layer_dry.Hl(), 0.1)
    new_sigma_x_dry = new_layer_dry.sigma_AEH(new_sigma_y_dry)
    new_sigma_y_wet = np.arange(0.0001, new_layer_wet.Hl(), 0.1)
    new_sigma_x_wet = new_layer_wet.sigma_AEH(new_sigma_y_wet) + new_layer_dry.sigma_AEH(new_layer_dry.Hl())
    new_x_sigma_all = new_sigma_x_dry.tolist()
    new_x_sigma_all.extend(new_sigma_x_wet.tolist())
    new_x_sigma_all.extend([0,0])
    new_y_sigma_all = new_sigma_y_dry.tolist()
    new_y_sigma_all.extend((new_sigma_y_wet + new_layer_dry.Hl()).tolist())
    new_y_sigma_all.extend([new_all_layer_Hl,0])

    ### Update the data
    source_wall.data = dict(x=new_x_wall, y=new_y_wall)
    source_earth.data = dict(x=new_x_earth, y=new_y_earth)
    source_ewt.data = dict(x=new_x_ewt, y=new_y_ewt)
    source_example.data = dict(x=new_x_sigma, y=new_y_sigma)
    sigma_data.data=dict(x=new_x_sigma_all, y=new_y_sigma_all)
    new_force = PolyArea(new_x_sigma_all,new_y_sigma_all) * np.cos(np.radians(omega))
    print(new_force)
    plot_sigma.y_range.start=0.99*new_example.Hl()
    plot_sigma.y_range.end=(new_example.Hl()-30)

    sigma_figure.y_range.start=0.99*(new_layer_dry.Hl()+new_layer_wet.Hl())
    sigma_figure.y_range.end= (new_layer_dry.Hl()+new_layer_wet.Hl()) - 30

    wall_label_data.data=dict(
                            x=[100,100,100,100],
                            y=[100-i*15 for i in range(4)],
                            names=['\u03B1h: {:.1f}g'.format(kh),
                                   '\u03B1v: {:.1f}g'.format(kv),
                                   '\u03C9: {:.0f}\u1d52'.format(omega),
                                   '\u03B2: {:.0f}\u1d52'.format(beta)])

    soil_label_data.data=dict(
                            x=[227,210,227,227,227],
                            y=[120*(H/10)-i*15 for i in range(5)],
                            names=['H: {:.1f} m'.format(H),
                                   'EWT: {:.1f} m'.format(ewt),
                                   '\u03C6: {:.0f}\u1d52'.format(phi),
                                   'c: {:.0f} kPa'.format(c),
                                   '\u03B3: {:.0f}\u1d52'.format(gamma)])


    # Update Mohr circle
    mohr_line_data.data = dict(
        x_fail = [0, 1.5*new_example.Ja(H)],
        y_fail = [c, c + (np.tan(np.radians(phi))*1.5*new_example.Ja(H))],
        x_conj = [0, 1.9*new_example.Ja(H)],
        y_conj = [0, (np.tan(np.radians(beta) + new_example.theta())
                     ) * 1.9 * new_example.Ja(H)],
        )

    mohr_circle_data.data = dict(
        y = [0],
        center = [new_example.Ja(H)],
        # Circle radius from equation 11
        radius = [(c * (1/np.tan(np.radians(phi)))
                   + new_example.Ja(H)) * np.sin(np.radians(phi))],
        )

    #mohr_depth.text="Depth (Zw): {:.1f} m".format(H)



    # Calculate intersection points for effective stress envelope
    x_line_inter, y_line_inter = line_circle_intersect(
                                    h=mohr_circle_data.data['center'][0],
                                    k=0,
                                    r=mohr_circle_data.data['radius'][0]*1.0000001,
                                    angle=phi,
                                    c=c)

    # Calculate intersection points for conjugate stress line
    x_conj_inter, y_conj_inter = line_circle_intersect(
                                    h=mohr_circle_data.data['center'][0],
                                    k=0,
                                    r=mohr_circle_data.data['radius'][0],
                                    angle=beta+np.degrees(new_example.theta()),
                                    c=0)

    # Calculate intersection points at zero (sigma1, sigma3)
    x_sigma_inter, y_sigma_inter = line_circle_intersect(
                                    h=mohr_circle_data.data['center'][0],
                                    k=0,
                                    r=mohr_circle_data.data['radius'][0],
                                    angle=0,
                                    c=0)


    # Store data in a ColumnDataSource
    intersect_data.data=dict(
                    x_line_inter = x_line_inter,
                    y_line_inter = y_line_inter,
                    x_conj_inter = x_conj_inter,
                    y_conj_inter = y_conj_inter,
                    x_sigma_inter = x_sigma_inter,
                    y_sigma_inter = y_sigma_inter
                    )

    # arc_data.data=dict(
    #                 x_phi = [-c/np.radians(phi)],
    #                 y_phi = [0],
    #                 phi_radius = [0.15 * 2.0 * new_example.Ja(H)],
    #                 phi_end_angle = [np.radians(phi)],
    #                 x_beta = [0],
    #                 y_beta = [0],
    #                 beta_radius = [0.16 * 2.0 * new_example.Ja(H)],
    #                 beta_end_angle = [np.radians(beta)+new_example.theta()]
    #                 )

    mohr_bold_label_data.data=dict(
                            x=[70,77.8,66,70,70],
                            y=[300-i*17 for i in range(5)],
                            names=['Zw: {:.1f} m'.format(H),
                                   '\u03D5: {:.0f}\u1d52'.format(phi),
                                   '\u03B2+\u03B8: {:.0f}\u1d52'.format(beta +
                                                np.degrees(new_example.theta())),
                                   "\u03C3'\u03B2: {:.0f} kPa".format(
                                        intersect_data.data['x_conj_inter'][1]),
                                   "\u03C3'\u03B8: {:.0f} kPa".format(
                                        intersect_data.data['x_conj_inter'][0]),
                                   ])

    mohr_sigma_label_data.data=dict(
                            x=intersect_data.data['x_conj_inter'],
                            y=intersect_data.data['y_conj_inter'],
                            names=["\u03C3'\u03B8",
                                   "\u03C3'\u03B2",
                                   ])


    ### IF clause for beta + theta < phi
    # UPDATE: IF clause for NO INTERSECTION
    #if beta + np.degrees(new_example.theta()) > phi:
    if np.isnan(intersect_data.data['x_conj_inter'][0]):
        # sigma_plot.fill_alpha = 0
        # sigma_plot.line_alpha = 0
        # arrow1.visible = False
        # load1.visible = False
        # h_load1.visible = False
        error.visible = True
        sigma_error.visible = True
        mohr_plot.background_fill_color = 'red'

    else:
        # sigma_plot.fill_alpha = 1
        # sigma_plot.line_alpha = 1
        # arrow1.visible = True
        # load1.visible = True
        # h_load1.visible = True
        error.visible = False
        sigma_error.visible = False
        mohr_plot.background_fill_color = None


    # Error for c>0 and ewt
    if (abs(ewt) < H) and (c > 0):
        error_c.visible = True
        plot.background_fill_color = 'red'
    else:
        error_c.visible = False
        plot.background_fill_color = None


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
    load1.x=0.3*new_example.sigma_AEH(new_example.H)
    load1.y=new_example.Hl()-new_example.Hp1()
    h_load1.text="@ {:.2f} m".format(new_example.Hp1())
    h_load1.x=0.3*new_example.sigma_AEH(new_example.H)
    h_load1.y=new_example.Hl()-new_example.Hp1()

    # New arrow
    new_load_height_top = new_all_layer_Hl - \
        (new_all_layer_Hl - new_layer_dry.Hzc()) / 3
    new_load_height_bot = (new_all_layer_Hl - new_layer_dry.Hzc()) / 3

    arrow_data.data=dict(
        x0= [max(new_x_sigma_all)],
        y0= [new_load_height_top],
        y1= [new_load_height_top]
    )

    arrow_load.x = 0.27 * max(new_x_sigma_all)
    arrow_load.y = new_load_height_top
    arrow_load.text = '{:.0f} kN'.format(new_force)
    arrow_height.x = 0.27 * max(new_x_sigma_all)
    arrow_height.y = new_load_height_top
    arrow_height.text='@ {:.2f} m'.format(new_load_height_bot)




# Callback function that updates mohr circle for zw
def update_mohr_zw(attr, old, new):
    omega = omega_slider.value
    beta = beta_slider.value
    phi = phi_slider.value
    H = H_value.value
    c = c_slider.value
    gamma = gamma_slider.value
    kh = kh_slider.value
    kv = kv_slider.value

    zwd = zw_slider.value

    new_example = sep(kh, kv, omega, beta, phi, gamma, c, H)

    # Update Mohr circle
    mohr_line_data.data = dict(
        x_fail = [0, 1.5*new_example.Ja(zwd)],
        y_fail = [c, c + (np.tan(np.radians(phi))*1.5*new_example.Ja(zwd))],
        x_conj = [0, 1.7*new_example.Ja(zwd)],
        y_conj = [0, (np.tan(np.radians(beta) + new_example.theta())
                     ) * 1.7 * new_example.Ja(zwd)],
        )

    mohr_circle_data.data = dict(
        y = [0],
        center = [new_example.Ja(zwd)],
        # Circle radius from equation 11
        radius = [(c * (1/np.tan(np.radians(phi)))
                   + new_example.Ja(zwd)) * np.sin(np.radians(phi))],
        )


    # Calculate intersection points for effective stress envelope
    x_line_inter, y_line_inter = line_circle_intersect(
                                    h=mohr_circle_data.data['center'][0],
                                    k=0,
                                    r=mohr_circle_data.data['radius'][0]*1.0000001,
                                    angle=phi,
                                    c=c)

    # Calculate intersection points for conjugate stress line
    x_conj_inter, y_conj_inter = line_circle_intersect(
                                    h=mohr_circle_data.data['center'][0],
                                    k=0,
                                    r=mohr_circle_data.data['radius'][0],
                                    angle=beta+np.degrees(new_example.theta()),
                                    c=0)

    # Calculate intersection points at zero (sigma1, sigma3)
    x_sigma_inter, y_sigma_inter = line_circle_intersect(
                                    h=mohr_circle_data.data['center'][0],
                                    k=0,
                                    r=mohr_circle_data.data['radius'][0],
                                    angle=0,
                                    c=0)


    # Store data in a ColumnDataSource
    intersect_data.data=dict(
                    x_line_inter = x_line_inter,
                    y_line_inter = y_line_inter,
                    x_conj_inter = x_conj_inter,
                    y_conj_inter = y_conj_inter,
                    x_sigma_inter = x_sigma_inter,
                    y_sigma_inter = y_sigma_inter
                    )

    mohr_bold_label_data.data=dict(
                            x=[70,77.8,66,70,70],
                            y=[300-i*17 for i in range(5)],
                            names=['Zw: {:.1f} m'.format(zwd),
                                   '\u03D5: {:.0f}\u1d52'.format(phi),
                                   '\u03B2+\u03B8: {:.0f}\u1d52'.format(beta +
                                                np.degrees(new_example.theta())),
                                   "\u03C3'\u03B2: {:.0f} kPa".format(
                                        intersect_data.data['x_conj_inter'][1]),
                                   "\u03C3'\u03B8: {:.0f} kPa".format(
                                        intersect_data.data['x_conj_inter'][0]),
                                   ])


    mohr_sigma_label_data.data=dict(
                            x=intersect_data.data['x_conj_inter'],
                            y=intersect_data.data['y_conj_inter'],
                            names=["\u03C3'\u03B8",
                                   "\u03C3'\u03B2",
                                   ])



    # UPDATE: IF clause for NO INTERSECTION
    if np.isnan(intersect_data.data['x_conj_inter'][0]):
        # sigma_plot.fill_alpha = 0
        # sigma_plot.line_alpha = 0
        # arrow1.visible = False
        # load1.visible = False
        # h_load1.visible = False
        # error.visible = True
        mohr_plot.background_fill_color = 'red'

    else:
        # sigma_plot.fill_alpha = 1
        # sigma_plot.line_alpha = 1
        # arrow1.visible = True
        # load1.visible = True
        # h_load1.visible = True
        # error.visible = False
        mohr_plot.background_fill_color = None


# Sliders
H_value = Slider(
                start=10,
                end=25,
                step=0.5,
                value=15,
                title='Wall height, H, (m)')
omega_slider = Slider(
                start=0,
                end=30,
                step=1,
                value=20,
                title='Wall inclination, \u03C9, (deg.)')
beta_slider = Slider(
                start=-30,
                end=30,
                step=1,
                value=15,
                title='Surface slope, \u03B2, (deg.)')
phi_slider = Slider(
                start=0,
                end=45,
                step=1,
                value=30,
                title='Ιnternal friction, \u03C6, (deg.)')
c_slider = Slider(
                start=0,
                end=100,
                step=5,
                value=20,
                title='Cohesion, c, (kPa)')
gamma_slider = Slider(
                start=16,
                end=25,
                step=1,
                value=23,
                title='Unit weight, \u03B3, (kN/m\u00B3)')
kh_slider = Slider(
                start=0,
                end=0.4,
                step=0.1,
                value=0.2,
                title='Horizontal, kh')
kv_slider = Slider(
                start=-0.4,
                end=0.4,
                step=0.1,
                value=-0.1,
                title='Vertical, kv')
ewt_slider = Slider(
                start=-H,
                end=0,
                step=0.1,
                value=-H,
                orientation="vertical",
                direction='rtl',
                show_value=False,
                height=340,
                width=25,
                callback_throttle=0,
                tooltips=False,
                bar_color='#e8f7f6')
zw_slider = Slider(
                start=0,
                end=H,
                step=0.1,
                value=H,
                orientation="vertical",
                direction='ltr',
                show_value=False,
                height=340,
                width=25,
                tooltips=False)


# Attach the callback to the 'value' property of slider
H_value.on_change('value', update_plot)
omega_slider.on_change('value', update_plot)
beta_slider.on_change('value', update_plot)
phi_slider.on_change('value', update_plot)
c_slider.on_change('value', update_plot)
gamma_slider.on_change('value', update_plot)
kh_slider.on_change('value', update_plot)
kv_slider.on_change('value', update_plot)
ewt_slider.on_change('value', update_plot)
zw_slider.on_change('value', update_mohr_zw)


page_header = Div(text=open(join(dirname(__file__), "page_header.html")).read(), width=1050)
page_footer = Div(text=open(join(dirname(__file__), "page_footer.html")).read(), width=1050)
wall_controls = [H_value,omega_slider,beta_slider]
wall_inputs = widgetbox(*wall_controls, width=200)

seismic_controls = [kh_slider,kv_slider]
seismic_inputs = widgetbox(*seismic_controls, width=180)

soil_controls = [phi_slider,c_slider,gamma_slider]
soil_inputs = widgetbox(*soil_controls, width=200)

# The layout function replaces the row and column functions
page_layout = layout([
                [page_header],
                [Div(text="<h3>Wall Properties:</h3>", width=200),
                 Div(text="<h3>Seismic Coefficients:</h3>", width=180),
                 Div(text="<h3>Soil Properties:</h3>", width=200)],
                [wall_inputs, seismic_inputs, soil_inputs],
                #[Div(text="<hr>", width=1050)],
                #[Div(text="<h3>Soil Properties:</h3>", width=260)],
                #[soil_inputs],
                [Div(text="<hr>", width=1050)],
                [Div(text="<h4></h4>", width=40),
                 Div(text="<h4>Retaining Wall and Backfill Properties</h4>",
                     width=300),
                 Div(text="<h4>EWT<br>&nbsp;(m)</h4>", width=95),
                 Div(text="<h4>Horizontal Pseudo-Static<br>"
                          "Lateral Earth Pressure</h4>",
                     width=225),
                 Div(text="<h4>Mohr's circle with failure envelopes at "
                          "depth, Zw,<br>from the top of wall surface</h4>",
                     width=350),
                 Div(text="<h4>Zw<br>(m)</h4>", width=70)],
                [plot, ewt_slider, Spacer(width=20),
                 sigma_figure, Spacer(width=20),
                 mohr_plot, zw_slider, mohr_plot2],
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
