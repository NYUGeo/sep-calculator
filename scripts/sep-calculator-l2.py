from os.path import dirname, join

import numpy as np
from analysis import sep, line_circle_intersect

from bokeh.client import push_session
from bokeh import models
from bokeh.plotting import figure
from bokeh.io import output_file, show, curdoc
from bokeh.models import (ColumnDataSource, Slider, Arrow, OpenHead, Label,
                            LabelSet, Div, Legend)
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models.glyphs import Patch
from bokeh.layouts import widgetbox, row, column, layout, Spacer


# Default input values
kh = 0.15
kv = 0
omega = 15
beta1 = 10
beta2 = 10

H1 = 6
phi1 = 30
gamma1 = 23
c1 = 0

H2 = 9
phi2 = 40
gamma2 = 20
c2 = 0

H = H1 + H2
ewt = -3
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

xC = xB + 4
yC = H

xD = xC + H * np.tan(np.radians(omega))
yD = 0


# Backfill Coordinates
# Layer 1
xE = xC + 15
yE = yC + 15 * np.tan(np.radians(beta1))

xM = xC + H1 * np.tan(np.radians(omega))
yM = H2

xN = xE
#yN = yM
yN = yM + (15 - H1 * np.tan(np.radians(omega))) * np.tan(np.radians(beta2))

# Layer 2
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
x_layer1 = [xC, xE, xN, xM]
y_layer1 = [yC, yE, yN, yM]
x_layer2 = [xM, xN, xF, xD]
y_layer2 = [yM, yN, yF, yD]
x_ewt = [xK, xL, xF, xD]
y_ewt = [yK, yL, yF, yD]


##############################################################################
###                                LAYERS                                  ###
##############################################################################

if abs(ewt) >= H:
    layer1_dry = sep(kh, kv, omega, beta1, phi1, gamma1, c1, H1)
    layer2_dry = sep(kh, kv, omega, beta2, phi2, gamma2, c2, H2)
    layer1_wet = sep(kh, kv, omega, beta1, phi1, gamma1-gamma_w, c1, 0.00000001)
    layer2_wet = sep(kh, kv, omega, beta2, phi2, gamma2-gamma_w, c2, 0.00000001)
elif H1 <= abs(ewt) < H:
    layer1_dry = sep(kh, kv, omega, beta1, phi1, gamma1, c1, H1)
    temp_var = abs(ewt)-H1
    layer2_dry = sep(kh, kv, omega, beta2, phi2, gamma2, c2,
                     temp_var if temp_var !=0 else 0.00000001)
    layer1_wet = sep(kh, kv, omega, beta1, phi1, gamma1-gamma_w, c1, 0.00000001)
    layer2_wet = sep(kh, kv, omega, beta2, phi2, gamma2-gamma_w, c2, H+ewt)
else:
    layer1_dry = sep(kh, kv, omega, beta1, phi1, gamma1, c1, abs(ewt))
    layer2_dry = sep(kh, kv, omega, beta2, phi2, gamma2, c2, 0.00000001)
    layer1_wet = sep(kh, kv, omega, beta1, phi1, gamma1-gamma_w, c1, H1+ewt)
    layer2_wet = sep(kh, kv, omega, beta2, phi2, gamma2-gamma_w, c2, H2)

all_layer_Hl = (layer1_dry.Hl() + layer2_dry.Hl()
                + layer1_wet.Hl() + layer2_wet.Hl())


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
###                               WALL PLOT                                ###
##############################################################################

# Define Data Sources
wall_data = ColumnDataSource(data=dict(x=x_wall,y=y_wall))
earth_data = ColumnDataSource(data=dict(x1=x_layer1,
                                        y1=y_layer1,
                                        x2=x_layer2,
                                        y2=y_layer2))
ewt_data = ColumnDataSource(data=dict(x=x_ewt,y=y_ewt))

# Add patches to figure
wall_plot = figure(
    # title="Retaining Wall and Backfill Properties",
    # title_location="above",
    x_axis_label='x-axis',
    y_axis_label='Height (m)',
    plot_width=375,
    plot_height=400,
    y_range=(0, 30),
    toolbar_location=None,
    background_fill_alpha=0.1)
wall_plot.xaxis.visible = False

wall = Patch(
    x='x',
    y='y',
    fill_color='#DCDDDE',
    line_color='black',
    line_width=1.5)
wall_plot.add_glyph(wall_data, wall)

layer1_patch = Patch(
    x='x1',
    y='y1',
    fill_color='#FFFF99',
    line_color='black',
    line_width=1.5)
wall_plot.add_glyph(earth_data, layer1_patch)
wall_plot.min_border_left = 50

layer2_patch = Patch(
    x='x2',
    y='y2',
    fill_color='#e5e589',
    line_color='black',
    line_width=1.5)
wall_plot.add_glyph(earth_data, layer2_patch)

water = Patch(
    x='x',
    y='y',
    fill_color='LightSeaGreen',
    fill_alpha=0.075,
    line_color='LightSeaGreen',
    line_width=1.5)
wall_plot.add_glyph(ewt_data, water)

wall_label_data = ColumnDataSource(data=dict(
    x=[100, 100, 100, 100, 100],
    y=[110 - i * 15 for i in range(5)],
    names=['WT: {:.1f} m'.format(ewt),
           'H: {:.1f} m'.format(H1+H2),
           '\u03B1h: {:.2f}g'.format(kh),
           '\u03B1v: {:.2f}g'.format(kv),
           '\u03C9: {:.0f}\u1d52'.format(omega)]))

wall_labels = LabelSet(
    x='x',
    y='y',
    x_units='screen',
    y_units='screen',
    text='names',
    text_font_size='9pt',
    text_color='black',
    text_font_style='bold',
    text_align='center',
    source=wall_label_data)
wall_plot.add_layout(wall_labels)

soil_label_data = ColumnDataSource(data=dict(
    x=[10, 14, 10, 10, 14, 14, 10, 14, 10, 10, 14, 14],
    y=[H2+H1/2+1, H2+H1/2+1, H2+H1/2, H2+H1/2-1, H2+H1/2, H2+H1/2-1,
       H2/2+1, H2/2+1, H2/2, H2/2-1, H2/2, H2/2-1],
    names=['LAYER 1',
           '\u03B2: {:.0f}\u1d52'.format(beta1),
           'H: {:.1f} m'.format(H1),
           '\u03C6: {:.0f}\u1d52'.format(phi1),
           'c: {:.0f} kPa'.format(c1),
           '\u03B3: {:.0f} kPa'.format(gamma1),
           'LAYER 2',
           '\u03B2: {:.0f}\u1d52'.format(beta2),
           'H: {:.1f} m'.format(H2),
           '\u03C6: {:.0f}\u1d52'.format(phi2),
           'c: {:.0f} kPa'.format(c2),
           '\u03B3: {:.0f} kPa'.format(gamma2)]))

soil_labels = LabelSet(
    x='x',
    y='y',
    #x_units='screen',
    #y_units='screen',
    text='names',
    text_font_size='9pt',
    text_color='black',
    text_font_style='bold',
    text_align='left',
    source=soil_label_data)
wall_plot.add_layout(soil_labels)


# error_c = Label(
#     x=150,
#     y=300,
#     x_units='screen',
#     y_units='screen',
#     render_mode='css',
#     text="Set c = 0 or drop EWT below wall",
#     text_font_style='bold',
#     text_align='center',
#     border_line_width=2,
#     text_font_size='12pt',
#     text_color='red')
# wall_plot.add_layout(error_c)
# error_c.visible = False





##############################################################################
###                              STRESS PLOT                               ###
##############################################################################

sigma_l1_y_dry = np.arange(0.0001, layer1_dry.Hl(), 0.1)
sigma_l1_y_wet = np.arange(0.0001, layer1_wet.Hl(), 0.1)
sigma_l1_x_dry = layer1_dry.sigma_AEH(sigma_l1_y_dry)
sigma_l1_x_wet = (layer1_wet.sigma_AEH(sigma_l1_y_wet)
                  + layer1_dry.sigma_AEH(layer1_dry.Hl()))
sigma_l2_y_dry = np.arange(0.0001, layer2_dry.Hl(), 0.1)
sigma_l2_y_wet = np.arange(0.0001, layer2_wet.Hl(), 0.1)
sigma_l2_x_dry = (layer2_dry.sigma_AEH(sigma_l2_y_dry)
                  + layer1_dry.sigma_AEH(
                            layer1_dry.Hl(),
                            user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001),
                            user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001)))
sigma_l2_x_wet = (layer2_wet.sigma_AEH(sigma_l2_y_wet)
                  + layer2_dry.sigma_AEH(layer2_dry.Hl())
                  + layer1_dry.sigma_AEH(
                            layer1_dry.Hl(),
                            user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001),
                            user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001))
                  + layer1_wet.sigma_AEH(
                            layer1_wet.Hl(),
                            user_Ka=layer2_wet.Ka(layer2_wet.H * 0.0001),
                            user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001)))

# sigma_l1_y_dry = np.arange(0.0001, layer1_dry.Hl(), 0.1)
# sigma_l1_x_dry = layer1_dry.sigma_AEH(sigma_l1_y_dry)
# sigma_l1_y_wet = np.arange(0.0001, layer1_wet.Hl(), 0.1)
# sigma_l1_x_wet = (layer1_wet.sigma_AEH(sigma_l1_y_wet)
#                   + layer1_dry.sigma_AEH(layer1_dry.Hl()))
# sigma_l2_y_dry = np.arange(0.0001, layer2_dry.Hl(), 0.1)
# sigma_l2_x_dry = (layer2_dry.sigma_AEH(sigma_l2_y_dry)
#                   + layer1_dry.sigma_AEH(layer1_dry.Hl()))
#
# sigma_l2_y_wet = np.arange(0.0001, layer2_wet.Hl(), 0.1)
# sigma_l2_x_wet = (layer2_wet.sigma_AEH(sigma_l2_y_wet)
#                   + layer2_dry.sigma_AEH(layer2_dry.Hl())
#                   + layer1_dry.sigma_AEH(layer1_dry.Hl())
#                   + layer1_wet.sigma_AEH(layer1_wet.Hl()))

y_sigma_all = sigma_l1_y_dry.tolist()
x_sigma_all = sigma_l1_x_dry.tolist()
y_sigma_all.extend((sigma_l1_y_wet + layer1_dry.Hl()).tolist())
x_sigma_all.extend(sigma_l1_x_wet.tolist())
y_sigma_all.extend((sigma_l2_y_dry + layer1_dry.Hl()).tolist())
x_sigma_all.extend(sigma_l2_x_dry.tolist())
y_sigma_all.extend((sigma_l2_y_wet + layer1_dry.Hl() + layer1_wet.Hl()
                    + layer2_dry.Hl()).tolist())
x_sigma_all.extend(sigma_l2_x_wet.tolist())
y_sigma_all.extend([all_layer_Hl])
x_sigma_all.extend([0])

sigma_data = ColumnDataSource(data=dict(x=x_sigma_all, y=y_sigma_all))

load_height_top = all_layer_Hl - all_layer_Hl / 3
load_height_bot = all_layer_Hl / 3


sigma_figure = figure(x_axis_label="\u03C3'\u1D00\u1D07\u029C (kPa)",  # sigma_AEH
                      y_axis_label="Depth Along Wall Length 'Zl' (m)",
                      y_range=(0.99 * all_layer_Hl, all_layer_Hl - 30),
                      plot_width=200,
                      plot_height=400,
                      toolbar_location=None,
                      toolbar_sticky=False)

sigma_patch = Patch(x='x', y='y', fill_color='#EEEEEE', line_color='black')
sigma_figure.add_glyph(sigma_data, sigma_patch)

force = PolyArea(x_sigma_all, y_sigma_all) * (np.cos(np.radians(omega))**2)

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
                     text='@ ~{:.2f} m'.format(load_height_bot),
                     text_font_style='bold',
                     border_line_width=2,
                     text_font_size='12pt',
                     text_color='red',
                     y_offset=-20)
sigma_figure.add_layout(arrow_height)
sigma_figure.min_border_left = 50


line_div_data = ColumnDataSource(data=dict(
    x1 = [0,layer1_dry.sigma_AEH(layer1_dry.Hl())+layer1_wet.sigma_AEH(layer1_wet.Hl())],
    y1 = [layer1_dry.Hl()+layer1_wet.Hl(),layer1_dry.Hl()+layer1_wet.Hl()],
    x2 = [0,layer1_dry.sigma_AEH(layer1_dry.Hl())+layer2_dry.sigma_AEH(layer2_dry.Hl())],
    y2 = [layer1_dry.Hl()+layer2_dry.Hl(),layer1_dry.Hl()+layer2_dry.Hl()]
))

sigma_figure.line(
    x='x1',
    y='y1',
    source=line_div_data,
    line_color='black',
    line_dash='dashed')

sigma_figure.line(
    x='x2',
    y='y2',
    source=line_div_data,
    line_color='LightSeaGreen',
    line_dash='dashed')


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
###                              MOHR CIRCLE                               ###
##############################################################################

# # Everything for Mohr works correctly ONLY for c=0
# if c == 0:
#     circle_center = (layer_dry.Ja(min(H,abs(ewt)))
#                     + layer_wet.Ja(max(0,H+ewt)))
# else:
#     circle_center = layer_dry.Ja(min(H,abs(ewt)))

l1_circle_center = layer1_dry.Ja(layer1_dry.H) + layer1_wet.Ja(layer1_wet.H)
l2_circle_center = (l1_circle_center
                    + layer2_dry.Ja(layer2_dry.H)
                    + layer2_wet.Ja(layer2_wet.H))


l1_mohr_plot = figure(
    x_axis_label="\u03C3' (kPa)",
    y_axis_label='\u03C4 (kPa)',
    plot_width=200,
    plot_height=200,
    toolbar_location=None,
    x_range=(0, 2.0 * l1_circle_center),
    y_range=(0, 2.03 * l1_circle_center),
    background_fill_alpha=0.1)

l2_mohr_plot = figure(
    x_axis_label="\u03C3' (kPa)",
    y_axis_label='\u03C4 (kPa)',
    plot_width=200,
    plot_height=200,
    toolbar_location=None,
    x_range=(0, 2.0 * l2_circle_center),
    y_range=(0, 2.03 * l2_circle_center),
    background_fill_alpha=0.1)

mohr_plot = figure(
    x_axis_label="\u03C3' (kPa)",
    y_axis_label='\u03C4 (kPa)',
    plot_width=400,
    plot_height=400,
    toolbar_location=None,
    x_range=(0, 2.0 * max(l1_circle_center,l2_circle_center)),
    y_range=(0, 2.03 * max(l1_circle_center,l2_circle_center)),
    background_fill_alpha=0.1)


mohr_line_data = ColumnDataSource(data=dict(
    x_fail1=[0, 1.5 * l1_circle_center],
    y_fail1=[c1, c1 + (np.tan(np.radians(phi1)) * 1.5 * l1_circle_center)],
    x_conj1=[0, 1.9 * l1_circle_center],
    y_conj1=[0, (np.tan(np.radians(beta1) + layer1_dry.theta())
                 ) * 1.9 * l1_circle_center],
    x_fail2=[0, 1.5 * l2_circle_center],
    y_fail2=[c2, c2 + (np.tan(np.radians(phi2)) * 1.5 * l2_circle_center)],
    x_conj2=[0, 1.9 * l2_circle_center],
    y_conj2=[0, (np.tan(np.radians(beta2) + layer2_dry.theta())
                 ) * 1.9 * l2_circle_center],
))

mohr_circle_data = ColumnDataSource(data=dict(
    y=[0],
    center1=[l1_circle_center],
    center2=[l2_circle_center],
    # Circle radius from equation 11
    radius1=[(c1 * (1 / np.tan(np.radians(phi1)))
              + l1_circle_center) * np.sin(np.radians(phi1))],
    radius2=[(c2 * (1 / np.tan(np.radians(phi2)))
              + l2_circle_center) * np.sin(np.radians(phi2))],
))

mohr_plot.line(
    x='x_fail1',
    y='y_fail1',
    source=mohr_line_data,
    line_dash='dashed',
    line_width=2,
    legend='Layer 1 Effective stress M-C failure envelope')

mohr_plot.line(
    x='x_conj1',
    y='y_conj1',
    source=mohr_line_data,
    line_dash='dashed',
    line_width=2,
    color='orange',
    legend='Layer 1 Conjugate stress line')

mohr_plot.line(
    x='x_fail2',
    y='y_fail2',
    source=mohr_line_data,
    line_width=2,
    legend='Layer 2 Effective stress M-C failure envelope')

mohr_plot.line(
    x='x_conj2',
    y='y_conj2',
    source=mohr_line_data,
    line_width=2,
    color='orange',
    legend='Layer 2 Conjugate stress line')

mohr_plot.circle(
    x='center1',
    y='y',
    radius='radius1',
    fill_color=None,
    line_dash='dashed',
    line_width=2,
    color='black',
    source=mohr_circle_data)

mohr_plot.circle(
    x='center2',
    y='y',
    radius='radius2',
    fill_color=None,
    line_width=2,
    color='black',
    source=mohr_circle_data)

mohr_plot.legend.location = "top_left"
mohr_plot.legend.label_text_font_size = "9pt"
mohr_plot.legend.spacing = -5
mohr_plot.legend.padding = 5


# Calculate intersection points for effective stress envelope
x_line1_inter, y_line1_inter = line_circle_intersect(
    h=mohr_circle_data.data['center1'][0],
    k=0,
    r=mohr_circle_data.data['radius1'][0],
    angle=phi1,
    c=c1)

x_line2_inter, y_line2_inter = line_circle_intersect(
    h=mohr_circle_data.data['center2'][0],
    k=0,
    r=mohr_circle_data.data['radius2'][0],
    angle=phi2,
    c=c2)

# Calculate intersection points for conjugate stress line
x_conj1_inter, y_conj1_inter = line_circle_intersect(
    h=mohr_circle_data.data['center1'][0],
    k=0,
    r=mohr_circle_data.data['radius1'][0],
    angle=beta1 + np.degrees(layer1_dry.theta()),
    c=0)

x_conj2_inter, y_conj2_inter = line_circle_intersect(
    h=mohr_circle_data.data['center2'][0],
    k=0,
    r=mohr_circle_data.data['radius2'][0],
    angle=beta2 + np.degrees(layer2_dry.theta()),
    c=0)


# Store in a ColumnDataSource
intersect_data = ColumnDataSource(data=dict(
                x_line1_inter = x_line1_inter,
                y_line1_inter = y_line1_inter,
                x_conj1_inter = x_conj1_inter,
                y_conj1_inter = y_conj1_inter,
                x_line2_inter = x_line2_inter,
                y_line2_inter = y_line2_inter,
                x_conj2_inter = x_conj2_inter,
                y_conj2_inter = y_conj2_inter,
                ))

mohr_plot.circle(x='x_conj1_inter',
                 y='y_conj1_inter',
                 line_color='orange',
                 line_width=2,
                 #line_dash='dashed',
                 fill_color='white',
                 size=7,
                 source=intersect_data)

mohr_plot.circle(x='x_conj2_inter',
                 y='y_conj2_inter',
                 line_color='orange',
                 line_width=2,
                 fill_color='white',
                 size=7,
                 source=intersect_data)

mohr_bold_label_data1 = ColumnDataSource(data=dict(
                        x=[75,72,79.8,68,72,72],
                        y=[285-i*15 for i in range(6)],
                        names=['LAYER 1',
                               'Zw: {:.1f} m'.format(H1),
                               '\u03D5: {:.0f}\u1d52'.format(phi1),
                               '\u03B2+\u03B8: {:.0f}\u1d52'.format(beta1 +
                                            np.degrees(layer1_dry.theta())),
                               "\u03C3'\u03B2: {:.0f} kPa".format(
                                    intersect_data.data['x_conj1_inter'][1]),
                               "\u03C3'\u03B8: {:.0f} kPa".format(
                                    intersect_data.data['x_conj1_inter'][0]),
                               ]))

mohr_bold_labels1 = LabelSet(x='x',
                       y='y',
                       x_units='screen',
                       y_units='screen',
                       text='names',
                       text_font_size='9pt',
                       text_color='black',
                       text_font_style='bold',
                       text_align='left',
                       background_fill_color='white',
                       source=mohr_bold_label_data1)
mohr_plot.add_layout(mohr_bold_labels1)

mohr_bold_label_data2 = ColumnDataSource(data=dict(
                        x=[175,172,179.8,168,172,172],
                        y=[285-i*15 for i in range(6)],
                        names=['LAYER 2',
                               'Zw: {:.1f} m'.format(H1+H2),
                               '\u03D5: {:.0f}\u1d52'.format(phi2),
                               '\u03B2+\u03B8: {:.0f}\u1d52'.format(beta2 +
                                            np.degrees(layer2_dry.theta())),
                               "\u03C3'\u03B2: {:.0f} kPa".format(
                                    intersect_data.data['x_conj2_inter'][1]),
                               "\u03C3'\u03B8: {:.0f} kPa".format(
                                    intersect_data.data['x_conj2_inter'][0]),
                               ]))

mohr_bold_labels2 = LabelSet(x='x',
                       y='y',
                       x_units='screen',
                       y_units='screen',
                       text='names',
                       text_font_size='9pt',
                       text_color='black',
                       text_font_style='bold',
                       text_align='left',
                       background_fill_color='white',
                       source=mohr_bold_label_data2)
mohr_plot.add_layout(mohr_bold_labels2)


beta_error = Label(
    x=170,
    y=100,
    x_units='screen',
    y_units='screen',
    render_mode='css',
    text="\u03B2 + \u03B8 > \u03D5",
    text_font_style='bold',
    text_align='center',
    border_line_width=2,
    text_font_size='14pt',
    text_color='red',
    background_fill_color='white')
mohr_plot.add_layout(beta_error)
beta_error.visible = False



# mohr_sigma_label_data = ColumnDataSource(data=dict(
#                         x=intersect_data.data['x_conj_inter'],
#                         y=intersect_data.data['y_conj_inter'],
#                         names=["\u03C3'\u03B8",
#                                "\u03C3'\u03B2",
#                                ]))
#
# mohr_sigma_labels = LabelSet(x='x',
#                        y='y',
#                        text='names',
#                        text_font_size='9pt',
#                        text_color='black',
#                        text_font_style='bold',
#                        text_align='left',
#                        background_fill_color='white',
#                        background_fill_alpha=0.6,
#                        x_offset=8,
#                        #y_offset=5,
#                        source=mohr_sigma_label_data)
# mohr_plot.add_layout(mohr_sigma_labels)
#
# mohr_plot.legend.location = "top_left"



##############################################################################
###                                 TABLE                                  ###
##############################################################################

columns = [
        TableColumn(field='loc',
                    title='Location'),
        TableColumn(field='zw',
                    title='Zw (m)',
                    formatter=models.NumberFormatter(format='0.00')),
        TableColumn(field='zl',
                    title='Zl (m)',
                    formatter=models.NumberFormatter(format='0.00')),
        TableColumn(field='z',
                    title='Z (m)',
                    formatter=models.NumberFormatter(format='0.00')),
        TableColumn(field='Ja',
                    title='J\u03B1',
                    formatter=models.NumberFormatter(format='0.00')),
        TableColumn(field='Ja_gz',
                    title='J\u03B1/\u03B3Z',
                    formatter=models.NumberFormatter(format='0.00')),
        TableColumn(field='a_a',
                    title='\u03B1\u2090 (deg)',
                    formatter=models.NumberFormatter(format='0.00')),
        TableColumn(field='Ka',
                    title='K\u03B1',
                    formatter=models.NumberFormatter(format='0.000')),
        TableColumn(field='s_a',
                    title="\u03C3'\u03B1 (kPa)",
                    formatter=models.NumberFormatter(format='0.00')),
        TableColumn(field='s_AEH',
                    title="\u03C3'\u1D00\u1D07\u029C (kPa)",
                    formatter=models.NumberFormatter(format='0.00'))
]

# Calculate before storing in ColumnDataSource
zw_ewt = abs(ewt)
zw_l1_top = 0.0001
zw_l1_bot = (layer1_dry.H + layer1_wet.H) * 0.9999
zw_l2_top = (layer1_dry.H + layer1_wet.H) + (layer2_dry.H + layer2_wet.H) * 0.0001
zw_l2_bot = (layer1_dry.H + layer1_wet.H) + (layer2_dry.H + layer2_wet.H)

zl_ewt = abs(ewt)/np.cos(np.radians(omega))
zl_l1_top = layer1_dry.zl(zw_l1_top) + layer1_wet.zl(zw_l1_top)
zl_l1_bot = (layer1_dry.zl(layer1_dry.H) + layer1_wet.zl(layer1_wet.H)) * 0.9999
zl_l2_top = (layer1_dry.zl(layer1_dry.H) + layer1_wet.zl(layer1_wet.H)
            + (layer2_dry.zl(layer2_dry.H) + layer2_wet.zl(layer2_wet.H)) * 0.0001)
zl_l2_bot = (layer1_dry.zl(layer1_dry.H) + layer1_wet.zl(layer1_wet.H)
            + (layer2_dry.zl(layer2_dry.H) + layer2_wet.zl(layer2_wet.H)))

z_ewt = (abs(ewt)
         * (np.cos(np.radians(beta1-omega))
           / (np.cos(np.radians(beta1))*np.cos(np.radians(omega)))))
z_l1_top = layer1_dry.z(zw_l1_top) + layer1_wet.z(zw_l1_top)
z_l1_bot = (layer1_dry.z(layer1_dry.H) + layer1_wet.z(layer1_wet.H)) * 0.9999
z_l2_top = (layer1_dry.z(layer1_dry.H) + layer1_wet.z(layer1_wet.H)
            + (layer2_dry.z(layer2_dry.H) + layer2_wet.z(layer2_wet.H)) * 0.0001)
z_l2_bot = (layer1_dry.z(layer1_dry.H) + layer1_wet.z(layer1_wet.H)
            + (layer2_dry.z(layer2_dry.H) + layer2_wet.z(layer2_wet.H)))

if abs(ewt) < H1:
    Ja_ewt = layer1_dry.Ja(layer1_dry.H)
else:
    Ja_ewt = layer1_dry.Ja(layer1_dry.H) + layer2_dry.Ja(layer2_dry.H)
    Ja_ewt0 = layer2_dry.Ja(layer2_dry.H)
Ja_l1_top = layer1_dry.Ja(zw_l1_top) + layer1_wet.Ja(zw_l1_top)
Ja_l1_bot = (layer1_dry.Ja(layer1_dry.H) + layer1_wet.Ja(layer1_wet.H)) * 0.9999
Ja_l2_top = (layer1_dry.Ja(layer1_dry.H) + layer1_wet.Ja(layer1_wet.H)
            + (layer2_dry.Ja(layer2_dry.H) + layer2_wet.Ja(layer2_wet.H)) * 0.0001)
Ja_l2_bot = (layer1_dry.Ja(layer1_dry.H) + layer1_wet.Ja(layer1_wet.H)
            + (layer2_dry.Ja(layer2_dry.H) + layer2_wet.Ja(layer2_wet.H)))
Ja_l2_top0 = (layer2_dry.Ja(layer2_dry.H) + layer2_wet.Ja(layer2_wet.H)) * 0.0001
Ja_l2_bot0 = layer2_dry.Ja(layer2_dry.H) + layer2_wet.Ja(layer2_wet.H)


if abs(ewt) < H1:
    Ja_gz_ewt = Ja_ewt / (gamma1 * layer1_dry.z(layer1_dry.H))
else:
    Ja_gz_ewt = Ja_ewt0 / (gamma2 * layer2_dry.z(layer2_dry.H))
Ja_gz_l1_top = Ja_l1_top / (gamma1 * layer1_dry.z(zw_l1_top)
                            + (gamma1-gamma_w) * layer1_wet.z(zw_l1_top))
Ja_gz_l1_bot = Ja_l1_bot / ((gamma1 * layer1_dry.z(layer1_dry.H)
                            + (gamma1-gamma_w) * layer1_wet.z(layer1_wet.H))
                            * 0.9999)
Ja_gz_l2_top = Ja_l2_top0 / ((gamma2 * layer2_dry.z(layer2_dry.H)
                            + (gamma2-gamma_w) * layer2_wet.z(layer2_wet.H))
                            * 0.0001)
Ja_gz_l2_bot = Ja_l2_bot0 / (gamma2 * layer2_dry.z(layer2_dry.H)
                            + (gamma2-gamma_w) * layer2_wet.z(layer2_wet.H))

if abs(ewt) < H1:
    a_a_ewt = layer1_dry.alpha_a(layer1_dry.H, degrees=True)
else:
    a_a_ewt = layer2_dry.alpha_a(layer2_dry.H, degrees=True)
a_a_l1_top = layer1_dry.alpha_a(layer1_dry.H, degrees=True)
a_a_l1_bot = layer1_dry.alpha_a(layer1_dry.H, degrees=True)
a_a_l2_top = layer2_dry.alpha_a(layer2_dry.H, degrees=True)
a_a_l2_bot = layer2_dry.alpha_a(layer2_dry.H, degrees=True)

if abs(ewt) < H1:
    Ka_ewt = layer1_dry.Ka(layer1_dry.H)
else:
    Ka_ewt = layer2_dry.Ka(layer2_dry.H)
Ka_l1_top = layer1_dry.Ka(layer1_dry.H)
Ka_l1_bot = layer1_dry.Ka(layer1_dry.H)
Ka_l2_top = layer2_dry.Ka(layer2_dry.H)
Ka_l2_bot = layer2_dry.Ka(layer2_dry.H)

s_a_l1_top = layer1_dry.sigma_a(zw_l1_top) + layer1_wet.sigma_a(zw_l1_top)
s_a_l1_bot = (layer1_dry.sigma_a(layer1_dry.H) + layer1_wet.sigma_a(layer1_wet.H)) * 0.9999
s_a_l2_top = (layer1_dry.sigma_a(layer1_dry.H,
                                 user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001))
              + layer1_wet.sigma_a(layer1_wet.H,
                                   user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001)))
s_a_l2_bot = (s_a_l2_top
              + layer2_dry.sigma_a(layer2_dry.H)
              + layer2_wet.sigma_a(layer2_wet.H))
if abs(ewt) < H1:
    s_a_ewt = layer1_dry.sigma_a(abs(ewt))
else:
    s_a_ewt = s_a_l2_top + layer2_dry.sigma_a(layer2_dry.H)

s_AEH_l1_top = layer1_dry.sigma_AEH(zw_l1_top) + layer1_wet.sigma_AEH(zw_l1_top)
s_AEH_l1_bot = (layer1_dry.sigma_AEH(layer1_dry.H) + layer1_wet.sigma_AEH(layer1_wet.H)) * 0.9999
s_AEH_l2_top = (layer1_dry.sigma_AEH(layer1_dry.H,
                                 user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001),
                                 user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001))
              + layer1_wet.sigma_AEH(layer1_wet.H,
                                   user_Ka=layer2_wet.Ka(layer2_wet.H * 0.0001),
                                   user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001)))
s_AEH_l2_bot = (s_AEH_l2_top
              + layer2_dry.sigma_AEH(layer2_dry.H)
              + layer2_wet.sigma_AEH(layer2_wet.H))
if abs(ewt) < H1:
    s_AEH_ewt = layer1_dry.sigma_AEH(abs(ewt))
else:
    s_AEH_ewt = s_AEH_l2_top + layer2_dry.sigma_AEH(layer2_dry.H)


# This one is for ewt within first layer
if 0 < abs(ewt) < H1:
    source_table = ColumnDataSource(data=dict(
            loc = ['Top of Layer 1',
                   '@ GW Level',
                   'Bottom Layer 1',
                   'Top of Layer 2',
                   'Bottom Layer 2'],
            zw  = [zw_l1_top, zw_ewt, zw_l1_bot, zw_l2_top, zw_l2_bot],
            zl  = [zl_l1_top, zl_ewt, zl_l1_bot, zl_l2_top, zl_l2_bot],
            z   = [z_l1_top, z_ewt, z_l1_bot, z_l2_top, z_l2_bot],
            Ja  = [Ja_l1_top, Ja_ewt, Ja_l1_bot, Ja_l2_top, Ja_l2_bot],
            Ja_gz=[Ja_gz_l1_top, Ja_gz_ewt, Ja_gz_l1_bot, Ja_gz_l2_top, Ja_gz_l2_bot],
            a_a = [a_a_l1_top, a_a_ewt, a_a_l1_bot, a_a_l2_top, a_a_l2_bot],
            Ka  = [Ka_l1_top, Ka_ewt, Ka_l1_bot, Ka_l2_top, Ka_l2_bot],
            s_a = [s_a_l1_top, s_a_ewt, s_a_l1_bot, s_a_l2_top, s_a_l2_bot],
            s_AEH=[s_AEH_l1_top, s_AEH_ewt, s_AEH_l1_bot, s_AEH_l2_top, s_AEH_l2_bot]
    ))

# This one for within second layer
elif H1 < abs(ewt) < H:
    source_table = ColumnDataSource(data=dict(
            loc = ['Top of Layer 1',
                   'Bottom Layer 1',
                   'Top of Layer 2',
                   '@ GW Level',
                   'Bottom Layer 2'],
            zw  = [zw_l1_top, zw_l1_bot, zw_l2_top, zw_ewt, zw_l2_bot],
            zl  = [zl_l1_top, zl_l1_bot, zl_l2_top, zl_ewt, zl_l2_bot],
            z   = [z_l1_top, z_l1_bot, z_l2_top, z_ewt, z_l2_bot],
            Ja  = [Ja_l1_top, Ja_l1_bot, Ja_l2_top, Ja_ewt, Ja_l2_bot],
            Ja_gz=[Ja_gz_l1_top, Ja_gz_l1_bot, Ja_gz_l2_top, Ja_gz_ewt, Ja_gz_l2_bot],
            a_a = [a_a_l1_top, a_a_l1_bot, a_a_l2_top, a_a_ewt, a_a_l2_bot],
            Ka  = [Ka_l1_top, Ka_l1_bot, Ka_l2_top, Ka_ewt, Ka_l2_bot],
            s_a = [s_a_l1_top, s_a_l1_bot, s_a_l2_top, s_a_ewt, s_a_l2_bot],
            s_AEH=[s_AEH_l1_top, s_AEH_l1_bot, s_AEH_l2_top, s_AEH_ewt, s_AEH_l2_bot]
    ))

# This one for all other cases, not showing EWT details
else:
    source_table = ColumnDataSource(data=dict(
            loc = ['Top of Layer 1',
                   'Bottom Layer 1',
                   'Top of Layer 2',
                   'Bottom Layer 2'],
            zw  = [zw_l1_top, zw_l1_bot, zw_l2_top, zw_l2_bot],
            zl  = [zl_l1_top, zl_l1_bot, zl_l2_top, zl_l2_bot],
            z   = [z_l1_top, z_l1_bot, z_l2_top, z_l2_bot],
            Ja  = [Ja_l1_top, Ja_l1_bot, Ja_l2_top, Ja_l2_bot],
            Ja_gz=[Ja_gz_l1_top, Ja_gz_l1_bot, Ja_gz_l2_top, Ja_gz_l2_bot],
            a_a = [a_a_l1_top, a_a_l1_bot, a_a_l2_top, a_a_l2_bot],
            Ka  = [Ka_l1_top, Ka_l1_bot, Ka_l2_top, Ka_l2_bot],
            s_a = [s_a_l1_top, s_a_l1_bot, s_a_l2_top, s_a_l2_bot],
            s_AEH=[s_AEH_l1_top, s_AEH_l1_bot, s_AEH_l2_top, s_AEH_l2_bot]
    ))

data_table = DataTable(source=source_table,
                       columns=columns,
                       width=1000,
                       height=175)


##############################################################################
###                               CALLBACK                                 ###
##############################################################################

# Callback function that updates all plots
def update_plot(attr, old, new):
    omega = omega_slider.value
    beta1 = beta1_slider.value
    beta2 = beta2_slider.value
    phi1 = phi1_slider.value
    if phi1 == 0:
        phi1 = 0.0000001
    phi2 = phi2_slider.value
    if phi2 == 0:
        phi2 = 0.0000001
    H1 = H1_slider.value
    H2 = H2_slider.value
    gamma1 = gamma1_slider.value
    gamma2 = gamma2_slider.value
    kh = kh_slider.value
    kv = kv_slider.value
    ewt = ewt_slider.value
    if ewt == 0:
        ewt = 0.0000001
    H = H1 + H2
    ewt_slider.start = -H
    if H >= 25:
        H1_slider.end = H1
        H2_slider.end = H2

    # NEW Retaining Wall Coordinates
    xB = H / np.tan(angleA)
    yB = H
    xC = xB + 4
    yC = H
    xD = xC + H * np.tan(np.radians(omega))

    # NEW Backfill Coordinates
    # Layer 1
    xE = xC + 15
    yE = yC + 15 * np.tan(np.radians(beta1))
    xM = xC + H1 * np.tan(np.radians(omega))
    yM = H2
    xN = xE
    yN = yM + (15 - H1 * np.tan(np.radians(omega))) * np.tan(np.radians(beta2))

    # Layer 2
    xF = xE

    # NEW EWT Coordinates
    xK = xC + abs(ewt) * np.tan(np.radians(omega))
    yK = H + ewt    # Note: ewt is negative
    xL = xE
    yL = H + ewt    # Note: ewt is negative

    # NEW List of coordinates
    x_wall = [xA, xB, xC, xD]
    y_wall = [yA, yB, yC, yD]
    x_layer1 = [xC, xE, xN, xM]
    y_layer1 = [yC, yE, yN, yM]
    x_layer2 = [xM, xN, xF, xD]
    y_layer2 = [yM, yN, yF, yD]
    x_ewt = [xK, xL, xF, xD]
    y_ewt = [yK, yL, yF, yD]

    # NEW Data Sources
    wall_data.data = dict(x=x_wall, y=y_wall)
    earth_data.data = dict(x1=x_layer1,
                           y1=y_layer1,
                           x2=x_layer2,
                           y2=y_layer2)
    ewt_data.data = dict(x=x_ewt, y=y_ewt)

    wall_label_data.data=dict(
        x=[100, 100, 100, 100, 100],
        y=[110 - i * 15 for i in range(5)],
        names=['WT: {:.1f} m'.format(ewt),
               'H: {:.1f} m'.format(H1+H2),
               '\u03B1h: {:.2f}g'.format(kh),
               '\u03B1v: {:.2f}g'.format(kv),
               '\u03C9: {:.0f}\u1d52'.format(omega)])

    soil_label_data.data=dict(
        x=[10, 14, 10, 10, 14, 14, 10, 14, 10, 10, 14, 14],
        y=[H2+H1/2+1, H2+H1/2+1, H2+H1/2, H2+H1/2-1, H2+H1/2, H2+H1/2-1,
           H2/2+1, H2/2+1, H2/2, H2/2-1, H2/2, H2/2-1],
        names=['LAYER 1',
               '\u03B2: {:.0f}\u1d52'.format(beta1),
               'H: {:.1f} m'.format(H1),
               '\u03C6: {:.0f}\u1d52'.format(phi1),
               'c: {:.0f} kPa'.format(c1),
               '\u03B3: {:.0f} kPa'.format(gamma1),
               'LAYER 2',
               '\u03B2: {:.0f}\u1d52'.format(beta2),
               'H: {:.1f} m'.format(H2),
               '\u03C6: {:.0f}\u1d52'.format(phi2),
               'c: {:.0f} kPa'.format(c2),
               '\u03B3: {:.0f} kPa'.format(gamma2)])



    # STRESS CALCULATIONS
    # New sep class calculations
    if abs(ewt) >= H:
        layer1_dry = sep(kh, kv, omega, beta1, phi1, gamma1, c1, H1)
        layer2_dry = sep(kh, kv, omega, beta2, phi2, gamma2, c2, H2)
        layer1_wet = sep(kh, kv, omega, beta1, phi1, gamma1-gamma_w, c1, 0.00000001)
        layer2_wet = sep(kh, kv, omega, beta2, phi2, gamma2-gamma_w, c2, 0.00000001)
    elif H1 <= abs(ewt) < H:
        layer1_dry = sep(kh, kv, omega, beta1, phi1, gamma1, c1, H1)
        temp_var = abs(ewt)-H1
        layer2_dry = sep(kh, kv, omega, beta2, phi2, gamma2, c2,
                         temp_var if temp_var !=0 else 0.00000001)
        layer1_wet = sep(kh, kv, omega, beta1, phi1, gamma1-gamma_w, c1, 0.00000001)
        layer2_wet = sep(kh, kv, omega, beta2, phi2, gamma2-gamma_w, c2, H+ewt)
    else:
        layer1_dry = sep(kh, kv, omega, beta1, phi1, gamma1, c1, abs(ewt))
        layer2_dry = sep(kh, kv, omega, beta2, phi2, gamma2, c2, 0.00000001)
        layer1_wet = sep(kh, kv, omega, beta1, phi1, gamma1-gamma_w, c1, H1+ewt)
        layer2_wet = sep(kh, kv, omega, beta2, phi2, gamma2-gamma_w, c2, H2)

    all_layer_Hl = (layer1_dry.Hl() + layer2_dry.Hl()
                    + layer1_wet.Hl() + layer2_wet.Hl())

    # New stress plot calculations
    sigma_l1_y_dry = np.arange(0.0001, layer1_dry.Hl(), 0.1)
    sigma_l1_y_wet = np.arange(0.0001, layer1_wet.Hl(), 0.1)
    sigma_l1_x_dry = layer1_dry.sigma_AEH(sigma_l1_y_dry)
    sigma_l1_x_wet = (layer1_wet.sigma_AEH(sigma_l1_y_wet)
                      + layer1_dry.sigma_AEH(layer1_dry.Hl()))
    sigma_l2_y_dry = np.arange(0.0001, layer2_dry.Hl(), 0.1)
    sigma_l2_y_wet = np.arange(0.0001, layer2_wet.Hl(), 0.1)
    sigma_l2_x_dry = (layer2_dry.sigma_AEH(sigma_l2_y_dry)
                      + layer1_dry.sigma_AEH(
                                layer1_dry.Hl(),
                                user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001),
                                user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001)))
    sigma_l2_x_wet = (layer2_wet.sigma_AEH(sigma_l2_y_wet)
                      + layer2_dry.sigma_AEH(layer2_dry.Hl())
                      + layer1_dry.sigma_AEH(
                                layer1_dry.Hl(),
                                user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001),
                                user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001))
                      + layer1_wet.sigma_AEH(
                                layer1_wet.Hl(),
                                user_Ka=layer2_wet.Ka(layer2_wet.H * 0.0001),
                                user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001)))


    y_sigma_all = sigma_l1_y_dry.tolist()
    x_sigma_all = sigma_l1_x_dry.tolist()
    y_sigma_all.extend((sigma_l1_y_wet + layer1_dry.Hl()).tolist())
    x_sigma_all.extend(sigma_l1_x_wet.tolist())
    y_sigma_all.extend((sigma_l2_y_dry + layer1_dry.Hl()).tolist())
    x_sigma_all.extend(sigma_l2_x_dry.tolist())
    y_sigma_all.extend((sigma_l2_y_wet + layer1_dry.Hl() + layer1_wet.Hl()
                        + layer2_dry.Hl()).tolist())
    x_sigma_all.extend(sigma_l2_x_wet.tolist())
    y_sigma_all.extend([all_layer_Hl])
    x_sigma_all.extend([0])

    sigma_data.data=dict(x=x_sigma_all, y=y_sigma_all)

    # Update sigma plot ranges
    sigma_figure.y_range.start = 0.99 * all_layer_Hl
    sigma_figure.y_range.end = all_layer_Hl - 30

    # Update arrow
    load_height_top = all_layer_Hl - all_layer_Hl / 3
    load_height_bot = all_layer_Hl / 3

    arrow_data.data=dict(
        x0 = [max(x_sigma_all)],
        y0 = [load_height_top],
        y1 = [load_height_top])

    force = PolyArea(x_sigma_all, y_sigma_all) * (np.cos(np.radians(omega))**2)

    arrow_load.x = 0.27 * max(x_sigma_all)
    arrow_load.y = load_height_top
    arrow_load.text = '{:.0f} kN'.format(force)
    arrow_height.x = 0.27 * max(x_sigma_all)
    arrow_height.y = load_height_top
    arrow_height.text='@ ~{:.2f} m'.format(load_height_bot)

    # Line div data
    line_div_data.data=dict(
        x1 = [0,layer1_dry.sigma_AEH(layer1_dry.Hl())+layer1_wet.sigma_AEH(layer1_wet.Hl())],
        y1 = [layer1_dry.Hl()+layer1_wet.Hl(),layer1_dry.Hl()+layer1_wet.Hl()],
        x2 = [0,layer1_dry.sigma_AEH(layer1_dry.Hl())+layer2_dry.sigma_AEH(layer2_dry.Hl())],
        y2 = [layer1_dry.Hl()+layer2_dry.Hl(),layer1_dry.Hl()+layer2_dry.Hl()]
    )


    # MOHR CIRCLE
    l1_circle_center = (layer1_dry.Ja(layer1_dry.H)
                        + layer1_wet.Ja(layer1_wet.H))
    l2_circle_center = (l1_circle_center
                        + layer2_dry.Ja(layer2_dry.H)
                        + layer2_wet.Ja(layer2_wet.H))

    mohr_plot.x_range.end = 2.0 * max(l1_circle_center,l2_circle_center)
    mohr_plot.y_range.end = 2.03 * max(l1_circle_center,l2_circle_center)

    mohr_line_data.data=dict(
        x_fail1=[0, 1.5 * l1_circle_center],
        y_fail1=[c1, c1 + (np.tan(np.radians(phi1)) * 1.5 * l1_circle_center)],
        x_conj1=[0, 1.9 * l1_circle_center],
        y_conj1=[0, (np.tan(np.radians(beta1) + layer1_dry.theta())
                     ) * 1.9 * l1_circle_center],
        x_fail2=[0, 1.5 * l2_circle_center],
        y_fail2=[c2, c2 + (np.tan(np.radians(phi2)) * 1.5 * l2_circle_center)],
        x_conj2=[0, 1.9 * l2_circle_center],
        y_conj2=[0, (np.tan(np.radians(beta2) + layer2_dry.theta())
                     ) * 1.9 * l2_circle_center],
    )

    mohr_circle_data.data=dict(
        y=[0],
        center1=[l1_circle_center],
        center2=[l2_circle_center],
        # Circle radius from equation 11
        radius1=[(c1 * (1 / np.tan(np.radians(phi1)))
                  + l1_circle_center) * np.sin(np.radians(phi1))],
        radius2=[(c2 * (1 / np.tan(np.radians(phi2)))
                  + l2_circle_center) * np.sin(np.radians(phi2))],
    )

    # Update line/circle intersection
    x_line1_inter, y_line1_inter = line_circle_intersect(
        h=mohr_circle_data.data['center1'][0],
        k=0,
        r=mohr_circle_data.data['radius1'][0],
        angle=phi1,
        c=c1)

    x_line2_inter, y_line2_inter = line_circle_intersect(
        h=mohr_circle_data.data['center2'][0],
        k=0,
        r=mohr_circle_data.data['radius2'][0],
        angle=phi2,
        c=c2)

    x_conj1_inter, y_conj1_inter = line_circle_intersect(
        h=mohr_circle_data.data['center1'][0],
        k=0,
        r=mohr_circle_data.data['radius1'][0],
        angle=beta1 + np.degrees(layer1_dry.theta()),
        c=0)

    x_conj2_inter, y_conj2_inter = line_circle_intersect(
        h=mohr_circle_data.data['center2'][0],
        k=0,
        r=mohr_circle_data.data['radius2'][0],
        angle=beta2 + np.degrees(layer2_dry.theta()),
        c=0)

    intersect_data.data=dict(
                    x_line1_inter = x_line1_inter,
                    y_line1_inter = y_line1_inter,
                    x_conj1_inter = x_conj1_inter,
                    y_conj1_inter = y_conj1_inter,
                    x_line2_inter = x_line2_inter,
                    y_line2_inter = y_line2_inter,
                    x_conj2_inter = x_conj2_inter,
                    y_conj2_inter = y_conj2_inter,
                    )

    # IF clause for NO INTERSECTION
    if (np.isnan(intersect_data.data['x_conj1_inter'][0]) or
            np.isnan(intersect_data.data['x_conj2_inter'][0])):
        sigma_error.visible = True
        beta_error.visible = True
        mohr_plot.background_fill_color = 'red'

    else:
        sigma_error.visible = False
        beta_error.visible = False
        mohr_plot.background_fill_color = None

    # Update Mohr plot label data

    mohr_bold_label_data1.data=dict(
                            x=[75,72,79.8,68,72,72],
                            y=[285-i*15 for i in range(6)],
                            names=['LAYER 1',
                                   'Zw: {:.1f} m'.format(H1),
                                   '\u03D5: {:.0f}\u1d52'.format(phi1),
                                   '\u03B2+\u03B8: {:.0f}\u1d52'.format(beta1 +
                                                np.degrees(layer1_dry.theta())),
                                   "\u03C3'\u03B2: {:.0f} kPa".format(
                                        intersect_data.data['x_conj1_inter'][1]),
                                   "\u03C3'\u03B8: {:.0f} kPa".format(
                                        intersect_data.data['x_conj1_inter'][0]),
                                   ])

    mohr_bold_label_data2.data=dict(
                            x=[175,172,179.8,168,172,172],
                            y=[285-i*15 for i in range(6)],
                            names=['LAYER 2',
                                   'Zw: {:.1f} m'.format(H1+H2),
                                   '\u03D5: {:.0f}\u1d52'.format(phi2),
                                   '\u03B2+\u03B8: {:.0f}\u1d52'.format(beta2 +
                                                np.degrees(layer2_dry.theta())),
                                   "\u03C3'\u03B2: {:.0f} kPa".format(
                                        intersect_data.data['x_conj2_inter'][1]),
                                   "\u03C3'\u03B8: {:.0f} kPa".format(
                                        intersect_data.data['x_conj2_inter'][0]),
                                   ])

    # UPDATE TABLE
        # Calculate before storing in ColumnDataSource
    zw_ewt = abs(ewt)
    zw_l1_top = 0.0001
    zw_l1_bot = (layer1_dry.H + layer1_wet.H) * 0.9999
    zw_l2_top = (layer1_dry.H + layer1_wet.H) + (layer2_dry.H + layer2_wet.H) * 0.0001
    zw_l2_bot = (layer1_dry.H + layer1_wet.H) + (layer2_dry.H + layer2_wet.H)

    zl_ewt = abs(ewt)/np.cos(np.radians(omega))
    zl_l1_top = layer1_dry.zl(zw_l1_top) + layer1_wet.zl(zw_l1_top)
    zl_l1_bot = (layer1_dry.zl(layer1_dry.H) + layer1_wet.zl(layer1_wet.H)) * 0.9999
    zl_l2_top = (layer1_dry.zl(layer1_dry.H) + layer1_wet.zl(layer1_wet.H)
                + (layer2_dry.zl(layer2_dry.H) + layer2_wet.zl(layer2_wet.H)) * 0.0001)
    zl_l2_bot = (layer1_dry.zl(layer1_dry.H) + layer1_wet.zl(layer1_wet.H)
                + (layer2_dry.zl(layer2_dry.H) + layer2_wet.zl(layer2_wet.H)))

    z_ewt = (abs(ewt)
             * (np.cos(np.radians(beta1-omega))
               / (np.cos(np.radians(beta1))*np.cos(np.radians(omega)))))
    z_l1_top = layer1_dry.z(zw_l1_top) + layer1_wet.z(zw_l1_top)
    z_l1_bot = (layer1_dry.z(layer1_dry.H) + layer1_wet.z(layer1_wet.H)) * 0.9999
    z_l2_top = (layer1_dry.z(layer1_dry.H) + layer1_wet.z(layer1_wet.H)
                + (layer2_dry.z(layer2_dry.H) + layer2_wet.z(layer2_wet.H)) * 0.0001)
    z_l2_bot = (layer1_dry.z(layer1_dry.H) + layer1_wet.z(layer1_wet.H)
                + (layer2_dry.z(layer2_dry.H) + layer2_wet.z(layer2_wet.H)))

    if abs(ewt) < H1:
        Ja_ewt = layer1_dry.Ja(layer1_dry.H)
    else:
        Ja_ewt = layer1_dry.Ja(layer1_dry.H) + layer2_dry.Ja(layer2_dry.H)
        Ja_ewt0 = layer2_dry.Ja(layer2_dry.H)
    Ja_l1_top = layer1_dry.Ja(zw_l1_top) + layer1_wet.Ja(zw_l1_top)
    Ja_l1_bot = (layer1_dry.Ja(layer1_dry.H) + layer1_wet.Ja(layer1_wet.H)) * 0.9999
    Ja_l2_top = (layer1_dry.Ja(layer1_dry.H) + layer1_wet.Ja(layer1_wet.H)
                + (layer2_dry.Ja(layer2_dry.H) + layer2_wet.Ja(layer2_wet.H)) * 0.0001)
    Ja_l2_bot = (layer1_dry.Ja(layer1_dry.H) + layer1_wet.Ja(layer1_wet.H)
                + (layer2_dry.Ja(layer2_dry.H) + layer2_wet.Ja(layer2_wet.H)))
    Ja_l2_top0 = (layer2_dry.Ja(layer2_dry.H) + layer2_wet.Ja(layer2_wet.H)) * 0.0001
    Ja_l2_bot0 = layer2_dry.Ja(layer2_dry.H) + layer2_wet.Ja(layer2_wet.H)

    if abs(ewt) < H1:
        Ja_gz_ewt = Ja_ewt / (gamma1 * layer1_dry.z(layer1_dry.H))
    else:
        Ja_gz_ewt = Ja_ewt0 / (gamma2 * layer2_dry.z(layer2_dry.H))
    Ja_gz_l1_top = Ja_l1_top / (gamma1 * layer1_dry.z(zw_l1_top)
                                + (gamma1-gamma_w) * layer1_wet.z(zw_l1_top))
    Ja_gz_l1_bot = Ja_l1_bot / ((gamma1 * layer1_dry.z(layer1_dry.H)
                                + (gamma1-gamma_w) * layer1_wet.z(layer1_wet.H))
                                * 0.9999)
    Ja_gz_l2_top = Ja_l2_top0 / ((gamma2 * layer2_dry.z(layer2_dry.H)
                                + (gamma2-gamma_w) * layer2_wet.z(layer2_wet.H))
                                * 0.0001)
    Ja_gz_l2_bot = Ja_l2_bot0 / (gamma2 * layer2_dry.z(layer2_dry.H)
                                + (gamma2-gamma_w) * layer2_wet.z(layer2_wet.H))

    if abs(ewt) < H1:
        a_a_ewt = layer1_dry.alpha_a(layer1_dry.H, degrees=True)
    else:
        a_a_ewt = layer2_dry.alpha_a(layer2_dry.H, degrees=True)
    a_a_l1_top = layer1_dry.alpha_a(layer1_dry.H, degrees=True)
    a_a_l1_bot = layer1_dry.alpha_a(layer1_dry.H, degrees=True)
    a_a_l2_top = layer2_dry.alpha_a(layer2_dry.H, degrees=True)
    a_a_l2_bot = layer2_dry.alpha_a(layer2_dry.H, degrees=True)

    if abs(ewt) < H1:
        Ka_ewt = layer1_dry.Ka(layer1_dry.H)
    else:
        Ka_ewt = layer2_dry.Ka(layer2_dry.H)
    Ka_l1_top = layer1_dry.Ka(layer1_dry.H)
    Ka_l1_bot = layer1_dry.Ka(layer1_dry.H)
    Ka_l2_top = layer2_dry.Ka(layer2_dry.H)
    Ka_l2_bot = layer2_dry.Ka(layer2_dry.H)

    s_a_l1_top = layer1_dry.sigma_a(zw_l1_top) + layer1_wet.sigma_a(zw_l1_top)
    s_a_l1_bot = (layer1_dry.sigma_a(layer1_dry.H) + layer1_wet.sigma_a(layer1_wet.H)) * 0.9999
    s_a_l2_top = (layer1_dry.sigma_a(layer1_dry.H,
                                     user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001))
                  + layer1_wet.sigma_a(layer1_wet.H,
                                       user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001)))
    s_a_l2_bot = (s_a_l2_top
                  + layer2_dry.sigma_a(layer2_dry.H)
                  + layer2_wet.sigma_a(layer2_wet.H))
    if abs(ewt) < H1:
        s_a_ewt = layer1_dry.sigma_a(abs(ewt))
    else:
        s_a_ewt = s_a_l2_top + layer2_dry.sigma_a(layer2_dry.H)

    s_AEH_l1_top = layer1_dry.sigma_AEH(zw_l1_top) + layer1_wet.sigma_AEH(zw_l1_top)
    s_AEH_l1_bot = (layer1_dry.sigma_AEH(layer1_dry.H) + layer1_wet.sigma_AEH(layer1_wet.H)) * 0.9999
    s_AEH_l2_top = (layer1_dry.sigma_AEH(layer1_dry.H,
                                     user_Ka=layer2_dry.Ka(layer2_dry.H * 0.0001),
                                     user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001))
                  + layer1_wet.sigma_AEH(layer1_wet.H,
                                       user_Ka=layer2_wet.Ka(layer2_wet.H * 0.0001),
                                       user_alpha_a=layer2_dry.alpha_a(layer2_dry.H * 0.0001)))
    s_AEH_l2_bot = (s_AEH_l2_top
                  + layer2_dry.sigma_AEH(layer2_dry.H)
                  + layer2_wet.sigma_AEH(layer2_wet.H))
    if abs(ewt) < H1:
        s_AEH_ewt = layer1_dry.sigma_AEH(abs(ewt))
    else:
        s_AEH_ewt = s_AEH_l2_top + layer2_dry.sigma_AEH(layer2_dry.H)


    # This one is for ewt within first layer
    if 0 < abs(ewt) < H1:
        source_table.data=dict(
                loc = ['Top of Layer 1',
                       '@ GW Level',
                       'Bottom Layer 1',
                       'Top of Layer 2',
                       'Bottom Layer 2'],
                zw  = [zw_l1_top, zw_ewt, zw_l1_bot, zw_l2_top, zw_l2_bot],
                zl  = [zl_l1_top, zl_ewt, zl_l1_bot, zl_l2_top, zl_l2_bot],
                z   = [z_l1_top, z_ewt, z_l1_bot, z_l2_top, z_l2_bot],
                Ja  = [Ja_l1_top, Ja_ewt, Ja_l1_bot, Ja_l2_top, Ja_l2_bot],
                Ja_gz=[Ja_gz_l1_top, Ja_gz_ewt, Ja_gz_l1_bot, Ja_gz_l2_top, Ja_gz_l2_bot],
                a_a = [a_a_l1_top, a_a_ewt, a_a_l1_bot, a_a_l2_top, a_a_l2_bot],
                Ka  = [Ka_l1_top, Ka_ewt, Ka_l1_bot, Ka_l2_top, Ka_l2_bot],
                s_a = [s_a_l1_top, s_a_ewt, s_a_l1_bot, s_a_l2_top, s_a_l2_bot],
                s_AEH=[s_AEH_l1_top, s_AEH_ewt, s_AEH_l1_bot, s_AEH_l2_top, s_AEH_l2_bot]
        )

    # This one for within second layer
    elif H1 < abs(ewt) < H:
        source_table.data=dict(
                loc = ['Top of Layer 1',
                       'Bottom Layer 1',
                       'Top of Layer 2',
                       '@ GW Level',
                       'Bottom Layer 2'],
                zw  = [zw_l1_top, zw_l1_bot, zw_l2_top, zw_ewt, zw_l2_bot],
                zl  = [zl_l1_top, zl_l1_bot, zl_l2_top, zl_ewt, zl_l2_bot],
                z   = [z_l1_top, z_l1_bot, z_l2_top, z_ewt, z_l2_bot],
                Ja  = [Ja_l1_top, Ja_l1_bot, Ja_l2_top, Ja_ewt, Ja_l2_bot],
                Ja_gz=[Ja_gz_l1_top, Ja_gz_l1_bot, Ja_gz_l2_top, Ja_gz_ewt, Ja_gz_l2_bot],
                a_a = [a_a_l1_top, a_a_l1_bot, a_a_l2_top, a_a_ewt, a_a_l2_bot],
                Ka  = [Ka_l1_top, Ka_l1_bot, Ka_l2_top, Ka_ewt, Ka_l2_bot],
                s_a = [s_a_l1_top, s_a_l1_bot, s_a_l2_top, s_a_ewt, s_a_l2_bot],
                s_AEH=[s_AEH_l1_top, s_AEH_l1_bot, s_AEH_l2_top, s_AEH_ewt, s_AEH_l2_bot]
        )

    # This one for all other cases, not showing EWT details
    else:
        source_table.data=dict(
                loc = ['Top of Layer 1',
                       'Bottom Layer 1',
                       'Top of Layer 2',
                       'Bottom Layer 2'],
                zw  = [zw_l1_top, zw_l1_bot, zw_l2_top, zw_l2_bot],
                zl  = [zl_l1_top, zl_l1_bot, zl_l2_top, zl_l2_bot],
                z   = [z_l1_top, z_l1_bot, z_l2_top, z_l2_bot],
                Ja  = [Ja_l1_top, Ja_l1_bot, Ja_l2_top, Ja_l2_bot],
                Ja_gz=[Ja_gz_l1_top, Ja_gz_l1_bot, Ja_gz_l2_top, Ja_gz_l2_bot],
                a_a = [a_a_l1_top, a_a_l1_bot, a_a_l2_top, a_a_l2_bot],
                Ka  = [Ka_l1_top, Ka_l1_bot, Ka_l2_top, Ka_l2_bot],
                s_a = [s_a_l1_top, s_a_l1_bot, s_a_l2_top, s_a_l2_bot],
                s_AEH=[s_AEH_l1_top, s_AEH_l1_bot, s_AEH_l2_top, s_AEH_l2_bot]
        )






# Sliders
H1_slider = Slider(
    start=1,
    end=24,
    step=0.5,
    value=6,
    title='Layer 1 height, H\u2081, (m)')
H2_slider = Slider(
    start=1,
    end=24,
    step=0.5,
    value=9,
    title='Layer 2 height, H\u2082, (m)')
omega_slider = Slider(
    start=0,
    end=30,
    step=1,
    value=15,
    title='Wall inclination, \u03C9, (deg.)')
beta1_slider = Slider(
    start=-30,
    end=30,
    step=1,
    value=beta1,
    title='Surface slope, \u03B2\u2081, (deg.)')
beta2_slider = Slider(
    start=-30,
    end=30,
    step=1,
    value=beta2,
    title='Interface slope, \u03B2\u2082, (deg.)')
phi1_slider = Slider(
    start=0,
    end=45,
    step=1,
    value=30,
    title='Layer 1 internal friction, \u03C6\u2081, (deg.)')
phi2_slider = Slider(
    start=0,
    end=45,
    step=1,
    value=40,
    title='Layer 2 internal friction, \u03C6\u2082, (deg.)')
gamma1_slider = Slider(
    start=16,
    end=25,
    step=1,
    value=23,
    title='Layer 1 unit weight, \u03B3\u2081, (kN/m\u00B3)')
gamma2_slider = Slider(
    start=16,
    end=25,
    step=1,
    value=20,
    title='Layer 2 unit weight, \u03B3\u2082, (kN/m\u00B3)')
kh_slider = Slider(
    start=0,
    end=0.4,
    step=0.05,
    value=0.15,
    title='Horizontal, kh')
kv_slider = Slider(
    start=-0.4,
    end=0.4,
    step=0.05,
    value=0,
    title='Vertical, kv')
ewt_slider = Slider(
    start=-H,
    end=0,
    step=0.1,
    value=ewt,
    orientation="vertical",
    direction='rtl',
    show_value=False,
    height=340,
    width=25,
    callback_throttle=0,
    tooltips=False,
    bar_color='#e8f7f6')



# Attach the callback to the 'value' property of slider
H1_slider.on_change('value', update_plot)
H2_slider.on_change('value', update_plot)
omega_slider.on_change('value', update_plot)
beta1_slider.on_change('value', update_plot)
beta2_slider.on_change('value', update_plot)
phi1_slider.on_change('value', update_plot)
phi2_slider.on_change('value', update_plot)
gamma1_slider.on_change('value', update_plot)
gamma2_slider.on_change('value', update_plot)
kh_slider.on_change('value', update_plot)
kv_slider.on_change('value', update_plot)
ewt_slider.on_change('value', update_plot)


page_header = Div(text=open(join(dirname(__file__), "page_header.html")).read(),
                  width=1000)
page_footer = Div(text=open(join(dirname(__file__), "page_footer.html")).read(),
                  width=1000)

wall_controls = [omega_slider,beta1_slider,beta2_slider]
wall_inputs = widgetbox(*wall_controls, width=200)

seismic_controls = [kh_slider,kv_slider]
seismic_inputs = widgetbox(*seismic_controls, width=180)

layer1_controls = [H1_slider,phi1_slider,gamma1_slider]
layer1_inputs = widgetbox(*layer1_controls, width=250)

layer2_controls = [H2_slider,phi2_slider,gamma2_slider]
layer2_inputs = widgetbox(*layer2_controls, width=250)

# The layout function replaces the row and column functions
page_layout = layout([
                [page_header],
                [Div(text="<h3>Wall Geometry:</h3>", width=200),
                 Div(text="<h3>Seismic Coefficients:</h3>", width=180),
                 Div(text="<h3>Soil Properties:</h3>", width=200)],
                [wall_inputs, seismic_inputs, layer1_inputs, layer2_inputs],
                [Div(text="<hr>", width=1000)],
                [Div(text="<h4></h4>", width=40),
                 Div(text="<h4>Retaining Wall and Backfill Geometry</h4>",
                     width=325),
                 Div(text="<h4>GWT<br>&nbsp;(m)</h4>", width=85),
                 Div(text="<h4>Horizontal Pseudo-Static<br>"
                          "Lateral Earth Pressure</h4>",
                     width=215),
                 Div(text="<h4>Mohr's circle with failure envelopes at "
                          "depth, Zw,<br>from the top of wall surface</h4>",
                     width=350)],
                [wall_plot, ewt_slider, Spacer(width=10),
                 sigma_figure, Spacer(width=10),mohr_plot],
                [Div(text="<h3>Calculated values at several vertical depths "
                          "from the top of wall surface, Zw</h3>",
                     width=600)],
                [data_table],
                [page_footer]
])


curdoc().add_root(page_layout)
curdoc().title = "SEP Calculator"

### test with:
### bokeh serve --show sep-calculator-l2.py

### run forever on server with:
### nohup bokeh serve sep-calculator-l2.py --allow-websocket-origin cue3.engineering.nyu.edu:5010 --host cue3.engineering.nyu.edu:5010 --port 5010
