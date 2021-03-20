#######  Hodgkin-Huxley Squid Giant Axon Action potential Simulation

from dearpygui.core import *
from dearpygui.simple import *
import numpy as np

add_additional_font("glacial_font.otf", 18)
################# Initialize global variables

Gna = 120.
Gk = 36.
Io = 0.04

dt = 0.025    # time step

#################### Hodgkin-Huxley Equations for rate constants

def alpha_m(v):
    return 0.1*(v+40)/(1-np.exp(-.1*(v+40)))

def beta_m(v):
    return 4.*np.exp(-.0556*(v+65))

def alpha_h(v):
    return 0.07*np.exp(-.05*(v+65))

def beta_h(v):
    return 1/(np.exp(-.1*(v+35)))

def alpha_n(v):
    return 0.01*(v+55)/(1-np.exp(-.1*(v+55)))

def beta_n(v):
    return 0.125*np.exp(-.125*(v+65.))


#################### Simulation based on given conductances for Na (Gna) and Potassium (Gk), and current injection (Io)

def solve(Gna, Gk, Io):
    dm = 0.
    dh = 0.
    m = 0.
    h = 0.
    n = 0.


    v = -65
    vlist = []
    Istep = 0

    m = alpha_m(v)/(alpha_m(v)+beta_m(v))
    h = alpha_h(v)/(alpha_h(v)+beta_h(v))
    n = alpha_n(v)/(alpha_n(v)+beta_n(v))

    ina = m*m*m*h*Gna*(v-55)
    ik = n*n*n*n*Gk*(v+77)
    el = -65 + ((ina + ik)/.2)

    inalist = []
    iklist = []

    for i in range(600):
        if i > 100:
            Istep = -1*Io
        dm = dt*( ((1.-m)*alpha_m(v)) - (m*beta_m(v))  )
        dh = dt*( ((1.-h)*alpha_h(v)) - (h*beta_h(v))  )
        dn = dt*( ((1.-n)*alpha_n(v)) - (n*beta_n(v))  )
        h += dh
        m += dm
        n += dn

        ina = m*m*m*h*Gna*(v-55)
        ik = n*n*n*n*Gk*(v+77)

        inalist.append(ina)
        iklist.append(ik)

        dv = -dt*( (0.2*(v-el)) + Istep + ina + ik)
        v += dv

        answer = ((v + 65)*1.75) - 65.

        vlist.append(answer)

    return vlist, inalist, iklist



######### Plot timecourse for potential and currents based on chosen conductances

def set_conductances(sender,data):        
    Gna = get_value("G-Na")
    Gk = get_value("G-K")
    Io = get_value("Io")
    potential, ina, ik = solve(Gna,Gk,Io)

    time = []
    for i in range(len(potential)):
        time.append(i*dt)

    add_line_series("Potential", "Potential", time, potential, weight=2, color=[255,255, 255, 100])
    add_line_series("Current", "I-Na", time, ina, weight=2, color=[0,255, 0, 100])
    add_line_series("Current", "I-K", time, ik, weight=2, color=[255,0, 0, 100])

def close(sender,data):             # Stop DPG
    stop_dearpygui()


######### Initialize potential array    

potential = solve(Gna, Gk, Io)


######################### DPG ##########################

with window("Action Potential", width=600, height=800):
    set_main_window_size(800,750)

    add_text("Hodgkin-Huxley Squid Giant Axon Action Potential Simulation")
    add_same_line(spacing=300)
    add_button('Close',callback=close)
    add_separator()
    add_spacing(count=2)



    add_slider_float("G-Na",default_value=120,max_value=200,callback=set_conductances)
    add_slider_float("G-K",default_value=40,max_value=100,callback=set_conductances)
    add_slider_float("Io",default_value=.04,max_value=20,callback=set_conductances)


    add_plot("Potential", height=300)
    set_plot_ylimits("Potential",ymin=-80,ymax=40)
    add_plot("Current", height=300)

start_dearpygui(primary_window="Action Potential")