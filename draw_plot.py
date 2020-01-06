import matplotlib
import plotly.tools as tls
import plotly.plotly as py
import plotly.graph_objs as go
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from tkinter import messagebox


def drawing_plots_by_call(calc_real):

    if len(calc_real) > 0:

        def parse(correct_realisations_sio2):
            vp_ar = []
            h_ar = []
            sio2_ar = []
            ro_ar = []

            for h in correct_realisations_sio2:
              #  print(h, correct_realisations_sio2[h])
                for single_real_arr in correct_realisations_sio2[h]:
                    for single_real in single_real_arr:
                        #print(h, single_real)
                        vp = single_real["efVp"]
                        vs = single_real["efVs"]
                        sio2 = single_real["sio2"]
                        ro = single_real["ro"]
                        vp_ar.append(vp)
                        h_ar.append(-h)
                        sio2_ar.append(sio2)
                        ro_ar.append(ro)
            drawing_plot(h_ar, vp_ar, sio2_ar, ro_ar)


        def drawing_plot(h, vp, sio2, ro):
            plt.plot(vp, h, 'ro', ro, h, 'bs')
            #plt.axis([5.2, 50, 0, 30])
            plt.ylabel('h')
            plt.xlabel("vp")
            y = h
            x = vp

            x1 = sio2
            x2 = ro

        #надо почитать про np и polyfit
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            plt.plot(x, p(x), "r--")

            #z1 = np.polyfit(x1, y, 1)
            #p1 = np.poly1d(z1)
            #plt.plot(x1, p1(x1), "r--")

            z2 = np.polyfit(x2, y, 1)
            p2 = np.poly1d(z2)
            plt.plot(x2, p2(x2), "r--")

            plt.show()
        parse(calc_real)
    else:
        messagebox.showwarning("Warning", "No correct realizations were founded")

def draw_profiles(calc_real):
    x = np.random.randn(500)
    data = [go.Histogram(x=x)]

    py.plot(data, filename='basic histogram')
