#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# Support module generated by PAGE version 4.14
# In conjunction with Tcl version 8.6
#    Aug 09, 2018 03:09:08 PM


import sys
import cStringIO
import numpy as np
from time import strftime
from datetime import datetime
from multiprocessing import Process

# First we should tell matplotlib that we are using a tkinter based GUI.
# Otherwise there might be backend conflicts which would cause the main
# window to freeze until matplotlib plots are closed. MARS import should
# happen after this setting too.
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import mars
from stl_tools import numpy2stl

import Adv_settings_support

try:
    from Tkinter import *
    import tkFileDialog
except ImportError:
    from tkinter import *

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True
    
def set_Tk_var():
    global N
    N = StringVar()
    N.set(128)
    global M
    M = StringVar()
    M.set(128)
    global dx
    dx = StringVar()
    dx.set(1)
    global dy
    dy = StringVar()
    dy.set(1)
    global attempts
    attempts = StringVar()
    attempts.set(1)
    global rmsheight
    rmsheight = StringVar()
    rmsheight.set(1)
    global skewness
    skewness = StringVar()
    skewness.set(0.5)
    global kurtosis
    kurtosis = StringVar()
    kurtosis.set(3.5)
    global acf_type
    acf_type = StringVar()
    acf_type.set("Exponential")
    global n
    n = StringVar()
    n.set(16)
    global m
    m = StringVar()
    m.set(16)
    global phi
    phi = StringVar()
    phi.set(0.0)

def update_box(box, stream, root):
    ### --- Function to add text to the output box in the GUI --- ###
    if len(stream.getvalue())>0:
            box.configure(state=NORMAL)
            time_stamp = datetime.now().strftime('%H:%M')
            new_text = "["+time_stamp + "] " + stream.getvalue()
            box.insert(END, new_text)
            box.configure(state=DISABLED)
            box.see(END)
            root.update_idletasks()
            stream.truncate(0)


def generate_surface(box, root, button_list):
    global N, M, dx, dy, attempts, rmsheight, skewness, kurtosis, acf_type
    global n, m, phi
    global hmap

    ### TO DO: Change this once you add support for advanced settings
    cutoff = 1e-5

    # Capture the stdoutput to a variable which then can be used to 
    # update the box that reports program updates.
    stdout_ = sys.stdout
    stream = cStringIO.StringIO()
    sys.stdout = stream

    # Try to convert all variables from the GUI to the respective data types
    try:
        N_int = int(N.get())
        M_int = int(M.get())
        dx_float = np.float64(dx.get())
        dy_float = np.float64(dy.get())
        n_int = int(n.get())
        m_int = int(m.get())
        cutoff_float = np.float64(cutoff)
        phi_float = np.float64(phi.get())
        attempts_int = int(attempts.get())
        rmsheight_float = np.float64(rmsheight.get())
        skewness_float = np.float64(skewness.get())
        kurtosis_float = np.float64(kurtosis.get())
        display_contour = bool(Adv_settings_support.display_contour.get())
        pl_residual = bool(Adv_settings_support.plot_residual.get())
        mean_float = np.float32(Adv_settings_support.johnson_mean.get())
        var_float = np.float32(Adv_settings_support.johnson_var.get())
        tolerance_float = np.float32(Adv_settings_support.krylov_tol.get())
        max_iter_int = int(Adv_settings_support.krylov_max_iter.get())
        krylov_method = Adv_settings_support.krylov_method.get()
    except:
        print ("ERROR: Could not convert input data to the correct type. Check for any non-number inputs.")
        update_box(box, stream, root)
        return

    # Button disabling at start of surface generation and cursor changing
    for entry in button_list:
        entry.configure(state=DISABLED)
    root.config(cursor="watch")


    # Repeat the surface generation steps however many times needed.
    for i in range(attempts_int):

        s= mars.surface(
                n_int, m_int, N_int, M_int, \
                cutoff_float, dx_float, \
                dy_float, phi_float)

        update_box(box, stream, root)

        print ("Attempting to generate a surface with Sk=" + skewness.get()\
               + " and Ku=" + kurtosis.get())
        update_box(box, stream, root)

        # Step 1: Specify ACF
        if (acf_type.get() == "Exponential"):
            acf= s.acf()
            if display_contour:
                s.plot_contour(acf)
        else:
            print ("ERROR: Unknown ACF type selected.")
            update_box(box, stream, root)
            restore_state(root, button_list)
            return

        # Step 2: Assemble & solve nonlinear system of equations
        guess= s.f0()
        alpha= s.krylov(method=krylov_method, plot=pl_residual,\
                        tolerance=tolerance_float, maxiter=max_iter_int)
        update_box(box, stream, root)
        if len(alpha) == 1 and alpha[0] == -1:
            sys.stdout = stdout_
            restore_state(root, button_list)
            return

        # Step 3: Generate a random number matrix
        rescaled_skew, rescaled_kurt = s.rescale(
                alpha, skewness_float, kurtosis_float)

        update_box(box, stream, root)

        if (np.isnan(rescaled_skew)) and (np.isnan(rescaled_kurt)):
            print("ERROR: The program can't generate the surface because of the input Skew and Kurtosis.")
            update_box(box,stream, root)
            sys.stdout = stdout_
            restore_state(root, button_list)
            return
        rand= s.johnson_eta(rescaled_skew, rescaled_kurt, mean=mean_float,\
                            var=var_float)
        update_box(box, stream, root)

        # Step 4: Generate the heightmap
        hmap= s.heightmap(alpha,rand)
        update_box(box, stream, root)



    # Restore stdout to initial state
    sys.stdout = stdout_
    restore_state(root, button_list)

    # Spawn a separate process that would plot the surface.
    # This way the main GUI would become interactive much quicker
    plot_process = Process(target = lambda: plot_surface(N_int, M_int, np.transpose(hmap)))
    plot_process.start()
    plot_process.join()


def plot_surface(N, M, hmap):
    # Plot the final surface with a colourbar
    plt.figure(2)
    plt.clf()
    plt.pcolormesh(range(N), range(M), hmap, cmap=plt.cm.RdYlBu_r)
    plt.title("Generated surface")
    plt.xlabel("Streamwise points")
    plt.ylabel("Spanwise points")
    plt.axis("tight")
    plt.colorbar()
    plt.show()


def restore_state(root, button_list):
    # Button enabling at end of surface generation and cursor resetting
    for entry in button_list:
        entry.configure(state=NORMAL)
    root.config(cursor="")

def save_as(self):
    global hmap, N, M, dx, dy
    self.f = tkFileDialog.asksaveasfilename(   
        filetypes = (("dat file", "*.dat"),    
                     ("CSV file", "*.csv"),
                     ("STL file", "*.stl")))
    file_ext = self.f[-4:]
    with open(self.f, 'w') as outputFile:
        if (file_ext) == ".dat":
            outputFile.write(hmap)
        elif (file_ext) == ".csv":
            N_int = int(N.get())
            M_int = int(M.get())
            dx_fl = np.float64(dx.get())
            dy_fl = np.float64(dy.get())
            outputFile.write("X,Y,height\n")
            for i in range(N_int):
                for j in range(M_int):
                    try:
                        line = str(i*dx_fl)+","+ str(j*dy_fl) +","+ str(hmap[i][j])+"\n"
                        outputFile.write(line)
                    except:
                        print (i, j)
        elif (file_ext) == ".stl":
            numpy2stl(hmap, self.f, solid=True, scale=5)


def init(top, gui, *args, **kwargs):
    global w, top_level, root
    w = gui
    top_level = top
    root = top

def destroy_window():
    # Function which closes the window.
    global top_level
    top_level.destroy()
    top_level = None

if __name__ == '__main__':
    import GUI
    GUI.vp_start_gui()


