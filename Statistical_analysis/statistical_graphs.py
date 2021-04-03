# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 00:50:43 2021

@author: claum
"""
# ===================================================================
# Preambolo - Inizio
# ===================================================================
# import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
# ===================================================================
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# ===================================================================
stat1  = PdfPages('figura1.pdf')
stat2  = PdfPages('figura2.pdf')
stat3  = PdfPages('figura3.pdf')
stat4  = PdfPages('figura4.pdf')
stat5  = PdfPages('figura5.pdf')
stat6  = PdfPages('figura6.pdf')
stat7  = PdfPages('figura7.pdf')
stat8  = PdfPages('figura8.pdf')
stat9  = PdfPages('figura9.pdf')
stat10 = PdfPages('figura10.pdf')
stat11 = PdfPages('figura11.pdf')
stat12 = PdfPages('figura12.pdf')
stat13 = PdfPages('figura13.pdf')
stat14 = PdfPages('figura14.pdf')
stat15 = PdfPages('figura15.pdf')
stat16 = PdfPages('figura16.pdf')
statEW = PdfPages('figura17.pdf')
# ===================================================================
# AIRCRAFT DATA
# ===================================================================
# MAXIMUM TAKE OFF WEIGHT
MTOW      = np.array([22800, 22930, 22850, 19950, 27987]).reshape((-1, 1))
MTOW1     = np.array([22800, 22930, 22850, 19950, 27987])
# MAXIMUM LANDING WEIGHT
MLW       = np.array([22350, 22930, 22000, 19500, 27442])
# BASIC OPERATING WEIGHT
BOC       = np.array([13450, 13995, 14200, 12520, 18418])
# USEFUL LOAD
PAYLOAD   = np.array([7550, 7000, 5500, 6080, 7356])
PAYLOAD1  = np.array([7550, 7000, 5500, 6080, 7356]).reshape((-1, 1))
# MAX TAKE OFF POWER
MTOP      = np.array([5500, 5500, 8304, 5000, 5071]) # [SHP]
FW        = np.array([5000, 5346, 5423, 4123, 5318]) # FUEL WEIGHT
EW        = MTOW - PAYLOAD1                          # EMPTY WEIGHT
print('-------------------')
print('EMPTY WEIGHT: \n', EW, '\n')
label     = ["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"]
# TIME TO CLIMB AT FL170
TTC_FL170 = np.array([17.50, 21.00, 15.00, 19.00, 22.62]) # [SHP]
# ===================================================================
# BOC VS MTOW GRAPH
# ===================================================================
# Y = b0 + b1*x + b2*x**2   b2=0 
model1   = LinearRegression().fit(MTOW, BOC)
r_sq1    = model1.score(MTOW, BOC)
print('-------------------')
print('Coefficient of determination:', r_sq1)
print('Intercept:', model1.intercept_)
print('Slope:', model1.coef_)
print('-------------------')
# ===================================================================
fig1  = plt.figure()
plt.scatter(MTOW, BOC)
for i, txt in enumerate(label):
    plt.annotate(txt, (MTOW[i], BOC[i]), size="xx-small")
plt.plot(MTOW, model1.intercept_+model1.coef_*MTOW, color="red")
plt.xlabel(r'\textsc{mtow} - $[kg]$')              # x-label to the axes.
plt.ylabel(r'\textsc{boc} - $[kg]$')     # y-label to the axes.
plt.title(r'Basic Operating Weight vs Maximum Take Off Weight') 
plt.grid(True, linestyle='-.')
plt.show()
stat1.savefig(fig1)
stat1.close()
# ===================================================================
# EW VS MTOW GRAPH
# ===================================================================
model2   = LinearRegression().fit(MTOW, PAYLOAD)
r_sq2    = model1.score(MTOW, PAYLOAD)
print('-------------------')
print('Coefficient of determination:', r_sq2)
print('Intercept:', model2.intercept_)
print('Slope:', model2.coef_)
print('-------------------')
# ===================================================================
fig2  = plt.figure()
plt.scatter(MTOW, PAYLOAD)
for i, txt in enumerate(label):
    plt.annotate(txt, (MTOW[i], PAYLOAD[i]), size="xx-small")
plt.plot(MTOW, model2.intercept_+model2.coef_*MTOW, color="red")
plt.xlabel(r'\textsc{mtow} - $[kg]$')              # x-label to the axes.
plt.ylabel(r'\textsc{Payload} - $[kg]$')     # y-label to the axes.
plt.title(r'Payload Weight vs Maximum Take Off Weight') 
plt.grid(True, linestyle='-.')
plt.show()
stat2.savefig(fig2)
stat2.close()
# ==================================================================
# # ==================================================================
# fig1  = plt.figure()
# plt.plot(x_c, eta_c, label="Codice BEMT", color='blue')
# plt.plot(x_x, eta_x, label="XROTOR", marker='o', linestyle='none', color='red')
# plt.plot(x_ex, eta_ex, label="Dati NACA", color='green')
# plt.xlabel(r'Advance ratio - $J$')              # x-label to the axes.
# plt.ylabel(r'Propeller efficiency - $\eta$')     # y-label to the axes.
# plt.title(r'Propeller efficiency - $\eta = \eta(J)$')    # Title to the axes.
# plt.legend(loc="upper right")
# # plt.legend(loc="upper right", prop={"size" : 6})
# plt.ylim(0, 1)
# plt.xlim(0, 1)
# plt.grid(True, linestyle='-.')
# plt.show()
# cr1.savefig(fig1)
# cr1.close()
# # ==================================================================
# ==================================================================
data1 = {r'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        r'MTOW':[22800, 22930, 22850, 19950, 27987]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data1)
ax1 = df.plot.bar(x=r'Aircraft', y=r'MTOW', rot=0, color=color, legend=False)
plt.title(r'\textsc{mtow} in kilograms') 
plt.grid(True, linestyle='-.')
plt.show()
ax1.figure.savefig('figura3.pdf')
# stat2.close()
# ==================================================================
data2 = {r'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        r'BOC':[13450, 13995, 14200, 12520, 18418]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data2)
ax2 = df.plot.bar(x=r'Aircraft', y=r'BOC', rot=0, color=color, legend=False)
plt.title(r'Basic Operating Weight in kilograms') 
plt.grid(True, linestyle='-.')
plt.show()
ax2.figure.savefig('figura4.pdf')
# stat3.savefig(fig3)
# stat3.close()
# ==================================================================
data3 = {r'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        r'Payload':[7550, 7000, 5500, 6080, 7356]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data3)
ax3 = df.plot.bar(x=r'Aircraft', y=r'Payload', rot=0, color=color, legend=False)
plt.title(r'Payload Weight in kilograms') 
plt.grid(True, linestyle='-.')
plt.show()
ax3.figure.savefig('figura5.pdf')
# stat4.savefig(fig4)
# stat4.close()
# ==================================================================
data4 = {'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        'fuel':[5000, 5346, 5423, 4123, 5318]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data4)
ax4 = df.plot.bar(x='Aircraft', y='fuel', rot=0, color=color, legend=False)
plt.title(r'Fuel Weight in kilograms') 
plt.grid(True, linestyle='-.')
plt.show()
ax4.figure.savefig('figura6.pdf')
# ==================================================================
# ==================================================================
# Maximum Take Off Power vs MTOW
# ==================================================================
# ==================================================================
model3   = LinearRegression().fit(MTOW, MTOP)
r_sq3    = model1.score(MTOW, MTOP)
print('Coefficient of determination:', r_sq3)
print('Intercept:', model3.intercept_)
print('Slope:', model3.coef_)
# ===================================================================
fig3  = plt.figure()
plt.scatter(MTOW, MTOP)
for i, txt in enumerate(label):
    plt.annotate(txt, (MTOW[i], MTOP[i]), size="xx-small")
plt.plot(MTOW, model3.intercept_+model3.coef_*MTOW, color="red")
plt.xlabel(r'\textsc{mtow} - $[kg]$')              # x-label to the axes.
plt.ylabel(r'Max Take Off Power - $[SHP]$')     # y-label to the axes.
# plt.xlim(0.23, 0.35)
plt.title(r'Maximum Take Off Power vs Maximum Take Off Weight') 
plt.grid(True, linestyle='-.')
plt.show()
stat7.savefig(fig3)
stat7.close()
# ==================================================================
data5 = {'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        'Surface':[61, 78.3, 55.7, 70, 64]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data5)
ax5 = df.plot.bar(x='Aircraft', y='Surface', rot=0, color=color, legend=False)
plt.title(r'Wing Surfaces in squared meters') 
plt.grid(True, linestyle='-.')
plt.show()
ax5.figure.savefig('figura8.pdf')
# ==================================================================
# ==================================================================
data6 = {'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        'Wing span':[27.05, 30.63, 24.76, 29, 28.40]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data6)
ax6 = df.plot.bar(x='Aircraft', y='Wing span', rot=0, color=color, legend=False)
plt.title(r'Wing span in meters') 
plt.grid(True, linestyle='-.')
plt.show()
ax6.figure.savefig('figura9.pdf')
# ==================================================================
# ==================================================================
data7 = {'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        'Aspect ratio':[11.995, 11.982, 11.00642, 12.0100, 11.600]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data7)
ax7 = df.plot.bar(x='Aircraft', y='Aspect ratio', rot=0, color=color, legend=False)
plt.title(r'Aspect ratio') 
plt.grid(True, linestyle='-.')
plt.show()
ax7.figure.savefig('figura10.pdf')
# ==================================================================
# ==================================================================
data8 = {'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        'Fuselage length':[27.166, 26.000, 27.300, 25.250, 32.800]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data8)
ax8 = df.plot.bar(x='Aircraft', y='Fuselage length', rot=0, color=color, legend=False)
plt.title(r'Fuselage length in meters') 
plt.grid(True, linestyle='-.')
plt.show()
ax8.figure.savefig('figura11.pdf')
# ==================================================================
# ==================================================================
data9 = {'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        'Range':[1403.816, 1320.476, 1527.900, 2780, 1877.928]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data9)
ax9 = df.plot.bar(x='Aircraft', y='Range', rot=0, color=color, legend=False)
plt.title(r'Range (\textsc{opt. fl} and max payload) in kilometers') 
plt.grid(True, linestyle='-.')
plt.show()
ax9.figure.savefig('figura12.pdf')
# ==================================================================
# ==================================================================
data10 = {'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        'MaxCruiseSpeed':[275, 271, 370, 250, 290]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data10)
ax10 = df.plot.bar(x='Aircraft', y='MaxCruiseSpeed', rot=0, color=color, legend=False)
plt.title(r'Max cruise speed in knots ($95\,\%$ \textsc{mtow}, \textsc{opt. fl})') 
plt.grid(True, linestyle='-.')
plt.show()
ax10.figure.savefig('figura13.pdf')
# ==================================================================
# ==================================================================
data11 = {'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        'LandingLength':[915, 1128, 982, 1017, 1290]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data11)
ax11 = df.plot.bar(x='Aircraft', y='LandingLength', rot=0, color=color, legend=False)
plt.title(r'Landing field length (\textsc{easa} Air Ops) in meters') 
plt.grid(True, linestyle='-.')
plt.show()
ax11.figure.savefig('figura14.pdf')
# ==================================================================
# ==================================================================
data12 = {'Aircraft':["ATR 72 - 600", "BAe ATP", "SAAB 2000", "Fokker 50", "Q - 400"],\
        'LandingSpeed':[113, 105, 116, 97, 100]}
color = ["blue", "green", "yellow", "red", "black"]
df = pd.DataFrame(data12)
ax12 = df.plot.bar(x='Aircraft', y='LandingSpeed', rot=0, color=color, legend=False)
plt.title(r'Landing speed in knots') 
plt.grid(True, linestyle='-.')
plt.show()
ax12.figure.savefig('figura15.pdf')
# ==================================================================
# ===================================================================
# BOC VS MTOW GRAPH
# ===================================================================
# Y = b0 + b1*x + b2*x**2   b2=0 
model4   = LinearRegression().fit(MTOW, TTC_FL170)
r_sq4    = model4.score(MTOW, TTC_FL170)
print('Coefficient of determination:', r_sq4)
print('Intercept:', model4.intercept_)
print('Slope:', model4.coef_)
# ===================================================================
fig4  = plt.figure()
plt.scatter(MTOW, TTC_FL170)
for i, txt in enumerate(label):
    plt.annotate(txt, (MTOW[i], TTC_FL170[i]), size="xx-small")
plt.plot(MTOW, model4.intercept_+model4.coef_*MTOW, color="red")
plt.xlabel(r'\textsc{mtow} - $[kg]$')              # x-label to the axes.
plt.ylabel(r'\textsc{time to climb} - $[min]$')     # y-label to the axes.
plt.title(r'Time to climb vs Maximum Take Off Weight') 
plt.grid(True, linestyle='-.')
plt.show()
stat16.savefig(fig4)
stat16.close()
# ===================================================================
# ===================================================================
# EMPTY WEIGHT vs MAXIMUM TAKE OFF WEIGHT
# ===================================================================
# DICIASSETTESIMA FIGURA
# ===================================================================
modelEW   = LinearRegression().fit(EW, MTOW1)
r_sqEW    = modelEW.score(EW, MTOW1)
print('-------------------')
print('Coefficient of determination:', r_sqEW)
print('Intercept:', modelEW.intercept_)
print('Slope:', modelEW.coef_)
print('-------------------')
# ===================================================================
figEW  = plt.figure()
plt.scatter(EW, MTOW1)
for i, txt in enumerate(label):
    plt.annotate(txt, (EW[i], MTOW[i]), size="xx-small")
plt.plot(EW, modelEW.intercept_+modelEW.coef_*(EW), color="red")
plt.xlabel(r'\textsc{Empty weight} - $[kg]$')              # x-label to the axes.
plt.ylabel(r'\textsc{Max Take Off Weight}- $[kg]$')     # y-label to the axes.
# plt.xlim(0.23, 0.35)
plt.title(r'Empty Weight vs Max Take Off Weight') 
plt.grid(True, linestyle='-.')
plt.show()
statEW.savefig(figEW)
statEW.close()
# ===================================================================
# ===================================================================
# USEFUL FUNCTIONS
# ===================================================================
def breguet_range(V, SFC, E, MTOW, FW):
    """
    A function to implement Breguet range equation

    Parameters
    ----------
    V    : FLOAT
        Flight velocity.
    SFC  : FLOAT
        Specific fuel consumption.
    E    : FLOAT
        Aerodynamic efficiency.
    MTOW : FLOAT
        Maximum Take Off Weight.
    FW   : FLOAT
        Fuel Weight fraction.

    Returns
    -------
    R : FLOAT
        Range.
    """
    return (V/SFC)*E*np.log((MTOW)/(MTOW-FW))
# ===================================================================
def breguet_range_inv(V, SFC, E, Range):
    """
    A function to implement Breguet range equation

    Parameters
    ----------
    V     : FLOAT
        Flight velocity.
    SFC   : FLOAT
        Specific fuel consumption.
    E     : FLOAT
        Aerodynamic efficiency.
    Range : FLOAT
        Range.

    Returns
    -------
    FW/MTOW : FLOAT
        Fuel fraction weight.
    """
    x = -(Range/V)*(SFC/E)
    return 1 - exp(x)