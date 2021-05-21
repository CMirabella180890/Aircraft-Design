import math
from sympy import *
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import json

# WING PLANFORM
# INPUT
S = 895.             #ft^2
AR = 12.5
#S1_ratio = 0.38
#S2_ratio = 1-S1_ratio
kink = 0.34
sweep_le2 = 2/(180/math.pi)        #dopo kink in rad
tap_rat = 5/9

#       WING CHARACTERISTICs
b_wing = sqrt(AR*S)                 #ft
b_kink = kink*b_wing/2            #ft
print('apertura alare {} ft'.format(b_wing))
print('apertura al kink {} ft'.format(b_kink))

A = np.array([[b_kink+(b_wing/2-b_kink)*(1+tap_rat)/2,0],[-tap_rat,1]], dtype=float)
b = np.array([S/2,0], dtype=float)
z = np.linalg.solve(A,b)
c_tip = z[1]
c_root = z[0]
print("corda all'estremità {} ft".format(c_tip));
print('corda alla radice {} ft'.format(c_root));
#taper_ratio = c_tip/c_root
#print("rapporto rastremazione {} ".format(taper_ratio))

Sa = c_root*b_kink/2
Sb = Sa + c_root*(b_wing/2-b_kink)/2
Sc = c_tip*(b_wing/2-b_kink)/2
print(Sa+Sb+Sc)

k1 = 2*Sa/S
k2 = 2*Sb/S
k3 = 2*Sc/S


print(k1)
print(k2)
print(k3)

# corda media aerodinamica
c_1 = c_root
c_2 = 2/3*c_root*(1+tap_rat+tap_rat**2)/(1+tap_rat)
S1 = c_root*b_kink*2;
S2 = (b_wing/2-b_kink)*c_root*(1+tap_rat)
c_mean = (S1*c_1+S2*c_2)/(S1+S2)
print(S1+S2)
print(c_mean)
Ac2=(c_tip-c_root)/(b_wing/2-b_kink)
y_c_mean = b_kink+(c_mean-c_root)/(Ac2)
xle_c_mean = (c_mean-c_root)/(Ac2)*math.tan(sweep_le2)
print(y_c_mean)
print(xle_c_mean)
x_tip_le = (b_wing/2-b_kink)*math.tan(sweep_le2)
x_tip_te = x_tip_le + c_tip

x_p1 = x_tip_le-4/(b_wing)*((b_wing/2+b_kink)*(x_tip_le)/2)
#print(x_p1)
x_p2 = +x_tip_te + 4/b_wing*((b_kink+b_wing/2)*math.fabs(c_root-x_tip_te)/2)
#print(x_p2)

#Seq = (math.fabs(x_p2-x_p1)+c_tip)*(b_wing/4)
#print(Seq)

c_root_eq=x_p2-x_p1

tap_rat_eq=c_tip/c_root_eq

sweep_eq = math.atan((x_tip_le-x_p1)/(b_wing/2))
#print('angolo di sweep al l.e. equivalente è {} gradi'.format(sweep_eq*180/math.pi))

sweep_eq_c4 = math.atan(math.tan(sweep_eq)-(1-tap_rat_eq)/(AR*(1+tap_rat_eq)))
print('angolo di sweep a c/4 equivalente è {} gradi'.format(sweep_eq_c4*180/math.pi))

Ac_eq = (c_tip-c_root_eq)/(b_wing/2)
y_c_mean_eq = (c_mean-c_root_eq)/Ac_eq
xle_c_mean_eq = y_c_mean_eq*math.tan(sweep_eq)+x_p1

# vettori per planform ala
x_wing = [0,0,x_tip_le,x_tip_te,c_root,c_root]
y_wing = [0,b_kink,b_wing/2,b_wing/2,b_kink,0]
x_wing_eq = [x_p1,x_tip_le,x_tip_te,x_p2]
y_wing_eq = [0,b_wing/2,b_wing/2,0]

thickness = 18*k1+18*k2+15*k3


# profilo medio aerodin

Naca615_0FL = np.genfromtxt( 'D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\polari profili\\Naca615_0fl.dat',dtype=float)
Naca615_250FL = np.genfromtxt( 'D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\polari profili\\Naca615_250fl.dat',dtype=float)
Naca618_0FL = np.genfromtxt( 'D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\polari profili\\Naca618_0fl.dat',dtype=float)
Naca618_250FL = np.genfromtxt( 'D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\polari profili\\Naca618_250fl.dat',dtype=float)
#print(Naca615_0FL[:,0])


alfa_0FL = Naca618_0FL[:,0]
Cl_mean_0FL = k1*Naca618_0FL[:,1] + k2*Naca618_0FL[:,1] + k3*Naca615_0FL[:,1]
Cm_c4_mean_0FL = k1*Naca618_0FL[:,3] + k2*Naca618_0FL[:,3] + k3*Naca615_0FL[:,3]
Cd_mean_0FL = k1*Naca618_0FL[:,2] + k2*Naca618_0FL[:,2] + k3*Naca615_0FL[:,2]

alfa_250FL = Naca618_250FL[:,0]
Cl_mean_250FL = k1*Naca618_250FL[:,1] + k2*Naca618_250FL[:,1] + k3*Naca615_250FL[:,1]
Cm_c4_mean_250FL = k1*Naca618_250FL[:,3] + k2*Naca618_250FL[:,3] + k3*Naca615_250FL[:,3]
Cd_mean_250FL = k1*Naca618_250FL[:,2] + k2*Naca618_250FL[:,2] + k3*Naca615_250FL[:,2]


# calcolo Cla del profilo medio 0FL e 250FL
Clalpha_mean = [None]*len(range(5,len(alfa_0FL)-30,1))
j = 0
for i in range(5,len(alfa_0FL)-30,1):
    #print(i)
    Clalpha_mean[j]=(Cl_mean_0FL[i]-Cl_mean_0FL[i-1])/(alfa_0FL[i]-alfa_0FL[i-1])
    j=j+1
    #print(Clalpha_mean[j-1])
Clalpha_mean_0FL = np.mean(Clalpha_mean)
print('il cl alpha del profilo medio in con 0 FL è {}'.format(Clalpha_mean_0FL))

Clalpha_mean2 = [None]*len(range(5,len(alfa_250FL)-30,1))
j = 0
for i in range(5,len(alfa_250FL)-30,1):
    #print(i)
    Clalpha_mean2[j] = (Cl_mean_250FL[i]-Cl_mean_250FL[i-1])/(alfa_250FL[i]-alfa_250FL[i-1])
    j=j+1
    #print(Clalpha_mean[j-1])
Clalpha_mean_250FL = np.mean(Clalpha_mean2)
print('il cl alpha del profilo medio in con 250 FL è {}'.format(Clalpha_mean_0FL))



alfa_0l_0FL = - Cl_mean_0FL[np.where(alfa_0FL == 0.)]/Clalpha_mean_0FL
print('alfa zero lift profilo medio 0 FL{}'.format(alfa_0l_0FL))

alfa_0l_250FL = - Cl_mean_250FL[np.where(alfa_0FL == 0.)]/Clalpha_mean_250FL
print('alfa zero lift profilo medio 250 FL {}'.format(alfa_0l_250FL))


#calcolo dei valori di fine linearita Cl* alfa*
tol = 0.05
x1 = alfa_0FL[np.where(alfa_0FL == 0.)]
x2 = alfa_0FL[-1]
step = (alfa_0FL[-1])/((alfa_0FL[-1]-alfa_0FL[np.where(alfa_0FL == 0.)])*100)
x = np.arange(x1,x2,step)
y = Clalpha_mean_0FL*(x-alfa_0l_0FL)

a = np.where(alfa_0FL == 0.)
a = int(a[0])
#print(a)
f2 = interp1d(alfa_0FL[a:], Cl_mean_0FL[a:], kind = 'cubic')
y2 = f2(x)

diff = 0.
i = 0
while diff < tol:
    diff = math.fabs(y[i]-y2[i])
    i=i+1
a_star_0FL = x[i]
Cl_star_OFL = y[i]

print('alfa^* 0FL {} Cl^* {}'.format(a_star_0FL,Cl_star_OFL))


# a* Cl* 250FL

x1 = alfa_250FL[np.where(alfa_250FL == 0.)]
x2 = alfa_250FL[-1]
step = (alfa_250FL[-1])/((alfa_250FL[-1]-alfa_250FL[np.where(alfa_250FL == 0.)])*100)
x = np.arange(x1,x2,step)
y = Clalpha_mean_250FL*(x-alfa_0l_250FL)

a = np.where(alfa_250FL == 0.)
a = int(a[0])
#print(a)
f2 = interp1d(alfa_250FL[a:], Cl_mean_250FL[a:], kind = 'cubic')
y2 = f2(x)

diff = 0.
i = 0
while diff < tol:
    diff = math.fabs(y[i]-y2[i])
    i=i+1
a_star_250FL = x[i]
Cl_star_25OFL = y[i]

print('alfa^* 250FL {} Cl^* {}'.format(a_star_250FL,Cl_star_25OFL))



# file output
#with open('D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\profili.json') as f:
    # python object to be appended
    #y = '{"c_root":c_root}'
# parsing JSON string:
    #z = json.loads(f)

# appending the data
   #z.update(y)

# the result is a JSON string:
    #print(json.dumps(z))

#f.close()
a_file = open("D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\profili.json", "r")
json_object = json.load(a_file)
a_file.close()

c_root_eq = str(c_root_eq)

y = {"spessore medio":str(thickness),"c_root":c_root,"c_tip":c_tip,"c_root_ala_eq":c_root_eq,"Sweept_ala_eq":sweep_eq,"Sweept_ala_eq_c4":sweep_eq_c4}
z={"corda_media_aerodin":str(c_mean),"Xle_tip":str(x_tip_le),"Cl_alpha_OFl":str(Clalpha_mean_0FL),\
     "Cl_alpha_250FL":str(Clalpha_mean_250FL),"alpha_0lift_0FL":str(alfa_0l_0FL),"alpha_0lift_250FL":str(alfa_0l_250FL),"alfa star 0FL":str(a_star_0FL),\
   "alfa star 250FL":str(a_star_250FL),"Cl star 0FL":str(Cl_star_OFL),"Cl star 250FL":str(Cl_star_25OFL)}
json_object.update(y)
json_object.update(z)
print(json_object)
with open("D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\output_2d.json", 'w') as outfile:
    json.dump(json_object, outfile)





# parte grafica
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


stat1 = PdfPages('D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\risultati immagini\\Cl_mean_0FL.pdf')
stat2 = PdfPages('D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\risultati immagini\\Cl_mean_250FL.pdf')
stat3 = PdfPages('D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\risultati immagini\\Cl_Cd_mean_0FL.pdf')
stat4 = PdfPages('D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\risultati immagini\\Cl_Cd_mean_250FL.pdf')
stat5 = PdfPages('D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\risultati immagini\\Cm_c4_mean_0FL.pdf')
stat6 = PdfPages('D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\risultati immagini\\Cm_c4_mean_250FL.pdf')
Wing = PdfPages('D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\risultati immagini\\wing_planform.pdf')
Wing_eq = PdfPages('D:\\AIRCRAFT DESIGN\\PROGETTO\\wing DABCM98\\risultati immagini\\wing_eq.pdf')

# cl mean airfoil 0 ft
figCl0 = plt.figure()
plt.plot(alfa_0FL, Cl_mean_0FL,color='k',markersize=12)
plt.xlabel(r'\textsc{$\alpha^\circ$}')
plt.ylabel(r'\textsc{$C_l$}')
plt.title(r" $C_l$ mean airfoil at $h = 0$ ft and $V = V_{stall}$ ")
plt.axhline(y=0,color='k',markersize=8)
plt.axvline(x=0,color='k',markersize=8)
plt.grid(True, linestyle='-.')
plt.show()
stat1.savefig(figCl0)
stat1.close()

# cl_cd mean airfoil 0 ft
figClCd0 = plt.figure()
plt.plot(Cd_mean_0FL*10000,Cl_mean_0FL,color='k',markersize=12)
plt.ylabel(r'\textsc{$C_l$}')
plt.xlabel(r'\textsc{$C_d$ drag count}')
plt.title(r" Mean airfoil drag polar at $h = 0$ ft and $V = V_{stall}$ ")
plt.axhline(y=0,color='k',markersize=8)
plt.axvline(x=0,color='k',markersize=8)
plt.grid(True, linestyle='-.')
plt.show()
stat3.savefig(figClCd0)
stat3.close()

# cl mean airfoil 25000ft
figCl250 = plt.figure()
plt.plot(alfa_250FL, Cl_mean_250FL,color='k',markersize=12)
plt.xlabel(r'\textsc{$\alpha^\circ$}')
plt.ylabel(r'\textsc{$C_l$}')
plt.title(r" $C_l$ mean airfoil at $h = 25000$ ft and $V = 1.1V_{cruise}$ ")
plt.axhline(y=0,color='k',markersize=8)
plt.axvline(x=0,color='k',markersize=8)
plt.grid(True, linestyle='-.')
plt.show()
stat2.savefig(figCl250)
stat2.close()

# cl_cd mean airfoil 25000 ft
figClCd250 = plt.figure()
plt.plot(Cd_mean_250FL*10000,Cl_mean_250FL,color='k',markersize=12)
plt.ylabel(r'\textsc{$C_l$}')
plt.xlabel(r'\textsc{$C_d$ drag count}')
plt.title(r" Mean airfoil drag polar at $h = 25000$ ft and $V = 1.1V_{cruise}$ ")
plt.axhline(y=0,color='k',markersize=8)
plt.axvline(x=0,color='k',markersize=8)
plt.grid(True, linestyle='-.')
plt.show()
stat4.savefig(figClCd250)
stat4.close()

# Cm c/4 mean airfoil a 0 ft
figCm0 = plt.figure()
plt.plot(Cl_mean_0FL,Cm_c4_mean_0FL,color='k',markersize=12)
plt.xlabel(r'\textsc{$C_l$}')
plt.ylabel(r'$C_{m,c/4}$')
plt.title(r'$C_{m,c/4}$ mean airfoil at $h = 25000$ ft and $V = V_{stall}$')
plt.axhline(y=0,color='k',markersize=8)
plt.axvline(x=0,color='k',markersize=8)
plt.grid(True, linestyle='-.')
plt.show()
stat5.savefig(figCm0)
stat5.close()

# Cm c/4 mean airfoil a 25000 ft
figCm250 = plt.figure()
plt.plot(Cl_mean_250FL,Cm_c4_mean_250FL,color='k',markersize=12)
plt.xlabel(r'\textsc{$C_l$}')
plt.ylabel(r'$C_{m,c/4}$')
plt.title(r'$C_{m,c/4}$ mean airfoil at $h = 25000$ ft and $V = 1.1V_{cruise}$')
plt.axhline(y=0,color='k',markersize=8)
plt.axvline(x=0,color='k',markersize=8)
plt.grid(True, linestyle='-.')
plt.show()
stat6.savefig(figCm250)
stat6.close()


#wing planform

wing = plt.figure()
plt.plot(y_wing,x_wing,color='k',markersize=12)
plt.xlabel(r'Y (ft)')
plt.ylabel(r'X (ft)')
plt.title(r'Wing Planform')
#plt.axhline(y=0,color='k',markersize=4)
#plt.axvline(x=0,color='k',markersize=4)
#plt.axvline(x = b_kink, ymin = 0., ymax = c_root,linewidth = 1, linestyle ="--",color ='black')
y_c_mean = float(y_c_mean)
xle_c_mean = float(xle_c_mean)
xte_c_mean = float(xle_c_mean+c_mean)
plt.vlines(x=y_c_mean, ymin=xle_c_mean, ymax=xte_c_mean, colors='black', linestyles='--',label=r'M.A.C')
plt.axis('equal')
xmin = float(0.)
xmax = float(b_wing/2+0.1*b_wing/2)
ymin = float(-0.1*c_root)
ymax = float(c_root+0.1*c_root)
plt.axis([xmin, xmax, ymin, ymax])
x_p1 = float(x_p1)
x_p2 = float(x_p2)
plt.ylim(max(x_wing),min(x_wing))
plt.legend(loc=0)
plt.show()
Wing.savefig(wing)
Wing.close()


wing_eq = plt.figure()
plt.plot(y_wing_eq,x_wing_eq,color='k',markersize=12)
plt.xlabel(r'Y (ft)')
plt.ylabel(r'X (ft)')
plt.title(r'Equivalent Wing Planform')
#plt.axhline(y=0,color='k',markersize=4)
#plt.axvline(x=0,color='k',markersize=4)
#plt.axvline(x = b_kink, ymin = 0., ymax = c_root,linewidth = 1, linestyle ="--",color ='black')
y_c_mean_eq = float(y_c_mean_eq)
xle_c_mean_eq = float(xle_c_mean_eq)
xte_c_mean_eq = float(xle_c_mean_eq+c_mean)
plt.vlines(x=y_c_mean_eq, ymin=xle_c_mean_eq, ymax=xte_c_mean_eq, colors='black', linestyles='--',label=r'M.A.C.')
plt.axis('equal')
xmin = float(0.)
xmax = float(b_wing/2+0.1*b_wing/2)
ymin = float(x_p1)
ymax = float(x_p2)
plt.axis([xmin, xmax, ymin, ymax])
plt.ylim(x_p2,x_p1)
plt.legend(loc=0)
plt.show()
Wing_eq.savefig(wing_eq)
Wing_eq.close()
