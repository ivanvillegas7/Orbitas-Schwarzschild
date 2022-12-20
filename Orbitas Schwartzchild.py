# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 09:28:52 2021

@author: Carlos Beltran
"""
import MyN_python_siempre_contigo

import numpy as np

import matplotlib.pyplot as plt

#

# Constantes del problema:
    
c = 1
mu = 1
R = 7 # La distancia inicial al centro medida por el observador lejano

# Consideramos w=[t,r,theta]. 

# Planteamos el sistema de ODEs de segundo orden:
def wddot(sigma, w, wdot):
    t = w[0]
    r = w[1]
    theta = w[2]
    tdot = wdot[0]
    rdot = wdot[1]
    thetadot = wdot[2]
    Factor = -1 + 2*mu/r
    tddot = 2*mu/r**2*rdot*tdot/Factor
    rddot = -Factor*(-mu/r**2*c**2*tdot**2+mu/(2*mu-r)**2*rdot**2+r*thetadot**2)
    thetaddot = -2*rdot*thetadot/r
    wddot = [tddot, rddot, thetaddot]
    return wddot
# Condiciones iniciales para la caida radial:

w0 = np.array([0, R, 0])
r0dot = c/5
theta0dot = 0 # Ha de cumplirse r0dot**2+R**2theta0dot**2<c**2
t0dot = np.sqrt((-c**2-1/(1-2*mu/R)*r0dot**2-R**2*theta0dot**2)/(-1+2*mu/R))/c
w0dot = np.array([t0dot, r0dot, theta0dot])
a = 0 # Tiempo inicial
b = 30# Tiempo final
N = 50000 # Numero de pasos
ss, W, Wdot = MyN_python_siempre_contigo.segundo_orden(wddot, a, b, w0, w0dot, N, verbose = 0)
t = W[:,0]
r = W[:,1]
theta = W[:,2]



plt.figure()
plt.plot(t, r, t, 0*t+2)
plt.xlim([0,max(t)])
plt.ylim([0,10])
plt.xlabel('Tiempo (medido por observador lejano)')
plt.ylabel('Distancia a la estrella (medida por observador lejano)')
plt.title('Objeto en caida libre en la metrica de Schwarzschild)')

#Añadido

from scipy import integrate

def integrando(r, v):
    
    return np.sqrt((1-2*mu/r)*c**2-v**2/(1-2*mu/r))

for i in range(len(r)):
    
    if (r[i] <= (R+1e-4) and r[i] >= (R-1e-4)) or r[i] == R:
        
        if i != 0:
            
            r_int = r[0:i]
            
            t_int = t[0:i]
            
            Wdot1_int = Wdot[0:i, 1]
            
            t_prop = integrate.simpson(integrando(r_int, Wdot1_int), t_int)/c
        
            print('')
            
            print(f'Aproximadamente, habrán pasado {t[i]:.3f} segundos para el observador externo (Schwarzschild) y {t_prop:.3f} para el objeto.')
    
            print('')
        
#%%
# Condiciones iniciales para la caida con velocidad angular inicial
# que lleva a la colision:
w0 = np.array([0, R, 0])
r0dot = c/5
theta0dot = 0 # Ha de cumplirse r0dot**2+R**2theta0dot**2<c**2
if r0dot**2+R**2*theta0dot**2>=c**2:
    print('Error: velocidad absurda.')
t0dot = np.sqrt((-c**2-1/(1-2*mu/R)*r0dot**2-R**2*theta0dot**2)/(-1+2*mu/R))/c
w0dot = np.array([t0dot, r0dot, theta0dot])
a = 0 # Tiempo inicial
b = 30# Tiempo final
N = 50000 # Numero de pasos
ss2, W2, Wdot2 = MyN_python_siempre_contigo.segundo_orden(wddot, a, b, w0, w0dot, N, verbose = 0)
t2 = W2[:,0]
r2 = W2[:,1]
theta2 = W2[:,2]
x=r2*np.cos(theta2)
y=r2*np.sin(theta2)
plt.figure()
plt.plot(x,y)
plt.plot(2*mu*np.cos(np.linspace(0,2*np.pi)),2*mu*np.sin(np.linspace(0,2*np.pi)))
indices = np.array([0, 500, 50000])
for i in indices:
    plt.plot(x[i],y[i],'o')
    cadena='t='+str(round(t2[i]))
    plt.text(x[i]+.2,y[i]+.2,cadena)
i=25000
plt.plot(x[i],y[i],'o')
cadena='t='+str(round(t2[i]))
plt.text(x[i]+.2,y[i]-.6,cadena)
    
#M=min(20,max(r2))*1.3
#xlim([-M,M])
#ylim([-M,M])
#plt.axes().set_aspect('equal')
plt.title('Órbita de colisión')

for i in range(len(r2)):
    
    if (r2[i] <= (R+1e-4) and r2[i] >= (R-1e-4)) or r2[i] == R:
        
        if i > 5000:
            
            r_int = r2[0:i]
            
            t_int = t2[0:i]
            
            Wdot2_int = Wdot2[0:i, 1]
            
            t_prop = integrate.simpson(integrando(r_int, Wdot2_int), t_int)/c
        
            print('')
            
            print(f'Aproximadamente, habrán pasado {t2[i]:.3f} segundos para el observador externo (Schwarzschild) y {t_prop:.3f} para el objeto.')
    
            print('')
                        
