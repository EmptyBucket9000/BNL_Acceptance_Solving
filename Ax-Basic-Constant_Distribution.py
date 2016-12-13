# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 11:17:19 2016

@author: Eric Schmidt
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.integrate import quad

save_plots = 1
through_contact = 1 # Set to 1 for positrons that passed through matter
fit = 2 # 1 for linear, 2 for quadratic

Ne = 2153
N = 1000
T = 1/10
TT = 3
t = np.linspace(0.0, N*T, N)
tf = np.linspace(0.0,1.0/(2.0*T),N/2)
w = ((6/5)*np.pi)*10**6 # (rad/s) Oscillating frequency
    
phi_x = -np.pi/8 # (rad) Phase-offset of x-position mean oscillations
phi_sigma = np.pi/2 # (rad) Phase-offset of beam distribution width oscillations
x_mid = .1 # (m) Middle of the acceptance function
x_0 = 0 # (m) Initial beam distribution mean
sigma_0 = 14.5*10**-3 # (m) Initial beam distribution width
D_0 = 80*10**-3 # (m) Maximum physical width of beam
A_x = 2*10**-3 # (m) Initial beam distribution mean oscillation amplitude
A_sigma = 5.5*10**-3 # (m) Initial beam distribution width oscillation amplitude

if through_contact == 1:
    if fit == 1:
        k_0 = 0.762043      # () Constant term
        k_1 = -0.5931      # () Linear term
    if fit == 2:
        k_0 = 0.760101      # () Constant term
        k_1 = -0.609356      # () Linear term
        k_2 = 7.97297     # () Quadratic term
else:
    if fit == 1:
        k_0 = 0.793666      # () Constant term
        k_1 = -1.277      # () Linear term
    if fit == 2:
        k_0 = 0.790605      # () Constant term
        k_1 = -1.31183      # () Linear term
        k_2 = 13.1419     # () Quadratic term

#phi_x = np.pi
phi_D = 0
#x_mid = 10
#x_0 = 10
#A_x = 2
#D_0 = 3
#A_D = 2
#k_1 = 0.05

def main():
#    xlinear(noplot=0)
#    sigmalinear(noplot=0)
#    combinedLinear()    
#    
#    xquad(noplot=0)
#    sigmaquad(noplot=0)
    combinedQuad()

def xlinear(noplot):

    # Initialize variables
    
    AD = np.zeros(len(t))

    x_B = x_0 + A_x*np.cos(w*t + phi_x)
    
    D = None
    D = D_0
    
    a = np.zeros(len(t))
    
    ul = x_B + D/2
    ll = x_B - D/2
    
    i = 0
    
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp1 = quad(integrandADlinear, ll[i], ul[i], args=(a[i],x_B[i]))
        AD[i] = temp1[0]
        
#        Nt = quad(integrandN, ll[i], ul[i], args=(a[i],x_B[i],D))
#        print(Nt[0])
        
        i = i + 1
    
    title = "Postion vs. Time"
    ylabel = "$x_m$"
    save_title = "CF_Linear_Position"
    stitle = "Constant Distribution: Linear Position"
    
    if noplot==0:
        plotSingle(AD,x_B,ylabel,title,save_title,stitle)
    
    return AD
    
def sigmalinear(noplot):
    
    AD = np.zeros(len(t))    
    a = np.zeros(len(t))

    D = (D_0-A_sigma) + (A_sigma)*np.cos(w*t + phi_D)
    
    x_B = None
    x_B = x_0
    
    ul = x_B + D/2
    ll = x_B - D/2
    
    i = 0
    
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp1 = quad(integrandADlinear, ll[i], ul[i], args=(a[i],x_B))
        AD[i] = temp1[0]
        
#        Nt = quad(integrandN, ll[i], ul[i], args=(a[i],x_B,D[i]))
#        print(Nt[0])
        
        i = i + 1
    
    title = "Sigma vs. Time"
    ylabel = "$\sigma$"
    save_title = "CF_Linear_Sigma"
    stitle = "Constant Distribution: Linear Sigma"
    
    if noplot==0:
        plotSingle(AD,D,ylabel,title,save_title,stitle)
    
    return AD
    
def combinedLinear():
    
    AD = np.zeros(len(t))        
    a = np.zeros(len(t))

    x_B = x_0 + A_x*np.cos(w*t + phi_x)
    D = (D_0-A_sigma) + (A_sigma)*np.cos(w*t + phi_D)
    
    ll = x_B - D/2
    ul = x_B + D/2
    
    i = 0
    
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp = quad(integrandADlinear, ll[i], ul[i], args=(a[i],x_B[i]))
        AD[i] = temp[0]
        
#        Nt = quad(integrandN, ll[i], ul[i], args=(a[i],x_B[i],D[i]))
#        print(Nt[0])
        
        i = i + 1
        
    datax = xlinear(1)
    datasig = sigmalinear(1)
    
    save_title = "CF_Linear_Combined"
    stitle = "Constant Distribution: Linear Combined"
    plotCombined(AD, save_title, datax, datasig, stitle)
    
def xquad(noplot):
    
    AD = np.zeros(len(t))    

    x_B = x_0 + A_x*np.cos(w*t + phi_x)
    
    D = None
    D = D_0
    
    ul = x_B + D/2
    ll = x_B - D/2
    
    a = np.zeros(len(t))
    
    i = 0
    
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp = quad(integrandADquad, ll[i], ul[i], args=(a[i],x_B[i]))
        AD[i] = temp[0]
        
#    return a*k_1*x**2
#    return a*((-(x - x_mid)**2 + (x_mid)**2)/((x_mid)**2))
#    return (a*(np.sqrt(2*np.pi*sigma**2)**(-1))*np.exp(-((x-x_B)**2)/(2*sigma**2)))*(k_1*x**2)
        
#        Nt = quad(integrandN, ll[i], ul[i], args=(a[i],x_B[i],D))
#        print(x_B[i])
#        print(ll[i])
#        print(ul[i])
#        print(a[i])
#        print(Nt[0])
#        print(D)
#        print(AD[i])
        
        i = i + 1
    
    title = "Postion vs. Time"
    ylabel = "$x_m$"
    save_title = "CF_Quad_Position"
    stitle = "Constant Distribution: Quadratic Position"
    
    if noplot==0:
        plotSingle(AD,x_B,ylabel,title,save_title,stitle)
    
    return AD
    
def sigmaquad(noplot):
    
    AD = np.zeros(len(t))    
    a = np.zeros(len(t))

    D = (D_0-A_sigma) + (A_sigma)*np.cos(w*t + phi_D)
    
    x_B = None
    x_B = x_0
    
    ul = x_B + D/2
    ll = x_B - D/2
    
    i = 0  
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp = quad(integrandADquad, ll[i], ul[i], args=(a[i],x_B))
        AD[i] = temp[0]
        
#        Nt = quad(integrandN, ll[i], ul[i], args=(a[i],x_B,D[i]))
#        print(x_B)
#        print(ll[i])
#        print(ul[i])
#        print(a[i])
#        print(Nt[0])
#        print(D[i])
#        print(AD[i])
        
        i = i + 1
    
    title = "Sigma vs. Time"
    ylabel = "$\sigma$"
    save_title = "CF_Quad_Sigma"
    stitle = "Constant Distribution: Quadratic Sigma"
    
    if noplot==0:
        plotSingle(AD,D,ylabel,title,save_title,stitle)
    
    return AD
    
def combinedQuad():
    
    AD = np.zeros(len(t))        
    a = np.zeros(len(t))

    x_B = x_0 + A_x*np.cos(w*t + phi_x)
    D = (D_0-A_sigma) + (A_sigma)*np.cos(w*t + phi_D)
    
    ll = x_B - D/2
    ul = x_B + D/2
    
    i = 0
    
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp = quad(integrandADquad, ll[i], ul[i], args=(a[i],x_B[i]))
        AD[i] = temp[0]
        
#        Nt = quad(integrandN, ll[i], ul[i], args=(a[i],x_B[i],D[i]))
#        print(Nt[0])
        
        i = i + 1
        
    datax = xquad(noplot = 1)
    datasig = sigmaquad(noplot = 1)
        
    save_title = "CF_Quad_Combined"
    stitle = "Constant Distribution: Quadratic Combined"
    plotCombined(AD,save_title,datax,datasig,stitle)
    
#def integranda(x, x_B, D):
#    return 1
    
def integrandN(x, a, x_B, D):
    return a
    
def integrandADlinear(x, a, x_B):
#    return a*x_B*k_0*x
    return a*(k_0 + k_1*x)
    
def integrandADquad(x, a, x_B):
#    return a*x_B*k_1*(x-x_mid)**2
    return a*(k_0 + k_1*x + k_2*x**2*100)
    
def plotSingle(AD,data,ylabel,title,save_title,stitle):

    n = 1

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plot1 = plt.subplot(2,1,1)
    plt.plot(t,data)
    plt.ylabel("%s"%(ylabel))
    plt.title("%s"%(title))
    
    plot2 = plt.subplot(2,1,2)
    plt.plot(t,AD)
    plt.ylabel("$N_D$")
    plt.xlabel("Time ($\mu$s)")
    plt.title("Particles Detected vs. Time")
    plt.tight_layout()
    plot1.locator_params(axis='y',nbins=5)
    plot2.locator_params(axis='y',nbins=5)
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    if save_plots == 1:
        plt.savefig('Plots/%s.png'%(save_title), bbox_inches='tight', dpi=300)

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plt.subplot(1,1,1)
    Af = fft(AD)
    plt.plot(tf[1:], 2/N * np.abs(Af[:N/2])[1:])
    plt.xlim(0,2)
    plt.xlabel("Frequency (MHz)")
    plt.title("DFT of Particles Detected")
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    if save_plots == 1:
        plt.savefig('Plots/%s-DFT.png'%(save_title), bbox_inches='tight',
                    dpi=300)
    
def plotCombined(AD,save_title,datax,datasig,stitle):
    n = 1

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plot1 = plt.subplot(3,1,1)
    plt.plot(t,datax)
    plt.ylabel("$N_D$")
    plt.title("Particles Detected from changing $x_m$")
    
    plot2 = plt.subplot(3,1,2)
    plt.plot(t,datasig)
    plt.ylabel("$N_D$")
    plt.title("Particles Detected from changing $\sigma$")
    
    plot3 = plt.subplot(3,1,3)
    plt.plot(t,AD)
    plt.ylabel("$N_D$")
    plt.xlabel("Time ($\mu$s)")
    plt.title("Particles Detected vs. Time")
    plt.tight_layout()
    plot1.locator_params(axis='y',nbins=5)
    plot2.locator_params(axis='y',nbins=5)
    plot3.locator_params(axis='y',nbins=5)
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    if save_plots == 1:
        plt.savefig('Plots/%s.png'%(save_title), bbox_inches='tight', dpi=300)

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plt.subplot(1,1,1)
    Af = fft(AD)
    plt.plot(tf[1:], 2/N * np.abs(Af[:N/2])[1:])
    plt.xlim(0,2)
#    plt.ylim(0,2)
    plt.xlabel("Frequency (MHz)")
    plt.title("DFT of Particles Detected")
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    if save_plots == 1:
        plt.savefig('Plots/%s-DFT.png'%(save_title), bbox_inches='tight', dpi=300)
    
main()
