# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 11:17:19 2016

@author: Eric Schmidt
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.integrate import quad


N = 2000
T = 1/1000
TT = 3
t = np.linspace(0.0, N*T, N)
tf = np.linspace(0.0,1.0/(2.0*T),N/2)
w = 2.0*1*np.pi
    
phi_x = np.pi
phi_D = 0
x_mid = 10
x_0 = 10
A_x = 2
D_0 = 3
A_D = 2
k_1 = 0.05

def main():
#    xlinear(noplot=0)
#    sigmalinear(noplot=0)
#    combinedLinear()    
#    
#    xquad(noplot=0)
#    sigmaquad(noplot=0)
    combinedQuad()

def xlinear(noplot):
    
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
        
        temp1 = quad(integrandADlinear, ll[i], ul[i], args=(a[i],x_B[i],D))
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

    D = D_0 + A_D*np.cos(w*t + phi_D)
    
    x_B = None
    x_B = x_0
    
    ul = x_B + D/2
    ll = x_B - D/2
    
    i = 0
    
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp1 = quad(integrandADlinear, ll[i], ul[i], args=(a[i],x_B,D[i]))
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
    D = D_0 + A_D*np.cos(w*t + phi_D)
    
    ll = x_B - D/2
    ul = x_B + D/2
    
    i = 0
    
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp = quad(integrandADlinear, ll[i], ul[i], args=(a[i],x_B[i],D[i]))
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
    D= D_0
    
    ul = x_B + D/2
    ll = x_B - D/2
    
    a = np.zeros(len(t))
    
    i = 0
    
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp = quad(integrandADquad, ll[i], ul[i], args=(a[i],x_B[i],D))
        AD[i] = temp[0]
        
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

    D = D_0 + A_D*np.cos(w*t + phi_D)
    
    x_B = None
    x_B = x_0
    
    ul = x_B + D/2
    
    a = np.zeros(len(t))
    ll = x_B - D/2
    
    i = 0  
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp2 = quad(integrandADquad, ll[i], ul[i], args=(a[i],x_B,D[i]))
        AD[i] = temp2[0]
        
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
    D = D_0 + A_D*np.cos(w*t + phi_D)
    
    ll = x_B - D/2
    ul = x_B + D/2
    
    i = 0
    
    while i < len(t):
        
        a[i] = N/(ul[i] - ll[i])
        
        temp = quad(integrandADquad, ll[i], ul[i], args=(a[i],x_B[i],D[i]))
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
    
def integrandADlinear(x, a, x_B, D):
    return a*(k_1*x)
    
def integrandADquad(x, a, x_B, D):
#    return a*(((-(x - x_mid)**2 + (x_mid)**2))/(x_mid)**2)
    return a*((-(x - x_mid)**2 + (x_mid)**2)/((x_mid)**2))
    
def plotSingle(AD,data,ylabel,title,save_title,stitle):

    n = 1

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plt.subplot(2,1,1)
    plt.plot(t,data)
    plt.ylabel("%s"%(ylabel))
    plt.title("%s"%(title))
    
    plt.subplot(2,1,2)
    plt.plot(t,AD)
    plt.ylabel("$N_D$")
    plt.xlabel("Time")
    plt.title("Particles Detected vs. Time")
    plt.tight_layout()
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    plt.savefig('Plots/%s.png'%(save_title), bbox_inches='tight', dpi=300)

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plt.subplot(1,1,1)
    Af = fft(AD)
    plt.plot(tf[1:], 2/N * np.abs(Af[:N/2])[1:])
    plt.xlim(0,10)
    plt.xlabel("Frequency")
    plt.title("DFT of Particles Detected")
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    plt.savefig('Plots/%s-DFT.png'%(save_title), bbox_inches='tight', dpi=300)
    
def plotCombined(AD,save_title,datax,datasig,stitle):

    n = 1

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plt.subplot(3,1,1)
    plt.plot(t,datax)
    plt.ylabel("$N_D$")
    plt.title("Particles Detected from changing $x_m$")
    
    plt.subplot(3,1,2)
    plt.plot(t,datasig)
    plt.ylabel("$N_D$")
    plt.title("Particles Detected from changing $\sigma$")
    
    plt.subplot(3,1,3)
    plt.plot(t,AD)
    plt.ylabel("$N_D$")
    plt.xlabel("Time")
    plt.title("Particles Detected vs. Time")
    plt.tight_layout()
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    plt.savefig('Plots/%s.png'%(save_title), bbox_inches='tight', dpi=300)

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plt.subplot(1,1,1)
    Af = fft(AD)
    plt.plot(tf[1:], 2/N * np.abs(Af[:N/2])[1:])
    plt.xlim(0,5)
#    plt.ylim(0,2)
    plt.xlabel("Frequency")
    plt.title("DFT of Particles Detected")
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    plt.savefig('Plots/%s-DFT.png'%(save_title), bbox_inches='tight', dpi=300)
    
main()