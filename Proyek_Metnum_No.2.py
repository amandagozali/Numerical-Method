# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 17:40:27 2023

@author: HP
"""

# NOMOR 2

import numpy as np
import pandas

# Diketahui pada soal
x0=0 #waktu awal
xn=0.1 #waktu akhir
C=10**(-3)
L=50*10**(-3)
R=4.7
dt=0.002
y0=15/4.7
z0=0

#Mendefinisikan PD Orde 2
def f(x,y,z):
    '''Fungsi f ini merupakan ruas kanan dari persamaan y'(x)=z(x),
    dengan y(x)=i(x) dan z(x)=i'(x)'''
    return z

def g(x,y,z):
    ''''Fungsi g ini merupakan ruas kanan dari persamaan z'(x)=(-R*y(x)-C*z(x)),
    dengan y(x)=i(x) dan z(x)=i'(x)'''
    return (-R*z-y/C)/L

''' Keterangan Variabel: 
    x = waktu (dalam detik sesuai dengan soal) 
    y = i 
    z = i (output)'''

# METODE_EULER
def euler(f,g,x0,xn,y0,z0,dt):
    # Membuat partisi dan perhitungan nilai iterasi i
    p = int((xn-x0)/dt)   # Menggunakan fungsi int agar nilainya positif
    Y = np.zeros(p+1)
    X = np.zeros(p+1)
    
    # Membuat tempat menyimpan hasil V yang melalui kapasitor
    V = np.zeros(p+1)
    
    #Nilai awal
    Y[0] = y0
    X[0] = x0
    V[0] = x0*y0/C
    x = x0
    y = y0
    z = z0
    
    for i in range(0,p):
        # Slope
        s =f(x,y,z)
        j =g(x,y,z)
        
        # Iterasi lanjutan
        yn = y + dt*s
        zn = z+dt*j
        
        # Persiapan iterasi lanjutan beserta hasil
        y = yn
        z = zn
        x += dt
        Y[i+1] = yn
        X[i+1] = x
        
        # Pencarian nilai v untuk iterasi i+1 atau untuk t=(i+1)*0,002
        V[i+1] = x*yn/C
    return X,Y,V


# METODE_HEUN
def heun(f,g,x0,xn,y0,z0,dt):
    # Membuat partisi dan perhitungan nilai iterasi i
    p = int((xn-x0)/dt) # Menggunakan fungsi int agar nilainya positif
    Y = np.zeros(p+1) 
    X = np.zeros(p+1)
    
    # Membuat tempat menyimpan hasil V yang melalui kapasitor
    V = np.zeros(p+1)
    
    # Nilai awal
    Y[0] = y0
    X[0] = x0
    V[0] = x0*y0/C
    x = x0
    y = y0
    z = z0
    
    for i in range(0,p):
        # Slope
        s = f(x,y,z) 
        j = g(x,y,z)
        
        # Iterasi lanjutan
        
        # Untuk PREDIKTOR
        yn =y + dt*s
        zn =z + dt*j
        # Untuk KOREKTOR
        y1 = y+dt/2*(f(x,y,z)+f(x+dt,yn,zn))
        z1 = z+dt/2*(g(x,y,z)+g(x+dt,yn,zn))
        
        # Persiapan iterasi lanjutan beserta hasil
        y = y1
        z = z1
        x += dt
        Y[i+1] = y1
        X[i+1] = x
        
        # Pencarian nilai v untuk iterasi i+1 atau untuk t=(i+1)*0,002
        V[i+1] = x*y1/C
    return X,Y,V
heun(f,g,0,0.1,15/4.7,0,0.002)

# METODE RUNGE-KUTTA ORDE 4
def rko4(f,g,x0,xn,y0,z0,dt):
    # Membuat partisi dan perhitungan nilai iterasi i
    p = int((xn-x0)/dt) # Menggunakan fungsi int agar nilainya positif
    Y = np.zeros(p+1)
    X = np.zeros(p+1)
    
    # Membuat tempat menyimpan hasil V yang melalui kapasitor
    V = np.zeros(p+1)
    
    #Nilai awal
    Y[0] = y0
    X[0] = x0
    V[0] = x0*y0/C
    x = x0
    y = y0
    z = z0
    
    for i in range(0,p):
        # Pencarian nilai k1-k4 dan l1-l4
        k1 =f(x,y,z)
        l1 =g(x,y,z)
       
        k2 =f((x+dt/2),(y+k1/2*dt),(z+l1/2*dt))
        l2 =g((x+dt/2),(y+k1/2*dt),(z+l1/2*dt))
        
        k3 =f((x+dt/2),(y+k2/2*dt),(z+l2/2*dt))
        l3 =g((x+dt/2),(y+k2/2*dt),(z+l2/2*dt))
        
        k4 =f((x+dt),(y+k3*dt),(z+l3*dt))
        l4 =g((x+dt),(y+k3*dt),(z+l3*dt))
        
        # Membuat faktor penambah di iterasi untuk fungsi f
        kdt = (k1+2*k2+2*k3+k4)/6*dt
        ldt = (l1+2*l2+2*l3+l4)/6*dt
        
        # Iterasi lanjutan
        yn = y+kdt
        zn = z+ldt
        
        # Persiapan iterasi lanjutan beserta hasil
        y = yn
        z = zn
        x += dt
        Y[i+1] = yn
        X[i+1] = x
        
        # Pencarian nilai v untuk iterasi i+1 atau untuk t=(i+1)*0,002
        V[i+1]=x*yn/C
    return X,Y,V
rko4(f,g,0,0.1,15/4.7,0,0.002)

print('Berikut merupakan pilihan metode pengerjaan yang ada')
bo=True
while bo==True:
    rules=input('\n1. Metode Euler \n2. Metode Heun \n3. Metode Runge Kutta Orde 4 \nPilihan Metode :')
    if rules=="1" or rules=="2" or rules=="3":
        bo=False
    else:
        print("Masukkan angka 1,2, atau 3 saja")
        bo=True
        
if rules=="1":
    a=euler(f,g,x0,xn,y0,z0,dt)
    d={'t':a[0],'i':a[1],'V_c':a[2]}
    table=pandas.DataFrame(data=d)
    print("-------------------------")
    print('''Berikut adalah nilai i dan V_c pada t tertentu yang 
          dihasilkan Metode Euler''')
    print("          ")
    print(table)
if rules=="2":
    a=heun(f,g,x0,xn,y0,z0,dt)
    d={'t':a[0],'i':a[1],'V_c':a[2]}
    table=pandas.DataFrame(data=d)
    print("-------------------------")
    print('''Berikut adalah nilai i dan V_c pada t tertentu yang 
          dihasilkan Metode Heun''')
    print("          ")
    print(table)
if rules=="3":
    a=rko4(f,g,x0,xn,y0,z0,dt)
    d={'t':a[0],'i':a[1],'V_c':a[2]}
    table=pandas.DataFrame(data=d)
    print("-------------------------")
    print('''Berikut adalah nilai i dan V_c pada t tertentu yang 
          dihasilkan Metode Runge Kutta Orde 4''')
    print("          ")
    print(table)
