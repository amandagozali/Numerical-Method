# -*- coding: utf-8 -*-
"""
Created on Tue Jan 3 17:45:57 2023

@author: Amanda Gozali, Lazarus Lie, Varani Clarissa Wedhasanti, Marnida Harnia
"""

import numpy as np
from numpy import *
import sys

#konstanta planck
h = 6.6256e-27

#kecepatan cahaya
c = 2.99792e+10

#konstanta boltzman
k = 1.3805e-16

#pi
pi = np.pi


def f(x,T):
    return (2*pi*h*c**2)/((x**5)*(np.exp((h*c)/(k*x*T)) - 1))

def ftakwajar(x,T):
    if x == 0 :
        return 0
    if x == sys.maxsize :
        return (2*pi*h*c**2)/(((1/x)**7)*(np.exp((h*c)/(k*(1/x)*T)) - 1))

def trapesium (a,b,n,T) :
    if a != 0 :
        p = (b - a) / n
        integration = f(a,T) + f(b,T)
        for i in range(1,n):
            c = a + i*p
            integration = integration + 2 * f(c,T)
        hasil = integration * p/2
        return hasil
    if a == 0 :
        if b != sys.maxsize :
            p = (b - a) / n
            integration = ftakwajar(a,T) + f(b,T)
            for i in range(1,n):
                c = a + i*p
                integration = integration + 2 * f(c,T)
            hasil = integration * p/2
            return hasil
        if b == sys.maxsize :
            p = (9223372036854 - a) / n
            integration = ftakwajar(a,T) + ftakwajar(b,T)
            for i in range(1,n):
                c = a + i*p
                integration = integration + 2 * f(c,T)
            hasil = integration * p/2
            return hasil

def simpson(a,b,n,T):
    if a != 0 :
        p = (b - a) / n
        integration = f(a,T) + f(b,T)
        for i in range(1,n):
            c = a + i*p
            if i%2 == 0:
                integration = integration + 2 * f(c,T)
            else:
                integration = integration + 4 * f(c,T)
        hasil = integration * p/3
        return hasil
    if a == 0 :
        if b != sys.maxsize:
            p = (b - a) / n
            integration = ftakwajar(a,T) + f(b,T)
            for i in range(1,n):
                c = a + i*p
                if i%2 == 0:
                    integration = integration + 2 * f(c,T)
                else:
                    integration = integration + 4 * f(c,T)
            hasil = integration * p/3
            return hasil 
        if b == sys.maxsize :
            p = (9223372036854 - a) / n
            integration = ftakwajar(a,T) + ftakwajar(b,T)
            for i in range(1,n):
                c = a + i*p
                if i%2 == 0:
                    integration = integration + 2 * f(c,T)
                else:
                    integration = integration + 4 * f(c,T)
            hasil = integration * p/3
            return hasil 
  
def Legendre(n,x):
	if (n==0):
		return x*0+1.0
	elif (n==1):
		return x
	else:
		return ((2.0*n-1.0)*x*Legendre(n-1,x)-(n-1)*Legendre(n-2,x))/n
    
def Turunan(n,x):
	if (n==0):
		return x*0
	elif (n==1):
		return x*0+1.0
	else:
		return (n/(x**2-1.0))*(x*Legendre(n,x)-Legendre(n-1,x))
    
def Akar(pangkat):
    roots=[]
    for i in range (1, int(pangkat)//2 + 1):
        x=cos(pi*(i-0.25)/(pangkat+0.5))
        iters = 0
        while iters < 1000 :
            dx=-Legendre(pangkat,x)/Turunan(pangkat,x)
            x=x+dx
            iters=iters+1
        roots.append(x)
    roots=array(roots)
    if pangkat %2==0:
        roots=concatenate((-1*roots,roots[::-1]))
    else:
        roots=concatenate((-1*roots,[0],roots[::-1]))
    return [roots]

def GaussLegendreWeights(pangkat) :
    W=[]
    [xis]=Akar(pangkat)
    W=2/((1-xis**2)*(Turunan(pangkat,xis)**2))
    return [W,xis]

def GaussLegendreQuadrature(f,pangkat,a,b,T):
    [Ws,xs]= GaussLegendreWeights(pangkat)
    if a != 0 :
        hasil=(b-a)*0.5*sum(Ws*f((b-a)*0.5*xs+(b+a)*0.5,T))
        return hasil
    if a == 0:
        if b != sys.maxsize :
            hasil=(b-a)*0.5*sum(Ws*f((b-a)*0.5*xs+(b+a)*0.5,T))
            return hasil
        if b == sys.maxsize :
            hasil=(9223372036854-a)*0.5*sum(Ws*f((9223372036854-a)*0.5*xs+(9223372036854+a)*0.5,T))
            return hasil

menu_options = {
    1: 'Aturan Trapesium',
    2: 'Aturan Simpson',
    3: 'Metode Gauss-Legendre',
    4: 'Keluar',
    }

def print_menu():
    for key in menu_options.keys():
        print (key, '--', menu_options[key] )
        
def option1():
    while(True):
        try:
            partisi = int(input('Masukan Partisi: '))
            break
        except:
            print('Salah Input. Tolong Masukan Angka Bilangan Bulat!')
    print('\n')
    print('Aturan Trapesium')
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [0,10]")
    print("Besar energi untuk kasus T =10 K")
    print(trapesium(0,10,partisi,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(trapesium(0,10,partisi,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(trapesium(0,10,partisi,1000))
    
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [100,110]")
    print("Besar energi untuk kasus T =10 K")
    print(trapesium(100,110,partisi,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(trapesium(100,110,partisi,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(trapesium(100,110,partisi,1000))
    
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [1000,1010]")
    print("Besar energi untuk kasus T =10 K")
    print(trapesium(1000,1010,partisi,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(trapesium(1000,1010,partisi,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(trapesium(1000,1010,partisi,1000))
    
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [0,infinity]")
    print("Besar energi untuk kasus T =10 K")
    print(trapesium(0,sys.maxsize,partisi,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(trapesium(0,sys.maxsize,partisi,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(trapesium(0,sys.maxsize,partisi,1000))
    print('--------------------------------------------------------')
    print("\n")
    
def option2():
    while (True) :
        try :
            partisi=int(input("Silahkan masukan berapa partisi yang diinginkan : "))
            if partisi %2 == 1 :
                raise ValueError
            break
        except ValueError : 
            print('Salah Input. Tolong Masukan Angka Bilangan Bulat Dan Genap!')
    
    print('\n')
    print('Aturan Gauss')
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [0,10]")
    print("Besar energi untuk kasus T =10 K")
    print(simpson(0,10,partisi,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(simpson(0,10,partisi,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(simpson(0,10,partisi,1000))
    
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [100,110]")
    print("Besar energi untuk kasus T =10 K")
    print(simpson(100,110,partisi,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(simpson(100,110,partisi,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(simpson(100,110,partisi,1000))
    
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [1000,1010]")
    print("Besar energi untuk kasus T =10 K")
    print(simpson(1000,1010,partisi,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(simpson(1000,1010,partisi,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(simpson(1000,1010,partisi,1000))
    
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [0,infinity]")
    print("Besar energi untuk kasus T =10 K")
    print(simpson(0,sys.maxsize,partisi,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(simpson(0,sys.maxsize,partisi,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(simpson(0,sys.maxsize,partisi,1000))
    print('--------------------------------------------------------')
    print("\n")
     
def option3():
    while(True):
        try:
            titik = int(input('Masukan Titik Interpolasi: '))
            if titik <= 1 :
                raise ValueError
            break
        except ValueError :
            print('Salah Input. Tolong Masukan Angka Bilangan Bulat Dan lebih dari 1!')
    pangkat = titik - 1
    [Ws,xs]= GaussLegendreWeights(pangkat)
    print('\n')
    print('Metode Gauss-Legendre')
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [0,10]")
    print("Besar energi untuk kasus T =10 K")
    print(GaussLegendreQuadrature(f,pangkat,0,10,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(GaussLegendreQuadrature(f,pangkat,0,10,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(GaussLegendreQuadrature(f,pangkat,0,10,1000))
    
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [100,110]")
    print("Besar energi untuk kasus T =10 K")
    print(GaussLegendreQuadrature(f,pangkat,100,110,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(GaussLegendreQuadrature(f,pangkat,100,110,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(GaussLegendreQuadrature(f,pangkat,100,110,100))
    
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [1000,1010]")
    print("Besar energi untuk kasus T =10 K")
    print(GaussLegendreQuadrature(f,pangkat,1000,1010,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(GaussLegendreQuadrature(f,pangkat,1000,1010,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(GaussLegendreQuadrature(f,pangkat,1000,1010,1000))
    
    print('--------------------------------------------------------')
    print("Untuk kasus interval panjang gelombang [0,infinity]")
    print("Besar energi untuk kasus T =10 K")
    print(GaussLegendreQuadrature(f,pangkat,0,sys.maxsize,10))
    
    print("Besar energi untuk kasus T =100 K")
    print(GaussLegendreQuadrature(f,pangkat,0,sys.maxsize,100))
    
    print("Besar energi untuk kasus T =1000 K")
    print(GaussLegendreQuadrature(f,pangkat,0,sys.maxsize,1000))
    print('--------------------------------------------------------')
    print("\n")
     
if __name__=='__main__':
    while(True):
        print_menu()
        option = ''
        try:
            option = int(input('Apakah pilihan mu: '))
        except:
            print('Salah Input. Tolong Masukan Angka!')
        if option == 1:
           option1()
        elif option == 2:
            option2()
        elif option == 3:
            option3()
        elif option == 4:
            print('Terima kasih')
            break
        else:
            print('Salah Input. Tolong masukan angka antara 1 sampai 4!')
            

            
            