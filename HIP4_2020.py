##### The Electric Field, and three different methods of calculating electric potential, of an L charge

##### Running Python 3.6.7 ("new" vpython)
### This should be identical to Python "Jupyter"
### "Analytical function" will not run in glowscript due to math functions like sqrt() and fabs(), but everything else should work
### Everything *should* run on older versions of vpython, but some syntax needs to be changed.
### Two examples:
# "from visual import *" is now "from vpython import *"
# "scene=display()" is now "scene=canvas()"

##### PH213 Jared Smith 30.April.2019

from __future__ import division, print_function
from math import *
from vpython import *
import numpy as np

print("Key: \n")
print("1. Yellow Sphere Height = Actual voltage given by analytical function")
print("---Should be 100% accurate at values for which it exists \n")
print("2. Transparent Rod Height = Voltage approximated by numerical integration of the voltage of a point charge")
print("---Dependent on dq, a small bit of charge \n")
print("3. Blue Sphere Height = Voltage approximated by numerical integration of E dot ds")
print("---Not just dependent on dq, but also dt, a small bit of path \n")
print("4. Blue Arrow = Vector pointing in the direction of ds at the moment when the E dot ds integral is evaluated")
print("---Options to slow down this animation exist in the code, so you can watch these vectors in more detail")
print("---There is also an option to toggle whether or not the vector is hidden after calculations finish \n")
print("5. Red, Yellow, & Green Arrows = The strength and direction of the electric field")
print("---In terms of strength: red>yellow>green \n")
print("Note: Options to print voltage exist in the code \n")

scene = canvas(title='The L Charge',
	width=600, height=600, range=0.2,center=vector(0.1,0.1,0))#range=20cm

vlabel=label(pixel_pos=True,pos=vector(6,16,0),text="",align="left")

k=8.99e9 #coulomb's constant
Q=50e-6 #total charge in C
N=20 #number of grid spheres
N1=20 #pieces of something--change this value to increase accuracy of dq
N2=N1-1 #1 less than pieces of something
M=N1+N2+2 #compound pieces of something (for dq)
L=15e-2 #total combined length in meters
dl=L/N1 #small piece of length
dq=Q/M #small piece of charge
dt=0.005 #not a bit of time--a bit of distance

#creates wire using M bits of dWire
wire=[]
for x in arange(0.0,L,dl):
	wire=wire+[sphere(pos=vector(x,0,0),radius=0.005,color=color.red)]
for y in arange(dl,L,(L-dl)/N2):
	wire=wire+[sphere(pos=vector(0,y,0),radius=0.005,color=color.red)]

#creates grid of exactly N*N little spheres
grid=[]
for x in arange(0.0,0.2,0.2/N):
	for y in arange(0.0,0.2,0.2/N):
		grid=grid+[sphere(pos=vector(x,y,0),radius=0.0001,visible=False)]

#loop for e-field and potential from dot product integral approximation
for i in grid:	
	E=vector(0,0,0) #initial E
	dP=0 #piece of potential
	PS=sphere(pos=i.pos,radius=0.001,color=color.blue) #Potential Sphere
	Pvec=arrow(pos=PS.pos,axis=vector(0,0,0),shaftwidth=0.0005,color=color.blue) #Direction of ds
	for p in wire:
		r=i.pos-p.pos
		if mag(r)>0:
			dE=norm(r)*k*dq/mag(r)**2
			E=E+dE
			#this is where the loop estimates potential
			###very low values of dt will yield the most accurate results, but will also make the calculation take much longer***
			t=mag(r) #initial distance
			while t<mag(r)+0.5:#the +c term should make kq/t^2 as close to 0 as possible for maximum accuracy
				#rate(5000) #for slowing down animation, if so desired
				dE=norm(r)*k*dq/t**2 #for this approximation to work, r=t
				ds=norm(r) #ds points in the radially outward direction
				Edr=dot(dE,ds) #note that this is equal to mag(dE) times dr
				dP=dP+Edr*dt #a small bit of potential
				t=t+dt
		PS.pos=i.pos+vector(0,0,dP)/1e8 #models potential as height of sphere
		Pvec.pos=PS.pos #positioned at sphere
		Pvec.axis=norm(r)/100 #pointing in direction of ds	
	#print("Potential at ",i.pos," =", PS.pos.z*1e8," Volts")  ###print volts
	Pvec.visible=False #comment this out to keep last position of ds displayed
	del Pvec #and possibly this too
	vlabel.text="Voltage in J/C", PS.pos.z*1e8
	vlabel.color=color.blue
	energy=arrow(pos=i.pos,axis=E/mag(E)/100,shaftwidth=0.0001,color=vector(mag(E)-2e14/mag(E),3e7/mag(E),0))	

#loop for potential using voltage from point charge
for i in grid:	
	rate(50) #for slowing down animation, if so desired
	V=0
	for p in wire:
		r=i.pos-p.pos
		if mag(r)>0:
			dV=k*dq/mag(r)
			V=V+dV
	if V>1e19:
		V=0
	potential=cylinder(pos=i.pos,axis=V*vector(0,0,1)/1e8,radius=0.001,opacity=0.4) #note that this is the same scale factor as before
	vlabel.text="Voltage in J/C",potential.axis.z*1e8
	vlabel.color=color.white
	#print("Potential at ",i.pos," =", V," Volts")	###print volts

#defines analytical function
def f(x,y):
	A=L-x
	B=L-y
	C=sqrt(y**2+A**2)
	D=sqrt(x**2+B**2)
	E=sqrt(x**2+y**2)
	if ((-x+E)*(-y+E)) == 0:
		return 0
	else:
		F=log(fabs(((A+C)*(B+D))/((-x+E)*(-y+E))))
		#F=log(fabs((A+C)/(-x+E)))+log(fabs((B+D)/(-y+E))) #alternate formula; more difficult to work with
		G=Q*k*F/(2*L)
		return G

#runs analytical function	
for i in grid:
	rate(50) #for slowing down animation, if so desired
	if f(i.pos.x,i.pos.y)>0:
		anvec=vector(i.pos.x,i.pos.y,f(i.pos.x,i.pos.y)/1e8)
		anshape=sphere(pos=anvec, radius=0.001, color=color.yellow)
		#print("Potential at ",i.pos," =", anvec.z*1e8," Volts") ###print volts
		vlabel.text="Voltage in J/C",anvec.z*1e8
		vlabel.color=color.yellow

vlabel.visible=False
