# Work done initially in matlab, refactored to python for function development

## Code for HW2

## Problem 2: Determine the 6 orbital elements.
# Find a, i, e, omega, w, theta
import numpy as np
import numpy.linalg as la
print("Problem 2 program output")
r = np.array([2,0.5,1])
v = np.array([0.5,0.5,-0.5])
# Calculate h
h = np.cross(r,v)
h0 = la.norm(h)
k = np.array([0,0,1])
n = np.cross(k,h)
n0 = la.norm(n)
#u = 3.986*10^5
u = 1
v0 = la.norm(v)
r0 = la.norm(r)
evec = ((1/u)*((v0**2)-(u/r0))*r)-(v*np.dot(r,v))
ecc = la.norm(evec)
p = (h0**2)/u
a = p/(1-ecc**2)
i_el = np.degrees(np.arccos(h[2]/h0))
omega = np.degrees(np.arccos(n[0]/n0))
w = np.degrees(np.arccos(np.dot(n,evec)/(n0*ecc)))
w2 = 360-w
theta = np.degrees(np.arccos(np.dot(evec,r)/(ecc*r0)))

## Problem 3: Getting r_ijk and v_ijk from 6 orbital elements. Using notation in Problem 3 for unit names.
# Clear previous problem variables
print("Problem 3 program output")
# Cannonical units
u = 1
a = 6.3920
ecc = 0.4880
# degrees for following units
i_el = 63.50 
RAAN = 96.40
omega = 246.00
theta = 18.00
# Note that reference frame is geocentric-eclipitic (geocentric equatorial)
p = a*(1-ecc**2)
r0 = p/(1+ecc*np.cos(np.radians(theta)))
# Using verical vector notation
r_pqw = np.array([r0*np.cos(np.radians(theta)), r0*np.sin(np.radians(theta)), 0])
v_pqw = np.sqrt(u/p)*np.array([-np.sin(np.radians(theta)), (ecc+np.cos(np.radians(theta))), 0])
R = np.array([[(np.cos(np.radians(RAAN))*np.cos(np.radians(omega)) - np.sin(np.radians(RAAN))*np.sin(np.radians(omega))*np.cos(np.radians(i_el))), 
    (-np.cos(np.radians(RAAN))*np.sin(np.radians(omega)) - np.sin(np.radians(RAAN))*np.cos(np.radians(omega))*np.cos(np.radians(i_el))), 
    np.sin(np.radians(RAAN))*np.sin(np.radians(i_el))], 
    [(np.sin(np.radians(RAAN))*np.cos(np.radians(omega)) + np.cos(np.radians(RAAN))*np.sin(np.radians(omega))*np.cos(np.radians(i_el))), 
    (-np.sin(np.radians(RAAN))*np.sin(np.radians(omega)) + np.cos(np.radians(RAAN))*np.cos(np.radians(omega))*np.cos(np.radians(i_el))), 
    -np.cos(np.radians(RAAN))*np.sin(np.radians(i_el))],
    [np.sin(np.radians(omega))*np.sin(np.radians(i_el)), np.cos(np.radians(omega))*np.sin(np.radians(i_el)), np.cos(np.radians(i_el))]])
r_ijk = np.dot(R,r_pqw)
v_ijk = np.dot(R,v_pqw)
print("r_ijk: ", r_ijk)
print("v_ijk: ", v_ijk)