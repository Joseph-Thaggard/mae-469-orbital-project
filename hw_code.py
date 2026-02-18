# Work done initially in matlab, refactored to python for function development

## Code for HW2

## Problem 2: Determine the 6 orbital elements.
# Find a, i, e, omega, w, theta
import numpy as np
print("Problem 2 program output")
r = [2,0.5,1]
v = [0.5,0.5,-0.5]
# Calculate h
h = cross(r,v)
h0 = norm(h)
k = [0,0,1]
n = cross(k,h)
n0 = norm(n)
#u = 3.986*10^5
u = 1
v0 = norm(v)
r0 = norm(r)
evec = ((1/u)*((v0**2)-(u/r0))*r)-(v*dot(r,v))
ecc = norm(evec)
p = (h0**2)/u
a = p/(1-ecc**2)
i_el = acosd(h[2]/h0)
omega = acosd(n[0]/n0)
w = acosd(dot(n,evec)/(n0*ecc))
w2 = 360-w
theta = acosd(dot(evec,r)/(ecc*r0))

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
r_pqw = [r0*np.cos(np.radians(theta)), r0*np.sin(np.radians(theta)), 0]
v_pqw = np.sqrt(u/p)*[-np.sin(np.radians(theta)), (ecc+np.cos(np.radians(theta))), 0]
R = [
    (np.cos(np.radians(RAAN))*np.cos(np.radians(omega)) - np.sin(np.radians(RAAN))*np.sin(np.radians(omega))*np.cos(np.radians(i_el))), ...
    (-np.cos(np.radians(RAAN))*np.sin(np.radians(omega)) - np.sin(np.radians(RAAN))*np.cos(np.radians(omega))*np.cos(np.radians(i_el))), ...
    np.sin(np.radians(RAAN))*np.sin(np.radians(i_el));
    
    (np.sin(np.radians(RAAN))*np.cos(np.radians(omega)) + np.cos(np.radians(RAAN))*np.sin(np.radians(omega))*np.cos(np.radians(i_el))), ...
    (-np.sin(np.radians(RAAN))*np.sin(np.radians(omega)) + np.cos(np.radians(RAAN))*np.cos(np.radians(omega))*np.cos(np.radians(i_el))), ...
    -np.cos(np.radians(RAAN))*np.sin(np.radians(i_el));
    
    np.sin(np.radians(omega))*np.sin(np.radians(i_el)), ...
    np.cos(np.radians(omega))*np.sin(np.radians(i_el)), ...
    np.cos(np.radians(i_el))
]
r_ijk = np.dot(R,r_pqw)
v_ijk = np.dot(R,v_pqw)