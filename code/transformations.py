#!/usr/bin/env python3.5

import numpy as np

MU = (3.986*10**5)*3600**2

def cartesian2spherical(coordinates):
    x,y,z = coordinates

    return np.arctan2(np.sqrt(x**2 + y**2), z), np.arctan2(y,x)

def spherical2cartesian(coordinates):
    lat,lon = coordinates
    lat = lat*np.pi/180
    lon = lon*np.pi/180
    xc = 6371*np.cos(lat)*np.cos(lon)
    yc = 6371*np.cos(lat)*np.sin(lon)
    zc = 6371*np.sin(lat)

    return np.array([xc,yc,zc]).T

def keplerian2cartesian(a,e,E,f):
    # positional vector in PQW frame from keplerian elements
    r = a*(1-e*np.cos(E)) # magnitude
    r_PQW = r*np.array([np.cos(f),
                        np.sin(f),
                        0])
    # velocity vector in PQW frame from keplerian elements
    v = (np.sqrt(MU*a)/r) # magnitude
    v_PQW = v*np.array([-np.sin(E),
                        np.sqrt(1-e**2)*np.cos(E),
                        0])

    return r_PQW, v_PQW # return states in orbital (PQW) frame

def cartesian2keplerian(r,v):
    """Compute conversion from cartesian state vector to keplerian orbital
    elements. Source for this conversion is found in the following reference:
    https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf"""
    # semi-major axis
    normr = np.linalg.norm(r)
    a = 1/((2/normr) - ((np.linalg.norm(v)**2)/MU))
    # orbital angular momentum
    h = np.cross(r,v)
    # inclination
    i = np.arccos(h[2]/np.linalg.norm(h))
    # eccentricity
    e_vec = np.cross(v,h)/MU - r/normr
    e = np.linalg.norm(e_vec)
    # n vector pointing towards ascending node
    n = np.array([-h[1],h[0],0])
    normn = np.linalg.norm(n)
    # right ascension
    raan = np.arccos(n[0]/normn)*(n[1]>=0)+(2*np.pi - np.arccos(n[0]/normn))*(n[1]<0)
    # argument of periapsis
    ap = np.arccos(np.dot(n,e_vec)/(normn*e))*(e_vec[2]>=0) + (2*np.pi -np.arccos(np.dot(n,e_vec)/(normn*e)))*(e_vec[2]<0)

    return a,e,i,raan,ap

def unit_vec(x):
    """Normalize vector x"""
    return x/np.linalg.norm(x)

def angle_vecs(v1,v2):
    """Compute angle between vectors v1 and v2 using dot product"""
    prod = np.dot(v1,v2)
    return np.arccos(prod)

def Rx(angle):
    """Rotation matrix around x axis"""
    return np.array([[1,0,0], [0,np.cos(angle),np.sin(angle)], [0,-np.sin(angle),np.cos(angle)]])

def Rz(angle):
    """Rotation matrix around z axis"""
    return np.array([[np.cos(angle),np.sin(angle),0], [-np.sin(angle),np.cos(angle),0], [0,0,1]])

def skew(x):
    """Generate skew matrix from vector x"""
    return np.array([[0,-x[2],x[1]], [x[2],0,-x[0]], [-x[1],x[0],0]])

def wrap2pi(angle):
    """Wrap angle between 0 and 2pi"""
    return angle + 2*np.pi*(angle<0)-2*np.pi*(angle>=2*np.pi) + 0*(not angle<0 and not angle>=2*np.pi)


