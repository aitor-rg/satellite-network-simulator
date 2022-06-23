#!/usr/bin/env python3

"""Simulation + 3DVisualizer class"""

import os
import sys
import math
import numpy as np
import transformations as tf
#import vpython as vp
import pyqtgraph.opengl as gl

from pyqtgraph.Qt import QtCore, QtGui
from objects import Satellite, Ground, Ship, Area

class Simulator:
    def __init__(self, args):
        os.system('cls' if os.name == 'nt' else 'clear')
#INIT PARAMETERS ------------------------------------------------------
        # put seed to the pseudorandom variables
        #np.random.seed(0)
        #print("[Initializing]: parameters..", end=' ', flush=True)
        # simulation parameters
        self.dt = args.dt                       #simulation step time
        self.k = 0                              #simulation step
        # 3Dvisualization True/False
        self.vis = (args.v=='y')
        # number of ground terminals (stations and ships) and satellites in the satnet
        self.num_areas   = len(args.data_areas)
        self.num_gs      = len(args.data_gs)
        self.num_ships   = len(args.data_ships)
        self.satnet_size = len(args.data_satnet)
        # initialize list of ground stations, ships and satellite objects
        self.areas  = [Area(args.data_areas[i]) for i in range(self.num_areas)]
        self.gs     = [Ground(args.data_gs[i], args.dt) for i in range(self.num_gs)]
        self.ships  = [Ship(args.data_ships[i], args.dt) for i in range(self.num_ships)]
        self.satnet = [Satellite(args.data_satnet[i], args) for i in range(self.satnet_size)]

        #print("T span:{}h".format(args.T), end=' / ', flush=True)
        #print("dt:{}".format(self.dt), end=' / ', flush=True)
        #print("satnet size:{}".format(self.satnet_size), end=' / ', flush=True)
        #print("n gs:{}".format(self.num_gs), end=' / ', flush=True)
        #print("n ships:{}".format(self.num_ships),end=' / ', flush=True)
        #print("Done", end='\n', flush=True)
        print("[Logging]: {}".format(args.log))

        #release memory
        del(args)

#INIT 3D TOOLS --------------------------------------------------------
        if self.vis:
            print("[Initializing]: 3Dvisualization tools..", end=' ', flush=True)
            self.init_graphics()
            print("Done", flush=True)

#SIMULATION LOOP----------------------------------------------------------
    def update_sim(self):
        # simulation time
        t = self.k*self.dt
        #if np.mod(self.k,100)==0: print('\rt:{:.2f}h'.format(self.k*self.dt), end='')

        # Terminals loop
        for ter in self.gs + self.ships:
            # TERMINAL COMMANDS
            ter.generate_message(t) if ter.is_generating_message() else None

        # Satellites loop
        for sat in self.satnet:
            # time propagation
            sat.propagate_GST()            #update ECI to ECEF matrix
            sat.propagate_satellite_pos()  #propagate pos in ECI and transf. to ECEF
            sat.propagate_satellite_att()  #propagate satellite attitude in PQW
            sat.propagate_antennae_att()   #update antennae atttitude in ECEF

            # COMMUNICATION COMMANDS
            sat.check_connection_to(self.gs+self.ships) #calculate signal budgets and connections
            if sat.mode == 'idle': sat.IDLE(self.gs+self.ships)
            elif sat.mode == 'tx': sat.TX(t)
            elif sat.mode == 'rx': sat.RX(self.gs+self.ships)

            # BROADCAST AREA COMMANDS

            # LOGS
            sat.logging()

        # set how often to update the visualization
        if self.vis: self.update_graphics() #if (np.mod(self.k,1)==0) else None
        self.k += 1

#GRAPHICS ------------------------------------------------------------
    def update_graphics(self):
        #update graphics for each satellite in satnet
        for s,sat in enumerate(self.satnet):
            #update satellite position
            self.satnet_scatter[s].setData(pos=sat.x_I[0:3])
            #update antennae orientation and position
            v_1 = sat.x_I[0:3] + sat.SRS3.xyz_I[0:3]
            v_2 = sat.x_I[0:3] + sat.SRS4.xyz_I[0:3]
            SRS3 = np.array([sat.x_I[0:3].transpose(), v_1.transpose()])
            SRS4 = np.array([sat.x_I[0:3].transpose(), v_2.transpose()])
            self.srs3_line[s].setData(pos=SRS3)
            self.srs4_line[s].setData(pos=SRS4)

            #update satellite trail every 50 steps
            if (np.mod(self.k,50)==0) and (self.trails_enabled):
                self.trail[s] = np.delete(np.vstack((self.trail[s], sat.x_I[0:3])), 0, axis=0)
                self.orbit_line[s].setData(pos=self.trail[s])

            #if there is connection link draw a line from satellite to ground station
            for t,term in enumerate(self.gs+self.ships):
                v_3 = tf.spherical2cartesian(term.coordinates)
                if term.name in sat.whitelist:
                    show_link_condition = (sat.connected_to[term.name]!=0)*(sat.mode!='idle')*\
                        ((sat.txrx["receiver"]==term.name) + (sat.txrx["sender"]==term.name))
                else:
                    show_link_condition = False
                #show_link_condition = sat.connected_to[term.name]
                line = np.array([sat.x_I[0:3].transpose(),
                                 v_3.transpose()])*show_link_condition
                self.link_line[s][t].setData(pos=line, color=self.color[sat.mode])

        for t,term in enumerate(self.gs + self.ships):
            req = term.is_requesting_allocation()
            color = self.color['y']*(not(req)) + self.color['r']*req
            self.term_scatter[t].setData(color = color)

    def init_graphics(self):
        self.app = QtGui.QApplication(sys.argv)
        self.viz = gl.GLViewWidget()
        self.viz.setWindowTitle('Constellation simulator')
        self.viz.setGeometry(0, 0, int(1920/2), int(1080/2))
        self.viz.show()
        self.viz.opts['distance'] = 30000
        self.color = {
            'r':(255,0,0,2),
            'g':(0,255,0,2),
            'b':(0,0,255,2),
            'y':(255,170,0,2),
            'w':(255,255,255,2),
            'tx':(0,255,0,2),
            'rx':(255,0,0,2),
            'idle':(0,255,0,2)
            }
        self.trails_enabled = False;
        if (self.trails_enabled):
            trail_length = 15
            self.trail = [0]*self.satnet_size

        #draw globe
        mesh = gl.MeshData.sphere(rows=20, cols=30, radius=6371)
        earth = gl.GLMeshItem(
            meshdata=mesh,
            smooth=True,
            color=(0, 0, 250, 0.2),
            shader="shaded",
            glOptions="opaque")
        self.viz.addItem(earth)
        #draw ECEF Reference
        x_axis = gl.GLLinePlotItem(pos=np.array([[0,0,0],[3000,0,0]]), color=(255, 0, 0,2))
        y_axis = gl.GLLinePlotItem(pos=np.array([[0,0,0],[0,3000,0]]), color=(0, 255, 0,2))
        z_axis = gl.GLLinePlotItem(pos=np.array([[0,0,0],[0,0,3000]]), color=(0, 0, 255,2))
        self.viz.addItem(x_axis)
        self.viz.addItem(y_axis)
        self.viz.addItem(z_axis)
        #draw ground stations and ships as points on the surface of the sphere
        self.term_scatter=[0]*(self.num_gs+self.num_ships)
        for t,term in enumerate(self.gs+self.ships):
            lat,lon = term.coordinates
            lat = lat*np.pi/180
            lon = lon*np.pi/180
            xc = 6371*np.cos(lat)*np.cos(lon)
            yc = 6371*np.cos(lat)*np.sin(lon)
            zc = 6371*np.sin(lat)
            term_color = self.color['r']*(term.type=="station")+\
                         self.color['y']*(term.type=="ship")
            self.term_scatter[t] = gl.GLScatterPlotItem(pos=np.array([xc,yc,zc]),
                                                   color=term_color,
                                                   size=3, glOptions='translucent')
            self.viz.addItem(self.term_scatter[t])

        # draw broadcast areas
        for area in self.areas:
            area.get_line_points()
            for p in range(len(area.midpoints)-1):
                point1 = area.midpoints[p,:]
                point2 = area.midpoints[p+1,:]
                r1 = np.array([6371*np.cos(point1[0])*np.cos(point1[1]),
                               6371*np.cos(point1[0])*np.sin(point1[1]),
                               6371*np.sin(point1[0])])
                r2 = np.array([6371*np.cos(point2[0])*np.cos(point2[1]),
                               6371*np.cos(point2[0])*np.sin(point2[1]),
                               6371*np.sin(point2[0])])
                borders = gl.GLLinePlotItem(pos=np.vstack([r1,r2]),
                                          color=self.color['r'],
                                          glOptions='translucent')
                self.viz.addItem(borders)


        #initialize satellite point plot
        self.satnet_scatter = [0]*self.satnet_size
        self.orbit_line = [0]*self.satnet_size
        self.srs3_line=[0]*self.satnet_size
        self.srs4_line=[0]*self.satnet_size
        self.link_line = []
        for s,sat in enumerate(self.satnet):
            self.satnet_scatter[s] = gl.GLScatterPlotItem(pos=sat.x_I[0:3],
                                                          color=self.color["w"],
                                                          glOptions='translucent')
            self.viz.addItem(self.satnet_scatter[s])

            if (self.trails_enabled):
                #init satellite trail
                self.trail[s] = np.zeros((trail_length,3))
                self.trail[s] = np.delete(np.vstack((self.trail[s], sat.x_I[0:3])), 0, axis=0)
                self.orbit_line[s] = gl.GLLinePlotItem(pos=self.trail[s],
                                                       color=self.color['w'],
                                                       glOptions='translucent')
                self.viz.addItem(self.orbit_line[s])

            #initialize satellite antennae plot
            v_1 = sat.x_I[0:3] + sat.SRS3.xyz_I[0:3]
            v_2 = sat.x_I[0:3] + sat.SRS4.xyz_I[0:3]
            SRS3 = np.array([sat.x_I[0:3].transpose(), v_1.transpose()])
            SRS4 = np.array([sat.x_I[0:3].transpose(), v_2.transpose()])
            self.srs3_line[s] = gl.GLLinePlotItem(pos=SRS3, color=(255, 0, 255,2),glOptions='translucent')
            self.srs4_line[s] = gl.GLLinePlotItem(pos=SRS4, color=(255, 0, 255,2),glOptions='translucent')
            self.viz.addItem(self.srs3_line[s])
            self.viz.addItem(self.srs4_line[s])

            #initialize connection link plot
            self.link_line.append([])
            for t,term in enumerate(self.gs+self.ships):
                v_3 = tf.spherical2cartesian(term.coordinates)
                #if there is connection draw a line from satellite to ground station
                if term.name in sat.whitelist:
                    show_link_condition = (sat.connected_to[term.name]!=0)*\
                        ((sat.txrx["receiver"]==term.name) + (sat.txrx["sender"]==term.name))
                else:
                    show_link_condition = False
                line = np.array([sat.x_I[0:3].transpose(),
                                 v_3.transpose()])*show_link_condition
                self.link_line[s].append(gl.GLLinePlotItem(pos=line,
                                                           color=(0, 50, 0,2),
                                                           glOptions='translucent'))
                self.viz.addItem(self.link_line[s][t])

#---------------------------------------------------------------------

    def start_without_graphics(self):
        #print("\n[Simulation]: started.. Press [Ctrl+C] to end", end='\n\n', flush=True)
        self.finished = False
        while not self.finished:
            try:
                #t1=time.time()
                self.update_sim()
                #print(time.time()-t1)
            except(KeyboardInterrupt,SystemExit):
                print("\n\n[Simulation]: exited successfully", flush=True)
                break

    def start_with_graphics(self):
        print("\n[Simulation]: started.. Close window to end", end='\n\n', flush=True)
        self.timer = QtCore.QTimer()
        self.timer.setInterval(1)
        self.timer.timeout.connect(self.update_sim)
        self.timer.start()
        QtGui.QApplication.instance().exec_()
