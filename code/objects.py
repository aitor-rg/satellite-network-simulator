#!/usr/bin/env python3
"""Satellite, Antenna and Terminal classes"""

import math
import numpy as np
import transformations as tf
#import matplotlib.pyplot as plt

from logger import Logger
from scipy import integrate
from scipy.linalg import expm


class Satellite:
    def __init__(self, data, args):
        #timestep size in hours
        self.dt = args.dt
        #satellite id
        self.name = data["id"]
        #satellite color
        self.marker = data["marker"]
        #satellite operative
        self.operative = data["operative"]
        #terminals admitted
        self.whitelist = data["grounds"]+data["ships"]
        #keplerian elements
        self.SMAJ = data["SMAJ"]
        self.ECCE = data["ECCE"]
        self.INCL = data["INCL"]*np.pi/180
        self.RAAN = data["RAAN"]*np.pi/180
        self.ARGP = data["ARGP"]*np.pi/180
        self.TRUE = data["TRUE"]*np.pi/180
        self.MANO = data["MANO"]*np.pi/180
        self.EANO = 2*np.arctan2(np.sqrt(1-self.ECCE)*np.sin(self.TRUE/2),
                                 np.sqrt(1+self.ECCE)*np.cos(self.TRUE/2))

        #Greenwich Siderial Time angle
        self.angleGST = 0 #rad

        #relevant constants
        J2 = 1.0826157E-3
        Re = 6378.125
        self.MU = (3.986*10**5)*3600**2
        self.C = 3*(J2*self.MU*Re**2)/2

        #earth angular rotation
        self.earth_omega = np.array([0,0,(7.29*10**-2)*3.6])

        #initial linear position, velocity, angular velocity and orientations
        self.set_initial_conditions()

        #antennae class
        self.SRS3 = Antenna(4,"gmsk")
        self.SRS4 = Antenna(4,"bpsk")
        #initial antennae orientation in ECI reference
        self.propagate_antennae_att()

        #modes: idle/tx/rx
        self.mode = "idle"
        self.txrx = {'sender':None, 'receiver':None, 'speed': None, 'prio': None,
                     'tstamp': None, 'size': None, 't': 0, 't_total': 0}

        #dictionaries
        self.buffer = {0:[], 1:[], 2:[]} #message buffer
        self.elev = {terminal:0 for terminal in self.whitelist} #current elevation
        self.connected_to = {terminal:0 for terminal in self.whitelist} #current connection

        self.logger = {l:Logger(self.name, l) for l in args.log.keys() if args.log[l]==True}

        #propagator: dynamics/kinematics
        #if args.prop == 'dyn':
        #if self.name == "Sternula-2":
            #self.propagate_satellite_pos = self.dynamics
        #else:
        self.propagate_satellite_pos = self.kinematics

    def set_initial_conditions(self):
        """Initial linear position, velocity and angular velocity
        conditions for the S/C. The angular position is not represented here
        since we keep track of the orientation with the rotation matrix R_os"""

        #transformation matrices
        self.R_ie = self.ECEF()     #ECI   --> ECEF
        self.R_oi = self.PQW()      #Orbit --> ECI
        self.R_so = self.SC()       #S/C   --> Orbit
        #initial S/C position and velocity in orbital frame from keplerian elements
        r_PQW,v_PQW = tf.keplerian2cartesian(self.SMAJ,self.ECCE,self.EANO,self.TRUE)
        #transformation from orbital to eci to ecef frame
        self.x_I = np.append(self.eci2ecef(self.orb2eci(r_PQW)),
                             self.eci2ecef(self.orb2eci(v_PQW)))
        #initial S/C angular velocity matching orbit periodicity
        #self.omega_O = self.set_angular_vel()
        #uncomment to set a randomized S/C angular velocity between [0, 2rpm]
        self.omega_O = self.set_rand_angular_vel()


    def set_rand_attitude(self):
        """Set a randomized initial orientation of the satellite wrt orbit"""
        #we employ Euler angles yaw-roll-yaw'
        Z1 = np.random.uniform(0,2*np.pi) #random yaw
        X1 = np.random.uniform(0,2*np.pi) #random roll
        Z2 = np.random.uniform(0,2*np.pi) #random yaw'
        #and rotation ZXZ'
        R_os = np.matmul(tf.Rz(Z2),np.matmul(tf.Rx(X1),tf.Rz(Z1)))
        return R_os

    def set_rand_angular_vel(self):
        """Set a fixed initial random angular velocity to the satellite
        characteristic from LEOP, i.e. <2rpm"""
        randdir = (np.random.rand(3,1)-0.5)*2
        randmod = np.random.uniform(0,2*np.pi*60)
        return np.array(randdir/np.linalg.norm(randdir)*randmod)

    def set_attitude():
        """Set a user-defined initial orientation of the satellite wrt orbit"""
        nadir_offset = 66 #SRS-4 offset from nadir-mode in degrees
        R_os = tf.Rz(self.TRUE - nadir_offset*np.pi/180)
        return R_os

    def set_angular_vel():
        """Set a user-defined initial angular velocity to the satellite wrt orbit
        frame. Now defined equal to the orbital period."""
        orbit_vel = np.sqrt(self.MU/self.SMAJ) #linear orbital velocity
        orbit_omega_O = orbit_vel/self.SMAJ #angular orbital velocity
        omega = np.array([0,0,orbit_omega_O]) #angular velocity in PQW frame
        return omega

    def propagate_GST(self):
        """Compute next GST angle"""
        #update GST angle
        self.angleGST = tf.wrap2pi(self.angleGST + self.dt*self.earth_omega[2])
        #update ECI to ECEF transformation matrix
        self.R_ie = self.ECEF()

    def propagate_satellite_att(self):
        """Next satellite orientation in Orbit reference"""
        self.R_so = np.matmul(expm(tf.skew(self.omega_O*self.dt)),self.R_so)

    def propagate_antennae_att(self):
        """Next antennae orientation in ECEF reference"""
        self.SRS3.xyz_I = self.eci2ecef(self.orb2eci(self.sat2orb(self.SRS3.xyz_S)))
        self.SRS4.xyz_I = self.eci2ecef(self.orb2eci(self.sat2orb(self.SRS4.xyz_S)))

    def kinematics(self):
        """Solve the Kepler equations for Elliptic Orbits using the Mean Anomaly M(t)
        and Eccentric Anomaly E(t) using the Newton-Raphson method"""
        self.MANO = tf.wrap2pi(self.MANO + self.dt*np.sqrt(self.MU/(self.SMAJ**3)))
        self.EANO = self.NewtonRaphson(self.MANO)
        self.TRUE = 2*np.arctan2(np.sqrt(1+self.ECCE)*np.sin(self.EANO/2),
                                 np.sqrt(1-self.ECCE)*np.cos(self.EANO/2))

        r_PQW,v_PQW = tf.keplerian2cartesian(self.SMAJ,self.ECCE,self.EANO,self.TRUE)
        self.x_I = np.append(self.eci2ecef(self.orb2eci(r_PQW)),
                             self.eci2ecef(self.orb2eci(v_PQW)))

    def NewtonRaphson(self, M):
        """Solve for E(t) equation E(t) - e*sin(E(t)) = M(t) iteratively. Use the
        tolerance variable TOL to twik the method accuracy"""
        E = M #Eccentric anomaly initialized as Mean anomaly, i.e. first guess.
        TOL = .001
        while(np.abs(E - self.ECCE*np.sin(E) - M) > TOL):
            E = E - 0.3*(E - self.ECCE*np.sin(E) - M)/(1 - self.ECCE*np.cos(E))
        return E

    def dynamics(self):
        """Next satellite position and velocity solving the ODE with order 4 R-K"""
        sol = integrate.odeint(self.odeECEF, self.x_I, np.array([0,self.dt]))
        self.x_I = sol[-1,:]

    def odeECEF(self, x, t):
        """Set of differential equations that describe the motion of the satellite
        under J2 perturbation expressed in the ECEF reference frame"""
        r = x[0:3]  #position in ECEF
        v = x[3:6]  #velocity in ECEF
        a = np.zeros(3)
        normr = np.linalg.norm(r)

        #coriolis
        a += -np.cross(2*self.earth_omega, v)
        #centrifugal
        a += -np.cross(self.earth_omega, np.cross(self.earth_omega, r))
        #gravity potential
        a += -r*self.MU/normr**3
        #J2
        a += (self.C/normr**7)*np.array([-r[0]**3 -r[0]*r[1]**2 +4*r[0]*r[2]**2,
                                         -r[1]*r[0]**2 -r[1]**3 +4*r[1]*r[2]**2,
                                         -3*r[2]*r[0]**2 -r[2]*r[1]**2 +2*r[2]**3])
        return np.append(v, a)

    def calculate_elevations(self, terminal):
        """Calculates the elevation angle over the horizon of the satellite with
        respect to the ground station/ship"""
        #lat,lon = tf.cartesian2spherical(self.x_I[0:3])

        v1 = self.pointing_to_terminal(terminal)
        v2 = tf.unit_vec(terminal.xyz_I)

        angle = tf.angle_vecs(v1,v2)
        #elevation is the complementary angle to 90deg
        self.elev[terminal.name] = (angle-np.pi/2)*(angle>=np.pi/2)+0*(angle<np.pi/2)

    def check_connection_to(self, terminals):
        """Computes the signal budget of the satellite antennae and determines
        if there is enough signal to stablish a connection link"""
        for terminal in terminals:
            if self.connectivity_conditions(terminal):
                #compute satellite elevation over terminal horizon
                self.calculate_elevations(terminal)
                #check elevation
                if self.elev[terminal.name]>0:
                    #unitary vector pointing from satellite to terminal
                    v_1 = self.pointing_to_terminal(terminal)
                    #unitary vectors of satellite antennae direction
                    v_a1 = tf.unit_vec(self.SRS3.xyz_I)
                    v_a2 = tf.unit_vec(self.SRS4.xyz_I)
                    #angle between satellite antenna and terminal center
                    beta_1 = tf.angle_vecs(v_1, v_a1)
                    beta_2 = tf.angle_vecs(v_1, v_a2)
                    #angle between ship horizon plane and satellite
                    alpha = self.elev[terminal.name]*(terminal.type=='ship')\
                            + (-1)*(terminal.type=='station')
                    #distance from satellite to terminal
                    dist = self.distance_to_terminal(terminal)
                    #compute signal budgets
                    GT = terminal.GT #terminal G/T noise
                    signal_1 = self.SRS3.signal_budget(alpha, beta_1, dist, GT)
                    signal_2 = self.SRS4.signal_budget(alpha, beta_2, dist, GT)
                    #check signal above threshold & terminal within apperture
                    link_1 = self.check_threshold(signal_1, self.SRS3.threshold)*self.check_fov(beta_1)
                    link_2 = self.check_threshold(signal_2, self.SRS4.threshold)*self.check_fov(beta_2)
                    #keep the highest value
                    self.connected_to[terminal.name] = link_1|link_2   #max(link_1, link_2)
                else:
                    #if satelite below terminal horizon
                    self.connected_to[terminal.name] = 0
            else:
                #if terminal is not whitelisted or inoperative
                self.connected_to[terminal.name] = 0

    def connectivity_conditions(self, terminal):
        #collect message destination names awaiting in buffers
        buffer = self.buffer[2]+self.buffer[1]+self.buffer[0]
        destinations = [message["receiver"] for message in buffer]
        #check if terminal is whitelisted and  operative
        cond0 = terminal.name in self.whitelist
        cond2 = terminal.name in destinations
        cond1 = terminal.is_requesting_allocation()
        cond3 = (terminal.name == self.txrx["receiver"])
        cond4 = (terminal.name == self.txrx["sender"])
        leop = (terminal.type == "station")

        return cond0&(cond1|cond2|cond3|cond4|leop)

    def check_threshold(self, signal, threshold):
        """Check signal budget above Es/N0 threshold of antenna"""
        return (signal>=threshold)

    def check_fov(self, angle):
        """Check that terminal is inside the field of view of the antenna"""
        return (angle<np.pi/2)

    def pointing_to_terminal(self, terminal):
        """Return unitary vector pointing from satellite to terminal"""
        return tf.unit_vec(terminal.xyz_I-self.x_I[0:3])

    def distance_to_terminal(self, terminal):
        """Return range from satellite to terminal"""
        return np.linalg.norm(terminal.xyz_I-self.x_I[0:3])

    def ready_to_transmit(self):
        """Generate a binary array that indicate messages in the buffer destined
        to terminals that are currently connected to this spacecraft. The array
        is ordered from high to low priority messages"""
        #array of messages in the buffer sorted from high to low priority
        buffer = self.buffer[2]+self.buffer[1]+self.buffer[0]
        #create the binary array of buffer
        return [self.connected_to[msg['receiver']] for msg in buffer]*(len(buffer)>0) + []

    def rx_param(self, message, allocation_time):
        """Uplink transmission parameters"""
        DC = 6 #number of dedicated channels
        self.txrx['speed']    = 600/DC
        self.txrx['sender']   = message['sender']
        self.txrx['receiver'] = message['receiver']
        self.txrx['prio']     = message['prio']
        self.txrx['size']     = message['size']
        self.txrx['tstamp']   = message['tstamp']
        self.txrx['t']        = 0
        self.txrx['t_total']  = (float(message['size'])*8)/(self.txrx['speed']*3600) + allocation_time

    def tx_param(self, message):
        """Downlink transmission parameters"""
        self.txrx['speed']    = 65 #downlink speed Kbit/s
        self.txrx['sender']   = message['sender']
        self.txrx['receiver'] = message['receiver']
        self.txrx['prio']     = message['prio']
        self.txrx['size']     = message['size']
        self.txrx['tstamp']   = message['tstamp']
        self.txrx['t']        = 0
        self.txrx['t_total']  = (float(message['size'])*8)/(self.txrx['speed']*3600)

    def IDLE(self, terminals):
        """Decide weather next operation is transmitting (TX) or receiving (RX) a
        private message. We prioritize TX before RX a message. Therefore first it
        is checked if any message can be transmitted from satellite. Otherwise,
        start listening to possible requests from ground station and ship terminals"""

        #first check if there is a connection stablished with any terminal
        connections = self.connected_to.values()
        if any(connections):
            #check if any candidate message to be sent to linked terminal
            waiting_messages = self.ready_to_transmit()
            if any(waiting_messages):
                #pick index i of first message in buffer waiting to be sent
                i = next(j[0] for j in enumerate(waiting_messages) if j[1]>0)
                #retrieve priority and index of the selected message in buffer based i
                prio, index = self.retrieve_priority_and_index(i)
                #change mode to transmit
                self.mode = "tx"
                self.tx_param(self.buffer[prio].pop(index))

            else:
                #create binary vector of connected terminals that are requesting allocation
                waiting_terminals = [self.connected_to[t.name]&t.is_requesting_allocation() for t in terminals]
                if any(waiting_terminals):
                    #pick index i of first terminal in bterm different than 0
                    i = next(j[0] for j in enumerate(waiting_terminals) if j[1]>0)
                    #take terminal's message
                    message = terminals[i].request_allocation()
                    if (message["sender"] in self.whitelist and message["receiver"] in self.whitelist):
                        allocation_size = 20 #Kbit
                        allocation_time = (allocation_size/3600)*(terminals[i].type=="ship")
                        self.rx_param(message,allocation_time)
                        #start receiving message
                        self.mode = "rx"
                    else:
                        #return message to terminal
                        self.return_message_to_terminal(terminals[i])

    def downlink_is(self, status, timedeliver):
        """Print downlink message status in terminal"""
        color = "\033[91m"*(status=='FAILED') + "\033[92m"*(status=='DELIVERED')
        time = "["+str(round(timedeliver - self.txrx['tstamp'],2)*60)+"min]"*(status=='DELIVERED') + ""
        print("\033[1;37m\u2193 LINK\033[0;37m\tSENDER: {}\tCARRIER: {}\tRECEIVER: {}".format(self.txrx['sender'], self.name ,self.txrx['receiver']), end=" ", flush=True)
        print(color+status,end=" ", flush=True)
        print(time+"\033[00m",flush=True)

    def uplink_is(self, status):
        """Print uplink message status in terminal"""
        color = "\033[91m"*(status=='FAILED') + "\033[92m"*(status=='RECEIVED')
        print("\033[1;30m\u2191 LINK\033[0;37m\tSENDER: {}\tCARRIER: {}\tRECEIVER: {}".format(self.txrx['sender'], self.name ,self.txrx['receiver']), end=" ", flush=True)
        print(color+status+"\033[00m",flush=True)

    def TX(self, t):
        """"""
        #get status of transmission
        status = self.transmission_status(self.txrx['receiver']) #return 0, 1 or 2
        if status == 0: #transmitting
            self.transmit_packet(self.txrx['receiver']) #to receiver
        elif status == 1: #delivered
            self.store_tx_delay(t)
            self.downlink_is('DELIVERED',t)
            self.mode = 'idle'
        else: #failed
            self.return_message_to_buffer()
            self.sort_buffer()
            self.downlink_is('FAILED',0)
            self.mode = 'idle'

    def RX(self, terminals):
        """"""
        #get status of transmission
        status = self.transmission_status(self.txrx['sender']) #return 0, 1 or 2
        if status == 0: #transmitting
            self.transmit_packet(self.txrx['sender']) #from sender
        elif status == 1: #received
            self.add_message_to_buffer()
            self.sort_buffer()
            self.uplink_is('RECEIVED')
            self.mode = "idle"
        else: #failed
            #get index of connected terminal
            i = next(j for j,terminal in enumerate(terminals) if terminal.name==self.txrx['sender'])
            self.return_message_to_terminal(terminals[i])
            self.uplink_is('FAILED')
            self.mode = "idle"

    def transmission_status(self, terminal):
        """Use binary conditions (b1,b2,b3) and the following equation to determine
        status of the transmission. Equation: 0*b1 + 1*b2 + 2*b1*b3, with the only
        possible cases shown below.
            not succeded(b1)    succeeded(b2)   failed(b3)
                0                   1               0           return 1
                0                   1               1           return 1
                1                   0               0           return 0
                1                   0               1           return 2
        return 0 = proceed with transmission
        return 1 = transmission succeeded
        return 2 = transmission failed"""
        success_condition = (self.txrx["t"] >= self.txrx["t_total"])
        proceed_condition = not(success_condition)
        fail_condition    = not(self.connected_to[terminal])
        return 0*proceed_condition + 1*success_condition + 2*proceed_condition*fail_condition

    def transmit_packet(self, terminal_name):
        """Transmit 1 bit of information and generate random error. If error produced
        the transmission of the bit will be regected and will be attempted to be resent"""
        #generate random bit error (bit error if <5)
        bit_error = np.random.uniform(0,100)
        #if biterror or link lost reattempt last transmission
        self.txrx["t"] = self.txrx["t"] + self.dt - self.dt*(bit_error<5)*(self.connected_to[terminal_name])
        #print("."*(bit_error>=5) + "b"*(bit_error<5),end="",flush=True)

    def add_message_to_buffer(self):
        """If message successfully transmited, add to corresponding satellite
        priority buffer"""
        message = self.restore_message()
        self.buffer[self.txrx['prio']].append(message)

    def restore_message(self):
        """If transmission failed then restore message in satellite buffer"""
        message = {'sender'   : self.txrx['sender'],
                   'receiver' : self.txrx['receiver'],
                   'size'     : self.txrx['size'],
                   'prio'     : self.txrx['prio'],
                   'tstamp'   : self.txrx['tstamp']}
        return message

    def return_message_to_buffer(self):
        message = self.restore_message()
        self.buffer[self.txrx['prio']].insert(0, message)

    def return_message_to_terminal(self, terminal):
        """If transmission failed then return message to terminal buffer"""
        message = self.restore_message()
        terminal.buffer[self.txrx['prio']].insert(0, message)

    def sort_buffer(self):
        """Sort messages in satellite buffer in order of earliest created"""
        def sort_criteria(message):
            return message['tstamp'] # time stamps
        self.buffer[self.txrx["prio"]].sort(key=sort_criteria)

    def retrieve_priority_and_index(self, index):
        prio = 2*(index<len(self.buffer[2])) +\
               1*(len(self.buffer[2])<=index)*\
                 (index<len(self.buffer[2] + self.buffer[1])) +\
               0*(len(self.buffer[2] + self.buffer[1])<=index)*\
                 (index<len(self.buffer[2] + self.buffer[1] + self.buffer[0]))
        index = (index)*(prio == 2) +\
                (index-len(self.buffer[2]))*(prio == 1) +\
                (index-len(self.buffer[2] + self.buffer[1]))*(prio == 0)
        return prio, index

    def store_tx_delay(self, tstamp):
        #self.delays.append([tstamp - self.txrx["tstamp0"], tstamp - self.txrx["t"]])
        #self.delays.append(tstamp - float(self.txrx["tstamp0"]))
        pass

    def add_terminal(self, sat_name, gs_name):
        self.whitelist[gs_name] = bool(self.whitelist[gs_name] + (self.name==sat_name))

    def rm_terminal(self, sat_name, gs_name):
        self.whitelist[gs_name] = bool(self.whitelist[gs_name]*(not self.name==sat_name))

    def logging(self):
        for l in self.logger.keys():
            if l == 'elevations':
                self.logger[l].log([self.elev['AAL'], self.elev['SVL'], self.elev['TRL']])
            if l == 'delays':
                pass

    def SC(self):
        """Initial Spacecraft reference frame relative to PQW frame. Describes
        Spacecraft orientation seen from the PQW frame."""
        R_os = self.set_rand_attitude() # set initial attitude to be random
        #R_os = self.set_attitude() # set intial attitude as desired
        return np.linalg.inv(R_os)

    def PQW(self):
        """Transformation matrix from ECI to Orbit (PQW) reference frames
        that is R_oi = Rz(-RAAN)Rx(-i)Rz(-ARGP)"""
        M3 = tf.Rz(-self.ARGP)
        M2 = tf.Rx(-self.INCL)
        M1 = tf.Rz(-self.RAAN)
        return np.matmul(M1,np.matmul(M2,M3))

    def ECEF(self):
        """Transformation matrix from ECI to ECEF reference frames. Update matrix
        each time step."""
        return tf.Rz(self.angleGST)

    def eci2ecef(self, v):
        return np.matmul(self.R_ie, v)

    def orb2eci(self, v):
        return np.matmul(self.R_oi, v)

    def sat2orb(self, v):
        return np.matmul(self.R_so, v)


class Antenna:
    def __init__(self, max_gain, modulation):
        #antenna tip position wrt ECI/ECEF reference
        self.xyz_I = np.array([0,0,0])
        #antenna tip position wrt satellite reference (vector fixed in SC)
        self.xyz_S = np.array([500,0,0])*(modulation=="gmsk") + np.array([-500,0,0])*(modulation=="bpsk")

        #bandwidth [dB]: gmsk->512kHz, bpsk->500kHz
        self.bandwidth = 57.09*(modulation=="gmsk") + 57*(modulation=="bpsk")
        #minimum Es/N0 for <1% PER [dB]
        self.threshold = 13.5*(modulation=="gmsk") + 3*(modulation=="bpsk")
        #dB antena taking into account the pfd limit [dBW]
        self.TX = -1*(max_gain==4) + (-3)*(max_gain==6)

        #sampled points on true radiation pattern diagram x: degrees, y: dBs
        x = np.array(([-11,27.85,-74.93,-150,90,-180,180,180,180,180,180,180,180]))*(max_gain==4) + np.array(([0,15,30,45,60,75,85,-15,-30,-45,-60,-75,-85]))*(max_gain==6)
        y = np.array(([4.73,1.6,1.6,-10,-5,-30,-30,-30,-30,-30,-30,-30,-30]))*(max_gain==4) + np.array(([6,5.5,4.5,2.5,-1,-6.5,-11,5.5,4.5,2.5,-1,-6.5,-11]))*(max_gain==6)

        #radiation pattern approx by polynomial fitting using previous points
        self.radiation_poly = np.poly1d(np.polyfit(x, y, 6))
        #noise density [dBW/Hz/K]
        self.boltzman = -228.6
        #0.5-pointing error, 0.5-polarization error, 0.2-rain loss, 0.3-random [dB]
        self.losses = 1.5

    def gain(self, angle):
        """Compute antenna gain with angle between satellite antenna and ground
        station antenna using the polynomial fit of the antenna radiation."""
        return self.radiation_poly(angle*180/np.pi)

    def path_loss(self, dist):
        """Compute path loss with distance from satellite to ground station"""
        return 20*np.log10(dist*1000)+39.44

    def signal_budget(self, alpha, beta, dist, gs_GT):
        """Compute signal budget Es/N0. The angle alpha only applies to ships"""
        return (self.TX + self.gain(alpha)*(alpha!=-1) + self.gain(beta) - self.bandwidth - self.path_loss(dist) + gs_GT - self.boltzman - self.losses)


class Terminal:
    def __init__(self, data, dt):
        #sampling time
        self.dt = dt
        #terminal data
        self.name = data["name"]                            #name of the station
        self.coordinates = np.array(data["coordinates"])    #lat and lon [degrees]
        self.operative = data["operative"]                  #operative:1, inoperative:0
        #xyz position of the terminal in ECEF reference
        self.xyz_I = tf.spherical2cartesian(self.coordinates)

    #def switch_on(self):
    #    self.operative=True

    #def switch_off(self):
    #    self.operative=False

    def is_generating_message(self):
        """Sample a random uniform variable between 0 and 1. If the sampled value
        is below or equal to the prob_generate_message and the terminal has at
        least one client, then return True. Otherwise, return False."""
        return (np.random.uniform(0,1)<=self.prob_generate_message)*(len(self.client_list)>0)

    def generate_message(self, tstamp):
        """ rclient: selects a client from client_list randomly
            rsize: random value between 1 and 3 determining slot size of message
            rtype: random value determining priority of msg
                0= priority 0 msg (60%)
                1= priority 1 msg (35%)
                2= priority 2 msg (5%)"""
        rtype = np.random.uniform(0,10)                         #choose random msg priority
        rsize = np.random.randint(1, self.max_message_size)     #choose random msg size
        rclient = np.random.randint(0, len(self.client_list))   #choose random client
        #associate priority to message
        priority = 0*(rtype<=self.prob_P0) +\
                   1*(rtype>self.prob_P0)*(rtype<=self.prob_P1) +\
                   2*(rtype>self.prob_P1)*(rtype<=10)
        #store message in dictionary
        #self.buffer[priority].append([self.client_list[rclient], rsize, tstamp])
        self.buffer[priority].append({'sender'  : self.name,
                                      'receiver': self.client_list[rclient],
                                      'size'    : rsize,
                                      'prio'    : priority,
                                      'tstamp'  : tstamp})

    def generate_predefined_message(self, tstamp):
        """ Generate top priority message to send to injected satellite in order
        to detumble it. Size of the message is 10KiB which is approximated by
        10KB, it is put to priority 2 and generated at tstamp=0"""
        # generate a predefined message. This is used for LEOP
        self.buffer[2].append([self.name, 10, tstamp])
        print(self.name+" generated the LEOP message")

    def is_requesting_allocation(self):
        return len(self.buffer[2]+self.buffer[1]+self.buffer[0])>0

    def request_allocation(self):
        """Return message requested to be sent to satellite"""
        #figure out what message to send according to priority
        prio = 2*(len(self.buffer[2])!=0) +\
               1*(len(self.buffer[2])==0)*(len(self.buffer[1])!=0)+\
               0*(len(self.buffer[2]+self.buffer[1])==0)*(len(self.buffer[0])!=0)
        return self.buffer[prio].pop(0)


class Ground(Terminal):
    def __init__(self, data, dt):
        Terminal.__init__(self, data, dt)
        self.type = "station"
        self.GT = data["G/T"] #G/T: performance metric [dB/K]

        #self.sat_list = data["satellites"]
        self.client_list = data["clients"]

        #mms
        self.prob_P0 = 6
        self.prob_P1 = 9.5
        """assuming 1msg/6h, we consider a poisson distribution with mean=1/6
        and approximate the probability of 1 message being generated (k=1) as the
        probability of more than 1 message happening (k>=1) since P(k>1) is very
        low for small means and dt's. Then,
            P(1msg)=P(>=1msg)=1 - P(<1msg)=1 - mean^k*exp(-mean*dt)/k!
        where k is the number of events happening in dt"""
        self.prob_generate_message = 1-np.exp((-1/6)*self.dt)

        self.max_message_size = 100 #KB
        self.buffer = {0:[], 1:[], 2:[]}


class Ship(Terminal):
    def __init__(self, data, dt):
        Terminal.__init__(self, data, dt)
        self.type = "ship"
        self.GT = data["G/T"] #G/T: performance metric [dB/K]

        #self.sat_list = data["satellites"]
        self.client_list = data["servers"]

        #mms
        self.prob_P0 = 6
        self.prob_P1 = 9.5
        """assuming 1msg/6h, we consider a poisson distribution with mean=1/6
        and approximate the probability of 1 message being generated (k=1) as the
        probability of more than 1 event happening (k>=1) since P(k>1) is very
        low for small means and dt's. Then,
            P(1msg)=P(>=1msg)=1 - P(<1msg)=1 - (mean^k)*exp(-mean*dt)/k!
        where k is the number of events happening in dt"""
        self.prob_generate_message = 1-np.exp((-1/6)*self.dt)

        self.max_message_size = 10 #KB
        self.buffer = {0:[], 1:[], 2:[]}

class Area():
    def __init__(self, data):
        self.points = np.array(data["points"])

    def get_line_points(self):
        # save middle points
        self.midpoints = np.zeros((1,2))
        # add initial point at end to close cycle
        self.points = np.vstack([self.points, self.points[0,:]])
        # get points between vertices to plot line on Earth
        for p in range(len(self.points)-1):
            theta1 = self.points[p,:]*np.pi/180 # Point A
            theta2 = self.points[p+1,:]*np.pi/180 # Point B
            self.midpoints = np.vstack([self.midpoints, np.linspace(theta1,theta2,25)]) # Points from A to B
        self.midpoints = self.midpoints[1:,:]
