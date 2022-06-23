#!/usr/bin/env python3

import sys
import argparse
import simulator
import filereader
import numpy as np


def main(args):
    #read terminals, constellation, and sea areas data from json files
    args.data_areas  = filereader.read_json(args.a)
    args.data_gs     = filereader.read_json(args.gs)
    args.data_ships  = filereader.read_json(args.sh)
    args.data_satnet = filereader.read_json(args.c)
    args.log = {
        "elevations": True,
        "link"      : False,
        "buffer"    : False,
        "delays"    : False,
        "COE"       : False
    }

    #initialize simulation
    sim = simulator.Simulator(args)
    #start simulation
    if args.v=='y': sim.start_with_graphics()
    else:           sim.start_without_graphics()

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='''
    .
    ''')
    #simulation setup parameters
    parser.add_argument('-dt', type=float, default=0.001,
                        help="step time of the simulation in hours. Default: 0.001h")
    parser.add_argument('-v', type=str, default='n',
                        help="visualize simulation <y/n>. Default: n")

    #terminals information
    parser.add_argument('-a', type=str, default="areas.json",
                        help="file [.json] name containing coordinates of convex area.")
    #terminals information
    parser.add_argument('-gs', type=str, default="gs.json",
                        help="file [.json] name containing ground station names, coordinates and G/T. See gs.json for an example.")
    parser.add_argument('-sh', type=str, default="ships.json",
                        help="file [.json] name containing ships names and coordinates. See ships.json for an example.")
    #constellation information
    parser.add_argument('-c', type=str, default='satellite.json',
                        help="file [.json] name containing keplerian elements and other relevant information of each satellite. See constellation.json for an example.")

    args = parser.parse_args()

    #check that arguments are correctly introduced
    if args.v!='y' and args.v!='n':
        sys.stderr.write('Argument [-v] must be either y or n\n')
        sys.exit(1)

    #call main function
    main(args)
