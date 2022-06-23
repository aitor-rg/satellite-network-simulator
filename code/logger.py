#!/usr/bin/env python3

import csv
import time
import datetime

class Logger:
    def __init__(self, name, table):
        #vectors with elevation angles over horizon in degrees for each gs
        fheader = ['Aalborg','Svalbard','Troll']
        self.csv_file = csv.writer(open("./logging/"+table+"_"+name+".csv","a"))
        self.csv_file.writerow(fheader)

    def log(self,values):
        self.csv_file.writerow(values)
