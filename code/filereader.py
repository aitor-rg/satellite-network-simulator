#!/usr/bin/env python3.5

import numpy as np

def read_json(file_name):
    """Return json file content in a list with dictionaries"""
    import json
    with open(file_name) as json_file:
        content = json.load(json_file)
        return [item for item in content]

def read_csv(file_name):
    """Return csv file content in a list"""
    import csv
    with open(file_name, 'r') as csv_file:
        content = csv.reader(csv_file)
        return [line for line in content]

def unpack_csv(content):
    """Turn list returned by read_csv into array for easy manipulation"""
    header = content.pop(0) # a header assumed in the csv file
    ncols = len(header) # number of columns of csv file based on header
    nrows = len(content) # number of rows in data

    # unpack into an array
    data = np.zeros((nrows,ncols))
    for row in range(nrows): # for each column of data do
        data[row,:] = content[row][:]
    return header, data

