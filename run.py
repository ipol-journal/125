#!/usr/bin/env python3

import subprocess
import argparse

# parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("nbscales", type=int)
ap.add_argument("stdfactor", type=float)
args = ap.parse_args()

def get_num_channels(filename):
    '''
    Reads the specified file and returns its number of channels
    '''
    f = open(filename)
    line = f.readline()
    f.close()
    num_channels = len(line.split()) / 2
    return num_channels

def nc_get_min_max(filename):
    '''
    Returns the minimum and maximum ranges for
    the X and Y axes
    '''
    x_min, x_max, y_min, y_max = 1E9, -1E9, 1E9, -1E9

    # Read data
    with open(filename, 'r') as f:
        EOF = False
        lines = []
        while not EOF:
            line = f.readline()
            s = line.split()
            EOF = (line == '')
            lineStrip = line.strip()
            isWhiteLine = (not EOF and lineStrip == '')
            isComment = lineStrip.startswith('#')
            if not (EOF or isWhiteLine or isComment):
                lines.append(s)

    # Guess number number of channels
    numBins = len(lines)
    if numBins > 0:
        numChannels = int(len(lines[0])/2)
    else:
        numChannels = 0

    # Read data values
    for b in range(numBins):
        line = lines[b]
        #
        for ch in range(numChannels):
            x_value = line[ch].strip().upper()
            y_value = line[ch+numChannels].strip().upper()
            #
            if x_value.upper() != 'X' and y_value.upper() != 'X':
                x = float(x_value)
                y = float(y_value)
                if x < x_min:
                    x_min = x
                if y < y_min:
                    y_min = y
                if x > x_max:
                    x_max = x
                if y > y_max:
                    y_max = y
    #
    return x_min, x_max, y_min, y_max


# Run Noise Clinic
subprocess.run(['MSD', 'input_0.png', 'denoised.png', 'diff.png', str(args.nbscales), str(args.stdfactor), '0'])
# The last parameter is "verbose".

# Get number of channels
num_channels = get_num_channels('denoised_noiseCurves_H_0.txt')

# Get noise curve bounds
x_min, x_max, y_min, y_max = 1E9, -1E9, 1E9, -1E9
for i in range(args.nbscales):
    for C in ["L", "H"]:
        filename = f'denoised_noiseCurves_{C}_{i}.txt'
        x_min_image, x_max_image, \
            y_min_image, y_max_image = \
            nc_get_min_max(filename)
        x_min = min(x_min, x_min_image)
        x_max = max(x_max, x_max_image)
        y_min = min(y_min, y_min_image)
        y_max = max(y_max, y_max_image)

# Generate noise curve figures
for i in range(args.nbscales):
    filename = f'denoised_noiseCurves_H_{i}.txt'
    for C in ["L", "H"]:
        estimation = f'denoised_noiseCurves_{C}_{i}.txt'
        figure = f'denoised_noiseCurves_{C}_{i}.png'
        subprocess.run(["writeNoiseCurve.sh", str(estimation), str(num_channels), str(x_min*0.9), str(x_max*1.1), \
            str(y_min*0.9), str(y_max*1.1), figure])