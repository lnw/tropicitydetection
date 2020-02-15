#!/usr/bin/env python3

import os
import signal
import sys
import subprocess
from argparse import ArgumentParser

sys.path.append('build/src')
import libtropicity as tr


def parseCommandline():
    parser = ArgumentParser(prog="tropicity.py", description="nothing yet", add_help=True)
    parser.add_argument("--input", help="the vti file", dest="inputfile", action="store", type=str, required=True)
    parser.add_argument("--grid", help="coordinates of the grid to be split", dest="grid", action="store", type=str, required=True)
    parser.add_argument("--weights", help="weights of the grid to be split", dest="weights", action="store", type=str, required=True)
    parser.add_argument("--bfielddir", help="{0,1,2,3,4,5} -> {+x, -x, +y, -y, +z, -z}", dest="bfielddir", action="store", type=int, required=True)
    argparse = parser.parse_args()
    return argparse


def main():
    args = parseCommandline()
    print(args)

    C = tr.cube(args.inputfile)
    C.splitgrid(args.grid, args.weights, args.bfielddir)


if __name__ == "__main__":
    main()

