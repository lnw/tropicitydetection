#!/usr/bin/env python3

import os
import signal
import sys
import subprocess
from argparse import ArgumentParser

import libtropicity as tr


def parseCommandline():
    parser = ArgumentParser(prog="tropicity.py", description="nothing yet", add_help=True)
    parser.add_argument("--input", help="the vti file", dest="inputfile", action="store", type=str, required=True)
    parser.add_argument("--output", help="a textfile", dest="outputfile", action="store", type=str, required=True)
    parser.add_argument("--bfielddir", help="{0,1,2,3,4,5} -> {+x, -x, +y, -y, +z, -z}", dest="bfielddir", action="store", type=int, required=True)
    parser.add_argument("--perpdir", help="vector perpendicular to the plane; {0,1,2} -> {x, y, z}", dest="perpdir", action="store", type=int, required=True)
    parser.add_argument("--perpcoord", help="position of the plane in the coord parpendicular to it", dest="perpcoord", action="store", type=float, required=True)
    argparse = parser.parse_args()
    return argparse


def main():
    args = parseCommandline()
    print(args)

    C = tr.cube(args.inputfile)
    tropplane = C.gettropplane(args.bfielddir, args.perpdir, args.perpcoord)
    C.writetropplane(args.outputfile, tropplane)


if __name__ == "__main__":
    main()

