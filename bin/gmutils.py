#!/usr/bin/env python

import sys
import csv

class Subintervals():
    infile = ""
    size = 10000

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-s":
                self.size = int(a)
                prev = ""
            elif a in ["-s"]:
                prev = a
            else:
                self.infile = a

    def run(self):
        with open(self.infile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                chrom = row[0]
                start = int(row[1])
                end   = int(row[2])
                a     = start
                b     = self.size
                while True:
                    if b > end:
                        sys.stdout.write("{}\t{}\t{}\n".format(chrom, a, end))
                        break
                    else:
                        sys.stdout.write("{}\t{}\t{}\n".format(chrom, a, b))
                        a = b
                        b += self.size

class bedToCoords():
    infile = ""
    column = 0

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-c":
                self.column = int(a) - 1
                prev = ""
            elif a in ["-c"]:
                prev = a
            else:
                self.infile = a
            
    def run(self):
        with open(self.infile) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                coords = row[self.column] + ":" + row[self.column+1] + "-" + row[self.column+2]
                outrow = row[:self.column] + [coords] + row[self.column+3:]
                sys.stdout.write("\t".join(outrow) + "\n")

class coordsToBed(bedToCoords):

    def run(self):
        with open(self.infile) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                coords = row[self.column]
                chrom, rest = coords.split(":")
                start, end = rest.split("-")
                outrow = row[:self.column] + [chrom, start, end] + row[self.column+3:]
                sys.stdout.write("\t".join(outrow) + "\n")

COMMANDS = {"subintervals": Subintervals,
            "bedtocoords": bedToCoords,
            "coordstobed": coordsToBed,
}

def main(args):
    cmd = args[0]
    cmdargs = args[1:]

    if cmd in COMMANDS:
        C = COMMANDS[cmd]()
        C.parseArgs(cmdargs)
        C.run()
                    
if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
