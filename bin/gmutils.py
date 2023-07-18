#!/usr/bin/env python

import sys
import csv
import numpy as np

class Adjust():
    infile = ""
    column = 0
    strandcol = None
    mode = "b"                  # or l, r
    left = 0
    right = 0
    upstream = 0
    downstream = 0
    _func = None

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-l":
                self.left = int(a)
                prev = ""
            elif prev == "-r":
                self.right = int(a)
                prev = ""
            elif prev == "-u":
                self.upstream = int(a)
                prev = ""
            elif prev == "-d":
                self.downstream = int(a)
                prev = ""
            elif prev == "-c":
                self.column = int(a)-1
                prev = ""
            elif prev == "-s":
                self.strandcol = int(a)-1
                prev = ""
            elif prev == "-m":
                self.mode = a
                prev = ""
            elif a in ["-l", "-r", "-u", "-d", "-c", "-s", "-m"]:
                prev = a
            else:
                self.infile = a

        if self.strandcol:
            if self.mode == 'l':
                self._func = self.sl
            elif self.mode == 'r':
                self._func = self.sr
            else:
                self._func = self.sb
        else:
            if self.mode == 'l':
                self._func = self.ul
            elif self.mode == 'r':
                self._func = self.ur
            else:
                self._func = self.ub

    def ub(self, a, b, row):
        return a + self.left, b + self.right

    def ul(self, a, b, row):
        return a + self.left, a + self.right

    def ur(self, a, b, row):
        return b + self.left, b + self.right

    def sb(self, a, b, row):
        strand = row[self.strandcol]
        if strand == "+":
            return a + self.upstream, b + self.downstream
        else:
            return a + self.downstream, b + self.upstream

    def sl(self, a, b, row):
        strand = row[self.strandcol]
        if strand == "+":
            return a + self.upstream, a + self.downstream
        else:
            return b + self.downstream, b + self.upstream

    def sr(self, a, b, row):
        strand = row[self.strandcol]
        if strand == "+":
            return b + self.upstream, b + self.downstream
        else:
            return a + self.downstream, a + self.upstream

    def run(self):
        with open(self.infile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                start = int(row[self.column+1])
                end = int(row[self.column+2])
                newstart, newend = self._func(start, end, row)
                newstart = max(newstart, 0)
                row[self.column+1] = str(newstart)
                row[self.column+2] = str(newend)
                sys.stdout.write("\t".join(row) + "\n")
        

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
                        a = b+1
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

class Correlation():
    infile = ""
    col1 = 0
    col2 = 1
    
    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-x":
                self.col1 = int(a) - 1
                prev = ""
            elif prev == "-y":
                self.col2 = int(a) - 1
                prev = ""
            elif a in ["-x", "-y"]:
                prev = a
            else:
                self.infile = a

    def run(self):
        xs = []
        ys = []
        with open(self.infile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                xs.append(float(row[self.col1]))
                ys.append(float(row[self.col2]))
        r = np.corrcoef(xs, ys)
        n = len(xs)
        sys.stdout.write("{}:{}:{}\t{}\t{}\t{}\t{}\n".format(self.infile, self.col1+1, self.col2+1, n, sum(xs)/n, sum(ys)/n, r[0,1]))

class Scale():
    column = 3
    factor = 10000000
    value = 1

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-c":
                self.column = int(a) - 1
                prev = ""
            elif prev == "-f":
                self.factor = int(a)
                prev = ""
            elif a in ["-c", "-f"]:
                prev = a
            else:
                self.value = int(a)

    def run(self):
        sf = self.factor / self.value
        c = csv.reader(sys.stdin, delimiter='\t')
        for row in c:
            row[self.column] = str(float(row[self.column]) * sf)
            sys.stdout.write("\t".join(row) + "\n")

COMMANDS = {"subintervals": Subintervals,
            "bedtocoords": bedToCoords,
            "coordstobed": coordsToBed,
            "correl": Correlation,
            "adjust": Adjust,
            "scale": Scale
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
