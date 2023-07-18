#!/usr/bin/env python

import sys
import csv

WORK_DIR = "_gmwork"

### Utils

def vval(name, bindings):
    return bindings[name] if name in bindings else name

def parse_arguments(words):
    args = {}
    while len(words) > 1:
        a = words.pop(0)
        b = words.pop(0)
        if a in args:
            if type(args[a]) == list:
                args[a].append(b)
            else:
                args[a] = [args[a], b]
        else:
            args[a] = b
    return args

def parse_number_list(s):
    """Parse a string of the form 1,2,3 into a list of numbers."""

    return [ int(x.strip()) for x in s.split(",") ]

### Commands 

class Command():
    name = ""
    file1 = ""
    destination = ""
    arguments = {}
    assignments = []
    work_dir = "_gmwork"

    def __init__(self, filenum):
        self.destination = self.work_dir + "/file-" + str(filenum)

    def getarg(self, name, default=None):
        return self.arguments[name] if name in self.arguments else default

    def write_prereq(self, out, *prereqs):
        out.write("{} : {}".format(self.destination, " ".join(prereqs)))
        for assg in self.assignments:
            out.write(" " + assg[1])
        out.write("\n")

    def write_assignments(self, out):
        for assg in self.assignments:
            out.write("\t{}=$$(head -1 {})\n".format(assg[0], assg[1]))

class Cmd_file(Command):
    name = "file"

    def example(self):
        return "a = file filename.txt"

    def parse(self, file1, bindings):
        self.destination = file1

    def write(self, out):
        pass

class Cmd_filter(Command):
    name = "filter"
    test = ""
    hdr = ""

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)
        self.test = self.getarg("test")
        self.hdr = self.getarg("header")

    def write(self, out):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        if self.hdr == "keep":
            h = "-P"
        elif self.hdr == "remove":
            h = "-p"
        else:
            h = ""
        out.write("\ttcalc.py {} \"filter {}\" {} > {}\n\n".format(h, self.test, self.file1, self.destination))

class Cmd_assoc(Command):
    name = "assoc"
    file2 = ""
    incol = "1"
    outcol = None
    missing = ""

    def example(self):
        return "a = assoc file.txt with file2.txt in 3 out 4 missing NA"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.file2 = vval(self.getarg("with"), bindings)
        self.destination = self.getarg("as", self.destination)
        self.incol = self.getarg("in") or "1"
        self.outcol = self.getarg("out")
        self.missing = self.getarg("missing")

    def write(self, out):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        oc = "-o " + self.outcol if self.outcol else "-w"
        mi = "-m " + self.missing if self.missing else "-x"
        out.write("\tcut -f 1 {} | assoc.py -i {} {} {} {} > {}\n\n".format(self.file1, self.incol, oc, mi, self.file2, self.destination))

class Cmd_intersect(Command):
    name = "intersect"
    file2 = ""

    def example(self):
        return "a = intersect file1.txt with file2.txt as file3.txt"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.file2 = vval(self.getarg("with"), bindings)
        self.destination = self.getarg("as", self.destination)

    def write(self, out):
        self.write_prereq(out, self.file1, self.file2)
        self.write_assignments(out)
        out.write("\tbedtools intersect -a {} -b {} | sort -k 1,1 -k 2,2n | bedtools merge -i - > {}\n\n".format(self.file1, self.file2, self.destination))

class Cmd_subtract(Command):
    name = "subtract"
    file2 = ""

    def example(self):
        return "a = subtract file2.txt from file1.txt"

    def parse(self, file1, bindings):
        self.file2 = vval(file1, bindings)
        self.file1 = vval(self.getarg("from"), bindings)
        self.destination = self.getarg("as", self.destination)

    def write(self, out):
        self.write_prereq(out, self.file1, self.file2)
        self.write_assignments(out)
        out.write("\tbedtools intersect -v -a {} -b {} | sort -k 1,1 -k 2,2n | bedtools merge -i - > {}\n\n".format(self.file1, self.file2, self.destination))

class Cmd_merge(Command):
    name = "merge"
    file2 = ""

    def example(self):
        return "a = merge file1.txt with file2.txt"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.file2 = vval(self.getarg("with"), bindings)
        self.destination = self.getarg("as", self.destination)

    def write(self, out):
        self.write_prereq(out, self.file1, self.file2)
        self.write_assignments(out)
        out.write("\tcat {} {} | sort -k 1,1 -k 2,2n | bedtools merge -i - > {}\n\n".format(self.file1, self.file2, self.destination))

class Cmd_cut(Command):
    name = "cut"
    columns = [1]

    def example(self):
        return "a = cut file1.txt columns 1,3,2,5"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)
        self.columns = self.getarg("cols") or self.getarg("columns")

    def write(self, out):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        out.write("\tkut -f {} {} > {}\n\n".format(self.columns, self.file1, self.destination))

class Cmd_paste(Command):
    name = "paste"
    otherfiles = []
    column = 1

    def example(self):
        return "a = paste file1.txt with file2.txt,file3.txt,file4.txt"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)
        self.columns = self.getarg("cols") or self.getarg("columns") or 1
        if "with" in self.arguments:
            self.otherfiles = [ vval(f, bindings) for f in self.arguments["with"].split(",") ]

    def write(self, out):
        self.write_prereq(out, self.file1, *self.otherfiles)
        self.write_assignments(out)
        out.write("\tpazte -o {} -F -c {} {} {}\n\n".format(self.destination, self.column, self.file1, " ".join(self.otherfiles)))

class Cmd_header(Command):
    name = "header"
    remove = 0
    colnames = ""

    def example(self):
        return "a = header b colnames name1,name2,name3"
    
    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)
        self.remove = self.getarg("remove", self.remove)
        if "colnames" in self.arguments:
            self.colnames = "\\t".join(self.arguments["colnames"].split(","))

    def write(self, out):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        if self.colnames:
            out.write('\techo -e "{}" > {}\n'.format(self.colnames, self.destination))
        if self.remove:
            out.write("\ttail -n +{} {} >> {}\n\n".format(self.remove, self.file1, self.destination))
        else:
            out.write("\tcat {} >> {}\n\n".format(self.file1, self.destination))

class Cmd_bedtocoords(Command):
    name = "bedtocoords"
    column = 1

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)
        self.columns = self.getarg("cols") or self.getarg("columns")
        
    def write(self, out):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        out.write("\tgmutils.py bedtocoords -c {} {} > {}\n\n".format(self.column, self.file1, self.destination))

class Cmd_subintervals(Command):
    name = "subintervals"
    size = 10000

    def example(self):
        return "a = subintervals file1.txt size 50000 as file2.txt"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)
        self.size = self.getarg("size") or self.size

    def write(self, out):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        out.write("\tgmutils.py subintervals -s {} {} > {}\n\n".format(self.size, self.file1, self.destination))

class Cmd_quantify(Command):
    name = "quantify"
    bedfile = ""
    norm = False

    def example(self):
        return "a = quantify file1.bam on file2.bed as file3.txt norm 1234567"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)
        self.norm = self.getarg("norm")
        self.bedfile = vval(self.getarg("on"), bindings)

    def write(self, out):
        self.write_prereq(out, self.file1, self.bedfile)
        self.write_assignments(out)
        if self.norm:
            out.write("\tsamtools bedcov {} {} | gmutils.py scale {} > {}\n\n".format(self.bedfile, self.file1, self.norm, self.destination))
        else:
            out.write("\tsamtools bedcov {} {} > {}\n\n".format(self.bedfile, self.file1, self.destination))

class Cmd_adjust(Command):
    name = "adjust"
    spec = ""

    def example(self):
        return "a = adjust file1.bed left 100 right 200"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)
        if "left" in self.arguments:
            self.spec += " -l " + self.arguments["left"]
        if "right" in self.arguments:
            self.spec += " -r " + self.arguments["right"]
        if "upstream" in self.arguments:
            self.spec += " -u " + self.arguments["upstream"]
        if "downstream" in self.arguments:
            self.spec += " -d " + self.arguments["downstream"]
        if "column" in self.arguments:
            self.spec += " -c " + self.arguments["column"]
        if "col" in self.arguments:
            self.spec += " -c " + self.arguments["col"]
        if "strandcol" in self.arguments:
            self.spec += " -s " + self.arguments["strandcol"]
        if "mode" in self.arguments:
            self.spec += " -m " + self.arguments["mode"]
        
    def write(self, out):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        out.write("\tgmutils.py adjust {} {} > {}\n\n".format(self.spec, self.file1, self.destination))

# Commands that extract values from files

class Cmd_nreads(Command):
    name = "nreads"
    
    def example(self):
        return "a = nreads file1.bam"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)

    def write(self, out):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        out.write("\tsamtools idxstats {} | tcalc.py -q \"do sum(C3)\" > {}\n\n".format(self.file1, self.destination))

class Cmd_nlines(Command):
    name = "nlines"
    
    def example(self):
        return "a = nlines file1.txt"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.destination = self.getarg("as", self.destination)

    def write(self, out):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        out.write("\tgrep -c ^ {} > {}\n\n".format(self.file1, self.destination))

## Operations on sets of objects

class Cmd_setmerge(Command):
    name = "setmerge"
    file2 = ""
    col1 = 1
    col2 = 1

    def example(self):
        return "a = setmerge file1 col1 3 with file2 col2 5"

    def parse(self, file1, bindings):
        self.file1 = vval(file1, bindings)
        self.file2 = vval(self.getarg("with"), bindings)
        self.destination = self.getarg("as", self.destination)
        self.col1 = self.getarg("col1", default="1")
        self.col2 = self.getarg("col2", default="1")

    def write(self, out):
        self.write_generic(out, "-u")

    def write_generic(self, out, op):
        self.write_prereq(out, self.file1)
        self.write_assignments(out)
        out.write("\tcolx.py -w {} {}:{} {}:{} > {}\n\n".format(op, self.file1, self.col1, self.file2, self.col2, self.destination))

class Cmd_setdifference(Cmd_setmerge):
    name = "setdifference"

    def example(self):
        return "a = setdifference file1 col1 3 with file2 col2 5"

    def write(self, out):
        self.write_generic(out, "-d")

class Cmd_setintersect(Cmd_setmerge):
    name = "setintersect"

    def example(self):
        return "a = setintersect file1 col1 3 with file2 col2 5"

    def write(self, out):
        self.write_generic(out, "-i")

COMMANDS = {"file": Cmd_file,
            "filter": Cmd_filter,
            "assoc" : Cmd_assoc,
            "intersect": Cmd_intersect,
            "subtract": Cmd_subtract,
            "merge": Cmd_merge,
            "cut": Cmd_cut,
            "paste": Cmd_paste,
            "header": Cmd_header,
            "subintervals": Cmd_subintervals,
            "quantify": Cmd_quantify,
            "bedtocoords": Cmd_bedtocoords,
            "adjust": Cmd_adjust,
            "nreads": Cmd_nreads,
            "nlines": Cmd_nlines,
            "setmerge": Cmd_setmerge,
            "setdifference": Cmd_setdifference,
            "setintersect": Cmd_setintersect,
}

### Main application class

class GenoMaker():
    scriptfile = None
    variables = {}
    commands = []
    targets = []
    filenum = 1

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-w":
                Command.work_dir = a
                prev = ""
            elif a in ["-w"]:
                prev = a
            else:
                self.scriptfile = a
        return True

    def parseScript(self, filename):
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line[0] == '#':
                    continue
                self.parseCommand(line)

    def parseAssignments(self, line):
        asslist = []
        p1 = line.find("{")
        p2 = line.find("}")
        if (p1 > 0) and (p2 > 0) and (p2 > p1):
            clean_line = line[:p1] + line[p2+1:]
            assignments = line[p1+1:p2].split(";")
            for assg in assignments:
                if "=" in assg:
                    parts = assg.split("=")
                    var = parts[0].strip()
                    source = vval(parts[1].strip(), self.variables)
                    asslist.append((var, source))
            return True, clean_line, asslist
        elif (p1 > 0) or (p2 > 0):
            return False, "", []
        else:
            return True, line, []

    def parseCommand(self, line):
        assignments = []
        good = True
        if "{" in line:
            good, line, assignments = self.parseAssignments(line)

        if good:
            words = line.split()
            if len(words) > 2 and words[0] == 'target':
                self.targets.append(words[1:])
            elif len(words) > 3 and words[1] == '=':
                cmdname = words[2]
                if cmdname in COMMANDS:
                    cmdclass = COMMANDS[cmdname]
                    cmd = cmdclass(self.filenum)
                    cmd.assignments = assignments
                    cmd.arguments = parse_arguments(words[4:])
                    cmd.parse(words[3], self.variables)
                    self.commands.append(cmd)
                    self.variables[words[0]] = cmd.destination
                    self.filenum += 1
                    #print(self.variables)
                else:
                    sys.stderr.write("Unknown command `{}'\n".format(cmdname))
        else:
            sys.stderr.write("Cannot parse line `{}'\n".format(line))

    def run(self):
        self.parseScript(self.scriptfile)
        sys.stdout.write("#!/usr/bin/env make -f\n\n")
        sys.stdout.write(".ONESHELL : \n\n")
        for cmd in self.commands:
            cmd.write(sys.stdout)
        for tg in self.targets:
            tgs = [vval(t, self.variables) for t in tg[1:]]
            sys.stdout.write("{} : {}\n\n".format(tg[0], " ".join(tgs)))

def main(args):
    GM = GenoMaker()
    if GM.parseArgs(args):
        GM.run()

if __name__ == "__main__":
  args = sys.argv[1:]
  main(args)


###

# a = file abc.bed
# b = file def.bed
# c = intersect b with c
# d = subtract c from b
# e = merge c with d
# x = filter a test "C3>0"
# quantify
# quartiles
# $p = 
