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
    work_dir = "_gmwork"

    def __init__(self, filenum):
        self.destination = self.work_dir + "/file-" + str(filenum)

class Cmd_file(Command):
    name = "file"

    def example(self):
        return "a = file filename.txt"

    def parse(self, file1, args, bindings):
        self.destination = file1

    def write(self, out):
        pass

class Cmd_filter(Command):
    name = "filter"
    column = 4
    test = ""

class Cmd_intersect(Command):
    name = "intersect"
    file2 = ""

    def example(self):
        return "a = intersect file1.txt with file2.txt as file3.txt"

    def parse(self, file1, args, bindings):
        self.file1 = vval(file1, bindings)
        if "with" in args:
            self.file2 = vval(args["with"], bindings)
        if "as" in args:
            self.destination = args["as"]

    def write(self, out):
        out.write("{} : {} {}\n".format(self.destination, self.file1, self.file2))
        out.write("\tbedtools intersect -a {} -b {} | sort -k 1,1 -k 2,2n | bedtools merge -i - > {}\n\n".format(self.file1, self.file2, self.destination))

class Cmd_subtract(Command):
    name = "subtract"
    file2 = ""

    def example(self):
        return "a = subtract file2.txt from file1.txt"

    def parse(self, file1, words, bindings):
        self.file2 = vval(file1, bindings)
        self.file1 = vval(words["from"], bindings)
        if "as" in words:
            self.destination = words["as"]

    def write(self, out):
        out.write("{} : {} {}\n".format(self.destination, self.file1, self.file2))
        out.write("\tbedtools intersect -v -a {} -b {} | sort -k 1,1 -k 2,2n | bedtools merge -i - > {}\n\n".format(self.file1, self.file2, self.destination))

class Cmd_merge(Command):
    name = "merge"
    file2 = ""

    def example(self):
        return "a = merge file1.txt with file2.txt"

    def parse(self, file1, words, bindings):
        self.file1 = vval(file1, bindings)
        self.file2 = vval(words["with"], bindings)
        if "as" in words:
            self.destination = words["as"]

    def write(self, out):
        out.write("{} : {} {}\n".format(self.destination, self.file1, self.file2))
        out.write("\tcat {} {} | sort -k 1,1 -k 2,2n | bedtools merge -i - > {}\n\n".format(self.file1, self.file2, self.destination))

class Cmd_cut(Command):
    name = "cut"
    columns = [1]

    def example(self):
        return "a = cut file1.txt columns 1,3,2,5"

    def parse(self, file1, args, bindings):
        self.file1 = vval(file1, bindings)
        if "as" in args:
            self.destination = args["as"]
        for k in ["cols", "columns"]:
            if k in args:
                self.columns = args[k]
                break

    def write(self, out):
        out.write("{} : {}\n".format(self.destination, self.file1))
        out.write("\tkut -f {} {} > {}\n\n".format(self.columns, self.file1, self.destination))

class Cmd_paste(Command):
    name = "paste"
    otherfiles = []

    def example(self):
        return "a = paste file1.txt with file2.txt,file3.txt,file4.txt"

    def parse(self, file1, args, bindings):
        self.file1 = vval(file1, bindings)
        if "as" in args:
            self.destination = args["as"]
        if "with" in args:
            self.otherfiles = [ vval(f, bindings) for f in args["with"].split(",") ]

    def write(self, out):
        out.write("{} : {} {}\n".format(self.destination, self.file1, " ".join(self.otherfiles)))
        out.write("\tpazte -o {} -F {} {}\n\n".format(self.destination, self.file1, " ".join(self.otherfiles)))

class Cmd_bedtocoords(Command):
    name = "bedtocoords"
    column = 1
    def parse(self, file1, args, bindings):
        self.file1 = vval(file1, bindings)
        if "as" in args:
            self.destination = args["as"]
        for k in ["col", "column"]:
            if k in args:
                self.column = args[k]
                break
        
    def write(self, out):
        out.write("{} : {}\n".format(self.destination, self.file1))
        out.write("\tgmutils.py bedtocoords -c {} {} > {}\n\n".format(self.column, self.file1, self.destination))

class Cmd_subintervals(Command):
    name = "subintervals"
    size = 10000

    def example(self):
        return "a = subintervals file1.txt size 50000 as file2.txt"

    def parse(self, file1, args, bindings):
        self.file1 = vval(file1, bindings)
        if "as" in args:
            self.destination = args["as"]
        if "size" in args:
            self.size = int(args["size"])

    def write(self, out):
        out.write("{} : {}\n".format(self.destination, self.file1))
        out.write("\tgmutils.py subintervals -s {} {} > {}\n\n".format(self.size, self.file1, self.destination))

class Cmd_quantify(Command):
    name = "quantify"
    bedfile = ""

    def example(self):
        return "a = quantify file1.bam on file2.bed as file3.txt"

    def parse(self, file1, args, bindings):
        self.file1 = vval(file1, bindings)
        if "as" in args:
            self.destination = args["as"]
        self.bedfile = vval(args["on"], bindings)

    def write(self, out):
        out.write("{} : {} {}\n".format(self.destination, self.file1, self.bedfile))
        out.write("\tsamtools bedcov {} {} > {}\n\n".format(self.bedfile, self.file1, self.destination))

COMMANDS = {"file": Cmd_file,
            "intersect": Cmd_intersect,
            "subtract": Cmd_subtract,
            "merge": Cmd_merge,
            "cut": Cmd_cut,
            "paste": Cmd_paste,
            "subintervals": Cmd_subintervals,
            "quantify": Cmd_quantify,
            "bedtocoords": Cmd_bedtocoords,
}

### Main application class

class GenoMaker():
    scriptfile = None
    variables = {}
    commands = []
    filenum = 1

    def parseArgs(args):
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

    def parseCommand(self, line):
        words = line.split()
        if len(words) > 3 and words[1] == '=':
            cmdname = words[2]
            if cmdname in COMMANDS:
                cmdclass = COMMANDS[cmdname]
                cmd = cmdclass(self.filenum)
                args = parse_arguments(words[4:])
                cmd.parse(words[3], args, self.variables)
                self.commands.append(cmd)
                self.variables[words[0]] = cmd.destination
                self.filenum += 1
                #print(self.variables)
            else:
                sys.stderr.write("Unknown command `{}'\n".format(cmdname))

    def run(self):
        self.parseScript(self.scriptfile)
        sys.stdout.write("#!/usr/bin/env make -f\n\n")
        sys.stdout.write(".ONESHELL : \n\n")
        for cmd in self.commands:
            cmd.write(sys.stdout)

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
