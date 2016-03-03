import subprocess
import time
from optparse import OptionParser
import random, string, os

# important global vars
TWOBITPROGRAMLOCATION="/srv/gs1/software/ucsc_tools/2.7.2/bin/x86_64/twoBitToFa"
TWOBITGENOMELOCATION="/srv/gs1/projects/snyder/collinmelton/bundle/2.3/b37/d5/hs37d5.2bit"
GFCLIENTLOCATION="/srv/gs1/software/blat/3.5/bin/x86_64/gfClient"

#
# To run you first need to start a blat server as follows:
# /srv/gs1/software/blat/3.5/bin/x86_64/gfServer start localhost 8888 /srv/gs1/projects/snyder/collinmelton/bundle/2.3/b37/d5/hs37d5.2bit&
# Then run python CheckForMappingErrors.py --I input file --O outputFile --L size of region to left of mut --R size of region to right of mut --P blat server port number --M fileType, T for merged file, F for not, default is F
# since reads are usually 100ish in length I usually run with --L and --R = 50
# default input file format is as follows: 
# chrom    pos    ref    var
# 1    18866399    A    G
# 1    18866407    G    A
# 1    18866403    G    A
# 1    18866411    G    A
# 1    228771106    C    G

# generates a random word
def randomword(length):
    return ''.join(random.choice(string.lowercase) for i in range(length))
    
# this writes a sequence to a fasta file
def writeSeqToFile(seq, filename):
    f=open(filename, "w")
    f.write(">tempseq\n"+seq)
    f.close()

# this function reads and returns the contents of a file
def readPSL(filename):
    f=open(filename, "r")
    result=f.read()
    f.close()
    return result

# this class stores the result of a blat alignment and has a method to 
# check for a variant at a particular position
class BlatResult:
    def __init__(self, blatData):
        self.blatData=blatData
        self.score=int(self.blatData["match"].strip())
        self.hasVariant=False
        self.hasReference=False
        self.negStrand = (self.blatData["strand"]=="-")
        self.starts=self.blatData["qstarts"].strip(",").split(",")
        self.blocksizes=self.blatData["block_sizes"].strip(",").split(",")
        self.querySeqs=self.blatData["queryseqs"].strip(",").split(",")
        self.refSeqs=self.blatData["refseqs"].strip(",").split(",")
        
    def switchStrand(self, seq, switch=True):
        if not switch: return seq
        if seq.upper() == "C": return "G"
        if seq.upper() == "G": return "C"
        if seq.upper() == "T": return "A"
        if seq.upper() == "A": return "T"
        return seq
        
    def checkVariantPosition(self, (pos, ref, var)):
        # go through blocks until you find the position of interest and 
        # determine if it is reference or not
        for i in range(len(self.blocksizes)):
            if self.starts[i]=="":break
            start=int(self.starts[i])
            end=start+int(self.blocksizes[i])
            if pos>=start and pos<end:
                # if negative strand switches base to positive strand
                queryallele=self.switchStrand(self.querySeqs[i][pos-start], self.negStrand)
                refallele=self.switchStrand(self.refSeqs[i][pos-start], self.negStrand)
                if refallele.upper()==var.upper():
                    self.hasVariant=True
                if refallele.upper()==ref.upper():
                    self.hasReference=True
                break

# this function parses the blat output
def parseBlat(data):
    header=["match","mis-match","rep_match","Ns","Q_gap_count","q_gap_bases","T_gap_count","T_gap_bases","strand","Q_name","Q_size","Q_start","Q_end","T_name","T_size","T_start","T_end","block_count","block_sizes","qstarts","tstarts","queryseqs", "refseqs"]
    lines=data.split("\n")
    next=False
    results=[]
    for line in lines:
        vals=line.strip().split("\t")
        if next and len(vals)>5:
            aResult={}
            for i in range(len(vals)):
                aResult[header[i]]=vals[i]
            #print aResult
            results.append(BlatResult(aResult))
        if "-------------" in line:
            next=True
    return results

# the function calls and returns the output of the blat gfClient for a particular query sequence
def callBlat(seq, tempFile, port="8888"):
    writeSeqToFile(seq, tempFile)
    blatData=subprocess.check_output([GFCLIENTLOCATION, "localhost", port, "/", tempFile, "stdout", "-out=pslx"])
    return(parseBlat(blatData))

# given to tuples (start1, end1) and (start2, end2) determines the length of 
# overlapping sequence between the two
def getIntervalOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

# this function takes a blat result and finds the largest matching region of all possible
# 100 base pair regions within the site
def getMax100merScore(variant, maxVarResult):
    varInts=map(lambda x: (x,x+100), range(1,len(variant)-100))
    resultIntervals=[]
    for i in range(len(maxVarResult.starts)):
        resultIntervals.append([int(maxVarResult.starts[i]), int(maxVarResult.starts[i])+int(maxVarResult.blocksizes[i])])
    maxScore=0
    for varInt in varInts:
        maxScore=max(maxScore, sum(map(lambda x: getIntervalOverlap(varInt, x), resultIntervals)))
    return maxScore

# finds the difference between the variant and reference seq and uses blat to identify the 
# maximal scoring alignment that maps with the var allele as reference
def getMappingErrors(variant, reference, tempFile, port="8888"):
    # find position in variant sequence that is the variant
    for i in range(len(reference)):
        if reference[i]!=variant[i]:
            dif=(i, reference[i], variant[i])
    # get blat results for the variant sequence
    blatResults=callBlat(variant, tempFile, port=port)
    # iterate through blat results to check for the given result if the variant is now reference
    for result in blatResults:
        result.checkVariantPosition(dif)
    # for each blat result determine the max score for a result where variant is now reference 
    maxVarScore=0
    maxRefScore=0
    maxVarResult=None
    for result in blatResults:
        if result.hasVariant:
            if result.score>maxVarScore:
                maxVarScore=result.score
                maxVarResult=result
        if result.hasReference:
            maxRefScore=max(maxRefScore, result.score)
    max100merScore=0
    if maxVarResult !=None:
        max100merScore=getMax100merScore(variant, maxVarResult)
    return maxVarScore, maxRefScore, max100merScore

# reads and returs the sequence in a fasta file
def readFasta(contents):
    lines=contents.split("\n")
    result="".join(map(lambda x:x.strip(), lines[1:]))
    return result

# calls the ucsc twoBitToFa tool to get genomic sequence
def getSequence(chrom, start, stop):
    query=TWOBITGENOMELOCATION+":"+chrom+":"+str(int(start)-1)+"-"+str(int(stop)-1)
    output=subprocess.check_output([TWOBITPROGRAMLOCATION,query, "stdout"])
    return readFasta(output)
    
# this function gets the variant and normal sequence
def getVariantAndNormal(v, left, right):
    chrom=v.chrom
    pos=v.pos
    ref=v.ref
    var=v.var
    leftSeq=getSequence(chrom, pos-left, pos)
    rightSeq=getSequence(chrom, pos+1, pos+1+right)
    refSeq=getSequence(chrom, pos, pos+1)
    v.blatRef=refSeq
    if ref=="":
        ref=refSeq
        v.ref=ref
    if var=="":
        var =random.choice(list(set(["A", "T", "C", "G"])-set([ref])))
        v.var=var
    return (leftSeq+var+rightSeq, leftSeq+ref+rightSeq)

# this class represents a variant and has fields to store the alternative/mismap alignment score
class variant:
    def __init__(self, chrom, pos, ref, var):
        self.chrom=chrom.strip("chr")
        self.pos=int(pos)
        self.ref=ref.upper()
        self.var=var.upper()
        self.blatRef=""
        self.maxVarScore=0
        self.maxRefScore=0
        self.max100merScore=0
        
    def toString(self):
        return("\t".join([self.chrom, str(self.pos), self.ref, self.blatRef, self.var, str(self.maxVarScore),str(self.maxRefScore), str(self.max100merScore)]))

# reads an input file
def readVariantFileFromPipeline(filename, mergeFile=True):
    f=open(filename, 'r')
    line=f.readline().strip()
    while line[0]=="#":
        line=f.readline().strip()
    line=f.readline().strip()
    header=line.split("\t")
    print header
    result=[]
    while line !="":
        vals=line.strip().split("\t")
        if mergeFile:
            if len(vals)>4:
                result.append(variant(vals[0], vals[1], vals[3], vals[4]))
            if len(vals)==2:
                chrom=vals[0].strip("chr")
                pos=int(vals[1])
                refSeq=getSequence(chrom, pos, pos+1)
                for var in (set(["A", "C", "T", "G"])-set([refSeq])):
                    result.append(variant(vals[0], vals[1], refSeq, var))
        else:
            if len(vals)>3:
                result.append(variant(vals[0], vals[1], vals[2], vals[3]))
            if len(vals)==2:
                chrom=vals[0].strip("chr")
                pos=int(vals[1])
                refSeq=getSequence(chrom, pos, pos+1)
                for var in (set(["A", "C", "T", "G"])-set([refSeq])):
                    result.append(variant(vals[0], vals[1], refSeq, var))
        line=f.readline()
    f.close()
    return result

# writes an output file
def writeOutput(variants, outputFile):
    f=open(outputFile, "w")
    f.write("\t".join(["Chrom", "Pos", "Ref", "BlatRef", "Var", "VarScore", "RefScore", "Max100merAlignment_CountingMismatches"]))
    for v in variants:
        f.write("\n"+v.toString())
    f.close()

# run mismap analysis on an input file, left and right are the amount of sequence to consider on either side of the variant
def run(inputFile, left, right, outputFile, tempFile, fileType, port="8888", verbose=True):
    if verbose: print "starting"
    variants=readVariantFileFromPipeline(inputFile, mergeFile=(fileType=="T")) 
    if verbose: print "read in file"
    if verbose: print "# of ", len(variants)
    i=0
    for v in variants:
        # get sequence reference sequence surrounding variant as well as 
        # reference sequence with variant instead of reference
        varSeq, refSeq = getVariantAndNormal(v, left, right)
        # get maximum score of variant sequence mapping elsewhere in genome where
        # called variant is actually the reference 
        (v.maxVarScore, v.maxRefScore, v.max100merScore)=getMappingErrors(varSeq, refSeq, tempFile, port=port)
        i+=1
        if verbose and i%100==0: print len(variants)-i
    writeOutput(variants, outputFile)

# function to get command line options 
def getOptions():
    parser = OptionParser()
    parser.add_option("--I", dest = "inputFile", help = "input file",
                      metavar = "FILE", type = "string", default = "/Users/cmelton/Downloads/3c975b68-4a6e-4c97-a8d1-b9a427d4957c.unannotated.varscan.mutect.merge")
    parser.add_option("--O", dest="outputFile", help = "output file is a tab delimited file of key and count",
                      metavar = "FILE", type = "string", default = "test.out")
    parser.add_option("--L", dest = "left", help = "size of region to left of mut"+
                      "concurrently", metavar = "INTEGER", default = 100, type = "int")
    parser.add_option("--R", dest = "right", help = "size of region to right of mut"+
                      "concurrently", metavar = "INTEGER", default = 100, type = "int")
    parser.add_option("--M", dest="fileType", help="T for merged file, F for not", default="F", type="string")
    parser.add_option("--P", dest="port", help="blat server port number", default="8888", type="string")
    parser.add_option("--V", dest="verbose", help="whether to print status of program as it is running, T/F", default="T", type="string")
    (options, args) = parser.parse_args()
    return options

# this is the main function
def runMain():
    # get options and defaults
    options = getOptions()
    print options.inputFile, options.outputFile, options.left+1-1, options.right
    run(options.inputFile, options.left, options.right, options.outputFile, tempFile, options.fileType, port=options.port, verbose=(options.verbos=="T"))
    
# generates a temp filename for fasta files needed for the blat search
tempFile=randomword(10)+".fa"
runMain()
# remove temp file if it exists
if os.path.isfile(tempFile):
    os.remove(tempFile)
