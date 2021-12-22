import sys
import os
from math import log
import subprocess as sp
import random
import time
from statistics import median
from random import randint
import argparse
import signal
from functools import partial

def receiveSignal(tempFiles, signalNumber, frame):
    return
    print('Received signal:', signalNumber)
    print('Cleaning tmp files')
    for f in tempFiles:
        if os.path.exists(f):
            print("removing", f, "...", end="")
            os.remove(f)
            print("removed")
    sys.exit()


#parse .gcnf instance, 
#returns a pair C,B where B contains the base (hard) clauses and C the other clauses
def parse(filename):
    C = []
    B = []
    with open(filename, "r") as f:
        lines = f.readlines()
        if filename[-5:] == ".gcnf":
            for line in lines[1:]:
                if line[0] in ["p","c"]: continue
                line = line.split(" ")
                cl = [int(i) for i in line[1:-1]]
                if len(cl) > 0:
                    if line[0] == "{0}":
                        B.append(cl)
                    else:
                        C.append(cl)
        else:
            for line in lines[1:]:
                if line[0] in ["p","c"]: continue
                line = line.split(" ")
                cl = [int(i) for i in line[:-1]]
                if len(cl) > 0:
                    C.append(cl)
    return C,B

def exportGCNF(soft, hard, filename):
    print("running export for ", filename)
    with open(filename, "w") as f:
        maxVar = max([max(abs(l) for l in cl) for cl in soft + hard])
        f.write("p gcnf {} {} {}\n".format(maxVar, len(soft + hard), len(soft)))
        for cl in hard:
            f.write("{0} " + " ".join([str(l) for l in cl]) + " 0\n")
        clid = 1
        for cl in soft:
            f.write("{" + str(clid)  + "} " + " ".join([str(l) for l in cl]) + " 0\n")
            clid += 1

#returns random Boolean value
def randomBool():
    return bool(random.getrandbits(1))

def run(cmd, timeout, ttl = 3):
    proc = sp.Popen([cmd], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    try:
        (out, err) = proc.communicate(timeout = int(timeout * 1.1) + 1)
        out = out.decode("utf-8")
    except sp.TimeoutExpired:
        proc.kill()
        try:
            (out, err) = proc.communicate()
            out = out.decode("utf-8")
        except ValueError:
            if ttl > 0:
                return run(cmd, timeout, ttl - 1)
            out = ""
    return out

class Counter:
    def __init__(self, filename, e, d):
        self.rid = randint(1,10000000)
        self.originalFilename = filename
        self.filename = filename
        self.C, self.B = parse(filename)
        self.trimFilename = filename
        self.autarkyTrim()
        self.MUSes = []
        self.counts = []
        self.dimension = len(self.C)
        self.XOR = None
        self.tresh = 1 + 9.84 * (1 + (e / (1 + e)))*(1 + 1/e)*(1 + 1/e)
        self.t = int(17 * log(3 / d,2));
        self.checks = 0
        self.unexXorFilename = "./tmp/unex_{}.cnf".format(self.rid)
        self.tmpFiles = [self.unexXorFilename]
        if self.trimFilename != self.filename:
            self.tmpFiles.append(self.trimFilename)
        self.QBF = "3QBF"
        self.findMAlgo = "hybridSearch"

    def initialThresholdCheck(self):
        if ".gcnf" in self.originalFilename: return self.tresh
        cmd = "timeout 10 ./unimus " + self.filename
        out = run(cmd, 10)
        lastMUS = ""
        for line in out.splitlines():
            if "Found MUS" in line: lastMUS = line
        if lastMUS == "": return self.tresh #no MUS found, we continue with amusic
        count = int(lastMUS.split("#")[1].split(",")[0])
        return count

    def autarkyTrim(self):
        if ".gcnf" in self.filename: return
        cmd = "timeout 3600 python3 autarky.py {}".format(self.filename)
        #print(cmd)
        out = run(cmd, 3600)
        if "autarky vars" in out:
            for line in out.splitlines():
                line = line.rstrip()
                if line[:2] == "v ":
                    autarky = [int(c) - 1 for c in line.split()[1:]]
        else: return

        imu = self.getImu()
        #print(len(self.C), len(autarky), len(set(autarky)))
        #print(autarky)

        C = [self.C[c] for c in sorted(set(autarky))]
        B = []
        if len(imu) > 0:
            B = [self.C[c] for c in imu]
        print("original size: {}, autarky: {}, IMU: {}".format(len(self.C), len(C), len(B)))
        self.C, self.B = C, B
        self.trimFilename = "./tmp/input_" + str(self.rid) + ".gcnf"
        exportGCNF(self.C, self.B, self.trimFilename) 

    def getImu(self):
        cmd = "timeout 3600 python3 gimu.py {}".format(self.filename)
        #print(cmd)
        out = run(cmd, 3600)
        if "imu size" in out and not "imu size: 0" in out:
            for line in out.splitlines():
                line = line.rstrip()
                if line[:2] == "v ":
                    return [int(c) - 1 for c in line.split()[1:]]
        else: return []

    #generate (and internaly store) a hash function from H_{xor}(dimension, dimension - 1) 
    def generateXOR(self):
        m = self.dimension - 1
        self.XOR = [[] for _ in range(m)]
        for i in range(m):
            for j in range(self.dimension):
                if randomBool():
                    self.XOR[i].append(j + 1)
            if (len(self.XOR[i]) > 0) and randomBool(): #this corresponds to a0
                self.XOR[i][0] *= -1

    #export the xor function in .cnf (xcnf) format
    #only individual clauses are exported (we do not export header)
    def exportXor(self, m):
        res = ""
        for i in range(m):
            if len(self.XOR[i]) > 0:
                res += "x " +  " ".join([str(l) for l in self.XOR[i]]) + " 0\n"
        return res

    #return the set (C \ N) (indices of the set)
    def complement(self, N):
        return [i for i in range(1, len(self.C) + 1) if i not in N]

    #find an unexplored MUS in the cell. If there is no such a MUS in the cell, returns []
    ## TODO: avoid external calling of gqbf.py, just integrate it
    def getMUS(self, m):
        self.checks += 1
        with open(self.unexXorFilename, "w") as f:
            f.write("p cnf 0 0\n")
            for MUS in self.MUSes:
                f.write(" ".join([str(-l) for l in MUS]) + " 0\n")
                f.write(" ".join([str(l) for l in self.complement(MUS)]) + " 0\n")
            #aternatively, we block only the MUS itself
            #for MUS in self.MUSes:
            #    f.write(" ".join([str(-l) for l in MUS]) + " ")
            #    f.write(" ".join([str(l) for l in self.complement(MUS)]) + " 0\n")
            f.write(self.exportXor(m))
        assert self.QBF in ["2QBF", "3QBF", "2QBFBeta"]
        cmd = "python3 2gqbf.py {} {}".format(self.trimFilename, self.unexXorFilename)
        if self.QBF == "3QBF":
            cmd = "python3 gqbf.py {} {}".format(self.trimFilename, self.unexXorFilename)
        if self.QBF == "2QBFBeta":
            cmd = "python3 2gqbfBeta.py {} {}".format(self.trimFilename, self.unexXorFilename)
        #print(cmd)
        proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.decode("utf-8")
        if not "SOLUTION" in out:
            return []
        reading = False
        for line in out.splitlines():
            if reading:
                #print(cmd)
                MUS = [int(l) for l in line.split(" ") if int(l) > 0]
                #print(MUS)
                return MUS
            if "SOLUTION" in line:
                reading = True
    # returns True if MUS is in the cell and False otherwise
    # this procedure needs a debug
    def isInCell(self, MUS, m):
        for i in range(m):
            satisfy = len(set(self.XOR[i][1:]).intersection(set(MUS))) % 2 == 1
            if self.XOR[i][0] > 0 and self.XOR[i][0] in MUS:
                satisfy = not satisfy
            if self.XOR[i][0] < 0 and self.XOR[i][0] not in MUS:
                satisfy = not satisfy
            if not satisfy:
                return False
        return True

    #Counts (and returns) the number of MUSes in the cell given by the m-th prefix of h
    def bsatXor(self, m, exploredMUSes, limit = None):
        if limit == None: limit = self.tresh
        print("start of bsatXor, MUSes:", len(self.MUSes), "m:", m)
        self.MUSes = []
        #assert len(exploredMUSes[m]) == 0
        exploredMUSes[m] = []
        for i in range(self.dimension -1, m, -1):
            if len(exploredMUSes[i]) > 0:
                self.MUSes = exploredMUSes[i][:]
                exploredMUSes[m] = exploredMUSes[i][:]
        print("---initial size", len(self.MUSes))
        while len(self.MUSes) < limit:
            MUS = self.getMUS(m)
            if len(MUS) == 0:
                return len(self.MUSes)
            self.MUSes.append(MUS)
            exploredMUSes[m].append(MUS)
        return len(self.MUSes)
  

    def findMPlainBinarySearch(self, mPrev):
        exploredMUSes = [[] for _ in range(self.dimension)]
        print("start of plain binary search")
        low = 0
        high = self.dimension - 1
        finalCount = -1
        finalM = -1
        m = int((low + high) / 2)
        while high - low > 1:
            count = self.bsatXor(m, exploredMUSes)
            print("m: {}, {}".format(m, count))
            if count >= self.tresh:
                low = m
            else:
                high = m
                if count >  finalCount:
                    finalCount = count
                    finalM = m
            m = int((low + high) / 2)
        if finalM < 0:
            return 0,0
        return finalCount * pow(2,finalM), finalM 

    def hybridSearch(self, mPrev):        
        exploredMUSes = [[] for _ in range(self.dimension)]
        print("start of hybrid search")
        #first identify the first cell with at least one MUS
        prevSafe = self.dimension - 1
        m = prevSafe
        count = self.bsatXor(m, exploredMUSes, 1)
        step = 1
        while True:
            while count == 0:
                prevSafe = m
                m -= step
                m = max(1, m)
                count = self.bsatXor(m, exploredMUSes, 1)
                step *= 2
            if prevSafe == m + 1: break        
            m = prevSafe
            step = 1
            count = 0
        
        finalCount = self.bsatXor(m, exploredMUSes)
        assert finalCount > 0
        if finalCount >= self.tresh:
            return 0,0
        while m > 0:
            m -= 1
            count = self.bsatXor(m, exploredMUSes)
            print("count", count, m)
            if count >= self.tresh: break
            finalCount = count
        return finalCount * pow(2,m + 1), m + 1 
                
    def mPrevLinearSearch(self, mPrev):       
        #if the last value of m was tight, i.e., the MUS count was close to the threshold, try first a bit larger cell to avoid 
        #exceeding the threshold
        if len(self.counts) > 0 and self.counts[-1] > 0.75*self.tresh:
            mPrev = min(mPrev + 1, self.dimension - 1) # try a bit larger cell first 
        exploredMUSes = [[] for _ in range(self.dimension)]
        low = 0
        high = self.dimension - 1
        finalCount = -1
        finalM = -1
        count = self.bsatXor(mPrev, exploredMUSes)
        if count >= self.tresh:
            low = mPrev
        else:
            high = mPrev
            finalCount = count
            finalM = mPrev
            print("first count: ", count, " with m:", mPrev)
            count = self.bsatXor(mPrev - 1, exploredMUSes)
            print("second count: ", count, " with m:", mPrev - 1)
            if count >= self.tresh:
                return finalCount * pow(2,finalM), finalM
            else:
                high = mPrev - 1
                finalCount = count
                finalM = mPrev - 1

        m = high
        if finalCount < 1: finalCouddnt = self.bsatXor(m, exploredMUSes)
        if finalCount >= self.tresh:
            return 0,0
        while m > 0:
            m -= 1
            count = self.bsatXor(m, exploredMUSes)
            print("count", count, m)
            if count >= self.tresh: break
            finalCount = count
        return finalCount * pow(2,m + 1), m + 1 

    def logSatSearch(self, mPrev):        
            low = 0
            high = self.dimension - 1
            finalCount = -1
            finalM = -1
            count = self.bsatXor(mPrev)
            if count >= self.tresh:
                low = mPrev
            else:
                high = mPrev
                finalCount = count
                finalM = mPrev
                print("first count: ", count, " with m:", mPrev)
                count = self.bsatXor(mPrev - 1)
                print("second count: ", count, " with m:", mPrev - 1)
                if count >= self.tresh:
                    return finalCount * pow(2,finalM), finalM
                else:
                    high = mPrev - 1
                    finalCount = count
                    finalM = mPrev - 1

            m = int((low + high) / 2)
            while high - low > 1:
                count = self.bsatXor(m)
                print("m: {}, {}".format(m, count))
                if count >= self.tresh:
                    low = m
                else:
                    high = m
                    if count >  finalCount:
                        finalCount = count
                        finalM = m
                m = int((low + high) / 2)
            if finalM < 0:
                return 0,0
            return finalCount * pow(2,finalM), finalM 

    def approxMC(self, mPrev):
        self.generateXOR()
        if mPrev == 0:
            if self.findMAlgo == "plainBinary":
                result = self.findMPlainBinarySearch(mPrev)
            else:
                result = self.hybridSearch(mPrev)
            print("FINDMCHECKS:", self.checks)
            return result
        else:
            return self.mPrevLinearSearch(mPrev)

    def run(self):
        MUSenumCount = self.initialThresholdCheck()
        print("initial MUS count", MUSenumCount)
        if MUSenumCount < self.tresh:
            print("a MUS enumerator identified within a timelimit of 10 seconds only {} MUSes, hence, this is either the exact MUS count or the enumeration is too expensive due to the hardness of the underlying SAT solver calls. Hence, we do not proceed with AMUSIC".format(MUSenumCount))
            return
        start = time.time()
        m = 0
        for iteration in range(self.t):
            print("iteration: " + str(iteration + 1))
            count, m = self.approxMC(m)
            if count > 0:
                self.counts.append(count)
                t = float(time.time() - start)
                print("# of MUSes in the last iteration: {}, average: {}, median: {}, m: {}, checks: {}, MUSes: {}, time: {}".format(count, sum(self.counts)/len(self.counts), median(self.counts), m, self.checks, len(self.MUSes), t))
            else:
                print("bsat failed")


def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

import sys
if __name__ == "__main__":
    parser = argparse.ArgumentParser("AMUSIC - Probabilistic Approximate Counter of Minimal Unsatisfiable Subsets") 
    parser.add_argument("--verbose", "-v", action="count", help = "Use the flag to increase the verbosity of the outputs. The flag can be used repeatedly.")
    parser.add_argument("--epsilon", "-e", type = float, help = "Set the epsilon parameter, i.e., controls the approximation factor of the algorithm. Allowed values: float (> 0). Default value is 0.8.", default = 0.8)
    parser.add_argument("--delta", "-d", type = restricted_float, help = "Set the delta parameter, i.e., controls the probabilistic guarantees of the algorithm. Allowed values: float (0-1). Default value is 0.2.", default = 0.2)
    parser.add_argument("--threshold", type = int, help = "Set manually the value of threshold. By default, the value of threshold is computed based on the epsilon parameter to guarantee the approximate guarantees that are required/set by epsilon. If you set threshold manually, you affect the guaranteed approximate factor of the algorithm.")
    parser.add_argument("--iterations", type = int, help = "Set manually the number of iterations the algorithm performs to find the MUS count estimate. By default, the number of iterations is determined by the value of the delta parameter (which controls the required probabilistic guarantees). By manually setting the number of iterations, you affect the probabilistic guarantees.")
    parser.add_argument("--qbf2", action="count", help = "Use the 2QBF encoding for finding an MUS in the cell instead of the default 3QBF encoding.")
    parser.add_argument("--qbf2beta", action="count", help = "Use the 2QBF encoding for finding an MUS in the cell instead of the default 3QBF encoding.")
    parser.add_argument("input_file", help = "A path to the input file. Either a .cnf or a .gcnf instance. See ./examples/")
    parser.add_argument("--findm-algo", default = "hybridSearch")
    args = parser.parse_args()

    counter = Counter(args.input_file, args.epsilon, args.delta)
    if args.threshold is not None:
        counter.tresh = args.threshold
    if args.iterations is not None:
        counter.t = args.iterations
    counter.QBF = "2QBF" if args.qbf2 is not None else "3QBF"
    if args.qbf2beta:
        counter.QBF = "2QBFBeta"
    counter.findMAlgo = args.findm_algo
    
    print("epsilon guarantee:", args.epsilon)
    print("delta guarantee:", args.delta)
    print("threshold", counter.tresh)
    print("iterations to complete:", counter.t)

    #clean temporal files in case of timeout or other kind of interruption
    signal.signal(signal.SIGHUP, partial(receiveSignal, counter.tmpFiles))
    signal.signal(signal.SIGINT, partial(receiveSignal, counter.tmpFiles))
    signal.signal(signal.SIGTERM, partial(receiveSignal, counter.tmpFiles))

counter.run()
