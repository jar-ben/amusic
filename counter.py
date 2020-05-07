import sys
from math import log
import subprocess as sp
#sys.path.insert(0, "/home/xbendik/usr/lib/lib/python3.7/site-packages")
import random
import time
from statistics import median
from random import randint

#parse .gcnf instance, 
#returns a pair C,B where B contains the base (hard) clauses and C the other clauses
def parse(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        assert lines[0][0] == "p"
        C = []
        B = []
        for line in lines[1:]:
            line = line.split(" ")
            cl = [int(i) for i in line[1:-1]]
            if len(cl) > 0:
                if line[0] == "{0}":
                    B.append(cl)
                else:
                    C.append(cl)
    return C,B

#returns random Boolean value
def randomBool():
    return bool(random.getrandbits(1))

class Counter:
    def __init__(self, filename, e, d):
        self.filename = filename
        self.C, self.B = parse(filename)
        self.MUSes = []
        self.dimension = len(self.C)
        self.XOR = None
        self.tresh = 1 + 9.84 * (1 + (e / (1 + e)))*(1 + 1/e)*(1 + 1/e)
        self.t = int(17 * log(3 / d,2));
        #self.tresh = 30 #lowered just for experimental reasons
        #self.t = 10 #lowered just for experimental reasons
        self.checks = 0
        self.rid = randint(1,10000000)

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
        unexXor = "unex_{}.cnf".format(self.rid)
        with open(unexXor, "w") as f:
            f.write("p cnf 0 0\n")
            for MUS in self.MUSes:
                f.write(" ".join([str(-l) for l in MUS]) + " 0\n")
                f.write(" ".join([str(l) for l in self.complement(MUS)]) + " 0\n")
            #aternatively, we block only the MUS itself
            #for MUS in self.MUSes:
            #    f.write(" ".join([str(-l) for l in MUS]) + " ")
            #    f.write(" ".join([str(l) for l in self.complement(MUS)]) + " 0\n")
            f.write(self.exportXor(m))
        cmd = "python3 gqbf.py {} {}".format(filename, unexXor)
        proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
        (out, err) = proc.communicate()
        out = out.decode("utf-8")
        if not "SOLUTION" in out:
            return []
        reading = False
        for line in out.splitlines():
            if reading:
                return [int(l) for l in line.split(" ") if int(l) > 0]
            if "SOLUTION" in line:
                reading = True

    # returns True if MUS is in the cell and False otherwise
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
    def bsatXor(self, m):
        found = 0
        self.MUSes = []
        for MUS in self.MUSes:
            if self.isInCell(MUS, m):
                found += 1
                if found >= self.tresh: return found
        while found < self.tresh:
            MUS = self.getMUS(m)
            if len(MUS) == 0:
                return found
            self.MUSes.append(MUS)
            found += 1
        return found

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
        return self.logSatSearch(mPrev)

    def run(self):
        start = time.time()
        counts = []
        m = int(self.dimension / 2)
        for iteration in range(self.t):
            print("iteration: " + str(iteration))
            count, m = self.approxMC(m)
            if count > 0:
                counts.append(count)
                t = float(time.time() - start)
                print("# of MUSes in the last iteration: {}, average: {}, median: {}, m: {}, checks: {}, MUSes: {}, time: {}".format(count, sum(counts)/len(counts), median(counts), m, self.checks, len(self.MUSes), t))
            else:
                print("bsat failed")


import sys
if __name__ == "__main__":
    assert len(sys.argv) == 2
    filename = sys.argv[1]
    counter = Counter(filename, 0.8, 0.2)
    print("thresh", counter.tresh)
    print("t:", counter.t)
    counter.run()
