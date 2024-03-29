import sys
from random import randint
import signal
from functools import partial
import os

def receiveSignal(tempFiles, signalNumber, frame):
    print('Received signal:', signalNumber)
    print('Cleaning tmp files')
    for f in tempFiles:
        if os.path.exists(f):
            os.remove(f)
            print("removed")
    sys.exit()

def prnt(cnf):
    result = "\n\n\n"
    for cl in cnf:
        result += " ".join([str(l) for l in cl]) + " 0\n"
    print(result + "\n\n\n")

def variables(cls):
    vrs = []
    for cl in cls:
        for l in cl:
            if l < 0: l *= -1
            vrs += [l]
    return list(set(vrs))

## input: l: a list of lists of numbers
## output: a maximum number in l
def maxVar(l):
    maximum = 0
    for i in l:
        maximum = max(maximum, max([abs(j) for j in i]))
    return maximum

def tseitinCube(cube, current):
    current += 1
    res = []
    for l in cube:
        res.append([-current,l])
    res.append([current] + [-l for l in cube])
    return res, current

def tseitinClause(clause, current):
    current += 1
    res = []
    for l in clause:
        res.append([current,-l])
    res.append([-current] + [l for l in clause])
    return res, current

## input: cnf: a cnf formula, given as a list of lists of literals
## input: current: available variable index (i.e. not representing any variable yet)
## output: cnf formula logacically equivalent to the input formula, however transformed via the tseitin transformation
## that is, several tseitin activation variables are created (starting with current) and the original cnf formula is now
## true if an only if the main tseitin activation literal holds. All the added activation literals are returned in the list "acts",
## the main activation literal is the last one in the list. 
def tseitinOnCnf(cnf, current):
    activators = []
    acts = []
    res = []
    for cl in cnf:
        tscl, current = tseitinClause(cl, current)
        activators.append(current)
        res += tscl
    current += 1
    for a in activators:
        res.append([a,-current])
    res.append([current] + [-1 * a for a in activators])    
    return res, current

def tseitinOnDnf(dnf, current):
    activators = []
    acts = []
    res = []
    for cl in dnf:
        tscl, current = tseitinCube(cl, current)
        activators.append(current)
        res += tscl
    current += 1    
    res.append([-current] + activators)
    for a in activators:
        res.append([-a,current])
    return res, current

## input: cl: a XOR represented as a set of literals
## output: tseitin representation of cl
def tseitinXOR(cl, current):
    if len(cl) == 1:
        return [cl], current
    fst = cl[0]
    cnf = []
    for snd in cl[1:]:
        c, current = tseitinOnCnf([[fst, snd],[-fst, -snd]], current)
        cnf += c
        fst = current
    return cnf, current

## input: two tseintin cnf formula, cn1 and cnf2, with main activators act1 and act2
## output: tseintin OR of cn1 and cn2 with main activator current
def orTwoTseitinCnf(cnf1, cnf2, act1, act2, current):
    current += 1
    return cnf1 + cnf2 + [[-current, act1, act2], [-act1, current], [-act2, current]], current

from pysat.card import *
def notClosestSubset(A, B, current):
    supset = []
    for i in range(len(A)):
        supset.append([-B[i], A[i]])

    Z = []
    ZActs = []
    for i in range(len(A)):
        cnf, current = tseitinOnCnf([[-A[i], -B[i]], [A[i], B[i]]], current)
        Z += cnf
        ZActs.append(current)
    al2 = CardEnc.atleast(lits=ZActs, bound=2, encoding=7, top_id = current)
    current = max(current, maxVar(al2.clauses))
    supsetT, supsetTAct = tseitinOnCnf(supset, current)
    current = supsetTAct
    al2T, al2TAct = tseitinOnCnf(al2.clauses, current)
    current = al2TAct
    
    #either is A a superset of B, or |B-A| > 1
    main, current = tseitinClause([supsetTAct, al2TAct], current)
    return main + supsetT + al2T + Z, current

## intput: C: set C = [c1, ..., cn] of constraints
## input: activators: set of activation literals [x1, ..., xn]
## output: tseitin cnf formula saying that C using the activators, i.e. \bigwedge_{i}(ci or -xi), is unsatisfiable
def unsat(C, B, activators, current):
    assert len(C) == len(activators)
    n = len(C)
    dnf = []
    for i in range(n):
        dnf.append([activators[i]] + [-l for l in C[i]])
    for b in B:
        dnf.append([-l for l in b])
    return tseitinOnDnf(dnf, current)

## intput: C: set C = [c1, ..., cn] of constraints
## input: activators: set of activation literals [x1, ..., xn]
## output: tseitin cnf formula saying that C using the activators, i.e. \bigwedge_{i}(ci or -xi), is satisfiable
def sat(C, B, activators, excluded = []):
    assert len(C) == len(activators)
    n = len(C)
    cnf = []
    for i in range(n):
        if i in excluded: continue
        cnf.append([-activators[i]] + [l for l in C[i]])
    cnf += B
    return cnf
    #return tseitinOnCnf(cnf, current)

def sign(l):
    if l > 0: return 1
    else: return -1

def parseUnex(filename, activators, current):
    with open(filename, "r") as f:
        lines = f.readlines()
        assert lines[0][0] == "p"
        C = []
        for line in lines[1:]:
            line = line.split(" ")
            if line[0] == "x":
                lits = [int(i) for i in line[1:-1]]
                lits = [sign(l) * activators[abs(l)-1] for l in lits]
                cl, current = tseitinXOR(lits, current)
                C += cl + [[current]]
            else:
                cl = [int(i) for i in line[:-1]]
                if len(cl) > 0:
                    C.append([sign(l) * activators[abs(l)-1] for l in cl])
        return C, current

def primeCls(C, offset):
    primeC = []
    for cl in C:
        primeCl = []
        for l in cl:
            if l > 0: primeCl.append(l + offset)
            else: primeCl.append(l - offset)
        primeC.append(primeCl)
    return primeC

import copy
def exMUS(constraints, unex):
    C,B  = parse(constraints)
    Vars = variables(C + B) ##these are in the 2nd quantificator level and are Exists
    n = len(C)
    m = max(Vars)

    primesC = []
    primesB = []
    for i in range(1, n + 1):
        primesC.append(primeCls(C, m * i))
        primesB.append(primeCls(B, m * i))

    primeVars = []
    for primes in primesC + primesB:
        primeVars += variables(primes)
    primeVars = list(set(primeVars))

    activators = [max(primeVars) + i + 1 for i in range(n)] ##these are in 1st quantificator level and are ForAll
    current = activators[-1]
    X = []

    #P is in the cell and unexplored
    tmp_current = current
    unexAndXor, current = parseUnex(unex, activators, current)
    unexCurrents = [i + 1 for i in range(tmp_current, current)]
    #prnt(unexAndXor)
    X += unexAndXor

    #print("activators", activators)
    #sol = [1208, -1209, 1210, 1211, 1212, 1213, 1214, 1215, 1216, -1217, -1218, -1219, -1220, -1221, -1222, -1223]
    #X += [[a] for a in sol]    

    #encode subsets of P
    for i in range(n):
        #(\neg p_i \vee \mathtt{sat(S_i)})
        satc = sat([primesC[i][j] + [-activators[i]] for j in range(n)], primesB[i], activators, excluded = [i])
        X += satc

    #X \implies sat(P)
    PsatCNF, PsatAct = tseitinOnCnf(sat(C, B, activators), current)
    current = PsatAct
    XCNF, XAct = tseitinOnCnf(X, current)    
    current = XAct
    main = PsatCNF + XCNF + [[-XAct, PsatAct]]
   
    currents = [i for i in range(activators[-1] + 1, current + 1) if i not in unexCurrents]

    result = "p cnf {} {}\n".format(maxVar(main), len(main))
    result += "a " + " ".join([str(i) for i in activators + primeVars + unexCurrents]) + " 0 \n"
    result += "e " + " ".join([str(i) for i in Vars + currents]) + " 0 \n"
    for cl in main:
        result += " ".join([str(l) for l in cl]) + " 0\n"

    
    return result, activators

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

import subprocess as sp

def compute(filename, activators):
    cmd = "./tools/caqe --qdo {}".format(filename)
    proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")

    N = []
    if "Satisfiable" in out:
        Cids = [0 for _ in range(activators[-1] + 1)]
        for i in range(len(activators)):
            act = activators[i]
            Cids[act] = i + 1
        reading = False
        N = []
        for line in out.splitlines():
            if reading:
                if len(line) == 0 or line[0] != "V": break
                cl = int(line.split(" ")[1])
                if abs(cl) < len(Cids):
                    N.append(sign(cl) * Cids[abs(cl)])
            if "s cnf" in line:
                reading = True
        print("SOLUTION")
        print(" ".join([str(n) for n in N]) + " 0")
    else:
        print("UNSATISFIABLE")

def computeCadet(filename, activators):
    cmd = "./tools/cadet --minimize {}".format(filename)
    proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")
    assert "SAT" in out
    if not "UNSAT" in out: 
        print("UNSATISFIABLE") #this means that there is no unexplored MUS in the cell
    else:        
        Cids = {}
        for i in range(len(activators)):
            Cids[activators[i]] = i + 1
        reading = False
        for line in out.splitlines():
            if reading:
                assert line[0] == "V"
                MUS = [Cids[int(l)] for l in line.split(" ")[1:] if int(l) in Cids]
            if "UNSAT" in line:
                reading = True
        print("SOLUTION")
        print(" ".join([str(n) for n in MUS]) + " 0")

def simplify2(filename, result):
    cmd = "./tools/qratpre+ --no-ble --no-qratu --no-qrate --print-formula {} > {}".format(filename, result)
    #cmd = "./tools/qratpre+ --print-formula {} > {}".format(filename, result)
    proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")
    print("simplified")

if __name__ == "__main__":
    assert len(sys.argv) > 2
    filename = sys.argv[1]
    unex = sys.argv[2]
    encoding, activators = exMUS(filename, unex)
    rid = randint(1,10000000)
    qdimacs = "./tmp/2exencoded_{}.qdimacs".format(rid)
    simpl = "./tmp/simplified_{}.qdimacs".format(rid)

    #clean up tmp files in case of timeout or other kind of interruption
    tmpFiles = [qdimacs, simpl]
    signal.signal(signal.SIGHUP, partial(receiveSignal, tmpFiles))
    signal.signal(signal.SIGINT, partial(receiveSignal, tmpFiles))
    signal.signal(signal.SIGTERM, partial(receiveSignal, tmpFiles))

    with open(qdimacs, "w") as f:
        f.write(encoding)
    print(qdimacs, simpl)
    #simplify2(qdimacs, simpl)
    simpl = qdimacs
    computeCadet(simpl, activators)
    #compute(simpl, activators)
    if os.path.exists(qdimacs): os.remove(qdimacs)
    if os.path.exists(simpl): os.remove(simpl)
