import sys
from math import log
import subprocess as sp
import random
from random import randint
import time
from statistics import median
import os

def is_number(n):
    is_number = True
    try:
        num = float(n)
    except ValueError:
        is_number = False
    return is_number

def isClause(line):
    sline = line.rstrip().split()
    return sline[-1] == "0" and set([is_number(l) for l in sline]) == set([True])

def parse(filename):
    with open(filename, "r") as f:
        C = []
        for line in f.readlines():
            if isClause(line.rstrip()):
                C.append([int(l) for l in line.rstrip().split()[:-1]])
    return C

def variables(cls):
    vrs = []
    for cl in cls:
        for l in cl:
            if l < 0: l *= -1
            vrs += [l]
    return list(set(vrs))

def renderCnf(cls):
    nvariables = max(variables(cls))
    clauses = len(cls)
    result = "p cnf {} {}\n".format(nvariables, clauses)
    for cl in cls:
        result += " ".join([str(l) for l in cl]) + " 0\n"
    return result

def renderWcnf(Hard, Soft):
    nvariables = max(variables(Hard + Soft))
    clauses = len(Hard) + len(Soft)
    hardWeight = len(Soft) + 1
    softWeight = 1

    result = "p wcnf {} {} {}\n".format(nvariables, clauses, hardWeight)
    for cl in Hard:
        result += str(hardWeight) + " " + " ".join([str(l) for l in cl]) + " 0\n"
    for cl in Soft:
        result += str(softWeight) + " " + " ".join([str(l) for l in cl]) + " 0\n"
    
    return result

def maxSat(Hard, Soft):
    satTimeout = 180
    wcnf = renderWcnf(Hard, Soft)
    file = "./tmp/testAutarky_{}.wcnf".format(randint(1,100000000))
    with open(file, "w") as f:
        f.write(wcnf)

    with open(file+".cnf", "w") as f:
        f.write(renderCnf(Hard+Soft))

    cmd = 'timeout {} ./tools/uwrmaxsat -m {}'.format(satTimeout, file)
    proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")         
    print(cmd)   
    model = []
    for line in out.splitlines():
        if line[:2] == "v ": 
            model = [int(l) for l in line.rstrip().split()[1:]]
   

    os.remove(file)
    os.remove(file+".cnf")
    return [x for x in model if x > 0]

def getAutarkyClauses(autarkyVars, C):
    av = set(autarkyVars)
    cls = []
    for i in range(len(C)):
        cl = C[i]
        cl = [abs(l) for l in cl]
        if len(av.intersection(set(cl))) == 0:
            cls.append(i)
    return cls

def exportAutarky(a, C, target):
    mVar = max(variables(C))
    dim = len(a)
    res = "p cnf {} {}\n".format(mVar, dim)
    for c in a:
        res += " ".join([str(l) for l in C[c]]) + " 0\n"
    open(target, "w").write(res)


#uses encoding based on the model 2 from Efficient Autarkies
def findAutarky(filename, target):
    C = parse(filename)
    X = variables(C)
    n = len(X)
    X1 = [n + i + 1 for i in range(n)]
    X0 = [2*n + i + 1 for i in range(n)]
    Xp = [3*n + i + 1 for i in range(n)]


    F = []
    #clauses 1 (equations (6) and (7))
    for c in C:
        for l in c:
            cl = []
            if l > 0: #rule 1a)
                cl.append(-X0[l - 1])
            else:
                cl.append(-X1[(-1 *l) - 1])
            for l2 in c:
                    if l2 > 0 and l2 != l:
                        cl.append(X1[l2 - 1])
                    elif l2 < 0 and l2 != l:
                        cl.append(X0[(-1 * l2) - 1])
            F.append(cl)

    #clauses 2 (equation 8)
    for x in X:
        F.append([-X1[x - 1], -X0[x - 1]])

    #clauses 3 (equation 9)
    for x in X:
        F.append([-Xp[x - 1], X1[x - 1], X0[x - 1]])
        F.append([Xp[x - 1], -X1[x - 1]])
        F.append([Xp[x - 1], -X0[x - 1]])

    #enforce unsat
    Soft = []
    for x in Xp:
        Soft.append([x])

    autarkyVars = [x - 3*n for x in maxSat(F, Soft)]
    autarkyClauses = getAutarkyClauses(autarkyVars, C)
    print("autarky vars size:", len(autarkyVars))
    print("autarky clauses size:", len(autarkyClauses))
    print("v " + " ".join([str(x + 1) for x in autarkyClauses]))
    #exportAutarky(autarkyClauses,C, target)

if __name__ == "__main__":
    assert len(sys.argv) > 1
    filename = sys.argv[1]
    target = sys.argv[2] if len(sys.argv) > 2 else "autarky.cnf"
    findAutarky(filename, target)
