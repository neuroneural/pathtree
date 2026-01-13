""" This module contains clingo interaction functions """
from __future__ import print_function

import subprocess
from subprocess import CalledProcessError
import sys, os
import numpy as np
import bfutils as bfu

THREADS='-t 20,split'
CLINGOPATH='/na/homes/splis/soft/tools/python-gringo/clingo-4.5.1-source/build/release/'
CAPSIZE=1000

def g2clingo_(g, file=sys.stdout):
    """ Save a graph to a file of grounded terms for clingo """
    n = len(g)
    print('node(1..'+str(n)+').', file=file)
    for v in g:
        for w in g[v]:
            if g[v][w] == 1: print('edgeu('+str(v)+','+str(w)+').', file=file)
            if g[v][w] == 2: print('confu('+str(v)+','+str(w)+').', file=file)
            if g[v][w] == 3:
                print('edgeu('+str(v)+','+str(w)+').', file=file)
                print('confu('+str(v)+','+str(w)+').', file=file)

def g2clingo(g, file=sys.stdout):
    """ Save a graph to a file of grounded terms for clingo """
    n = len(g)
    print('node(1..'+str(n)+').', file=file)
    for v in g:
        for w in g[v]:
            if (0,1) in g[v][w]: print('edgeu('+v+','+w+').', file=file)
            if (2,0) in g[v][w]: print('confu('+v+','+w+').', file=file)

def c2edgepairs(clist):
    return [x[6:-1].split(',') for x in clist]
def nodenum(edgepairs):
    nodes = 0
    for e in edgepairs:
        nodes = np.max([nodes, np.max(map(int,e))])
    return nodes
def edgepairs2g(edgepairs):
    n = nodenum(edgepairs)
    g = {str(x+1):{} for x in range(n)}
    for e in edgepairs:
        g[e[0]][e[1]] = set([(0,1)])
    return g

def filterAnswers(slist):
    alist = []
    anAnswer = False
    for e in slist:
        if anAnswer:
            alist.append(e.split(' '))
            anAnswer=False
        if  e[:6] == "Answer":
            anAnswer = True
    return alist

def clingo(g,
           timeout=0,
           threads=THREADS,
           capsize=CAPSIZE,
           graphfile='gu.pl',
           ufile='drawu.pl',
           program='supersample.pl',
           cpath=CLINGOPATH,
           nameid=''):

    cmdline = cpath+'clingo '+threads+' --time-limit='+str(timeout)+' -n '+str(capsize)+' '+cpath+graphfile+' '+cpath+ufile+' '+cpath+program
    with open(cpath+graphfile,'w') as f:
        g2clingo(g,file=f)
    try:
        p = subprocess.Popen(cmdline.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        pass
    p_status = p.wait()
    (output, err) = p.communicate()

    os.remove(cpath+graphfile)
    return clingo2g(output)

def a2edgetuple(answer):
    edges = [x for x in answer if 'edge1' in x]
    u = [x for x in answer if x[0]=='u'][0]
    return edges,u

def clingo2g(output):
    s = set()
    answers = filterAnswers(output.split('\n'))
    answers = [a2edgetuple(x) for x in answers]
    l = [(c2edgepairs(x[0]),x[1]) for x in answers]
    l = [(bfu.g2num(edgepairs2g(x[0])),int(x[1][2:-1])) for x in l]
    return l
