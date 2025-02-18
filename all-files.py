import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import bfutils as bfu
import graphkit as gk
import traversal as trv
import linear_model as lm
import numpy as np
import unknownrate as ur
import pc
import zickle as zkl
import simpleloops as sls
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import bfutils as bfu
import graphkit as gk
import traversal as trv
import linear_model as lm
import numpy as np
import unknownrate as ur
import pc
import zickle as zkl
import simpleloops as sls
import functools
import sys
from collections import Counter
sys.path.append('tools')

# local package
from bfutils import *

nodes = 5 # number of nodes in the graph
PNUM = 75 # number of processes to use

template = "{0:9}{1:9}{2:9}{3:9}{4:9}{5:10}"

def wrapper_c(fold):
    return icompat(fold, nodes)
def wrapper_l(fold):
    return ilength(fold, nodes)
def wrapper_list(fold):
    return iall(fold, nodes)

def wrapper_u(fold, steps):
    return cc_all(fold, nodes, steps)

def make_rect(l):
    max_seq = max(map(len,l))
    for e in l:
        e += [e[-1]] * (max_seq - len(e))
    return l

def wrapper_unique(fold):
    counter = 0
    for i in results[number]:
        if i not in unique_appeared_graph[counter]:
            unique_appeared_graph[counter].append(i)
        counter += 1
    return 1

def wrapper_non_eliminate(fold):
	return tuple(results[fold][counter])
resultsmp=True
results=[]
pool=Pool(processes=PNUM)

print template.format(*('        u', ' all_uniq', '   unique', '  seen_Gu', ' converge',  ' uconverge'))
cumset = set()
clen = 0
for s in xrange(100):

    results=pool.map(functools.partial(wrapper_u, steps=s),
                     range(2**(nodes**2)))
    
    converged = len([e for e in results if e[1]==[]])
    notconverged = len(results) - converged

    if notconverged < 100 and notconverged > 0:
        survivors = [e for e in results if e[1]]
    if notconverged == 0: break
    results = filter(lambda (x): x[1] != [], results)

    r = set(results)
    d = r.difference(cumset)
    cumset = cumset.union(r)

    cl = 2**(nodes**2) - len(results) - clen
    clen += cl

    print template.format(*(s, len(r), len(d), len(cumset),
                            converged, notconverged))

pool.close()
pool.join()
import sys, os
sys.path.append('./tools/')
import traversal, bfutils, graphkit
import unknownrate as ur
from multiprocessing import Pool,Process, Queue, cpu_count, current_process, active_children
import functools
import zickle as zkl
import time, socket
import scipy

KEY='rasl_il_u2'
UMAX = 2
INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 1000 # stop traversing after growing equivalence class tothis size
REPEATS = 100
if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=60
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'mars':
    PNUM=20
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'hooke':
    PNUM=21
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'saturn':
    PNUM=12
    PNUM=max((1,PNUM/INPNUM))
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count())
    PNUM=max((1,PNUM/INPNUM))
print 'processes: ',PNUM, INPNUM

def multiprocess(argslist, ncpu):
    total = len(argslist)
    done = 0
    result_queue = Queue()
    jobs = []

    def ra_wrapper_(fold, n=10, k=10):
        scipy.random.seed()
        l = {}
        while True:
            try:
                g = bfutils.ringmore(n,k) # random ring of given density
                gs= bfutils.call_undersamples(g)
                for u in range(1,min([len(gs),UMAX])):
                    g2 = bfutils.undersample(g,u)
                    print fold,': ',traversal.density(g),':',
                    startTime = int(round(time.time() * 1000))
                    s = ur.liteqclass(g2, verbose=False, capsize=CAPSIZE)
                    endTime = int(round(time.time() * 1000))
                    print len(s), u
                    l[u] = {'eq':s,'ms':endTime-startTime}
            except MemoryError:
                print 'memory error... retrying'
                continue
            break
        result_queue.put( {'gt':g,'solutions':l} )

    while argslist != [] and done<10 :
        if len(active_children()) < ncpu:
            p = Process(target=ra_wrapper_,args=(argslist.pop(),))
            jobs.append(p)
            p.start()
            done+=1
            print "\r",float(done)/total,"%",
    #get results here
    res = [result_queue.get() for p in jobs]
    print res

def ra_wrapper(fold, n=10, k=10):
    scipy.random.seed()
    l = {}
    while True:
        try:
            g = bfutils.ringmore(n,k) # random ring of given density
            gs= bfutils.call_undersamples(g)
            for u in range(1,min([len(gs),UMAX])):
                g2 = bfutils.undersample(g,u)
                print fold,': ',traversal.density(g),':', traversal.density(g2),':',
                startTime = int(round(time.time() * 1000))
                #s = ur.liteqclass(g2, verbose=False, capsize=CAPSIZE)
                s = ur.eqclass(g2)
                endTime = int(round(time.time() * 1000))
                print len(s), u
                l[u] = {'eq':s,'ms':endTime-startTime}
        except MemoryError:
            print 'memory error... retrying'
            continue
        break
    return {'gt':g,'solutions':l}

def ra_wrapper_preset(fold, glist=[]):
    scipy.random.seed()
    l = {}
    while True:
        try:
            g = glist[fold]
            gs= bfutils.call_undersamples(g)
            for u in range(1,min([len(gs),UMAX])):
                g2 = bfutils.undersample(g,u)
                print fold,': ',traversal.density(g),':',
                startTime = int(round(time.time() * 1000))
                s = ur.liteqclass(g2, verbose=False, capsize=CAPSIZE)
                endTime = int(round(time.time() * 1000))
                print len(s), u
                l[u] = {'eq':s,'ms':endTime-startTime}
        except MemoryError:
            print 'memory error... retrying'
            continue
        break
    return {'gt':g,'solutions':l}

def killall(l):
    for e in l:
        e.join(timeout=0.001)
        if not e.is_alive():
            #print 'first result'
            for p in l:
                if p != e:
                    #print 'terminating ', p.name
                    p.terminate()
                    p.join()
                else:
                    p.join()
            return True
    return False

def fan_wrapper(fold,n=10,k=10):
    scipy.random.seed()
    curr_proc=current_process()
    curr_proc.daemon=False
    output = Queue()
    while True:
        try:
            g = bfutils.ringmore(n,k)
            gdens = traversal.density(g)
            g2 = bfutils.increment_u(g,g)
            #g2 = bfutils.undersample(g,2)
            def inside_wrapper():
                scipy.random.seed()
                try:
                    startTime = int(round(time.time() * 1000))
                    #s = traversal.v2g22g1(g2, capsize=CAPSIZE)
                    s = traversal.backtrack_more2(g2, rate=2, capsize=CAPSIZE)
                    endTime = int(round(time.time() * 1000))
                    print "{:2}: {:8} : {:4}  {:10} seconds".\
                        format(fold, round(gdens,3), len(s),
                               round((endTime-startTime)/1000.,3))
                    output.put({'gt':g,'eq':s,'ms':endTime-startTime})
                except MemoryError:
                    print 'memory error...'
		    raise
            pl = [Process(target=inside_wrapper) for x in range(INPNUM)]
            for e in pl: e.start()
            while True:
                if killall(pl): break
            r = output.get()
        except MemoryError:
            print 'memory error... retrying'
            for p in pl:
                p.terminate()
                p.join()
            continue
        break
    for p in pl: p.join()
    return r

densities = {5: [0.2],
             6: [0.2, .25, .3],
             7: [0.2, .25, .3],
             8: [0.15, 0.2, 0.25, 0.3],
             9: [.2],
             10:[.2],
             15:[0.2],
             20:[0.1],# 0.15, 0.2, 0.25, 0.3],
             25:[0.1],
             30:[0.1],
             35:[0.1],
             40:[0.1],
             50:[0.05, 0.1],
             60:[0.05, 0.1]}

for nodes in [5]:
    z = {}
    #pool=Pool(processes=PNUM)
    for dens in densities[nodes]:
        print "{:2}: {:8} : {:10} : {:10}  {:10}".format('id', 'densityi(G)', 'density(H)', 'eq class', 'time')
        e = bfutils.dens2edgenum(dens, n=nodes)

        eqclasses = []
        for x in pool.imap(functools.partial(ra_wrapper, n=nodes, k=e),
                           range(REPEATS)):
            eqclasses.append(x)
            z[dens] = eqclasses
            zkl.save(eqclasses,
                     socket.gethostname().split('.')[0]+\
                     '_nodes_'+str(nodes)+'_density_'+\
                     str(dens)+'_'+KEY+'_.zkl')

        print ''
        print '----'
        print ''
        #pool.close()
        #pool.join()
        #zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_'+KEY+'_.zkl')
import seaborn as sb
import zickle as zkl
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append('./tools/')
densities = [.15, 0.2, 0.25, 0.3, 0.35]

def gettimes(d):
    t = [x['ms'] for x in d]
    time  = map(lambda x: x/1000./60., t)
    return time

l = ['leibnitz_nodes_15_density_0.1_newp_.zkl',
     'leibnitz_nodes_20_density_0.1_newp_.zkl',
     'leibnitz_nodes_25_density_0.1_newp_.zkl',
     'leibnitz_nodes_30_density_0.1_newp_.zkl',
     'leibnitz_nodes_35_density_0.1_newp_.zkl']

alltimes_new = []
for fname in l:
    d = zkl.load(fname)
    alltimes_new.append(gettimes(d))

shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

plt.figure(figsize=[10,2])


g = sb.boxplot(alltimes_new,names=map(lambda x: str(int(x*100))+"",
                                      densities),
               widths=wds, color="Reds",fliersize=fliersz,
               linewidth=lwd,
               **{'positions':np.arange(len(densities))+shift,
                  'label':'MSL'})

# plt.plot(np.arange(len(densities))-shift,
#         map(np.median,alltimes_old), 'ro-', lw=0.5, mec='k')
# plt.plot(np.arange(len(densities))+shift,
#         map(np.median,alltimes_new), 'bo-', lw=0.5, mec='k')
g.figure.get_axes()[0].set_yscale('log')
plt.xlabel('number of nodes in a graph')
plt.ylabel('computation time (minutes)')
#plt.title('100 6 node graphs per density\n$G_2 \\rightarrow G_1$',
#          multialignment='center')
plt.subplots_adjust(right=0.99, left=0.2)
plt.legend(loc=0)
plt.show()
import pandas as pd
import pylab as pl
import seaborn as sb
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np

SBDIR = '~/soft/src/dev/tools/stackedBarGraph/'
GFDIR = '/na/home/splis/soft/src/dev/craft/gunfolds/tools/'
import sys, os

sys.path.append(os.path.expanduser(SBDIR))
sys.path.append(os.path.expanduser(GFDIR))

import zickle as zkl
from stackedBarGraph import StackedBarGrapher
SBG = StackedBarGrapher()

def gettimes(d):
    t = [x['ms'] for x in d]
    time  = map(lambda x: x/1000./60., t)
    return time

l = [(0.15, 'leibnitz_nodes_15_density_0.1_newp_.zkl'),
     (0.20, 'leibnitz_nodes_20_density_0.1_newp_.zkl'),
     (0.25, 'leibnitz_nodes_25_density_0.1_newp_.zkl'),
     (0.30, 'leibnitz_nodes_30_density_0.1_newp_.zkl'),
     (0.35, 'leibnitz_nodes_35_density_0.1_newp_.zkl')]

fig = pl.figure(figsize=[10,3])
#Read in data & create total column

d = zkl.load("hooke_nodes_6_g32g1_.zkl")#hooke_nodes_35_newp_.zkl")
densities = [.15, .20 ,.25, .30, .35]
d = {}
for fname in l:
    d[fname[0]] = zkl.load(fname[1])

def get_counts(d):
    eqc = [len(x['eq']) for x in d]
    keys = np.sort(np.unique(eqc))
    c = {}
    for k in keys:
        c[k] = len(np.where(eqc == k)[0])
    return c

# unique size
usz = set()
dc = {}
for u in densities:
    dc[u] = get_counts(d[u])
    for v in dc[u]:
        usz.add(v)

for u in densities:
    for c in usz:
        if not c in dc[u]:
            dc[u][c] = 0

A = []
for u in densities:
    A.append([dc[u][x] for x in np.sort(dc[u].keys())])
    
#print A
#A = np.array(A)
pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.color_palette("Paired",len(usz)))
#pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.dark_palette("#5178C7",len(usz)))
#pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.blend_palette(["mediumseagreen", "ghostwhite", "#4168B7"],len(usz)))
scalarMap = mpl.cm.ScalarMappable(norm = lambda x: x/np.double(len(usz)),
                                  cmap=pp)

d_widths = [.5]*len(densities)
d_labels = map(lambda x: str(int(x*100))+"%",densities)
#u = np.sort(list(usz))
d_colors = [scalarMap.to_rgba(i) for i in range(len(A[0]))]
#d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', '#fee090']

#ax = fig.add_subplot(211)
ax = plt.subplot2grid((3,1), (0, 0), rowspan=2)

SBG.stackedBarPlot(ax,
                   A,
                   d_colors,
                   xLabels=d_labels,
                   yTicks=3,
                   widths=d_widths,
                   gap = 0.005,
                   scale=False
)

for i in range(len(A)):
    Ai = [x for x in A[i] if x>0]
    y = [x/2.0 for x in Ai]
    for j in range(len(Ai)):
        if j>0:
            yy = y[j]+np.sum(Ai[0:j])
        else:
            yy = y[j]        
        pl.text(0.5*i-0.02,yy-1.2,str(Ai[j]),fontsize=12,zorder=10)
    
# #Set general plot properties
# sns.set_style("white")
# sns.set_context({"figure.figsize": (24, 10)})

# for i in np.sort(list(usz))[::-1]:
#     y = [100-dc[u][i] for u in np.sort(dc.keys())]
#     bottom_plot=sns.barplot(x=np.asarray(densities)*100, y=y)
#    #                         color=scalarMap.to_rgba(i))
#     #y = (sbd[i+1]-sbd[i])/2.+sbd[i]scala
#     #for j in range(len(sbd.Density)):
#     #    pl.text(j-0.1,y[j],'1',fontsize=16,zorder=i)

# #Optional code - Make plot look nicer
sns.despine(left=True)
# #Set fonts to consistent 16pt size
ax.set(xticklabels="")
for item in ([ax.xaxis.label, ax.yaxis.label] +
             #ax.get_xticklabels() +
             ax.get_yticklabels()):
    item.set_fontsize(12)



alltimes_new = []
for fname in l:
    dp = zkl.load(fname[1])
    alltimes_new.append(gettimes(dp))

shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

ax = plt.subplot2grid((3,1), (2, 0))
g = sb.boxplot(alltimes_new,names=map(lambda x: str(int(x*100))+"",
                                      densities),
               widths=wds, color="Reds",fliersize=fliersz,
               linewidth=lwd,
               **{'positions':np.arange(len(densities))+shift,
                  'label':'MSL'})

# plt.plot(np.arange(len(densities))-shift,
#         map(np.median,alltimes_old), 'ro-', lw=0.5, mec='k')
# plt.plot(np.arange(len(densities))+shift,
#         map(np.median,alltimes_new), 'bo-', lw=0.5, mec='k')
g.figure.get_axes()[1].set_yscale('log')
plt.xlabel('number of nodes in a graph')
plt.ylabel('computation time\n(minutes)')
#plt.title('100 6 node graphs per density\n$G_2 \\rightarrow G_1$',
#          multialignment='center')
#plt.subplots_adjust(right=0.99, left=0.2)
plt.legend(loc=0)
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() +
             ax.get_yticklabels()):
    item.set_fontsize(12)

pl.subplots_adjust(bottom=0.1,hspace=0.01,top=0.98)
# plt.show()

pl.show()
import seaborn as sb
import zickle as zkl
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append('./tools/')
densities = [.15, 0.2, 0.25, 0.3, 0.35]

def gettimes(d):
    t = [x['ms'] for x in d]
    time  = map(lambda x: x/1000./60., t)
    return time

def getalltimes(data, densities=densities):
    alltimes = []
    for dens in densities:        
        alltimes.append(gettimes(data[dens]))
    return alltimes

def timesfromfile(fname):
    d = zkl.load(fname)
    return getalltimes(d)

#alltimes_old = timesfromfile("hooke_nodes_8_old_nomem.zkl")
#alltimes_new = timesfromfile("hooke_nodes_8_newp_.zkl")

#alltimes_old = timesfromfile("hooke_nodes_10_newp_.zkl")
#alltimes_old = timesfromfile("hooke_nodes_6_g32g1_.zkl")
#alltimes_new = timesfromfile("leibnitz_nodes_35_newp_.zkl")


shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

plt.figure(figsize=[10,2])

# g = sb.boxplot(alltimes_old,names=map(lambda x: str(int(x*100))+"%",
#                                      densities),
#               widths=wds, color="Reds", fliersize=fliersz, linewidth=lwd,
#               **{'positions':np.arange(len(densities))-shift,
#                  'label':'naive approach'})

g = sb.boxplot(alltimes_new,names=map(lambda x: str(int(x*100))+"",
                                      densities),
               widths=wds, color="Blues",fliersize=fliersz,
               linewidth=lwd,
               **{'positions':np.arange(len(densities))+shift,
                  'label':'MSL'})

# plt.plot(np.arange(len(densities))-shift,
#         map(np.median,alltimes_old), 'ro-', lw=0.5, mec='k')
# plt.plot(np.arange(len(densities))+shift,
#         map(np.median,alltimes_new), 'bo-', lw=0.5, mec='k')
g.figure.get_axes()[0].set_yscale('log')
plt.xlabel('density (% of 36 total possible edges)')
plt.ylabel('computation time (minutes)')
plt.title('100 6 node graphs per density\n$G_2 \\rightarrow G_1$',
          multialignment='center')
plt.subplots_adjust(right=0.99, left=0.2)
plt.legend(loc=0)
plt.show()
import pandas as pd
import pylab as pl
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np

import sys
sys.path.append('/na/home/splis/soft/src/dev/craft/gunfolds/tools/')
sys.path.append('/na/homes/splis/soft/src/dev/tools/stackedBarGraph/')
import zickle as zkl
from stackedBarGraph import StackedBarGrapher
SBG = StackedBarGrapher()

fig = pl.figure(figsize=[10,1.3])
#Read in data & create total column

d = zkl.load("hooke_nodes_6_g32g1_.zkl")#hooke_nodes_35_newp_.zkl")
densities = np.sort(d.keys())

def get_counts(d):
    eqc = [len(x['eq']) for x in d]
    keys = np.sort(np.unique(eqc))
    c = {}
    for k in keys:
        c[k] = len(np.where(eqc == k)[0])
    return c

# unique size
usz = set()
dc = {}
for u in densities:
    dc[u] = get_counts(d[u])
    for v in dc[u]:
        usz.add(v)

for u in densities:
    for c in usz:
        if not c in dc[u]:
            dc[u][c] = 0

A = []
for u in densities:
    A.append([dc[u][x] for x in np.sort(dc[u].keys())])
    
#print A
#A = np.array(A)
pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.color_palette("Paired",len(usz)))
#pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.dark_palette("#5178C7",len(usz)))
#pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.blend_palette(["mediumseagreen", "ghostwhite", "#4168B7"],len(usz)))
scalarMap = mpl.cm.ScalarMappable(norm = lambda x: x/np.double(len(usz)),
                                  cmap=pp)

d_widths = [.5]*len(densities)
d_labels = map(lambda x: str(int(x*100))+"%",densities)
#u = np.sort(list(usz))
d_colors = [scalarMap.to_rgba(i) for i in range(len(A[0]))]
#d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', '#fee090']

ax = fig.add_subplot(111)
SBG.stackedBarPlot(ax,
                   A,
                   d_colors,
                   xLabels=d_labels,
                   yTicks=3,
                   widths=d_widths,
                   gap = 0.005,
                   scale=False
)

for i in range(len(A)):
    Ai = [x for x in A[i] if x>0]
    y = [x/2.0 for x in Ai]
    for j in range(len(Ai)):
        if j>0:
            yy = y[j]+np.sum(Ai[0:j])
        else:
            yy = y[j]        
        pl.text(0.5*i-0.02,yy-1.2,str(Ai[j]),fontsize=12,zorder=10)
    
# #Set general plot properties
# sns.set_style("white")
# sns.set_context({"figure.figsize": (24, 10)})

# for i in np.sort(list(usz))[::-1]:
#     y = [100-dc[u][i] for u in np.sort(dc.keys())]
#     bottom_plot=sns.barplot(x=np.asarray(densities)*100, y=y)
#    #                         color=scalarMap.to_rgba(i))
#     #y = (sbd[i+1]-sbd[i])/2.+sbd[i]scala
#     #for j in range(len(sbd.Density)):
#     #    pl.text(j-0.1,y[j],'1',fontsize=16,zorder=i)

# #Optional code - Make plot look nicer
sns.despine(left=True)
# #Set fonts to consistent 16pt size
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(12)
pl.subplots_adjust(bottom=0.2)
# plt.show()
pl.show()
import sys, os
sys.path.append('./tools/')
import traversal, bfutils, graphkit
import unknownrate as ur
from multiprocessing import Pool,Process, Queue, cpu_count, current_process
import functools
import zickle as zkl
import time, socket
import scipy
import clingo as cg

UMAX = 1
INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 1000 # stop traversing after growing equivalence class tothis size
REPEATS = 100
if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=60
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'mars':
    PNUM=12
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'hooke':
    PNUM=21
    PNUM=max((1,PNUM/INPNUM))
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count())
    PNUM=max((1,PNUM/INPNUM))
print 'processes: ',PNUM, INPNUM

def wrapper_rate_agnostic(fold, n=10, k=10):
    scipy.random.seed()
    l = {}
    while True:
        try:
            g = bfutils.ringmore(n,k) # random ring of given density
            gs= bfutils.call_undersamples(g)
            for u in range(1,min([len(gs),UMAX])):
                g2 = bfutils.undersample(g,u)
                print fold,': ',traversal.density(g),':',
                startTime = int(round(time.time() * 1000))
                s = ur.iteqclass(g2, verbose=False)
                endTime = int(round(time.time() * 1000))
                print len(s)
                l[u] = {'eq':s,'ms':endTime-startTime}
        except MemoryError:
            print 'memory error... retrying'
            continue
        break
    return {'gt':g,'solutions':l}

def killall(l):
    for e in l:
        e.join(timeout=0.001)
        if not e.is_alive():
            #print 'first result'
            for p in l:
                if p != e:
                    #print 'terminating ', p.name
                    p.terminate()
                    p.join()
                else:
                    p.join()
            return True
    return False

def fan_wrapper(fold,n=10,k=10):
    scipy.random.seed()
    curr_proc=current_process()
    curr_proc.daemon=False
    output = Queue()
    while True:
        try:
            g = bfutils.ringmore(n,k)
            gdens = traversal.density(g)
            g2 = bfutils.increment_u(g,g)
            #g2 = bfutils.undersample(g,2)
            def inside_wrapper():
                scipy.random.seed()
                try:
                    startTime = int(round(time.time() * 1000))
                    s = traversal.v2g22g1(g2, capsize=CAPSIZE)
                    #s = traversal.backtrack_more2(g2, rate=2, capsize=CAPSIZE)
                    endTime = int(round(time.time() * 1000))
                    print "{:2}: {:8} : {:4}  {:10} seconds".\
                        format(fold, round(gdens,3), len(s),
                               round((endTime-startTime)/1000.,3))
                    output.put({'gt':g,'eq':s,'ms':endTime-startTime})
                except MemoryError:
                    print 'memory error...'
		    raise
            pl = [Process(target=inside_wrapper) for x in range(INPNUM)]
            for e in pl: e.start()
            while True:
                if killall(pl): break
            r = output.get()
        except MemoryError:
            print 'memory error... retrying'
            for p in pl:
                p.terminate()
                p.join()
            continue
        break
    for p in pl: p.join()
    return r

densities = {6: [0.2, 0.25, 0.3, 0.35],
             8: [0.3],
             10:[0.1],# 0.15, 0.2, 0.25, 0.3],
             15:[0.25, 0.3],
             20:[0.1],# 0.15, 0.2, 0.25, 0.3],
             25:[0.1],
             30:[0.1],
             35:[0.1],
             40:[0.1],
             50:[0.05, 0.1],
             60:[0.05, 0.1]}

for nodes in [15]:
    z = {}
    pool=Pool(processes=PNUM)
    for dens in densities[nodes]:
        print "{:2}: {:8} : {:10}  {:10}".format('id', 'density', 'eq class', 'time')
        e = bfutils.dens2edgenum(dens, n=nodes)
        eqclasses = pool.map(functools.partial(fan_wrapper, n=nodes, k=e),
                             range(REPEATS))
        z[dens] = eqclasses
        zkl.save(z[dens],
                 socket.gethostname().split('.')[0]+\
                     '_nodes_'+str(nodes)+'_density_'+str(dens)+'_newp_.zkl')
        print ''
        print '----'
        print ''
    pool.close()
    pool.join()
    zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_newp_.zkl')
# This is a use-case of the tools in the tools directory. The example defines a graph and shows how to generate a figure that shows the graph at different undersampling rates. Running the file in python (python dbnplot.py) generates a figure in figures folder: shipfig.pdf


# system packages
import os, sys
import numpy as np
from random import random
sys.path.append('tools')
import zickle as zkl
# local packages
import dbn2latex as d2l
import bfutils as bfu
from bfutils import jason2graph

def ring(n):
    g = {}
    g['1'] = {'1': {(0,1)}, '2': {(0,1)}}
    for i in range(2,n):
        g[str(i)] = {str(i+1): {(0,1)}}
    g[str(n)] = {'1': {(0,1)}}
    return g

def listplot(fname, mname='JJ', stl='', width=5, R=2):
    l = zkl.load(fname)
    l = l[:17*10]
    y = min(width,len(l))
    x = np.int(np.ceil(len(l)/float(y)))
    d2l.matrix_list(l,y,x, R=R, w_gap=1, h_gap=2, mname=mname, stl=stl)


g = {'1': {'2': set([(0, 1)]), '7': set([(0, 1)])},
     '10': {'1': set([(0, 1)]), '5': set([(0, 1)]), '9': set([(0, 1)])},
     '2': {'3': set([(0, 1)]),
           '4': set([(0, 1)]),
           '6': set([(0, 1)]),
           '7': set([(0, 1)])},
     '3': {'4': set([(0, 1)])},
     '4': {'1': set([(0, 1)]), '4': set([(0, 1)]), '5': set([(0, 1)])},
     '5': {'10': set([(0, 1)]),
           '5': set([(0, 1)]),
           '6': set([(0, 1)]),
           '8': set([(0, 1)]),
           '9': set([(0, 1)])},
     '6': {'2': set([(0, 1)]), '7': set([(0, 1)])},
     '7': {'8': set([(0, 1)])},
     '8': {'4': set([(0, 1)]),
           '7': set([(0, 1)]),
           '8': set([(0, 1)]),
           '9': set([(0, 1)])},
     '9': {'1': set([(0, 1)]),
           '10': set([(0, 1)]),
           '2': set([(0, 1)]),
           '6': set([(0, 1)]),
           '7': set([(0, 1)])}}



# generation of the output
g = {'1': {'2':set([(0,1)])},
     '2': {'3':set([(0,1)]), '1':set([(0,1)])},
     '3': {'4':set([(0,1)])},
     '4': {'1':set([(0,1)])},
}


#d2l.matrix_unfold(l[0],2,1,R=5, w_gap=1, h_gap=2, mname='TT1')

n = 20
dens = 0.07

for i in range(10):
    print i
    g = bfu.ringmore(n,bfu.dens2edgenum(dens,n))
    gl = bfu.all_undersamples(g)
    zkl.save(gl, 'list1.zkl')

    # output file
    foo = open('figures/shipfig_figure.tex', 'wb')
    sys.stdout = foo
    listplot('list1.zkl', width=17, R=6)
    sys.stdout = sys.__stdout__              # remember to reset sys.stdout!
    foo.flush()
    foo.close()
    PPP = os.getcwd()
    os.chdir('figures')
    os.system('pdflatex --shell-escape shipfig.tex 2>&1 > /dev/null')
    os.system('mv shipfig.pdf /tmp/all_graphs_'+str(n)+'nodes_density_'+str(dens)+'_'+str(i)+'.pdf')
    os.chdir(PPP)
import signal
import pprint
import time, socket
import numpy as np
import scipy
import functools, itertools
import progressbar as pb
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))
import copy

from multiprocessing import Pool,Process, Queue, cpu_count, current_process
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed
import linear_model as lm
import traversal as trv
import bfutils as bfu
import graphkit as gk
import zickle as zkl
import pc
import pylab as plt
import unknownrate as ur
NOISE_STD = '0.1'
DEPTH=2
URATE=3
DIST='beta'
BURNIN=100
SAMPLESIZE=2000
PARALLEL=True
POSTFIX='_rasl_u'+str(URATE)
EST = 'svar'
INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 100 # stop traversing after growing equivalence class tothis size
REPEATS = 20
if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=30
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'mars':
    PNUM=21
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'saturn':
    PNUM=12
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'hooke':
    PNUM=22
    PNUM=max((1,PNUM/INPNUM))
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count())
    PNUM=max((1,PNUM/INPNUM))
print 'processes: ',PNUM, INPNUM

def timeout(func, args=(), kwargs={}, timeout_duration=1, default=None):
    import signal

    class TimeoutError(Exception):
        pass

    def handler(signum, frame):
        raise TimeoutError()

    # set the timeout handler
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout_duration)
    try:
        result = func(*args, **kwargs)
    except TimeoutError as exc:
        result = default
    finally:
        signal.alarm(0)

    return result

def hamming_neighbors(v, step):
    l = []
    for e in itertools.combinations(range(len(v)),step):
        b = copy.copy(v)
        for i in e: b[i] = int(not b[i])
        l.append(b)
    return l

def find_nearest_reachable(g2, max_depth=4):
	s = ur.liteqclass(g2, capsize=CAPSIZE, verbose=False)	
	if s: return s
	step = 1
	n = len(g2)
	v = bfu.g2vec(g2)
	while True:
		l = hamming_neighbors(v,step)
		pbar = ProgressBar(widgets=['neighbors checked @ step '+str(step)+': ', Percentage(), ' '], maxval=len(l)).start()
		c = 0
		for e in l:
			g = bfu.vec2g(e,n)
			if not gk.scc_unreachable(g):
				s = ur.liteqclass(g, capsize=CAPSIZE, verbose=False)
			else:
				s = set()
			if s: return s
			pbar.update(c)
			c += 1
		pbar.finish()
		if step > max_depth:
			return set()
		step += 1


def wrapper(fold,n=10,dens=0.1, urate=URATE):
    scipy.random.seed()
    rate = urate

    r = None
    s = set()
    counter = 0
    while not s:
        scipy.random.seed()
        sst = 0.9
        r = None
        while not r:
            r = lm.getAring(n, dens, sst, False, dist=DIST)
            print sst,
            sys.stdout.flush()
            if sst < 0.03:
                sst -= 0.001
            else:
                sst -= 0.01
            if sst < 0: sst = 0.02
        #pprint.pprint(r['transition'].round(2),width=200)
        #d = zkl.load('leibnitz_nodes_'+str(n)+'_OCE_model_.zkl')
        #r = d[dens][fold]
        g = r['graph']
        true_g2 = bfu.undersample(g, rate-1)
        data = lm.drawsamplesLG(r['transition'], samples=BURNIN+SAMPLESIZE*2,
                                nstd=np.double(NOISE_STD))
        data = data[:,BURNIN:]
        if np.max(data) > 1000.:
            pprint.pprint(r['transition'].round(2),width=200)
            #raise ValueError
        startTime = int(round(time.time() * 1000))
        if EST=='pc':
            g2 = pc.dpc(data[:,::rate], pval=0.0001)
        elif EST=='svar':
            g2 = lm.data2graph(data[:,::rate])
        if trv.density(g2) < 0.7:
            print gk.OCE(g2,true_g2)
            #s = examine_bidirected_flips(g2, depth=DEPTH)
            s = find_nearest_reachable(g2, max_depth=1)
            #s = trv.v2g22g1(g2, capsize=CAPSIZE, verbose=False)
            #s = trv.edge_backtrack2g1_directed(g2, capsize=CAPSIZE)
            #s = timeout(trv.v2g22g1,
            #s = timeout(trv.edge_backtrack2g1_directed,
            #            args=(g2,CAPSIZE),
            #            timeout_duration=1000, default=set())
            print 'o',
            sys.stdout.flush()
            if -1 in s: s=set()
        endTime = int(round(time.time() * 1000))
        #if counter > 3:
        #    print 'not found'
        #    return None
        counter += 1
    print ''
    oce = [gk.OCE(bfu.num2CG(x,n),g) for x in s]
    cum_oce = [sum(x['directed'])+sum(x['bidirected']) for x in oce]
    idx = np.argmin(cum_oce)
    print "{:2}: {:8} : {:4}  {:10} seconds".\
          format(fold, round(dens,3), cum_oce[idx],
                 round((endTime-startTime)/1000.,3))
    #np.set_printoptions(formatter={'float': lambda x: format(x, '6.3f')+", "})
    #pprint.pprint(r['transition'].round(2))
    #np.set_printoptions()

    return {'gt':r,
            'eq':s,
            'OCE':oce[idx],
            'tries_till_found': counter,
            'estimate': g2,
            'graphs_tried': counter,
            'strength':sst+0.01,
            'ms':endTime-startTime}


def wrapgen(fold,n=10,dens=0.1):
    scipy.random.seed()
    rate = 2

    s = set()
    sst = 0.06
    r = None
    while not r:
        r = timeout(lm.getAring, args=(n, dens, sst, False),
                    timeout_duration=3)
        print sst,
        if sst < 0.03:
            sst -= 0.002
        else:
            sst -= 0.01
        if sst < 0: break
    print 'model '+str(fold)+' found \n'+str(r['transition'].round(2))
    sys.stdout.flush()
    return r

densities = {5: [0.25, 0.3, 0.35],
			 6: [0.25, 0.3, 0.35],
             8: [.15, .2, 0.25, 0.3],
             10:[0.3],
             15:[0.1],
             20:[0.1],
             25:[0.1],
             30:[0.1],
             35:[0.1]}

wrp = wrapper

for nodes in [6]:
    z = {}
    pool=Pool(processes=PNUM)
    for dens in densities[nodes]:
        print "{:2}: {:8} : {:10}  {:10}".format('id', 'density', 'OCE', 'time')

        if PARALLEL:
            errors = pool.map(functools.partial(wrp, n=nodes,
                                                dens=dens),
                              range(REPEATS))
            print 'done'
        else:
            errors = []
            for i in range(REPEATS):
                errors.append(wrp(i,n=nodes,dens=dens))
        print 'computed'
        z[dens] = errors
        zkl.save(z[dens],
                 socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_samples_'+str(SAMPLESIZE)+'_density_'+str(dens)+'_noise_'+NOISE_STD+'_OCE_b_'+EST+'_'+DIST+POSTFIX+'.zkl')
        print ''
        print '----'
        print ''
    pool.close()
    pool.join()
    zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_samples_'+str(SAMPLESIZE)+'_noise_'+NOISE_STD+'_OCE_b_'+EST+'_'+DIST+POSTFIX+'.zkl')
import sys
sys.path.append('./tools/')
import seaborn as sb
import zickle as zkl
import graphkit as gk
import bfutils as bfu
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


#

uplimit = 0.25

#

#d05 = zkl.load('leibnitz_nodes_5_samples_2000_noise_0.1_OCE_b_svar_beta_rasl_more.zkl')
#d05 = zkl.load('oranos_nodes_5_samples_2000_noise_0.1_OCE_b_svar_beta_rasl_more.zkl')
#d05 = zkl.load('leibnitz_nodes_6_samples_2000_noise_0.1_OCE_b_svar_beta_rasl.zkl')
d05 = zkl.load('leibnitz_nodes_6_samples_2000_noise_0.1_OCE_b_svar_beta_rasl.zkl')

def estOE(d):
    gt= d['gt']['graph']
    gt=bfu.undersample(gt,1)
    e = gk.OCE(d['estimate'],gt)
    N = np.double(len(gk.edgelist(gt))) +\
        np.double(len(gk.bedgelist(gt)))
    return (e['directed'][0]+e['bidirected'][0])/N

def estCOE(d):
    gt= d['gt']['graph']
    gt=bfu.undersample(gt,1)
    e = gk.OCE(d['estimate'],gt)
    n = len(gt)
    N = np.double(n**2+(n-1)**2/2.0\
                  -len(gk.edgelist(gt))
                  -len(gk.bedgelist(gt)))
    return (e['directed'][1]+e['bidirected'][1])/N

d = d05#d08_01
density = np.sort(d.keys())
n = len(d[density[0]][0]['gt']['graph'])
OE = [[gk.oerror(x) for x in d[dd]] for dd in density]
COE = [[gk.cerror(x) for x in d[dd]] for dd in density]

eOE = [[estOE(x) for x in d[dd]] for dd in density]
eCOE = [[estCOE(x) for x in d[dd]] for dd in density]

samplesnum = 20
denscol = [0.25]*samplesnum+[0.30]*samplesnum+[0.35]*samplesnum
OE = OE[0]+OE[1]+OE[2]
eOE = eOE[0]+eOE[1]+eOE[2]
COE = COE[0]+COE[1]+COE[2]
eCOE = eCOE[0]+eCOE[1]+eCOE[2]


OE = pd.DataFrame(np.asarray([denscol+denscol, OE+eOE, pd.Categorical(['RASL']*samplesnum*3+['SVAR']*samplesnum*3)]).T, columns=['density', 'time', 'OE'])

COE = pd.DataFrame(np.asarray([denscol+denscol, COE+eCOE, pd.Categorical(['RASL']*samplesnum*3+['SVAR']*samplesnum*3)]).T, columns=['density', 'time', 'COE'])


shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

plt.figure(figsize=[14,6])

plt.subplot(121)
ax = sb.boxplot(x="density", y="time", hue="OE",
                data=OE,
                palette="Set2",
                linewidth=lwd,
                width=wds,
                fliersize=fliersz)
sb.stripplot(x="density", y="time", hue="OE",
                data=OE,
                palette='Set2',
                size=4, jitter=True, edgecolor="gray")
#ax.figure.get_axes()[0].set_yscale('log')
plt.xticks([0,1,2],('25%','30%','35%'))

plt.ylim([-0.02,uplimit])
plt.xlabel('density (% of '+str(n**2)+' total possible edges)')
plt.ylabel('Edge omission error')
plt.title(str(samplesnum)+' '+str(n)+'-node graphs per density',
          multialignment='center')
plt.legend(loc=0)


plt.subplot(122)
ax = sb.boxplot(x="density", y="time", hue="COE",
                data=COE,
                palette="Set2",
                linewidth=lwd,
                width=wds,
                fliersize=fliersz)
sb.stripplot(x="density", y="time", hue="COE",
                data=COE,
                palette='Set2',
                linewidth=lwd,
                size=4, jitter=True, edgecolor="gray")

#ax.figure.get_axes()[0].set_yscale('log')
plt.xticks([0,1,2],('25%','30%','35%'))

plt.ylim([-0.02,uplimit])
plt.xlabel('density (% of '+str(n**2)+' total possible edges)')
plt.ylabel('Edge comission error')
plt.title(str(samplesnum)+' '+str(n)+'-node graphs per density',
          multialignment='center')
plt.legend(loc=2)

sb.set_context('poster')
plt.savefig('/tmp/RASL_simulation.svgz',transparent=False, bbox_inches='tight', pad_inches=0)
plt.show()
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
import sys

sys.path.append('./tools/')
from pathtree import PathTree
from ortools.constraint_solver import pywrapcp
from matplotlib.cbook import flatten
from functools import wraps
import numpy as np
import bisect
from sortedcontainers import SortedDict
import ipdb


class SolutionNotFoundInTime(Exception):
    pass


def ptloopnum(pt):
    """
    Given a PathTree object returns the number of loops in it
    :param pt: PathTree object
    :return: number of loops (n)
    """

    def ptn(pt, n=0):
        for e in pt.loopset:
            if type(e) is int:
                n += 1
                continue
            n += ptn(e, n=1)
        return n

    return ptn(pt)


def ptnodenum(pt):
    """
    Given a PathTree object returns the number of latents that comprise it
    :param pt: PathTree object
    :return: number of nodes (n)
    """
    n = pt.preset - 1

    def ptn(pt, n=0):
        for e in pt.loopset:
            if type(e) is int:
                n += e - 1
                continue
            n += ptn(e, n=1)
        return n

    return n + ptn(pt)


def ptelement(pt, w):
    """
    An element generated by a PathTree with a given weight setting
    :param pt: PathTree
    :param w: a list of weights
    :return: an integer
    """
    n = pt.preset

    def sumloops(pt, w):
        n = 0
        ls = list(pt.loopset)
        for i in range(len(ls)):
            if type(ls[i]) is int:
                n += w[i] * ls[i]
                continue
            n += w[i][0] * ls[i].preset \
                 + min(1, w[i][0]) * sumloops(ls[i], w[i][1])
        return n

    return n + sumloops(pt, w)


def weights_pt(pt, weights):
    c = [0]

    def crawl(pt, w, c):
        wl = []
        for e in pt.loopset:
            if type(e) is int:
                wl.append(w[c[0]])
                c[0] += 1
                continue
            ww = w[c[0]]
            c[0] += 1
            wl.append([ww, crawl(e, w, c)])
        return wl

    return crawl(pt, weights, c)


def extraloops_pt(pt, loops):  # loops are tuples (loop, weight)
    c = [0]

    def crawl(pt, l, c):
        first = [l[c[0]]]
        wl = []
        for e in pt.loopset:
            c[0] += 1
            if type(e) is int:
                wl.append(l[c[0]])
                continue
            wl.append(crawl(e, l, c))
        return first + [wl]

    return crawl(pt, loops, c)


def ptelement_extraloop(pt, w, eloops):
    """
    An element generated by a PathTree with a given weight setting and extra loops on each level
    :param pt: PathTree
    :param w: a list of list of weights
    :param eloops: a list of tuples with lengths of extra loops and their weights
    :return: an integer
    """
    n = pt.preset + eloops[0][0] * eloops[0][1]

    def sumloops(pt, w, lps):
        ls = list(pt.loopset)
        n = 0
        for i in range(len(ls)):
            if type(ls[i]) is int:
                n += w[i] * ls[i] + min(1, w[i]) * lps[i][0] * lps[i][1]
                continue
            n += w[i][0] * ls[i].preset \
                 + min(1, w[i][0]) * (lps[i][0][0] * lps[i][0][1] + sumloops(ls[i], w[i][1], lps[i][1]))
        return n

    return n + sumloops(pt, w, eloops[1])


def isptelement_el(el, pt, w, eloops):
    return el == ptelement_extraloop(pt, w, eloops)


def isptsubset_el(elist, pt, w, eloops):
    for i in range(elist[-1]):
        if isptelement_el(i, pt, w, eloops):
            if not i in elist:
                return False
    return True


def isrightpt(el, elist, pt, w, eloops):
    for i in range(elist[-1]):
        if isptelement_el(i, pt, w, eloops):
            if not i in elist:
                return False
        if i == el and not isptelement_el(i, pt, w, eloops):
            return False
    return True


def ptelements(pt, seqlen=100, verbose=False, maxloop=100):
    """
    Generate first `seqlen` elements from a pathtree
    :param pt: a path tree object from pathtree.py
    :param seqlen: number of elements to generate in ascending order
    :param verbose: whether to print debugging information
    :return: a list of elements
    """
    solver = pywrapcp.Solver("pt-elements")

    # declare variables
    weights = []
    N = ptloopnum(pt)
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))

    # declare constraints
    # solver.Add()

    # run the solver
    solution = solver.Assignment()
    solution.Add(weights)
    db = solver.Phase(weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    num_solutions = 0
    els = set()
    while solver.NextSolution():
        w = [x.Value() for x in weights]
        num_solutions += 1
        els.add(ptelement(pt, w))
        if len(els) == seqlen:
            break
    solver.EndSearch()

    # output solutions
    if verbose:
        print "num_solutions:", num_solutions
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()

    return list(els)


def isptelement(pt, element, verbose=False, maxloop=100):
    """
    Check if an integer element is in the weight set represented by the path tree
    :param pt: a path tree object from pathtree.py
    :param element: an integer to check for presence in the weight
    :param verbose: whether to print debugging information
    :return: True or False
    """
    solver = pywrapcp.Solver("isptelement")

    # declare variables
    weights = []
    N = ptloopnum(pt)
    if not N:
        return element == pt.preset
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))

    wpt = weights_pt(pt, weights)

    # declare constraints
    solver.Add(element == ptelement(pt, wpt))

    # run the solver
    solution = solver.Assignment()
    solution.Add(weights)
    db = solver.Phase(weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    solution_exists = False
    while solver.NextSolution():
        solution_exists = True
        break
    solver.EndSearch()

    # output solutions
    if verbose:
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()

    return solution_exists


def loops_and_weights(solver, loops, weights):
    """
    Add constraints to solver that make sure loops are not generated if subtree is not active due to a zero weight upstream
    :param solver:
    :param loops:
    :param weights:
    :return:
    """

    def recurse(s, l, w):
        for ww, ll in zip(w, l):
            if type(ww) is list:
                for e in flatten(ll):
                    s.Add((ww[0] == 0) <= (e == 0))
                recurse(s, ll[1:], ww[1:])
            else:
                for e in flatten(ll):
                    s.Add((ww == 0) <= (e == 0))

    recurse(solver, loops[1], weights)


def eloops_simplify(eloops):
    l = []
    for e in eloops:
        if type(e) is list:
            l.append(eloops_simplify(e))
        else:
            l.append(int(e[0].Value()))
    return l


def ptaugmented(pt, eloops):
    def augment(pt, ls):
        pre = pt.preset
        loop = pt.loopset
        s = set()
        if ls[0]:
            s.add(ls[0])
        for l, el in zip(loop, ls[1]):
            if type(l) is int:
                if not el:
                    s.add(l)
                else:
                    s.add(PathTree({el}, pre=l))
                continue
            s.add(augment(l, el))

        return PathTree(s, pre=pre)

    t = augment(pt, eloops)

    return t


def ptsubset(pt, elist):
    for i in range(elist[-1]):
        if isptelement(pt, i) and not i in elist:
            return False
    return True


def smallest_pt(ptlist):
    if ptlist:
        idx = np.argsort(map(ptnodenum, ptlist))
        sol = ptlist[idx[0]]
    else:
        sol = None
    return sol


def pairprint(pt1, pt2, k=40):
    print np.c_[pt2seq(pt1, k), pt2seq(pt2, k)]


def etesteq(pt1, pt2, k=100):
    a1 = np.asarray(pt2seq(pt1, k))
    a2 = np.asarray(pt2seq(pt2, k))
    return np.sum(a1 - a2) == 0


def keeptreegrow(pt, e, seq, cutoff=10, cap=1000):
    t = None
    while t is None:
        t = growtree(pt, e, seq, cutoff=cutoff)
        cutoff += 10
        if cutoff > cap:
            raise SolutionNotFoundInTime("Cannot keep the tree growing")
    return t


def add_element(d, pt):
    """
    Add a PathTree to dictionary d such that it is either appended to the list or added anew
    Args:
        d: a dictionary
        pt: a PathTree

    Returns:

    """
    key = ptnodenum(pt)
    if key in d:
        d[key].append(pt)
    else:
        d[key] = pt


def del_element(d, pt, key=None):
    """
    Delete a PathTree from dictionary d such that it is either removed from the list or the list that only contains one element is removed
    Args:
        d: a dictionary
        pt: a PathTree

    Returns:

    """
    if key is None:
        key = ptnodenum(pt)
    if len(d[key]) == 1:
        del d[key]
    else:
        d[key].remove(pt)


def swap_elements(d, pt1, pt2, key=None):
    del_element(d, pt1, key=key)
    add_element(d, pt2)


def seq2pt(seq, verbose=False, cutoff=100):
    if not seq:
        return None

    pt = PathTree({}, pre=seq[0])
    pts = SortedDict()  # PathTrees
    pts[ptnodenum(pt)] = [pt]

    for e in seq[1:]:
        e_is_in = False
        for key in pts:
            for pt in pts[key]:
                if verbose:
                    print e
                try:
                    newpt = keeptreegrow(pt, e, seq, cutoff=cutoff)
                    swap_elements(pts, pt, newpt, key=key)
                    e_is_in = True
                    break
                except SolutionNotFoundInTime:
                    continue
        if not e_is_in:
            newpt = PathTree({}, pre=e)
            add_element(d, newpt)

    return pt


def growtree(pt, element, ref_elements, verbose=False, maxloop=100, cutoff=100):
    """
    Add a loop with the minimal length to a path tree to enable it to generate a given element and still be a subset of a given list
    :param pt: a path tree object from pathtree.py
    :param element: an integer to check for presence in the weight
    :param ref_elements: a (finite) list that should be a superset of numbers generated by the new path tree, for numbers smaller than tosubset[-1]
    :param verbose: whether to print debugging information
    :return: a PathTree augmented with a new loop
    """
    solver = pywrapcp.Solver("loop_an_element")

    # PathTree already can generate that number. Just to foolproof
    if isptelement(pt, element):
        return pt

    # declare variables
    weights = []  # weights denoting how many times a loop is active (marginalized)
    loops = []  # extra loops that can be potentially added
    lweights = []  # weights for the extra loops (marginalized out in the end)
    ltuples = []  # tuple list to hold loops and weights together

    N = ptloopnum(pt)  # number of loops in the PathTree
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))

    for i in range(N + 1):
        w = solver.IntVar(0, maxloop, "lw[%04i]" % i)
        l = solver.IntVar(0, maxloop, "l[%04i]" % i)
        lweights.append(w)  # loop related weight
        loops.append(l)
        ltuples.append((l, w))

    eloops = extraloops_pt(pt, ltuples)
    ws = weights_pt(pt, weights)

    # declare constraints
    solver.Add(solver.MemberCt(ptelement_extraloop(pt, ws, eloops), ref_elements))
    solver.Add(element == ptelement_extraloop(pt, ws, eloops))  # make sure the element can be generated
    solver.Add(solver.Count(loops, 0, len(loops) - 1))  # only one loop is on
    solver.Add(solver.Count(lweights, 0, len(lweights) - 1))  # only one loop is weighted
    for i in range(len(lweights)):
        solver.Add((lweights[i] == 0) <= (loops[i] == 0))  # if a loop has weight zero then it can't be active
        # solver.Add(lweights[i] >= loops[i])
    loops_and_weights(solver, eloops, ws)  # if a subtree is off (weight zero) no need to add loops

    # run the solver
    solution = solver.Assignment()
    solution.Add(loops)
    db = solver.Phase(loops + lweights + weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    numsol = 0
    pts = []
    while solver.NextSolution():
        # print numsol,
        new_pt = ptaugmented(pt, eloops_simplify(eloops))
        if verbose:
            print "trying PathTree: ", new_pt
        if ptsubset(new_pt, ref_elements):
            pts.append(new_pt)
            if verbose:
                print "OK PathTree: ", pts[-1]
        numsol += 1
        if numsol >= cutoff:
            break
    solver.EndSearch()

    # output solutions
    if verbose:
        print "solutions:", numsol
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()
        print "for ", element, "solutions found ", numsol

    return smallest_pt(pts)


def pt2seq(pt, num):
    if not pt.loopset:
        return [pt.preset]
    i = 0
    s = set()
    while len(s) < num:
        if isptelement(pt, i, maxloop=10 * num):
            s.add(i)
        i += 1
    l = list(s)
    l.sort()
    return l


def s2spt(s):  # convert edge set to pt
    ss = set()
    for e in s:
        if type(e) is int:
            ss.add(PathTree({0}, pre={e}))
            continue
        ss.add(e)
    return ss


def spt_elements(spt, num):
    """
    Generate numbers from a set of PathTrees
    :param spt: set of PathTrees
    :param num: number of elements (from the first) to generate
    :return: list of num numbers
    """
    i = 0
    s = set()
    while len(s) < num:
        if issptelement(spt, i):
            s.add(i)
        i += 1
    return list(s)


def issptelement(spt, element):
    a = False
    for pt in s2spt(spt):
        a = a or isptelement(pt, element)
    return a
import sys

sys.path.append('./tools/')


def osumnum(s, num):
    return set(num + x for x in s)


def osumset(s1, s2):
    s = set()
    for e in s1:
        s.update(osumnum(s2, e))
    return s


class PathTree:
    def __init__(self, lset, pre=0):
        self.loopset = lset
        self.preset = pre

    def __add__(self, other):
        if type(other) is int:
            return PathTree(self.loopset, pre=self.preset + other)
        if type(other) is set:
            return PathTree(self.loopset, pre=osumnum(other, self.preset))
        return PathTree(self.loopset.union(other.loopset), pre=self.preset + other.preset)

    def __radd__(self, other):
        return self.__add__(other)

    def __str__(self):
        s = '('
        comma = False
        if not self.preset == 0:
            s += str(self.preset) + ' '
            comma = True
        if not self.loopset == {0}:
            if comma:
                s += ', '
            s += '<' + ', '.join(map(str, self.loopset)) + '>'
        s += ')'
        return s
import sys

sys.path.append('./tools/')
import traversal, bfutils
import numpy as np
from ortools.constraint_solver import pywrapcp

U = 2 # undersampling rate: 1 means no undersampling
N = 10 # number of nodes
k = 25 # number of extra edges

solver = pywrapcp.Solver("MSL")

# generate a random graph and undersample
g = bfutils.ringmore(N,k)
gdens = traversal.density(g)
g2 = bfutils.undersample(g,U-1)

# undersampled edges
dedgeu = {}
bedgeu = {}
for i in range(N):
    for j in range(N):
        dedgeu[(i,j)] = 0
        bedgeu[(i,j)] = 0
        v = str(i+1)
        w = str(j+1)
        if w in g2[v]:
            if (0,1) in g2[v][w]: dedgeu[(i,j)] = 1
            if (2,0) in g2[v][w]: bedgeu[(i,j)] = 1
        

# declare variables
edges = []
for i in range(N):
    e = []
    for j in range(N):
        e.append(solver.IntVar(0,1,"%i -> %i" % (i,j)))
    edges.append(e)

# path constraint
def apath(i,j,u, e=edges):
    n = len(e)
    if u <= 1: return [e[i][k]*e[k][j] for k in range(n)]
    l = []
    for k in range(n):
        for z in range(n):
            l.extend(map(lambda x: e[i][k]*x*e[z][j], apath(k,z,u-1)))
    return l                

# directed path constraint
def dcons(i, j, u, e=edges, s=solver):
    return s.Sum(apath(i,j,u,e=e))
# bidirected edge constraint (a balanced trek constraint)
def bcons(i, j, u, e=edges, s=solver):
    n = len(e)
    l = [e[k][i]*e[k][j] for k in range(n)]
    for ui in range(1,u):
        for k in range(n):
            l.extend([x*y for x in apath(k,i,ui,e=e) for y in apath(k,j,ui,e=e)])
    return s.Sum(l)
    
# declare constraints
for i in range(N):
    for j in range(N):
        # directed edge constraints
        de = dcons(i,j,U-1) 
        if dedgeu[(i,j)] == 1:
            solver.Add( 0 < de )
        else:
            solver.Add( 0 == de)

# bidirected edge constraints
for i in range(N):
    for j in range(i,N):
        if j == i: continue        
        be = bcons(i,j,U-1) #solver.Sum([edges[k][i]*edges[k][j] for k in range(N)])
        if bedgeu[(i,j)] == 1:
            solver.Add( 0 < be )
        else:
            solver.Add( 0 == be)


# run the solver
solution = solver.Assignment()
solution.Add([edges[i][j] for i in range(N) for j in range(N)])
collector = solver.AllSolutionCollector(solution)
solver.Solve(solver.Phase([edges[i][j] for i in range(N) for j in range(N)],
                          solver.CHOOSE_FIRST_UNBOUND,
                          solver.ASSIGN_MIN_VALUE),
                          [collector])
num_solutions = collector.SolutionCount()

# output solutions
print "num_solutions:", num_solutions
print "failures:", solver.Failures()
print "branches:", solver.Branches()
print "WallTime:", solver.WallTime()

if num_solutions > 0 and num_solutions < 5:
    for s in range(num_solutions):
        qval = [collector.Value(s, edges[i][j]) for i in range(N) for j in range(N)]
        for i in range(len(qval)):
            if qval[i]:
                e = np.unravel_index(i,[N,N])
                print e[0],"->",e[1]
        print
        print
import sys, os
import sys, os

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import graphkit as gk
import numpy as np
from random import shuffle
import ecj
from copy import deepcopy
from pathtree import PathTree, osumset
import ipdb


def g2lg(g):
    """
    Convert a data structure encoding the MSL-type graph into a structure encoding latents graph
    :return: a graph with integer indices and sets as weights

    Args:
        g (MSL-graph): 
    """
    edge_type = {(0, 1): 1, (2, 0): 2}
    edge_weight = {(0, 1): 1, (2, 0): 0}
    lg = {int(e): {int(c): {edge_type[w]: {edge_weight[w]}
                            for w in g[e][c]}
                   for c in g[e]}
          for e in g}
    return fix_selfloops(lg)


def fix_selfloops(g):
    for v in g:
        if v in g[v]:
            g[v][v] = {1: {PathTree({1})}}
    return g


def gtranspose(G):  # Transpose (rev. edges of) G
    GT = {u: {} for u in G}
    for u in G:
        for v in G[u]:
            if 1 in G[u][v]:
                GT[v][u] = {1: {1}}  # Add all reverse edges
    return GT


def parents(g, N):  # an inefficient way to do it
    t = gtranspose(g)  # Transpose the graph
    return {n: g[n][N][1]
            for n in t[N] if n != N}  # Return all children


def children(g, N):
    return {n: g[N][n][1]
            for n in g[N] if n != N}


def remove_node(g, N):
    del g[N]
    for V in g:
        if N in g[V]:
            del g[V][N]


def iterate_ws(ws):
    starts = []
    for e in ws:
        if type(e) is int:
            starts.append(e)
            continue
        if type(e) is set:
            starts.extend(list(e))
            continue
        starts.append(e.preset)
        starts.append([e.preset, iterate_ws(e.loopset)])
    return starts


def iterate_pt(pt):  # iterate over a path tree
    starts = [pt.preset]
    for e in pt.loopset:
        if type(e) is int:
            starts.append(e)
            continue
        if type(e) is set:
            starts.extend(list(e))
            continue
        starts.append(e.preset)
        starts.append([e.preset, iterate_ws(e.loopset)])
    return starts


def merge_weightsets(ab, ah, hb, hh):
    ws = osumset(ah, hb)
    if hh:
        ws = osumset(ws, hh)
    if not type(ws) is set:
        ws = {ws}
    return ws.union(ab)


def hide_node(g, H):
    """
    Removes a node from a graph taking care of the edges and weights as if the node was not observed
    :param g: input graph
    :param H: node to hide
    :return: the graph
    """

    gg = deepcopy(g)

    if not H in g:
        raise KeyError
    ch = children(g, H)
    pa = parents(g, H)
    if H in g[H]:
        sl = g[H][H][1]
    else:
        sl = set()
    remove_node(gg, H)

    for p in pa:
        for c in ch:
            if c in gg[p]:
                ab = gg[p][c][1]  # self loop
            else:
                ab = set()
            w = merge_weightsets(ab, pa[p], ch[c], sl)  # new weight set
            if c == p:  # a new loop is forming
                w = {PathTree(w)}
            gg[p][c] = {1: w}

    return gg


def degrees(nodes, g):
    return [len(parents(g, v)) + len(children(g, v)) for v in nodes]


def sortbydegree(nodes, g):
    idx = np.argsort(degrees(nodes, g))
    return list(np.asarray(nodes)[idx])


def hide_nodes(g, nodelist, dosort=True):
    nodeset = set()  # make sure not to delete a node twice
    if dosort: nodelist = sortbydegree(nodelist, g)
    gg = deepcopy(g)
    for n in nodelist:
        if n in nodeset: continue
        gg = hide_node(gg, n)
        nodeset.add(n)
    return gg


def hide_random(g, ratio):
    """
    Hire random modes in the `ratio` proportion from graph g
    :param g: input graph
    :param ratio: what percentage of nodes to hide
    :return: the graph with hidden variables
    """
    nodes = g.keys()
    shuffle(nodes)
    return hide_nodes(g, nodes[:int(len(g) * ratio)])


def print_ws(ws):
    print '{',
    for e in ws:
        print e, ', ',
    print '}'


def test_osumnum():
    assert osumnum(set(range(5)), 1) == set(range(1, 5 + 1))


def testcase(n):
    g1 = {1: {2: {1: {1}}, 4: {1: {1}}},
          2: {3: {1: {1}}, 7: {1: {1}}},
          3: {},
          4: {5: {1: {1}}},
          5: {3: {1: {1}}, 6: {1: {1}}},
          6: {5: {1: {1}}},
          7: {8: {1: {1}}},
          8: {2: {1: {1}}}}

    g2 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 6: {1: {1}}},
          3: {4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {},
          6: {7: {1: {1}}},
          7: {3: {1: {1}}}}

    g3 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 6: {1: {1}}},
          3: {4: {1: {1}}, 8: {1: {1}}},
          4: {5: {1: {1}}},
          5: {},
          6: {7: {1: {1}}},
          7: {2: {1: {1}}},
          8: {9: {1: {1}}},
          9: {10: {1: {1}}},
          10: {11: {1: {1}}},
          11: {3: {1: {1}}}}

    g4 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {}}

    g5 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    g6 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}, 6: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    g7 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}, 6: {1: {1}}},
          7: {8: {1: {1}}, 7: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    g8 = {1: {2: {1: {1}}, 5: {1: {1}}},
          2: {3: {1: {1}}, 2: {1: {1}}},
          3: {4: {1: {1}}},
          4: {8: {1: {1}}},
          5: {6: {1: {1}}},
          6: {7: {1: {1}}},
          7: {4: {1: {1}}},
          8: {9: {1: {1}}},
          9: {9: {1: {1}}, 10: {1: {1}}},
          10: {}}

    cases = [g1, g2, g3, g4, g5, g6, g7, g8]

    return fix_selfloops(cases[n])
import graph_tool as gt
from graph_tool import draw as gtd
import numpy as np

def lg2gt(g):
    gr = gt.Graph()
    vlabel = gr.new_vertex_property("string")
    verts = {}
    edges = {}
    for v in g:
        verts[v] = gr.add_vertex()
        vlabel[verts[v]] = str(v)
    gr.vertex_properties["label"] = vlabel
    for v in g:
        for w in g[v]:
            edges[(v,w)] = gr.add_edge(verts[v], verts[w])
    return gr

def plotg(g, layout='sfdp', pos=True):
    gg = lg2gt(g)
    if not pos:
        if layout=='fr':
            pos = gtd.fruchterman_reingold_layout(gg)
        else:
            pos = gtd.sfdp_layout(gg)
    else:
        pos = gg.new_vertex_property("vector<double>")
        n = gg.num_vertices()
        s = 2.0*np.pi/n
        for v in range(gg.num_vertices()):
            idx = int(gg.vertex_properties['label'][gg.vertex(v)]) - 1
            pos[gg.vertex(v)] = (n * np.cos(s * idx),
                                 n * np.sin(s * idx))

    gtd.graph_draw(gg, pos,
               vertex_text=gg.vertex_properties['label'],
               vertex_font_size=32,
               edge_pen_width=1,
               edge_marker_size=15,
               vertex_pen_width=1,
               vertex_fill_color=[0.62109375,
                                  0.875     ,
                                  0.23828125,
                                  1])
# tools to construct (random) graphs
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import ecj
import bfutils as bfu
import traversal as trv

import random as std_random
import numpy as np
import scipy
import networkx as nx
import jgraph as igr
from numpy.random import randint
from comparison import nx2graph

def edgelist(g): # directed
    '''
    return a list of tuples for edges of g
    '''
    l = []
    for n in g:
        l.extend([(n,e) for e in g[n] if (0,1) in g[n][e]])
    return l

def edgenumber(g):
    return sum([sum([len(g[y][x]) for x in g[y]]) for y in g])

def iedgelist(g): # directed iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    for v in g:
        for w in g[v]:
            if (0,1) in g[v][w]: yield (v,w)
def inedgelist(g): # missing directed iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    n = len(g)
    for v in g:
        for i in xrange(1,n+1):
            w = str(i)
            if not w in g[v]:
                yield (v,w)
            elif not (0,1) in g[v][w]:
                yield (v,w)
def ibedgelist(g): # directed iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    for v in g:
       for w in g[v]:
           if (2,0) in g[v][w]: yield (v,w)
def inbedgelist(g): # missing bidirected iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    for v in g:
       for w in g:
           if v!=w:
               if not w in g[v]:
                   yield (v,w)
               elif not (2,0) in g[v][w]:
                   yield (v,w)

def bedgelist(g): # bidirected edge list with flips
    l = []
    for n in g:
        l.extend([tuple(sorted((n,e))) for e in g[n] if (2,0) in g[n][e]])
    l = list(set(l))
    l = l + map(lambda x: (x[1],x[0]), l)
    return l

def rnd_edges(n):
    """ generate a random uniformly distributed mask
    """
    rnum = std_random.getrandbits(n**2)
    l = list(bin(rnum)[2:])
    l = ['0' for i in range(0,n**2 - len(l))] + l
    return l

def list2dbn(l):
    """ convert list of edge presences/absences (0,1) to a DBN graph
    """
    n = scipy.sqrt(len(l))
    l = scipy.reshape(map(int, l), [n, n])
    G = ecj.adj2DBN(l)
    return G

def list2CG(l):
    """ convert list of edge presences/absences (0,1) to a compressed
    graph (CG) representation
    """
    n = scipy.sqrt(len(l))
    l = scipy.reshape(map(int, l), [n, n])
    G = bfu.adj2graph(l)
    return G

def rnd_dbn(n): return list2dbn(rnd_edges(n))
def rnd_cg(n):  return list2CG(rnd_edges(n))

def rnd_adj(n, maxindegree=5):
    l = scipy.zeros([n,n])
    for u in range(0,n):
        cap = scipy.random.randint(min([n,maxindegree+1]))
        idx = scipy.random.randint(n,size=cap)
        l[u, idx] = 1
    return l

def sp_rnd_edges(n, maxindegree=5):
    '''
    a sparse set of random edges
    '''
    l = rnd_adj(n, maxindegree=maxdegree)
    return scipy.reshape(l, n**2)

def sp_rnd_dbn(n, maxindegree=3):
    '''
    a sparse random DBN graph
    '''
    l = sp_rnd_edges(n)
    return list2dbn(l)

def emptyG(n):
    A = [[0 for j in range(n)] for i in range(n)]
    return bfu.adj2graph(np.asarray(A))

def fullG(n):
    A = [[1 for j in range(n)] for i in range(n)]
    return bfu.adj2graph(np.asarray(A))


def CG2uCG(cg):
    """
    convert to an undirected graph
    """
    G = {}
    for u in cg:
        G[u] = cg[u].copy()
    for u in cg:
        for v in cg[u]:
            G[v][u] = cg[u][v]
    return G

def connected(cg):
    n = len(cg)
    return sum(1 for _ in ecj.traverse(CG2uCG(cg),'1')) == n

def sp_rnd_CG(n, maxindegree=3, force_connected=False):
    l = sp_rnd_edges(n, maxindegree=maxindegree)
    cg = list2CG(l)
    if force_connected:
        while not connected(cg):
            cg = list2CG(sp_rnd_edges(n, maxindegree=maxindegree))
    return cg

def CG2adj(G):
    n = len(G)
    A = [[0 for i in range(0,n)] for j in range(0,n)]
    for v in G:
        if G[v]:
            directed = [w for w in G[v] if (0,1) in G[v][w]]
            for w in directed:
                A[int(w)-1][int(v)-1] = 1
    A = np.double(np.asarray(A))
    return A

def g2ig(g):
    """
    Converts our graph represenataion to an igraph for plotting
    """
    t = scipy.where(CG2adj(g)==1)
    l = zip(t[0],t[1])
    ig = igr.Graph(l,directed=True)
    ig.vs["name"] = scipy.sort([u for u in g])
    ig.vs["label"] = ig.vs["name"]
    return ig

def superclique(n):
    g = {}
    for i in range(n):
        g[str(i+1)] = {str(j+1):set([(0,1),(2,0)])
                       for j in range(n) if j!=i}
        g[str(i+1)][str(i+1)] = set([(0,1)])
    return g

def complement(g):
    n = len(g)
    sq = superclique(n)
    for v in g:
        for w in g[v]:
            sq[v][w].difference_update(g[v][w])
            if not sq[v][w]: sq[v].pop(w)
    return sq

def gtranspose(G):                      # Transpose (rev. edges of) G
    GT = {u:{} for u in G}
    for u in G:
        for v in G[u]:
            if (0,1) in G[u][v]:
                GT[v][u] = set([(0,1)])        # Add all reverse edges
    return GT

def scale_free(n, alpha=0.7, beta=0.25,
               delta_in=0.2, delta_out=0.2):
    g = nx.scale_free_graph(n, alpha=alpha,
                            beta=beta,
                            delta_in=delta_in, delta_out=delta_out)
    g = nx2graph(g)
    g = gtranspose(g)
    addAring(g)
    return g

def ring(n, permute=False):
    g = {}
    names = [str(x+1) for x in range(n)]
    if permute: names = np.random.permutation(names) 
    for i in range(n-1):
        g[names[i]] = {names[i+1]: set([(0,1)])}
    g[names[n-1]] = {names[0]: set([(0,1)])}
    return g

def addAring(g):
    for i in range(1,len(g)):
        if str(i+1) in g[str(i)]:
            g[str(i)][str(i+1)].add((0,1))
        else:
            g[str(i)][str(i+1)] = set([(0,1)])
    if '1' in g[str(len(g))]:
        g[str(len(g))]['1'].add((0,1))
    else:
        g[str(len(g))]['1'] = set([(0,1)])

def upairs(n,k):
    '''
    n unique nonsequential pairs
    '''
    s = set()
    for p in randint(n, size=(3*k, 2)):
        if p[1]-p[0] == 1: continue
        s.add(tuple(p))
    l = [e for e in s]
    return l[:k]

def ringarcs(g,n):
    for edge in upairs(len(g),n):
        g[str(edge[0]+1)][str(edge[1]+1)] = set([(0,1)])
    return g
def ringmore(n,m, permute=False):
    return ringarcs(ring(n,permute=permute),m)

def digonly(H):
    """returns a subgraph of H contatining all directed edges of H

    Arguments:
    - `H`: undersampled graph
    """
    g = {n:{} for n in H}
    for v in g:
        g[v] = {w:set([(0,1)]) for w in H[v] if not H[v][w] == set([(2,0)])}
    return g

# Justin's ternary representation: 1 = directed edge; 2 = bidirected; 3 = both
def justin2graph(g):
    r = {}
    d = {1: set([(0,1)]),
         2: set([(2,0)]),
         3: set([(0,1),(2,0)]) }
    for head in g:
        r[head] = {}
        for tail in g[head]:
            r[head][tail] = d[g[head][tail]]
    return r

def graph2justin(g):
    r = {}
    for head in g:
        r[head] = {}
        for tail in g[head]:
            if g[head][tail] == set([(0,1)]):
                r[head][tail] = 1
            elif g[head][tail] == set([(2,0)]):
                r[head][tail] = 2
            elif g[head][tail] == set([(0,1),(2,0)]):
                r[head][tail] = 3
    return r

def OCE(g1,g2):
    '''
    omission/commision error of g1 referenced to g2
    '''
    s1 = set(edgelist(g1))
    s2 = set(edgelist(g2))
    omitted = len(s2 - s1)
    comitted = len(s1 - s2)

    s1 = set(bedgelist(g1))
    s2 = set(bedgelist(g2))
    bomitted = len(s2 - s1)
    bcomitted = len(s1 - s2)

    return {'directed': (omitted, comitted),
            'bidirected': (bomitted, bcomitted)}

def clean_leaf_nodes(g):
    for v in g: g[v] = {w:g[v][w] for w in g[v] if g[v][w]}

def cerror(d):
    return d['OCE']['directed'][1]/np.double(len(d['gt']['graph'])**2-len(edgelist(d['gt']['graph'])))

def oerror(d):
    return d['OCE']['directed'][0]/np.double(len(edgelist(d['gt']['graph'])))

def bidirected_no_fork(g):
    be = bedgelist(g)
    T = gtranspose(g)
    for e in be:
        if not set(T[e[0]].keys())&set(T[e[1]].keys()):
            return True
    return False

def fork_mismatch(g):
    be = bedgelist(g)
    benum = len(be)/2
    forknum = 0
    for v in g:
        fn = len([w for w in g[v] if (0,1) in g[v][w]])
        forknum += fn*(fn-1)/2.
    if benum < len(g)*(len(g)-1)/2.:
        return (forknum-benum) > benum
    else:
        return False

def no_parents(g):
    T = gtranspose(g)
    for n in T:
        if not T[n]: return True
    return False

def no_children(g):
    for n in g:
        if not g[n]: return True
    return False

def scc_unreachable(g):
    if bidirected_no_fork(g): return True
    if no_parents(g): return True
    if no_children(g): return True
    return False

# unlike functions from traversal package these do no checking
def addanedge(g,e): g[e[0]][e[1]] =  set([(0,1)])
def delanedge(g,e): g[e[0]].pop(e[1], None)
def addedges(g,es):
    for e in es: addanedge(g,e)
def deledges(g,es):
    for e in es: delanedge(g,e)

def checkequality(H,G_test, au = None):
    if not au:
        allundersamples = bfu.call_undersamples(G_test)
    else:
        allundersamples = au
    for graph in allundersamples:
        if graph == H: return True
    return False

def isdedgesubset(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                if not (0,1) in g2[n][h]:
                    return False
            else:
                    return False
    return True

def isedgesubset(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                #if not (0,1) in g2[n][h]:
                if not g2star[n][h].issubset(g2[n][h]):
                    return False
            else:
                    return False
    return True

def isedgesubset_(g,H):
    '''
    check if g edges are a subset of those of H
    '''
    for e in inbedgelist(H):
        if e[1] in g[e[0]] and (2,0) in g[e[0]][e[1]]: return False
    for e in inedgelist(H):
        if e[1] in g[e[0]] and (0,1) in g[e[0]][e[1]]: return False
    return True

def checkconflict(H,G_test, au = None):
    if not au:
        allundersamples = bfu.call_undersamples(G_test)
    else:
        allundersamples = au
    for graph in allundersamples:
        if isedgesubset(graph,H): return False
    return True

def checkconflict_(Hnum, G_test, au = None):
    if not au:
        allundersamples = bfu.call_undersamples(G_test)
    else:
        allundersamples = au
    #Hnum = bfu.ug2num(H)
    for graph in allundersamples:
        gnum = bfu.ug2num(graph)
        if gnum[0]&Hnum[0] == gnum[0] and gnum[1]&Hnum[1] == gnum[1]:
            return False
    return True
import zickle as zkl
alloops = zkl.load('allloops.zkl')
import scipy
from numpy.lib.arraysetops import intersect1d
from itertools import combinations
from testgraphs import *

def walk(G, s, S=set()):
    P, Q = dict(), set()
    P[s] = None
    Q.add(s)
    while Q:
        u = Q.pop()
        for v in G[u].difference(P,S):
            Q.add(v)
            P[v] = u
    return P

def traverse(G, s, qtype=set):
    S, Q = set(), qtype()
    Q.add(s)
    while Q:
        u = Q.pop()
        if u in S: continue
        S.add(u)
        for v in G[u]:
            Q.add(v)
        yield u

def dfs_topsort(G):
    S, res = set(), []
    def recurse(u):
        if u in S: return
        S.add(u)
        for v in G[u]:
            recurse(v)
        res.append(u)
    for u in G:
        recurse(u)
    res.reverse()
    return res

def tr(G):                      # Transpose (rev. edges of) G
    GT = {}
    for u in G: GT[u] = set()   # Get all the nodes in there
    for u in G:
        for v in G[u]:
            GT[v].add(u)        # Add all reverse edges
    return GT

def scc(G):                   # Kosaraju's algorithm
    GT = tr(G)                # Get the transposed graph
    sccs, seen = [], set()
    for u in dfs_topsort(G):   # DFS starting points
        if u in seen: continue # Ignore covered nodes
        C = walk(GT, u, seen)  # Don't go "backward" (seen)
        seen.update(C)         # We've now seen C
        sccs.append(C)         # Another SCC found
    return sccs


TT = {
    'A': ['B', 'C'],
    'B': ['D','E'],
    'C': [],
    'D': [],
    'E': []
}
def bfs_print_tree(tree,r):
    """
    A modified single list solution
    """
    Q = []
    idx = 1
    print str(idx)+': '+r
    Q.extend(tree[r])
    while Q:
        idx += 1
        print str(idx)+':',
        for u in range(0,len(Q)):
            e = Q.pop(0)
            print e,
            Q.extend(tree[e])
        print ''

def bfs_dict(tree,r):
    """
    Bob's suggested dictionary based solution
    """
    D = {}
    idx = 1
    D[idx] = [r]
    while D[idx]:
        idx += 1
        D[idx] = []
        for u in D[idx-1]: D[idx].extend(tree[u])
    D.pop(idx) # the last dictionary element is empty - must go
    for idx in D: print str(idx)+': '+' '.join(D[idx])


def cloneBfree(G):
    D = {}
    for v in G:
        D[v] = {}
        for u in G[v]:            
            if not (len(G[v][u].intersection(set([(0,1)]))) == 0):
                D[v][u] = set([(0,1)])
    return D

def clrbi(G):    
    for v in G:
        d = []
        for u in G[v]:
            try:
                G[v][u].remove((edge_type['bidirected'],0))
                if len(G[v][u]) == 0:
                    d.append(u)
            except KeyError:
                pass
        for e in d:
            G[v].pop(e)

def ecj(G,s,sccs=set()):           # elementary circuit by Johnson
    blocked = {v:False for v in G} # unblock all
    B = {v:[] for v in G}
    stack = []
    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            B[u].remove(w)
            if blocked[w]: unblock(w)
    def circuit(v, stack):
        f = False
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                f = True
                #print stack
                sccs.add(len(stack))
            elif not blocked[u]:
                if circuit(u, stack):
                    f = True
        if f:
            unblock(v)
        else:
            for w in G[v]:
                if v not in B[w]:
                    B[w].append(v)
        stack.pop()
        return f
    circuit(s,stack)
    return sccs

def ecj_loops(G,s,sl=set()):           # elementary circuit by Johnson
    blocked = {v:False for v in G} # unblock all
    B = {v:[] for v in G}
    stack = []
    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            B[u].remove(w)
            if blocked[w]: unblock(w)
    def circuit(v, stack):
        f = False
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                f = True
                #print scipy.sort(stack)
                sl.add(tuple(scipy.sort(stack)))
            elif not blocked[u]:
                if circuit(u, stack):
                    f = True
        if f:
            unblock(v)
        else:
            for w in G[v]:
                if v not in B[w]:
                    B[w].append(v)
        stack.pop()
        return f
    circuit(s,stack)
    return sl

def gcd(a, b):
    while b != 0: a, b = b, a%b
    return a

def listgcd(l):
    if len(l)>0:
        return gcd(l[0],listgcd(l[1:]))
    else:
        return 0
def lcm(a, b): return a*b/gcd(a,b)
def chmatch(n,m,delta):
    m,n = scipy.sort([n,m])
    sq = scipy.mod(range(n,lcm(n,m)+1,n),m)
    return scipy.mod(delta,m) in sq

def reachable(s, G, g):
    S, Q = set(), []
    Q.append(s)
    while Q:
        u = Q.pop()
        if u in S: continue
        if g in G[u]: return True
        S.add(u)
        Q.extend(G[u])
    return False

def allpaths(G, s, g, S=[]):
    if S is None: S = []
    S.append(s)
    if s == g:
        print S
    else:
        for u in G[s]:
            if u in S: continue
            allpaths(G,u,g,S)
    S.remove(s)

def lz_ecj(G,s,sccs=set()):           # elementary circuit by Johnson
    blocked = {v:False for v in G}    # unblock all
    B = {v:set() for v in G}
    stack = []
    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            if blocked[w]: unblock(w)
        B[u].clear()
    def circuit(v, stack):
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                print 'bottom'
                unblock(v)
                yield len(stack)
            elif not blocked[u]:
                print 'recurse'
                for x in circuit(u, stack):
                    unblock(v)
                    yield x
            else:
                print 'unmet'
                for w in G[v]: B[w].add(v)
        stack.pop()
#    circuit(s,stack)
    for v in circuit(s,stack): yield v

def iterate_allpaths(G, s, g, d=0, S=[], c=True):
    if S is None: S = []
    S.append(s)
    d += 1
    if s == g:
        if c:
            yield d-1
        else:
            yield list(S)
    else:
        for u in G[s]:
            if u in S: continue
            for v in iterate_allpaths(G,u,g,d,S,c):
                yield v
    S.remove(s)

def iddfs(G, s, g): # iterative depening DFS paths
    yielded = set()
    def recurse(G, s, g, d, S=None):
        if s not in yielded:
            yielded.add(s)
        if d == 0: return
        if S is None: S = []
        S.append(s)
        if s == g:
            yield list(S)
        else:
            for u in G[s]:
                if u in S: continue
                for v in recurse(G, u, g, d-1, S):
                    yield v
        S.remove(s)
    n = len(G)
    for d in range(n):
        #if len(yielded) == n: break
        for u in recurse(G, s, g, d):
            yield u

def reached_at_step(G, s, d):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    yielded = set()
    def recurse(G, s, d, B=None):
        if d == 0:
            if s not in yielded: # this avoids yielding duplicates
                yielded.add(s)
                yield s
            return
        if B is None: B = [] # black - backed out of this path
        for u in G[s]:
            #if u in B: continue
            if G[s][u] == (edge_type['bidirected'],0): continue
            for v in recurse(G, u, d-1, B):
                yield v
        B.append(s)
    for u in recurse(G, s, d):
            yield u

def d_trek(h,G,a,b,d):
    """
    Does there exist a trek with head h connecting a and b in d steps.
    """
    return set([a,b]).issubset(reached_at_step(G,h,d))

def d_biegde(G,a,b,d):
    """
    Do  a and  b  become connected  by  a bidirectional  edge after  d
    undersamples
    """
    for i in range(1,d+1):
        for u in G:
            if d_trek(u,G,a,b,i):
                return True
    return False

def undersample(G,d,bid=True):
    """
    
    """
    N = {}
    for u in G: 
        N.update({u:{v:set([(0,1)]) for v in reached_at_step(G,u,d+1)}})
    if bid:
        items = G.keys()
        for i in range(len(items)):
            for j in range(i+1, len(items)):
                u,v = items[i],items[j]
                if d_biegde(G,u,v,d):
                    try:
                        N[u][v].add((edge_type['bidirected'],0))
                    except KeyError:
                        N[u].update({v:set([(edge_type['bidirected'],0)])})
                    try:
                        N[v][u].add((edge_type['bidirected'],0))
                    except KeyError:
                        N[v].update({u:set([(edge_type['bidirected'],0)])})
    return N

import itertools as itt
def exist_equal_paths(h,G,a,b):
    Sa, Sb, Dff = set(), set(), set()
    ag=iterate_allpaths(G,h,a,0,[],True)
    bg=iterate_allpaths(G,h,b,0,[],True)
    for v in itt.izip_longest(ag,bg):
        print v
        Sa.add(v[0])
        Sb.add(v[1])
        if v[0] in Sb or v[1] in Sa: return True
    return False

# checks if there  exist exact length paths from the  head node to the
# nodes at question  by iterative deepining to avoid  oing through all
# paths
def iexist_equal_paths(h,G,a,b):
    Sa, Sb = set(), set()
    Pa, Pb = [], []
    ag=iddfs(G,h,a)
    bg=iddfs(G,h,b)
    for v in itt.izip(ag,bg):
        print v
        Sa.add(len(v[0])); Pa.append(v[0])
        Sb.add(len(v[1])); Pb.append(v[1])
        if len(v[0]) in Sb or len(v[1]) in Sa: return True
    return False

# check  if two  unequal  length  paths can  be  compensated by  their
# elementary cycles

def has_unit_cycle(G,path):
    for v in path:
        if v in G[v]: return True
    return False
def ecj_compat(G,p1,p2):
    n = len(p1)
    m = len(p2)
    p2, p1 = [[p1,p2][i] for i in scipy.argsort([n,m])]
    m,n = scipy.sort([n,m])
    delta = n - m
    if not delta: return True # redundant check for 0
    if has_unit_cycle(G,p2): return True # equivalent
    # if the shorter path does not have cycles they are not compatible
    # if the  longer path does not  have cycles: check  if the shorter
    #                                            path    has    cycles
    #                                            divisible by delta

    # otherwise start checking
    print p1, p2, n, m, delta


def wc(n):
    n = n*3
    a={str(v):set([str(v+1),str(v+2),str(v+3)]) for v in range(1,n,3)}
    b={str(v):set([str(v+2)]) for v in range(2,n,3)}
    c={str(v):set([str(v+1)]) for v in range(3,n,3)}
    a.update(b)
    a.update(c)
    a.update({str(n):set()})
    return a

# Frobenius number from here: http://cgi.gladman.plus.com/wp/?page_id=563 
def residue_table(a):
  n = [0] + [None] * (a[0] - 1)
  for i in range(1, len(a)):
    d = gcd(a[0], a[i])
    for r in range(d):
      try:
        nn = min(n[q] for q in range(r, a[0], d) if n[q] != None)
      except:
        continue
      if nn != None: 
        for c in range(a[0] // d):
          nn += a[i]
          p = nn % a[0]
          nn = min(nn, n[p]) if n[p] != None else nn
          n[p] = nn
  return n
 
def frobenius_number(a):
  return max(residue_table(sorted(a))) - min(a)
 
def isSclique(G):
    n = len(G)
    for v in G:
        if sum([(0,1) in G[v][w] for w in G[v]]) < n: return False
        if sum([(2,0) in G[v][w] for w in G[v]]) < n-1: return False
    return True

# Jianyu does not use bidirected edges
def isJclique(G):
    return (sum([len(G[w].keys()) for w in G]) == len(G)**2)

def directed_inc(G,D):
    G_un = {}
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in [el for el in D[v] if (0,1) in D[v][el]]:
            for e in G[w]:
                G_un[v][e] = set([(0,1)])
    return G_un
def bidirected_inc(G,D):
    # bidirected edges
    for w in D:
        # transfer old bidirected edges
        l = [e for e in D[w] if (2,0) in D[w][e]]
        for p in l:
            try: 
                 G[w][p].add((2,0))
            except KeyError: 
                 G[w][p] = set([(2,0)])
        # new bidirected dges
        l = [e for e in D[w] if (0,1) in D[w][e]]
        for p in list(combinations(l, 2)):
            try: 
                 G[p[0]][p[1]].add((2,0))
            except KeyError: 
                 G[p[0]][p[1]] = set([(2,0)])
            try: 
                G[p[1]][p[0]].add((2,0))
            except KeyError: 
                G[p[1]][p[0]] = set([(2,0)])
    return G
def increment_u(G_star, G_u):
    # directed edges
    G_un = directed_inc(G_star,G_u)
    # bidirected edges
    G_un = bidirected_inc(G_un,G_u)
    return G_un

def sample_graph(graph_g,steps=5):
    graph_g_list = [graph_g]
    for i in range(0,steps):
        g = increment_u(graph_g,graph_g_list[-1])
        graph_g_list.append(g)
    return graph_g_list
#BFS implementation of subgraph and supergraph
#Gu to G1 algorithm
from itertools import combinations, permutations
from functools import wraps
import copy
import time
import sys,os
import numpy as np
import ipdb
import operator
from scipy.misc import comb
import math
import gmpy as gmp
import gmpy as gmp
from scipy.misc import comb
import zickle as zkl
import simpleloops as sls
import math
import load_loops
from matplotlib.cbook import flatten
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed

TOOLSPATH='./tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))


circp = zkl.load('circular_p.zkl')
alloops = load_loops.alloops

import pprint
import bfutils as bfu
import traversal as trv
import graphkit as gk
import comparison as cmp
import simpleloops as sl

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = trv.gsig(args[0])         # Signature: just the g
        #s = tool.signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def prune_conflicts(H, g, elist):
    """checks if adding an edge from the list to graph g causes a
    conflict with respect to H and if it does removes the edge
    from the list

    Arguments:
    - `H`: the undersampled graph
    - `g`: a graph under construction
    - `elist`: list of edges to check
    """
    l  = []
    for e in elist:
        gk.addanedge(g,e)
        if not bfu.call_u_conflicts(g, H): l.append(e)
        gk.delanedge(g,e)
    return l

def eqclass(H):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    g = {n:{} for n in H}
    s = set()

    @memo
    def addedges(g,H,edges):
        if edges:
            nedges = prune_conflicts(H, g, edges)
            n = len(nedges)

            if n == 0: return None

            for i in range(n):
                gk.addanedge(g,nedges[i])
                if bfu.call_u_equals(g, H): s.add(bfu.g2num(g))
                addedges(g,H,nedges[:i]+nedges[i+1:])
                gk.delanedge(g,nedges[i])
    edges = gk.edgelist(gk.complement(g))
    addedges(g,H,edges)

    return s-set([None])

# these two functions come fromt his answer:
# http://stackoverflow.com/a/12174125
def set_bit(value, bit): return value | (1<<bit)
def clear_bit(value, bit): return value & ~(1<<bit)
def e2num(e,n): return (1<<(n*n +n - int(e[0],10)*n - int(e[1],10)))
def le2num(elist,n):
    num = 0
    for e in elist:
        num |= e2num(e,n)
    return num
def ekey2e(ekey,n):
    idx = np.unravel_index(n*n - ekey .bit_length() - 1 + 1,(n,n))
    idx = tuple([x+1 for x in idx])
    return ('%i %i'%idx).split(' ')

def cacheconflicts(num, cache):
    """Given a number representation of a graph and an iterable of
    conflicting subgraphs return True if the graph conflicts with any
    of them and false otherwise

    Arguments:
    - `num`: the number representation of a graph
    - `cache`: an iterable of number representations of conflicting
      graphs
    """
    conflict = False
    for c in cache:
        if num & c == c:
            return True
    return False

class nobar:
    def update(self,c): return None
    def finish(self): return None

def start_progress_bar(iter, n, verbose = True):
    if verbose:
        pbar = ProgressBar(widgets=['%3s' % str(iter) +
                                '%10s' % str(n)+' ',
                                Bar('-'), ' '],
                        maxval=n).start()
    else:
        pbar = nobar()
    return pbar

def add2set_loop(ds, H, cp, ccf, iter=1, verbose=True,
                 capsize=100, currsize=0):
    n = len(H)
    n2 = n*n +n
    dsr = {}
    s = set()
    ss = set()

    pbar = start_progress_bar(iter, len(ds), verbose = verbose)

    c = 0

    for gnum in ds:
        c += 1
        pbar.update(c)
        gset = set()
        eset = set()
        for sloop in ds[gnum]:
            if sloop & gnum == sloop: continue
            num = sloop | gnum
            if sloop in ccf and skip_conflictors(num,ccf[sloop]):continue
            if not num in s:
                g = bfu.num2CG(num, n)
                if not bfu.call_u_conflicts(g, H):
                    s.add(num)
                    gset.add((num,sloop))
                    eset.add(sloop)
                    if bfu.call_u_equals(g, H):
                        ss.add(num)
                        if capsize <= len(ss)+currsize: return dsr, ss

        #for gn,e in gset:
        #   if e in cp:
        #       dsr[gn] = eset - cp[e] - set([e])
        #   else:
        #       dsr[gn] = eset - set([e])
        for gn in gset: dsr[gn[0]] = eset - set([gn[1]])
    pbar.finish()
    return dsr, ss

def add2set_(ds, H, cp, ccf, iter=1, verbose=True, capsize=100):
    n = len(H)
    n2 = n*n +n
    dsr = {}
    s = set()
    ss = set()

    pbar = start_progress_bar(iter, len(ds), verbose = verbose)

    c = 0

    for gnum in ds:
        g = bfu.num2CG(gnum, n)
        c += 1
        pbar.update(c)
        glist = []
        elist = []
        eset = set()
        for e in ds[gnum]:
            if not e[1] in g[e[0]]:
                gk.addanedge(g,e)
                num = bfu.g2num(g)
                ekey = (1<<(n2 - int(e[0],10)*n - int(e[1],10)))
                if ekey in ccf and skip_conflictors(num,ccf[ekey]):
                    gk.delanedge(g,e)
                    continue
                if not num in s:
                    s.add(num)
                    if not bfu.call_u_conflicts(g, H):
                        #cf, gl2 = bfu.call_u_conflicts2(g, H)
                        #if not cf:
                        glist.append((num,ekey))
                        elist.append(e)
                        eset.add(ekey)
                        if bfu.call_u_equals(g, H):
                            ss.add(num)
                            #if bfu.call_u_equals2(g, gl2, H): ss.add(num)
                        if capsize <= len(ss): break
                gk.delanedge(g,e)

        for gn,e in glist:
            if e in cp:
                dsr[gn] = [ekey2e(k,n) for k in eset - cp[e]]
            else:
                dsr[gn] = elist
        if capsize <= len(ss): return dsr, ss

    pbar.finish()
    return dsr, ss

def skip_conflictors(gnum, ccf):
    pss = False
    for xx in ccf:
        if xx&gnum == xx:
            pss = True
            break
    return pss

def bconflictor(e,H):
    n = len(H)
    s = set()
    for v in H:
        s.add(e2num((v,e[0]),n)|e2num((v,e[1]),n))
    return s

def conflictor(e,H):
    n = len(H)
    def pairs(e):
        ekey = e2num(e,n)
        return [ekey|e2num((e[0],e[0]),n),
                ekey|e2num((e[1],e[1]),n)]

    def trios(e,H):
        s = set()
        for v in H:
            if not v in e:
                s.add(e2num((e[0],v),n)|
                      e2num((v,e[1]),n)|
                      e2num((e[1],e[1]),n))
                s.add(e2num((e[0],e[0]),n)|
                      e2num((e[0],v),n)|
                      e2num((v,e[1]),n))
                s.add(e2num((e[0],v),n)|
                      e2num((v,v),n)|
                      e2num((v,e[1]),n))
        return s

    return trios(e,H).union(pairs(e))

def conflictor_set(H):
    s = set()
    for x in gk.inedgelist(H):  s = s | conflictor(x,H)
    for x in gk.inbedgelist(H): s = s | bconflictor(x,H)
    return s

def conflictors(H):
    s = conflictor_set(H)
    ds = {}
    num = reduce(operator.or_,s)
    for i in xrange(gmp.bit_length(num)):
        if num & 1<<i:
            ds[1<<i] = [x for x in s if x&(1<<i)]
    return ds

def may_be_true_selfloop(n,H):
    for v in H[n]:
        if v == n: continue
        if (0,1) in H[n][v] and not ((2,0) in H[n][v]): return False
    return True

def issingleloop(num):
    bl = gmp.bit_length(num)
    idx = [1 for i in xrange(bl) if num & (1<<i)]
    return len(idx) == 1

def nonbarren(H):
    for v in H:
        if H[v]: return v
    return False

def prune_loops(loops, H):
    l = []
    n = len(H)
    for loop in loops:
        g = bfu.num2CG(loop, n)
        x = [k for k in g if g[k]]
        if len(x) == 1:
            s = reduce(lambda x, s: s.union(x),
                       [H[x[0]][w] for w in H[x[0]]])
            if not (2,0) in s: continue
        if not bfu.call_u_conflicts_d(g, H): l.append(loop)        
    return l

def lconflictors(H, sloops=None):
    if not sloops: sloops = prune_loops(allsloops(len(H)),H)
    s = conflictor_set(H)
    ds = {}
    num = reduce(operator.or_,s)
    for i in xrange(gmp.bit_length(num)):
        if num & 1<<i:
            cset = [x for x in s if x&(1<<i)]
            for sloop in sloops:
                if sloop & 1<<i:
                    ds.setdefault(sloop,[]).extend(cset)
    return ds

def confpairs(H):
    n = len(H)
    g = {n:{} for n in H}
    d = {}

    edges = gk.edgelist(gk.complement(g))
    edges = prune_conflicts(H, g, edges)

    for p in combinations(edges,2):
        gk.addedges(g,p)
        if bfu.call_u_conflicts(g, H):
            n1 = e2num(p[0],n)
            n2 = e2num(p[1],n)
            d.setdefault(n1,set()).add(n2)
            d.setdefault(n2,set()).add(n1)
        gk.deledges(g,p)

    return d

def lconfpairs(H, cap=10, sloops=None):
    n = len(H)
    d = {}
    if not sloops: sloops = prune_loops(allsloops(len(H)),H)
    c = 0
    for p in combinations(sloops,2):
        g = bfu.num2CG(p[0]|p[1], n)
        if bfu.call_u_conflicts(g, H):
            d.setdefault(p[0],set()).add(p[1])
            d.setdefault(p[1],set()).add(p[0])
        if c >= cap: break
        c +=1
    return d


def iteqclass(H, verbose=True, capsize=100):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    if cmp.isSclique(H):
        print 'not running on superclique'
        return None
    g = {n:{} for n in H}
    s = set()
    Hnum = bfu.ug2num(H)
    if Hnum[1]==0: s.add(Hnum[0])

    cp = confpairs(H)
    ccf = conflictors(H)

    edges = gk.edgelist(gk.complement(g))
    ds = {bfu.g2num(g): edges}

    if verbose: print '%3s'%'i'+'%10s'%' graphs'
    for i in range(len(H)**2):
        ds, ss = add2set_(ds, H, cp, ccf, iter=i,
                            verbose=verbose,
                            capsize=capsize)
        s = s | ss
        if capsize <= len(ss): break
        if not ds: break

    return s

def liteqclass(H, verbose=True, capsize=100, asl=None):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set([-1])
    g = {n:{} for n in H}
    s = set()



    if asl:
        sloops = asl
    else:
        sloops = prune_loops(allsloops(len(H)),H)

    cp  = []#lconfpairs(H, sloops=sloops)
    ccf = lconflictors(H, sloops=sloops)
    ds = {0: sloops}

    if verbose: print '%3s'%'i'+'%10s'%' graphs'
    i=0
    while ds:
        ds, ss = add2set_loop(ds, H, cp, ccf, iter=i,
                              verbose=verbose,
                              capsize=capsize,
                              currsize=len(s))
        s = s | ss
        i += 1
        if capsize <= len(s): break

    return s

def edgemask(gl,H, cds):
    """given a list of encoded graphs and observed undersampled graph
    H returns a matrix with -1 on diagonal, 0 at the conflicting graph
    combination and encoded graph at non-conflicted
    positions. Furthermore, returns a set of graphs that are in the
    equivalence class of H

    Arguments:
    - `gl`: list of integer encoded graphs
    - `H`: the observed undersampled graph
    """
    n = len(H)
    nl= len(gl)
    s = set()
    mask = np.zeros((nl,nl),'int')
    np.fill_diagonal(mask,-1)

    for i in xrange(nl):
        for j in xrange(i+1,nl):

            if gl[i] & gl[j]: continue
            if skip_conflict(gl[i], gl[j], cds): continue

            gnum = gl[i] | gl[j]
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                if bfu.call_u_equals(g, H): s.add(gnum)
                mask[i,j] = gnum
                mask[j,i] = gnum
    return mask, s

def ledgemask(gl,H, cds):
    """given a list of encoded graphs and observed undersampled graph
    H returns a matrix with -1 on diagonal, 0 at the conflicting graph
    combination and encoded graph at non-conflicted
    positions. Furthermore, returns a set of graphs that are in the
    equivalence class of H

    Arguments:
    - `gl`: list of integer encoded graphs
    - `H`: the observed undersampled graph
    """
    n = len(H)
    nl= len(gl)
    s = set()
    mask = np.zeros((nl,nl),'int')
    np.fill_diagonal(mask,-1)

    for i in xrange(nl):
        for j in xrange(i+1,nl):

            if gl[i] & gl[j]: continue
            gnum = gl[i] | gl[j]
            if skip_conflictors(gnum, cds): continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                if bfu.call_u_equals(g, H): s.add(gnum)
                mask[i,j] = gnum
                mask[j,i] = gnum
    return mask, s

def edgeds(mask):
    """construct an edge dictionary from the mask matrix

    Arguments:
    - `mask`:
    """
    ds = {}
    nl = mask.shape[0]
    idx = np.triu_indices(nl,1)
    for i,j in zip(idx[0], idx[1]):
        if mask[i,j]:
            ds[(i,j)] = set()
            conf = set([i,j])
            conf = conf.union(np.where(mask[i,:]==0)[0])
            conf = conf.union(np.where(mask[j,:]==0)[0])
            for k,m in zip(idx[0], idx[1]):
                if not mask[k,m]: continue
                if k in conf: continue
                if m in conf: continue
                if not (k,m) in ds: ds[(i,j)].add(mask[k,m])
            if not ds[(i,j)]: ds.pop((i,j))
    return ds

def edgedsg(mask):
    """construct an edge dictionary from the mask matrix

    Arguments:
    - `mask`:
    """
    ds = {}
    nl = mask.shape[0]
    idx = np.triu_indices(nl,1)
    for i,j in zip(idx[0], idx[1]):
        if mask[i,j]:
            ds[mask[i,j]] = set()
            conf = set([i,j])
            conf = conf.union(np.where(mask[i,:]==0)[0])
            conf = conf.union(np.where(mask[j,:]==0)[0])
            for k,m in zip(idx[0], idx[1]):
                if not mask[k,m]: continue
                if k in conf: continue
                if m in conf: continue
                if not (k,m) in ds: ds[mask[i,j]].add(mask[k,m])
            if not ds[mask[i,j]]: ds.pop(mask[i,j])
    return ds


def quadlister(glist, H, cds):
    n = len(H)
    s = set()
    cache = {}

    def edgemask(gl, H, cds):
        nl= len(gl)
        ss = set()
        mask = np.zeros((nl,nl),'int')
        np.fill_diagonal(mask,-1)
        idx = np.triu_indices(nl,1)
        for i,j in zip(idx[0], idx[1]):
            if gl[i] & gl[j]:
                mask[i,j] = -1
                mask[j,i] = -1
                continue
            if skip_conflict(gl[i], gl[j], cds):
                gnum = gl[i] | gl[j]
                cache[gnum] = False
                continue

            gnum = gl[i] | gl[j]
            if gnum in cache:
                if cache[gnum]:
                    mask[i,j] = gnum
                    mask[j,i] = gnum
            else:
                cache[gnum] = False
                g = bfu.num2CG(gnum,n)
                if not bfu.call_u_conflicts(g, H):
                    if bfu.call_u_equals(g, H): ss.add(gnum)
                    mask[i,j] = gnum
                    mask[j,i] = gnum
                    cache[gnum] = True
        return mask, ss


    def quadmerger(gl, H, cds):
        mask, ss = edgemask(gl, H, cds)
        ds = edgeds(mask)
        #ipdb.set_trace()
        return [[mask[x]]+list(ds[x]) for x in ds], ss

    l = []
    for gl in glist:
        ll, ss = quadmerger(gl, H, cds)
        l.extend(ll)
        s = s|ss

    return l, s


def dceqc(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cds = confpairs(H)

    glist =  [2**np.arange(n**2)]
    i = 1
    #for i in range(int(np.log2(n**2))):
    while glist != []:
        print i, np.max(map(len,glist)), len(glist)
        glist, ss = quadlister(glist, H, cds)
        s = s|ss
        i += 1
    return s

def quadmerge(gl, H, cds):
    n = len(H)
    l = set()
    s = set()
    mask, ss = edgemask(gl, H, cds)
    s = s | ss
    ds = edgeds(mask)

    #pp = pprint.PrettyPrinter(indent=1)
    #pp.pprint(ds)

    for idx in ds:
        for gn in ds[idx]:
            if mask[idx]&gn: continue
            if skip_conflict(mask[idx], gn, cds): continue
            gnum = mask[idx] | gn
            if gnum in l or gnum in ss: continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                l.add(gnum)
                if bfu.call_u_equals(g, H): s.add(gnum)

    return list(l), s

def skip_conflict(g1, g2, ds):
    pss = False
    for ekey in ds:
        if (g1 & ekey) == ekey:
            if ekey in ds and cacheconflicts(g2,ds[ekey]):
                pss = True
                break
    return pss

def edgemask2(gl,H, cds):
    n = len(H)
    nl= len(gl)
    s = set()
    o = set()
    mask = np.zeros((nl,nl),'int')
    np.fill_diagonal(mask,-1)
    for i in xrange(nl):
        for j in xrange(i+1,nl):
            if gl[i] & gl[j]: continue
            gnum = gl[i] | gl[j]
            if skip_conflictors(gnum, cds): continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                if bfu.call_u_equals(g, H): s.add(gnum)
                mask[i,j] = gnum
                mask[j,i] = gnum
            elif bfu.overshoot(g, H): o.add(gnum)
    return mask, s, o # mask, found eqc members, overshoots

def patchmerge(ds, H, cds):
    n = len(H)
    l = set()
    s = set()
    o = set()
    for gkey in ds:
        for num in ds[gkey]:
            if gkey & num: continue
            gnum = gkey | num
            if gnum is s: continue
            if skip_conflictors(gnum, cds): continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                l.add(gnum)
                if bfu.call_u_equals(g, H): s.add(gnum)
            elif not gnum in o and bfu.overshoot(g, H): o.add(gnum)
    return l, s, o

def quadmerge2(gl, H, cds):
    n = len(H)

    mask, s, o = edgemask2(gl, H, cds)
    #ipdb.set_trace()
    ds = edgedsg(mask)
    l, ss, oo = patchmerge(ds, H, cds)

    o = o | oo
    s = s | ss

    print 'overshoots: ', len(o)

    return list(l), s

def quadmerge21(gl, H, cds):
    n = len(H)
    l = set()

    mask, ss, o = edgemask2(gl, H, cds)
    idx = np.triu_indices(mask.shape[0], 1)
    print len(o)
    for i in range(len(idx[0])):
        if mask[idx[0][i],idx[1][i]]: l.add(mask[idx[0][i],idx[1][i]])

    return list(l), ss

def dceqclass2(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cp = confpairs(H)
    confs = conflictor_set(H)
    ccf = conflictors(H)

    def prune_loops(gl, H):
        l = []
        for e in gl:
            if e[0] == e[1] and not (e[1] in H[e[0]] and (1,0) in H[e[0]][e[1]]): continue
            l.append(e)
        return l
    edges = gk.edgelist(gk.complement(bfu.num2CG(0,n)))
    edges = prune_loops(edges, H)
    glist = map(lambda x: e2num(x,n),edges)

    #glist =  list(2**np.arange(n**2))
    i = 0
    while glist != []:
        print 2**i, len(glist)
        glist_prev = glist
        glist, ss = quadmerge21(glist, H, confs)
        s = s|ss
        i += 1


    ds = {x: edges for x in glist_prev}

    for j in range(i, len(H)**2):
        ds, ss = add2set_(ds, H, cp, ccf, iter=j, verbose=True)
        s = s | ss
        if not ds: break

    return s


def dceqclass(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cds = confpairs(H)

    glist =  [0]+list(2**np.arange(n**2))
    i = 1
    while glist != []:
        print i, len(glist)
        glist, ss = quadmerge(glist, H, cds)
        s = s|ss
        i += 1
    return s

def ldceqclass(H,asl=None):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cds = lconfpairs(H)
    if asl:
        sloops = asl
    else:
        sloops = prune_loops(allsloops(len(H)),H)

    glist =  sloops
    i = 1
    while glist != []:
        print i, len(glist)
        glist, ss = lquadmerge(glist, H, cds)
        s = s|ss
        i += 1
    return s

def lquadmerge(gl, H, cds):
    n = len(H)
    l = set()
    s = set()
    mask, ss = ledgemask(gl, H, cds)
    s = s | ss
    ds = edgeds(mask)

    #pp = pprint.PrettyPrinter(indent=1)
    #pp.pprint(ds)

    for idx in ds:
        for gn in ds[idx]:
            if mask[idx]&gn: continue
            if skip_conflictors(mask[idx], gn, cds): continue
            gnum = mask[idx] | gn
            if gnum in l or gnum in ss: continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                l.add(gnum)
                if bfu.call_u_equals(g, H): s.add(gnum)

    return list(l), s

def quadmerge_(glist, H, ds):
    n = len(H)
    gl = set()
    ss = set()
    conflicts = set()
    for gi in combinations(glist, 2):
        if gi[0] & gi[1]: continue
        #if skip_conflict(gi[0], gi[1], ds): continue
        gnum = gi[0] | gi[1]
        if gnum in conflicts: continue
        if skip_conflictors(gnum, ds):
            conflicts.add(gnum)
            continue
        if gnum in gl: continue
        g = bfu.num2CG(gnum,n)
        if not bfu.call_u_conflicts(g, H):
            gl.add(gnum)
            if bfu.call_u_equals(g, H): ss.add(gnum)
        else:
            conflicts.add(gnum)
    return gl, ss

def ecmerge(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return None
    n = len(H)
    s = set()
    ds = confpairs(H)
    ccf = conflictors(H)
    cset = set()
    for e in ccf:
        cset = cset.union(ccf[e])

    glist =  np.r_[[0],2**np.arange(n**2)]
    #glist =  2**np.arange(n**2)

    #glist, ss = quadmerge(glist,H)

    for i in range(int(2*np.log2(n))):
        print i, len(glist)
        glist, ss = quadmerge_(glist,H, cset)
        s = s | ss
    return s

def getrates(g,H):
    n = len(H)
    au = bfu.call_undersamples(g)
    return list(np.where(map(lambda x: x == H, au))[0])

def withrates(s,H):
    n = len(H)
    d = {g:set() for g in s}
    for g in s:
        d[g] = getrates(bfu.num2CG(g,n),H)
    return d

def add2set(gset, elist, H):
    n = len(H)

    s = set()
    ss = set()

    eremove = {e: True for e in elist}

    for gnum in gset:
        g = bfu.num2CG(gnum, n)
        for e in elist:
            if not e[1] in g[e[0]]:
                gk.addanedge(g,e)
                num = bfu.g2num(g)
                if not num in s:
                    au = bfu.call_undersamples(g)
                    if not gk.checkconflict(H, g, au=au):
                        eremove[e] = False
                        s.add(num)
                        if gk.checkequality(H, g, au=au): ss.add(num)
                gk.delanedge(g,e)

    for e in eremove:
        if eremove[e]: elist.remove(e)

    return s, ss, elist

def eqclass_list(H):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    g = {n:{} for n in H}
    s = set()

    edges = gk.edgelist(gk.complement(g))
    #edges = prune_conflicts(H, g, edges)

    gset = set([bfu.g2num(g)])
    for i in range(len(H)**2):
        print i
        gset, ss, edges = add2set(gset, edges, H)
        s = s | ss
        if not edges: break

    return s

def loop2graph(l,n):
    g = {str(i):{} for i in range(1,n+1)}
    for i in range(len(l)-1):
        g[l[i]][l[i+1]] = set([(0,1)])
    g[l[-1]][l[0]] = set([(0,1)])
    return g

def set_loop(loop, graph):
    for i in range(0,len(loop)-1):
        graph[loop[i]][loop[i+1]] = set([(0,1)])
    graph[loop[-1]][loop[0]] = set([(0,1)])

def rotate(l,n): return l[n:] + l[:n]
def get_perm(loop1, loop2, n=None):
    if not n: n = len(loop1)
    basel = [str(i) for i in xrange(1,n+1)]
    diff1 = set(basel) - set(loop1)
    diff2 = set(basel) - set(loop2)
    if loop1[0] in loop2:
        l2 = rotate(loop2, loop2.index(loop1[0]))
    else:
        l2 = loop2
    mp = {}
    for x,y in zip(loop1+list(diff1),l2+list(diff2)):
        mp[x] = y
    return mp

def permute(g, perm):
    gn = {x:{} for x in g}
    for e in g:
        gn[perm[e]] = {perm[x]:g[e][x] for x in g[e]}
    return gn

def permuteAset(s, perm):
    n = len(perm)
    ns = set()
    for e in s:
        ns.add(bfu.g2num(permute(bfu.num2CG(e,n),perm)))
    return ns

def noverlap_loops(loops):
    d = {}
    for l in loops:
        el = []
        for k in loops:
            if not set(l) & set(k):
                el.append(tuple(k))
                #d.setdefault(tuple(l),set()).add(tuple(k))
        d[tuple(l)] = noverlap_loops(el)
    return d

def loop_combinations(loops):
    s = set()
    d = noverlap_loops(loops)
    def dfs_traverse(d, gs):
        if d:
            for e in d:
                dfs_traverse(d[e], gs|set([e]))
        else:
            s.add(frozenset(gs))
    for e in d:
        dfs_traverse(d[e],set([e]))
    return list(s)

def sorted_loops(g):
    l = [x for x in sl.simple_loops(g,0)]
    s = {}
    for e in l:
        s.setdefault(len(e),[]).append(e)
    return s

def loopgroups(g, n=None):
    d = sorted_loops(g)
    if n:
        return loop_combinations(d[n])
    else:
        l=[]
        for key in d:
            l.append(loop_combinations(d[key]))
        return l

def count_loops(n):
    s = 0
    for i in range(1,n+1):
        s += comb(n,i) * math.factorial(i-1)
    return s

def perm_cyclic(l): return [tuple(l[i:]+l[:i]) for i in range(len(l))]
def hashloop(l):
    t = [int(x) for x in l]
    idx = np.argmin(t)
    return tuple(l[idx:]+l[:idx])
def perm_circular_slow2(l):
    s = [tuple(l)]
    c = {}
    c[hashloop(l)] = True
    for e in permutations(l):
        if not hashloop(e) in c:
            s.append(e)
            c[hashloop(e)] = True
    return s

def perm_circular_slow(l):
    s = [tuple(l)]
    c = set(perm_cyclic(l))
    for e in permutations(l):
        if not e in c:
            s.append(e)
            c = c | set(perm_cyclic(e))
    return s
def perm_circular(l, cp=circp):
    r = []
    n = len(l)
    for e in cp[n]:
        r.append([l[i] for i in e])
    return r

def gen_loops(n):
    l = [str(i) for i in range(1,n+1)]
    s = []
    for i in range(1,n+1):
        for e in combinations(l,i):
            s.extend(perm_circular(e))
    return s

def allsloops(n, asl = alloops):
    if asl: return asl[n]
    s = []
    l = gen_loops(n)
    for e in l:
        s.append(bfu.g2num(loop2graph(e,n)))
    return s

def reverse(H, verbose=True, capsize=1000):
    n = len(H)
    s = set()

    g      = gk.superclique(n)
    sloops = set(allsloops(n))

    ds = {bfu.g2num(g): sloops}

    if verbose: print '%3s'%'i'+'%10s'%' graphs'
    i=0
    while ds:
        ds, ss = del_loop(ds, H, iter=i,
                          verbose=verbose,
                          capsize=capsize)
        s = s | ss
        i += 1
        if capsize <= len(s): break

    return s

# ----------------------

def build_loop_step(ds, loop, n, iter=1):
    n2 = n*n +n
    dsr = {}
    s = set()
    ss = set()

    pbar = start_progress_bar(iter, len(ds))

    c = 0

    for gnum in ds:
        c += 1
        pbar.update(c)
        gset = set()
        eset = set()
        for sloop in ds[gnum]:
            num = sloop | gnum
            if not num in s:
                g = bfu.num2CG(num, n)
                s.add(num)
                if bfu.forms_loop(g, loop):
                    ss.add(num)
                else:
                    gset.add((num,sloop))
                    eset.add(sloop)

        for gn,e in gset:
            dsr[gn] = eset - set([e])

    pbar.finish()
    return dsr, ss

def forward_loop_match(loop, n):
    """start with an empty graph and keep adding simple loops until
    the loop is generated at some undersampling rate

    Arguments:
    - `loop`: binary encoding of the loop
    - `n`: number of nodes in the graph
    """
    s = set()
    sloops = allsloops(n)
    ds = {0: sloops}

    i=0
    while ds:
        ds, ss = build_loop_step(ds, loop, n, iter=i)
        s = s | ss
        i += 1

    return s

def delAloop(g, loop):
    n = len(g)
    l = []
    l = [bfu.g2num(ur.loop2graph(s,n)) for s in sls.simple_loops(g,0)]
    l = [num for num in l if not num == loop ]
    print loop, ': ',  l
    return bfu.num2CG(reduce(operator.or_, l),n)

def reverse_loop_match(g,loop):
    """start with a graph and keep removing loops while the loop is still matched

    Arguments:
    - `g`: graph that generates the loop
    - `loop`: the reference loop
    """
    s = set()
    n = len(g)

    def prune(g):
        numh = bfu.g2num(g)
        cannotprune = True
        for l in sls.simple_loops(gk.digonly(g),0):
            gg = delAloop(g,bfu.g2num(loop2graph(l,n)))
            if bfu.forms_loop(gg, loop):
                cannotprune = False
                prune(gg)
        if cannotprune:
            print 'one'
            s.add(g)

    prune(g)
    return s

def reverse_edge_match(g,loop):
    """start with a graph and keep removing loops while the loop is still matched

    Arguments:
    - `g`: graph that generates the loop
    - `loop`: the reference loop
    """
    s = set()
    n = len(g)

    def prune(g):
        numh = bfu.g2num(g)
        cannotprune = True
        for l in gk.edgelist(gk.digonly(g)):
            gk.delanedge(g,l)
            if bfu.forms_loop(g, loop):
                cannotprune = False
                prune(g)
            gk.addanedge(g,l)
        if cannotprune: s.add(bfu.g2num(g))

    prune(g)
    return s

def matchAloop(loop, n):
    """returns a set of minimal graphs that generate this loop

    Arguments:
    - `loop`: binary encoding of the loop
    - `n`: number of nodes in the graph
    """
    s = set()
    l = forward_loop_match(loop,n)
    print len(l)
    for g in l:
        s = s | reverse_edge_match(bfu.num2CG(g,n),loop)

    return s

# ----------------------

def del_loop(ds, H, iter=0, verbose=True, capsize=1000):
    n = len(H)

    dsr = {}
    s = set()
    ss = set()
    print iter,
    for gnum in ds:
        gset = []
        s = set()
        for sloop in ds[gnum]:
            rset = ds[gnum] - set([sloop])
            num = reduce(operator.or_, rset)
            if not num in s:
                g = bfu.num2CG(num, n)
                if bfu.overshoot(g, H):
                    s.add(num)
                    gset.append((num,rset))

        if gset == []:
            print '.',
            ss.add(gnum)

        for gn in gset: dsr[gn[0]] = gn[1]
    print ''
    return dsr, ss

def main():
    g = bfu.ringmore(6,1);
    H = bfu.undersample(g,1);
    ss = liteqclass(H)
    print ss

if __name__ == "__main__":
    main()
from networkx import strongly_connected_components
from functools import wraps
import scipy
import numpy as np
import itertools, copy, time
import random
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

# local packages
import bfutils as bfu
import graphkit as gk
import comparison
import ecj


def isedgesubsetD(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if (0,1) in g2star[n][h]:
                if h in g2[n]:
                    if not (0,1) in g2[n][h]:
                        return False
                else:
                    return False
    return True

def isedgesubset(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                #if not (0,1) in g2[n][h]:
                if not g2star[n][h].issubset(g2[n][h]):
                    return False
            else:
                    return False
    return True

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def purgepath(path, l):
    for i in range(1,len(path)-1):
        l.remove((path[i],path[i+1]))

def next_or_none(it):
    try:
        n = it.next()
    except StopIteration:
        return None
    return n

def try_till_d_path(g,d,gt,order=None):
    k = []
    i = 1
    while not k:
        if order:
            k = [x for x in length_d_paths(g,order(i),d)]
        else:
            k = [x for x in length_d_paths(g,str(i),d)]
        i += 1
        if i > len(g): return []

    ld = []
    for i in range(min(10, len(k))):
        ld.append(len(checkApath(['2']+k[i], gt)))
    idx = np.argmin(ld)
    return k[0]

def try_till_path(g, gt):
    d = len(g)-1
    gx = comparison.graph2nx(g)
    sccl = [x for x in strongly_connected_components(gx)]
    # take the largest
    ln = [len(x) for x in sccl]
    idx = np.argsort(ln)
    d = len(sccl[idx[-1]])-1
    sccl = [sccl[i] for i in idx[::-1]]
    order = [item for sublist in sccl for item in sublist]
    k = []
    while not k:
        if d < 5: return []
        k = try_till_d_path(g,d,gt)
        d -= 1
    return k

def gpurgepath(g,path):
    for i in range(1,len(path)-1):
        del g[path[i]][path[i+1]]

def forks(n,c,el,bl,doempty=lambda x,y: True):
    '''
    INPUT:
    n - single node string
    c - mutable list of children of node n (will be changed as a side effect)
    el - list of edges available for construction (side effect change)
    bl - list of bidirected edges

    OUTPUT:
    l - list of forks
    '''
    l = []
    r = set()
    for p in [x for x in itertools.combinations(c,2)]:
        if doempty(p, bl) and (n,p[0]) in el and (n,p[1]) in el:
            l.append((n,)+p)
            el.remove((n,p[0]))
            el.remove((n,p[1]))
            r.add(p[0])
            r.add(p[1])
    for e in r: c.remove(e)
    return l

def childrenedges(n,c,el,bl):
    '''
    INPUT:
    n - single node string
    c - mutable list of children of node n (will be changed as a side effect)
    el - list of edges available for construction (side effect change)
    bl - list of bidirected edges

    OUTPUT:
    l - list of
    '''
    l = []
    r = set()
    for ch in c:
        if (n,ch) in el:
            l.append((n,ch))
            el.remove((n,ch))
            r.add(ch)
    for e in r: c.remove(e)
    return l

# an empty fork is a for without the bidirected edge
def make_emptyforks(n,c,el,bl):
    return forks(n,c,el,bl,doempty=lambda x,y: not x in y)
def make_fullforks(n,c,el,bl):
    return forks(n,c,el,bl,doempty=lambda x,y: x in y)

def make_longpaths(g, el):
    '''
    el - list of edges that is modified as a side effect
    '''
    l = []
    gc = copy.deepcopy(g)
    for i in range(16):
        k = try_till_path(gc,g)
        if len(k) < 5: break
        if k:
            l.append(('2',)+tuple(k))
            purgepath(l[-1],el)
            gpurgepath(gc,l[-1])
        else:
            break
    return l

def make_allforks_and_rest(g,el,bl, dofullforks=True):
    '''
    el - list of edges that is modified as a side effect
    '''
    l = []
    r = []
    nodes = [n for n in g]
    random.shuffle(nodes)
    for n in nodes:

        c = [e for e in g[n] if (0,1) in g[n][e]]# all children
        if len(c) == 1:
            if (n,c[0]) in el:
                r.append((n,c[0]))
                el.remove((n,c[0]))
        elif len(c) > 1:
            l.extend(make_emptyforks(n,c,el,bl))
            if dofullforks: l.extend(make_fullforks(n,c,el,bl))
            r.extend(childrenedges(n,c,el,bl))
    return l, r

def vedgelist(g, pathtoo=False):
    """ Return a list of tuples for edges of g and forks
    a superugly organically grown function that badly needs refactoring
    """
    l = []
    el = gk.edgelist(g)
    bl = gk.bedgelist(g)

    if pathtoo: l.extend(make_longpaths(g,el))
    l2,r = make_allforks_and_rest(g,el,bl,dofullforks=True)
    l.extend(l2)

    A, singles = makechains(r)

    if singles:
        B, singles = makesinks(singles)
    else:
        B, singles = [], []

    l = longpaths_pick(l)+threedges_pick(l) + A + B + singles
    return l

def twoedges_pick(l):  return [e for e in l if len(e)==2]
def threedges_pick(l): return [e for e in l if len(e)==3]
def longpaths_pick(l): return [e for e in l if len(e)>3 and e[0]=='2']
def makechains(l):
    """ Greedily construct 2 edge chains from edge list
    """
    ends = {e[1]:e for e in l}
    starts = {e[0]:e for e in l}
    r = []
    singles = []
    while l:

        e = l.pop()
        if e[1] in starts and e[0] != e[1] and starts[e[1]] in l:
            r.append(('0', e[0],)+starts[e[1]])
            l.remove(starts[e[1]])
        elif e[0] in ends and e[0] != e[1] and ends[e[0]] in l:
            r.append(('0',)+ends[e[0]]+(e[1],))
            l.remove(ends[e[0]])
        else:
            singles.append(e)
    return r, singles
def makesink(es): return ('1', es[0][0],) + es[1]
def makesinks(l):
    """ Greedily construct 2 edge sinks ( b->a<-c ) from edge list
    """
    sinks = {}
    for e in l:
        if e[1] in sinks:
            sinks[e[1]].append(e)
        else:
            sinks[e[1]] = [e]
    r = []
    singles = []
    for e in sinks:
        if len(sinks[e])>1:
            for es in chunks(sinks[e],2):
                if len(es)==2:
                    r.append(makesink(es))
                else:
                    singles.append(es[0])
        else:
            singles.append(sinks[e][0])
    return r, singles

def selfloop(n,g):
    return n in g[n]

def selfloops(l,g):
    return reduce(lambda x,y: x and y, map(lambda x: selfloop(x,g), l))

def checkbedges(v,bel,g2):
    r = []
    for e in bel:
        if e == tuple(v[1:]) and not selfloops(e, g2):
            r.append(e)
        if e == (v[2],v[1]) and not selfloops(e, g2):
            r.append(e)
    for e in r: bel.remove(e)
    return bel

def checkedge(e, g2):
    if e[0] == e[1]:
        l = [n for n in g2 if n in g2[n]]
        s = set()
        for v in g2[e[0]]: s=s.union(g2[e[0]][v])
        if not (2,0) in s:
            l.remove(e[0])
        return l
    else:
        return g2.keys()

def single_nodes(v,g2):
    """ Returns a list of singleton nodes allowed for merging with v
    """
    l = [(n,n) for n in g2 if not n in v and len(g2[n])>1]
    return l

def checkvedge(v, g2):
    """ Nodes to check to merge the virtual nodes of v ( b<-a->c )
    """
    l = gk.bedgelist(g2)
    if (v[1],v[2]) in l:
        l = single_nodes(v,g2) + checkbedges(v,l,g2)
        for n in v:
            if n in g2[n]: l.append((n,n))
    else:
        l = checkbedges(v,l,g2)
    return list(set(l))

def checkAedge(v, g2):
    """ Nodes to check to merge the virtual nodes of A ( b->a<-c )
    """
    l = []
    # try all pairs but the sources
    for pair in itertools.combinations(g2,2):
        #if pair == (v[1],v[2]): continue
        #if pair == (v[2],v[1]): continue
        l.append(pair)
        l.append(pair[::-1])
    for n in g2:
        l.append((n,n))
    return l

def checkcedge(c, g2):
    """ Nodes to check to merge the virtual nodes of c ( a->b->c )
    """
    l = gk.edgelist(g2)
    return list(set(l))

def checkApath(p, g2):
    sl = [x for x in g2 if selfloop(x,g2)]
    d = len(p) - 2
    l = []
    for n in g2:
        l.extend([tuple(x) for x in length_d_loopy_paths(g2, n, d, p[1:])])
    #k = prunepaths_1D(g2, p, l)
    return l


def isedge(v):  return len(v) == 2 # a->b
def isvedge(v): return len(v) == 3 # b<-a->c
def isCedge(v): return len(v) == 4 and v[0] == '0' # a->b->c
def isAedge(v): return len(v) == 4 and v[0] == '1'# a->c<-b
def isApath(v):  return len(v) >= 4 and v[0] == '2'# a->b->...->z

def checker(n,ee):
    g = bfu.ringmore(n,ee)
    g2 = bfu.increment(g)
    d = checkable(g2)
    t = [len(d[x]) for x in d]
    r = []
    n = len(g2)
    ee= len(gk.edgelist(g2))
    for i in range(1,len(t)):
        r.append(sum(np.log10(t[:i])) - ee*np.log10(n))
    return r

def checkerDS(n,ee):
    g = bfu.ringmore(n,ee)
    g2 = bfu.increment(g)
    gg = checkable(g2)
    d,p,idx = conformanceDS(g2,gg,gg.keys())
    t = [len(x) for x in p]
    r = []
    n = len(g2)
    ee= len(gk.edgelist(g2))
    for i in range(1,len(t)):
        r.append(sum(np.log10(t[:i])) - ee*np.log10(n))
    return r

def fordens(n,denslist, repeats=100):
    rl={}
    for d in denslist:
        ee = bfu.dens2edgenum(d,n)
        l=[checker(n,ee)[-1] for i in range(repeats)]
        rl[d] = (round(scipy.mean(l),3),round(scipy.std(l),3))
    return rl

def fordensDS(n,denslist, repeats=100):
    rl={}
    for d in denslist:
        print d
        ee = bfu.dens2edgenum(d,n)
        l=[checkerDS(n,ee)[-1] for i in range(repeats)]
        rl[d] = (round(scipy.mean(l),3),round(scipy.std(l),3))
    return rl

def checkable(g2):
    d = {}
    g = cloneempty(g2)
    vlist = vedgelist(g2,pathtoo=False)
    for v in vlist:
        if isvedge(v):
            d[v] = checkvedge(v,g2)
        elif isCedge(v):
            d[v] = checkcedge(v,g2)
        elif isAedge(v):
            d[v] = checkAedge(v,g2)
        elif isApath(v):
            d[v] = checkApath(v,g2)
        else:
            d[v] = checkedge(v,g2)

    # # check if some of the otherwise permissible nodes still fail
    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge),
         (addaAedge,delaAedge),
         (addapath,delapath)]
    c = [ok2add2edges,
         ok2addavedge,
         ok2addacedge,
         ok2addaAedge,
         ok2addapath]

    for e in d:
        adder, remover = f[edge_function_idx(e)]
        checks_ok = c[edge_function_idx(e)]
        for n in d[e]:
            if not checks_ok(e,n,g,g2):
                d[e].remove(n)

    return d

def inorder_check2(e1, e2, j1, j2, g2, f=[], c=[]):
    g = cloneempty(g2) # the graph to be used for checking

    if f==[]:
        f = [(add2edges,del2edges,mask2edges),
             (addavedge,delavedge,maskavedge),
             (addacedge,delacedge,maskaCedge),
             (addaAedge,delaAedge,maskaAedge),
             (addapath,delapath,maskapath)]

    if c==[]:
        c = [ok2add2edges,
             ok2addavedge,
             ok2addacedge,
             ok2addaAedge,
             ok2addapath]

    adder, remover, masker = f[edge_function_idx(e1)]
    checks_ok = c[edge_function_idx(e2)]

    d = {}
    s1 = set()
    s2 = set()
    for c1 in j1: # for each connector
        mask = adder(g,e1,c1)
        d[c1] = set()
        for c2 in j2:
            if checks_ok(e2,c2,g,g2):
               d[c1].add(c2)
               s2.add(c2)
        remover(g,e1,c1,mask)
        if d[c1]: s1.add(c1)
    return d,s1,s2

def check3(e1, e2, e3, j1, j2, j3, g2, f=[], c=[]):
    g = cloneempty(g2) # the graph to be used for checking
    if f==[]:
        f = [(add2edges,del2edges,mask2edges),
             (addavedge,delavedge,maskavedge),
             (addacedge,delacedge,maskaCedge),
             (addaAedge,delaAedge,maskaAedge),
             (addapath,delapath,maskapath)]
    if c==[]:
        c = [ok2add2edges,
             ok2addavedge,
             ok2addacedge,
             ok2addaAedge,
             ok2addapath]

    adder1, remover1, masker1 = f[edge_function_idx(e1)]
    adder2, remover2, masker2 = f[edge_function_idx(e2)]

    checks_ok2 = c[edge_function_idx(e2)]
    checks_ok3 = c[edge_function_idx(e3)]

    d = {}
    s1 = set()
    s2 = set()
    s3 = set()

    for c1 in j1: # for each connector
        mask1 = adder1(g,e1,c1)
        append_set1 = False
        for c2 in j2:
            append_set2 = False
            if checks_ok2(e2,c2,g,g2):
                mask2 = adder2(g,e2,c2)
                for c3 in j3:
                    if checks_ok3(e3,c3,g,g2):
                        append_set1 = append_set2 = True
                        s3.add(c3)
                remover2(g,e2,c2,mask2)
            if append_set2: s2.add(c2)
        if append_set1: s1.add(c1)
        remover1(g,e1,c1,mask1)
    return s1, s2, s3

def del_empty(d):
    l = [e for e in d]
    for e in l:
        if d[e]==set(): del d[e]
    return d
def inorder_checks(g2, gg):
    #idx = np.argsort([len(gg[x]) for x in gg.keys()])
    #ee = [gg.keys()[i] for i in idx] # to preserve the order
    ee = [e for e in gg] # to preserve the order
    #cds = conformanceDS(g2, ee)
    #oo = new_order(g2, ee, repeats=100, cds=None)
    #ee = oo[0]
    random.shuffle(ee)
    d = {} # new datastructure
    d[ee[0]] = {('0'):gg[ee[0]]}
    for i in range(len(ee)-1):
        d[ee[i+1]] = del_empty(inorder_check2(ee[i], ee[i+1],
                                              gg[ee[i]], gg[ee[i+1]], g2)[0])
    return ee, d

def cloneempty(g): return {n:{} for n in g} # return a graph with no edges

def ok2addanedge1(s, e, g, g2,rate=1):
    """
    s - start,
    e - end
    """
    # directed edges
    # self-loop
    if s == e and not e in g2[s]: return False
    for u in g: # Pa(s) -> e
        if s in g[u] and not (e in g2[u] and (0,1) in g2[u][e]):
            return False
    for u in g[e]: # s -> Ch(e)
        if not (u in g2[s] and (0,1) in g2[s][u]):return False
    # bidirected edges
    for u in g[s]: # e <-> Ch(s)
        if u!=e and not (u in g2[e] and (2,0) in g2[e][u]):return False
    return True

def ok2addanedge2(s, e, g, g2, rate=1):
    mask = addanedge(g,(s,e))
    value = bfu.undersample(g,rate) == g2
    delanedge(g,(s,e),mask)
    return value

def ok2addanedge_sub(s, e, g, g2, rate=1):
    mask = addanedge(g,(s,e))
    value = isedgesubset(bfu.undersample(g,rate),g2)
    delanedge(g,(s,e),mask)
    return value

def ok2addanedge(s, e, g, g2, rate=1):
    f = [ok2addanedge1, ok2addanedge2]
    return f[min([1,rate-1])](s,e,g,g2,rate=rate)

def ok2addanedge_(s, e, g, g2, rate=1):
    f = [ok2addanedge1, ok2addanedge_sub]
    return f[min([1,rate-1])](s,e,g,g2,rate=rate)


def ok2add2edges(e,p,g,g2): return edge_increment_ok(e[0],p,e[1],g,g2)

def maskanedge(g,e): return [e[1] in g[e[0]]]
def mask2edges(g,e,p): return [p in g[e[0]], e[1] in g[p]]
def maskavedge(g,e,p):
    return [p[0] in g[e[0]], p[1] in g[e[0]],
            e[1] in g[p[0]], e[2] in g[p[1]]]
def maskaAedge(g,e,p):
    return [p[0] in g[e[1]], p[1] in g[e[2]],
            e[3] in g[p[0]], e[3] in g[p[1]]]
def maskaCedge(g,e,p):
    return [p[0] in g[e[1]], e[2] in g[p[0]],
            p[1] in g[e[2]], e[3] in g[p[1]]]
def maskapath(g,e,p):
    mask = []
    for i in range(len(p)):
        mask.append(p[i] in g[e[i+1]])
        mask.append(e[i+2] in g[p[i]])
    return mask
def maskaVpath(g,e,p):
    mask = []
    mask.extend([p[0] in g[e[0]], e[1] in g[p[-1]]])
    for i in range(1,len(p)):
        mask.append(p[i] in g[p[i-1]])
    return mask

def addanedge(g,e):
    '''
    add edge e[0] -> e[1] to g
    '''
    mask = maskanedge(g,e)
    g[e[0]][e[1]] =  set([(0,1)])
    return mask
def delanedge(g,e,mask):
    '''
    delete edge e[0] -> e[1] from g if it was not there before
    '''
    if not mask[0]: g[e[0]].pop(e[1], None)


def add2edges(g,e,p):
    '''
    break edge e[0] -> e[1] into two pieces
    e[0] -> p and p -> e[1]
    and add them to g
    '''
    mask = mask2edges(g,e,p)
    g[e[0]][p] = g[p][e[1]] = set([(0,1)])
    return mask

def del2edges(g,e,p,mask):
    '''
    restore the graph as it was before adding e[0]->p and p->e[1]
    '''
    if not mask[0]: g[e[0]].pop(p, None)
    if not mask[1]: g[p].pop(e[1], None)

def ok2addavedge(e,p,g,g2):
    if p[1] == e[0]:
        if p[0] != p[1] and p[0] != e[2] and not (e[2] in g2[p[0]] and (2,0) in g2[p[0]][e[2]]):
                return False
        if p[0] == p[1] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False
        if p[0] == e[1] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False

    if p[0] == e[0]:
        if p[0] != p[1] and p[1] != e[1] and not (e[1] in g2[p[1]] and (2,0) in g2[p[1]][e[1]]):
                return False
        if p[0] == p[1] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False
        if p[1] == e[2] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False

    if p[0] == e[1] and p[1] == e[2] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False
    if p[0] == e[2] and not (e[1] in g2[p[1]] and (0,1) in g2[p[1]][e[1]]):
            return False
    if p[1] == e[1] and not (e[2] in g2[p[0]] and (0,1) in g2[p[0]][e[2]]):
            return False
    if p[0] == p[1] == e[0] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False

    if not edge_increment_ok(e[0],p[0],e[1],g,g2): return False
    if not edge_increment_ok(e[0],p[1],e[2],g,g2): return False

    return  True

def addavedge(g,v,b):
    mask = maskavedge(g,v,b)
    g[v[0]][b[0]] = g[v[0]][b[1]] = g[b[0]][v[1]] = g[b[1]][v[2]] = set([(0,1)])
    return mask

def delavedge(g,v,b,mask):
    if not mask[0]: g[v[0]].pop(b[0], None)
    if not mask[1]: g[v[0]].pop(b[1], None)
    if not mask[2]: g[b[0]].pop(v[1], None)
    if not mask[3]: g[b[1]].pop(v[2], None)

def ok2addaAedge(e,p,g,g2):
    if p[1] == e[1] and not (p[0] in g2[e[2]] and (0,1) in g2[e[2]][p[0]]): return False
    if p[0] == e[2] and not (p[1] in g2[e[1]] and (0,1) in g2[e[1]][p[1]]): return False

    if not edge_increment_ok(e[1],p[0],e[3],g,g2): return False
    if not edge_increment_ok(e[2],p[1],e[3],g,g2): return False

    return True

def addaAedge(g,v,b):
    mask = maskaAedge(g,v,b)
    g[v[1]][b[0]] = g[v[2]][b[1]] = g[b[0]][v[3]] = g[b[1]][v[3]] = set([(0,1)])
    return mask

def delaAedge(g,v,b,mask):
    if not mask[0]: g[v[1]].pop(b[0], None)
    if not mask[1]: g[v[2]].pop(b[1], None)
    if not mask[2]: g[b[0]].pop(v[3], None)
    if not mask[3]: g[b[1]].pop(v[3], None)

def cleanedges(e,p,g, mask):
    i = 0
    for m in mask:
        if not m[0]: g[e[i+1]].pop(p[i], None)
        if not m[1]: g[p[i]].pop(e[i+2], None)
        i += 1
def cleanVedges(g, e,p, mask):

    if mask:
        if not mask[0]: g[e[0]].pop(p[0], None)
        if not mask[1]: g[p[-1]].pop(e[1], None)

        i = 0
        for m in mask[2:]:
            if not m: g[p[i]].pop(p[i+1], None)
            i += 1

def ok2addapath(e,p,g,g2):
    mask = []
    for i in range(len(p)):
        if not edge_increment_ok(e[i+1],p[i],e[i+2],g,g2):
            cleanedges(e,p,g,mask)
            return False
        mask.append(add2edges(g,(e[i+1],e[i+2]),p[i]))
    cleanedges(e,p,g,mask)
    return True

def ok2addaVpath(e,p,g,g2,rate=2):
    mask = addaVpath(g,e,p)
    if not isedgesubset(bfu.undersample(g,rate), g2):
        cleanVedges(g,e,p,mask)
        return False
    cleanVedges(g,e,p,mask)
    return True

def ok2addapath1(e,p,g,g2):
    for i in range(len(p)):
        if not edge_increment_ok(e[i+1],p[i],e[i+2],g,g2):
            return False
    return True

def addapath(g,v,b):

    mask = maskapath(g,v,b)
    s = set([(0,1)])
    for i in range(len(b)):
        g[v[i+1]][b[i]] = g[b[i]][v[i+2]] = s

    return mask

def addaVpath(g,v,b):
    mask = maskaVpath(g,v,b)
    s = set([(0,1)])
    l = [v[0]] + list(b) + [v[1]]
    for i in range(len(l)-1):
        g[l[i]][l[i+1]] = s
    return mask

def delaVpath(g, v, b, mask):
    cleanVedges(g, v, b, mask)

def delapath(g, v, b, mask):
    for i in range(len(b)):
        if not mask[2*i]: g[v[i+1]].pop(b[i], None)
        if not mask[2*i+1]:g[b[i]].pop(v[i+2], None)

def prunepaths_1D(g2, path, conn):
    c = []
    g = cloneempty(g2)
    for p in conn:
        mask = addapath(g,path,p)
        if isedgesubset(bfu.increment(g), g2): c.append(tuple(p))
        delapath(g,path,p,mask)
    return c

def ok2addacedge(e,p,g,g2):

    if p[0] == p[1]:
        if not e[2] in g2[e[2]]: return False
        if not p[0] in g2[p[0]]: return False
        if not (e[3] in g2[e[1]] and (0,1) in g2[e[1]][e[3]]): return False

    if not edge_increment_ok(e[1],p[0],e[2],g,g2): return False
    if not edge_increment_ok(e[2],p[1],e[3],g,g2): return False

    return True

def addacedge(g,v,b): # chain
    mask = maskaCedge(g,v,b)
    g[v[1]][b[0]] = g[v[2]][b[1]] = g[b[0]][v[2]] = g[b[1]][v[3]] = set([(0,1)])
    return mask

def delacedge(g,v,b,mask):
    if not mask[0]: g[v[1]].pop(b[0], None)
    if not mask[1]: g[b[0]].pop(v[2], None)
    if not mask[2]: g[v[2]].pop(b[1], None)
    if not mask[3]: g[b[1]].pop(v[3], None)

def rotate(l): return l[1:] + l[:1] # rotate a list
def density(g): return len(gk.edgelist(g))/np.double(len(g)**2)
def udensity(g): return (len(gk.edgelist(g))+len(gk.bedgelist(g))/2.)/np.double(len(g)**2 + len(g)*(len(g)-1)/2.)

def esig(l,n):
    '''
    turns edge list into a hash string
    '''
    z = len(str(n))
    n = map(lambda x: ''.join(map(lambda y: y.zfill(z),x)), l)
    n.sort()
    n = ''.join(n[::-1])
    return int('1'+n)

def gsig(g):
    '''
    turns input graph g into a hash string using edges
    '''
    return bfu.g2num(g)

def signature(g, edges): return (gsig(g),esig(edges,len(g)))

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def memo1(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = gsig(args[0])             # Signature: g
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def memo2(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return set()#cache[s]               # Return the cached solution
    return wrap

def eqsearch(g2, rate=1):
    '''Find  all  g  that are also in  the equivalence
    class with respect to g2 and the rate.
    '''

    s = set()
    noop = set()

    @memo1
    def addnodes(g,g2,edges):
        if edges:
            masks  = []
            for e in edges:
                if ok2addanedge_(e[0],e[1],g,g2,rate=rate):
                    masks.append(True)
                else:
                    masks.append(False)
            nedges = [edges[i] for i in range(len(edges)) if masks[i]]
            n = len(nedges)
            if n:
                for i in range(n):
                    mask = addanedge(g,nedges[i])
                    if bfu.undersample(g,rate) == g2: s.add(bfu.g2num(g))
                    addnodes(g,g2,nedges[:i]+nedges[i+1:])
                    delanedge(g,nedges[i],mask)
                return s
            else:
                return noop
        else:
            return noop


    g = cloneempty(g2)
    edges = gk.edgelist(gk.complement(g))
    addnodes(g,g2,edges)
    return s


def supergraphs_in_eq(g, g2, rate=1):
    '''Find  all supergraphs of g  that are also in  the same equivalence
    class with respect to g2 and the rate.
    Currently works only for bfu.undersample by 1
    '''
    if bfu.undersample(g,rate) != g2:
        raise ValueError('g is not in equivalence class of g2')

    s = set()

    def addnodes(g,g2,edges):
        if edges:
            masks  = []
            for e in edges:
                if ok2addanedge(e[0],e[1],g,g2,rate=rate):
                    masks.append(True)
                else:
                    masks.append(False)
            nedges = [edges[i] for i in range(len(edges)) if masks[i]]
            n = len(nedges)
            if n:
                for i in range(n):
                    mask = addanedge(g,nedges[i])
                    s.add(bfu.g2num(g))
                    addnodes(g,g2,nedges[:i]+nedges[i+1:])
                    delanedge(g,nedges[i],mask)

    edges = gk.edgelist(gk.complement(g))
    addnodes(g,g2,edges)
    return s

def edge_backtrack2g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    single_cache = {}

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            e = edges.pop()
            ln = [n for n in g2]
            for n in ln:
                if (n,e) in single_cache: continue
                mask = add2edges(g,e,n)
                if isedgesubset(bfu.increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and bfu.increment(r)==g2:
                        s.add(bfu.g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements in eqclass')
                del2edges(g,e,n,mask)
            edges.append(e)
        else:
            return g
    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:
            mask = add2edges(g,e,n)
            if not isedgesubset(bfu.increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def edge_backtrack2g1_directed(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    single_cache = {}

    def edgeset(g):
        return set(gk.edgelist(g))
    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            e = edges.pop()
            ln = [n for n in g2]
            for n in ln:
                if (n,e) in single_cache: continue
                mask = add2edges(g,e,n)
                if isedgesubsetD(bfu.increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and edgeset(bfu.increment(r))==edgeset(g2):
                        s.add(bfu.g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements in eqclass')
                del2edges(g,e,n,mask)
            edges.append(e)
        else:
            return g
    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:
            mask = add2edges(g,e,n)
            if not isedgesubsetD(bfu.increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def g22g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    single_cache = {}

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            if bfu.increment(g) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                return g
            e = edges[0]
            for n in g2:

                if (n,e) in single_cache: continue
                if not edge_increment_ok(e[0],n,e[1],g,g2): continue

                mask = add2edges(g,e,n)
                r = nodesearch(g,g2,edges[1:],s)
                del2edges(g,e,n,mask)

        elif bfu.increment(g)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements in eqclass')
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:

            mask = add2edges(g,e,n)
            if not isedgesubset(bfu.increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def backtrack_more(g2, rate=1, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    single_cache = {}
    if rate == 1:
        ln = [n for n in g2]
    else:
        ln = []
        for x in itertools.combinations_with_replacement(g2.keys(),rate):
            ln.extend(itertools.permutations(x,rate))
        ln = set(ln)

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            if bfu.undersample(g,rate) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                return g
            e = edges[0]
            for n in ln:

                if (n,e) in single_cache: continue
                if not ok2addaVpath(e, n, g, g2, rate=rate): continue

                mask = addaVpath(g,e,n)
                r = nodesearch(g,g2,edges[1:],s)
                delaVpath(g,e,n,mask)

        elif bfu.undersample(g,rate)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements in eqclass')
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in ln:

            mask = addaVpath(g,e,n)
            if not isedgesubset(bfu.undersample(g,rate), g2):
                single_cache[(n,e)] = False
            delaVpath(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def backtrackup2u(H,umax=2):
    s = set()
    for i in xrange(1,umax+1):
        s = s | backtrack_more(H,rate=i)
    return s

def vg22g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge),
         (addaAedge,delaAedge),
         (addapath,delapath)]
    c = [ok2add2edges,
         ok2addavedge,
         ok2addacedge,
         ok2addaAedge,
         ok2addapath]
    @memo2 # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            #key, checklist = edges.popitem()
            key = random.choice(edges.keys())
            checklist = edges.pop(key)
            adder, remover = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]
            for n in checklist:
                mask = adder(g,key,n)
                if isedgesubset(bfu.increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and bfu.increment(r)==g2:
                        s.add(bfu.g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements')
                remover(g,key,n,mask)
            edges[key] = checklist
        else:
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    chlist = checkable(g2)
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g,g2,chlist,s)
    except ValueError:
        s.add(0)
    return s

def edge_function_idx(edge):
    return min(4,len(edge))-2+min(max(3,len(edge))-3,1)*int(edge[0])

def v2g22g1(g2, capsize=None, verbose=True):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    f = [(add2edges,del2edges,mask2edges),
         (addavedge,delavedge,maskavedge),
         (addacedge,delacedge,maskaCedge),
         (addaAedge,delaAedge,maskaAedge),
         (addapath,delapath,maskapath)]
    c = [ok2add2edges,
         ok2addavedge,
         ok2addacedge,
         ok2addaAedge,
         ok2addapath]

    def predictive_check(g,g2,pool,checks_ok, key):
        s = set()
        for u in pool:
            if not checks_ok(key,u,g,g2): continue
            s.add(u)
        return s

    @memo2 # memoize the search
    def nodesearch(g, g2, order, inlist, s, cds, pool, pc):
        if order:
            if bfu.increment(g) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                s.update(supergraphs_in_eq(g, g2))
                return g

            key = order[0]
            if pc:
                tocheck = [x for x in pc if x in cds[len(inlist)-1][inlist[0]]]
            else:
                tocheck = cds[len(inlist)-1][inlist[0]]

            if len(order) > 1:
                kk = order[1]
                pc = predictive_check(g,g2,pool[len(inlist)],
                                      c[edge_function_idx(kk)],kk)
            else:
                pc = set()

            adder, remover, masker = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]

            for n in tocheck:
                if not checks_ok(key,n,g,g2): continue
                masked = np.prod(masker(g,key,n))
                if masked:
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                else:
                    mask = adder(g,key,n)
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                    remover(g,key,n,mask)

        elif bfu.increment(g)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements')
            return g

    @memo2 # memoize the search
    def nodesearch0(g, g2, order, inlist, s, cds):

        if order:
            key = order.pop(0)
            tocheck = cds[len(inlist)-1][inlist[0]]

            adder, remover, masker = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]

            if len(tocheck) > 1:
                for n in tocheck:
                    if not checks_ok(key,n,g,g2): continue
                    mask = masker(g,key,n)
                    if not np.prod(mask):
                        mask = adder(g,key,n)
                        r = nodesearch0(g,g2,order, [n]+inlist, s, cds)
                        if r and bfu.increment(r)==g2:
                            s.add(bfu.g2num(r))
                            if capsize and len(s)>capsize:
                                raise ValueError('Too many elements')
                        remover(g,key,n,mask)
                    else:
                        r = nodesearch0(g,g2,order, [n]+inlist, s, cds)
                        if r and bfu.increment(r)==g2:
                            s.add(bfu.g2num(r))
                            if capsize and len(s)>capsize:
                                raise ValueError('Too many elements')
            elif tocheck:
                (n,) = tocheck
                mask = adder(g,key,n)
                r = nodesearch0(g,g2, order, [n]+inlist, s, cds)
                if r and bfu.increment(r) == g2:
                    s.add(bfu.g2num(r))
                    if capsize and len(s)>capsize:
                        raise ValueError('Too many elements')
                remover(g,key,n,mask)

            order.insert(0,key)

        else:
            return g

    # find all directed g1's not conflicting with g2

    startTime = int(round(time.time() * 1000))
    gg = checkable(g2)

    idx = np.argsort([len(gg[x]) for x in gg.keys()])
    keys = [gg.keys()[i] for i in idx]

    cds, order, idx = conformanceDS(g2, gg, keys)
    endTime = int(round(time.time() * 1000))
    if verbose:
        print "precomputed in {:10} seconds".format(round((endTime-startTime)/1000.,3))
    if 0 in [len(x) for x in order]:
        return set()
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g, g2, [keys[i] for i in idx], ['0'], s, cds, order, set())
        #nodesearch0(g, g2, [gg.keys()[i] for i in idx], ['0'], s, cds)
    except ValueError, e:
        print e
        s.add(0)
    return s

def backtrack_more2(g2, rate=2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    f = [(addaVpath,delaVpath,maskaVpath)]
    c = [ok2addaVpath]

    def predictive_check(g,g2,pool,checks_ok, key):
        s = set()
        for u in pool:
            if not checks_ok(key,u,g,g2,rate=rate): continue
            s.add(u)
        return s

    @memo2 # memoize the search
    def nodesearch(g, g2, order, inlist, s, cds, pool, pc):
        if order:
            if bfu.undersample(g,rate) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                s.update(supergraphs_in_eq(g, g2, rate=rate))
                return g

            key = order[0]
            if pc:
                tocheck = [x for x in pc if x in cds[len(inlist)-1][inlist[0]]]
            else:
                tocheck = cds[len(inlist)-1][inlist[0]]

            if len(order) > 1:
                kk = order[1]
                pc = predictive_check(g,g2,pool[len(inlist)],
                                      c[edge_function_idx(kk)],kk)
            else:
                pc = set()

            adder, remover, masker = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]

            for n in tocheck:
                if not checks_ok(key,n,g,g2,rate=rate): continue
                masked = np.prod(masker(g,key,n))
                if masked:
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                else:
                    mask = adder(g,key,n)
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                    remover(g,key,n,mask)

        elif bfu.undersample(g,rate)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements')
            return g

    # find all directed g1's not conflicting with g2

    startTime = int(round(time.time() * 1000))
    ln = [x for x in itertools.permutations(g2.keys(),rate)] + \
         [(n,n) for n in g2]
    gg = {x:ln for x in gk.edgelist(g2)}
    keys = gg.keys()
    cds, order, idx = conformanceDS(g2, gg, gg.keys(), f=f, c=c)
    endTime = int(round(time.time() * 1000))
    print "precomputed in {:10} seconds".format(round((endTime-startTime)/1000.,3))
    if 0 in [len(x) for x in order]:
        return set()
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g, g2, [keys[i] for i in idx], ['0'], s, cds, order, set())
    except ValueError, e:
        print e
        s.add(0)
    return s


def unionpool(idx, cds):
    s = set()
    for u in cds[idx]:
        for v in cds[idx][u]:
            s = s.union(cds[idx][u][v])
    return s

def conformanceDS(g2, gg, order, f=[], c=[]):
    CDS = {}
    pool = {}

    CDS[0] = set(gg[order[0]])
    pool = [set(gg[order[i]]) for i in range(len(order))]

    for x in itertools.combinations(range(len(order)),2):

        d, s_i1, s_i2 = inorder_check2(order[x[0]], order[x[1]],
                                       pool[x[0]], pool[x[1]],
                                       g2, f=f, c=c)

        pool[x[0]] = pool[x[0]].intersection(s_i1)
        pool[x[1]] = pool[x[1]].intersection(s_i2)

        d = del_empty(d)
        if not x[1] in CDS:
            CDS[x[1]] = {}
            CDS[x[1]][x[0]] = d
        else:
            CDS[x[1]][x[0]] = d
    if density(g2) > 0.35:
        itr3 = [x for x in itertools.combinations(range(len(order)),3)]
        for x in random.sample(itr3, min(10,np.int(scipy.misc.comb(len(order),3)))):
            s1, s2, s3 = check3(order[x[0]], order[x[1]], order[x[2]],
                                pool[x[0]], pool[x[1]], pool[x[2]],
                                g2, f=f, c=c)

            pool[x[0]] = pool[x[0]].intersection(s1)
            pool[x[1]] = pool[x[1]].intersection(s2)
            pool[x[2]] = pool[x[2]].intersection(s3)

    return prune_sort_CDS(CDS, pool)

def prune_modify_CDS(cds, pool):
    ds = {}
    ds[0]={}
    ds[0]['0'] = pool[0]
    for i in range(1,len(pool)):
        ds[i] = {}
        for j in cds[i].keys():
            for e in pool[i-1].intersection(cds[i][j].keys()):
                ds[i][e] = pool[i].intersection(cds[i][j][e])
    return ds, pool, range(len(pool))

def prune_sort_CDS(cds, pool):

    idx = np.argsort([len(x) for x in pool])
    p = [pool[i] for i in idx]

    ds = {}
    ds[0]={}
    ds[0]['0'] = pool[idx[0]]

    for i in range(1,len(idx)):
        ds[i] = {}
        for j in range(i):
            if idx[j] > idx[i]:
                dd = invertCDSelement(cds[idx[j]][idx[i]])
            else:
                dd = cds[idx[i]][idx[j]]
            for e in pool[idx[i-1]].intersection(dd.keys()):
                ds[i][e] = pool[idx[i]].intersection(dd[e])

    return ds, p, idx

def invertCDSelement(d_i):
    d = {}
    for e in d_i:
        for v in d_i[e]:
            if v in d:
                d[v].add(e)
            else:
                d[v]=set([e])
    return d

def conformant(cds, inlist):

    if inlist[len(inlist)-2] in cds[len(inlist)-1][0]:
        s = cds[len(inlist)-1][0][inlist[len(inlist)-2]]
    else:
        return set()
    for i in range(1,len(inlist)-1):
        if inlist[len(inlist)-i-2] in cds[len(inlist)-1][i]:
            s = s.intersection(cds[len(inlist)-1][i][inlist[len(inlist)-i-2]])
        else:
            return set()
    return s

def length_d_paths(G, s, d):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    def recurse(G, s, d, path=[]):
        if d == 0:
            yield path
            return

        for u in G[s]:
            if G[s][u] == set([(2,0)]) or u in path: continue
            for v in recurse(G, u, d-1, path+[u]):
                yield v

    for u in recurse(G, s, d, [s]):
            yield u

def edge_increment_ok(s,m,e,g,g2):
    """
    s - start,
    m - middle,
    e - end
    """
    # bidirected edges
    for u in g[s]:
        if u!=m and not (m in g2[u] and (2,0) in g2[u][m]):return False

    # directed edges
    if s == e:
        if not (m in g2[m] and (0,1) in g2[m][m]):return False
        if not (s in g2[s] and (0,1) in g2[s][s]):return False
    for u in g[m]:
        if not (u in g2[s] and (0,1) in g2[s][u]):return False
        # bidirected edges
        if u!=e and not (e in g2[u] and (2,0) in g2[u][e]):return False
    for u in g[e]:
        if not (u in g2[m] and (0,1) in g2[m][u]):return False

    for u in g:
        if s in g[u] and not (m in g2[u] and (0,1) in g2[u][m]):
            return False
        if m in g[u] and not (e in g2[u] and (0,1) in g2[u][e]):
            return False

    return True


def length_d_loopy_paths(G, s, dt, p):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    g1 = cloneempty(G)

    def recurse(g, g2, s, d, path=[]):

        if edge_increment_ok(p[-d-2],s,p[-d-1],g,g2):

            if d == 0:
                yield path
                return

            mask = add2edges(g,(p[-d-2],p[-d-1]),s)
            for u in g2[s]:
                if g2[s][u] == set([(2,0)]): continue
                for v in recurse(g, g2, u, d-1, path+[u]):
                    yield v
            del2edges(g,(p[-d-2],p[-d-1]),s, mask)

    for u in recurse(g1, G, s, dt-1, [s]):
            yield u
# system packages
import StringIO
import scipy
import sys
import os
import igraph
import numpy as np

# local packages
import zickle
import ecj
from ecj import undersample
import testgraphs
reload(testgraphs)
from testgraphs import *

colors = zickle.load('data/colors.zkl') # load the colors

def graph2dict(g):
    D = {}
    for v in range(0,len(g.vs)):
        D[g.vs[v]["label"]] = {}
        for u in g.neighbors(v,mode="OUT"):
            D[g.vs[v]["label"]][g.vs[u]["label"]] = set([(0,1)])
    return D

def paintSCC(g, cm):
    nameidx = {}
    D = graph2dict(g)
    scc = ecj.scc(D)
    for i in range(0,len(scc)):
        for v in scc[i]:
            g.vs[g.vs['label'].index(v)]["color"] = cm[i]
fstr = "%.5f" # precision to use
def dict2graph(D):
    A = scipy.zeros([len(D),len(D)])
    nodes = []
    indices = {}
    c = 0
    for v in D:
        nodes.append(v)
    idx = np.argsort([int(v) for v in nodes])
    nodes = [nodes[i] for i in idx]
    for v in nodes:
        indices[v] = c
        c +=1
    for v in nodes:
        for u in D[v]:
            if  not ((2,0) in D[v][u]):
                A[indices[v],indices[u]] = 1
    g = igraph.Graph.Adjacency(A.tolist())
    g.vs["label"] = nodes
    return g

def unroll(G,steps):
    N = {}
    for i in range(0,steps):
        N.update({v+str(i):set([u+str(i+1) for u in G[v]]) for v in G})
    N.update({v+str(steps):set() for v in G})
    return N

def unroll_undersample(G,steps):
    # does not provide isochronal bidirectional edges
    N = {}
    steps += 2
    U = unroll(G,steps)
    nodes = G.keys()
    for v in G:
        N.update({v:set([nodes[k] for k in scipy.where([ecj.reachable(v+'0',U,u+str(steps-1)) for u in G])[0]])})
    return N

def matrix_start(mname='generic',w_gap=0.45,h_gap=0.5,stl=''):
    print "\matrix ("+mname+") [matrix of nodes, row sep="\
        +str(h_gap)+"cm,column sep="+str(w_gap)+"cm"+stl+"]"
    print "{"
def matrix_end():
    print "};"

def matrix_grid(G,s,mname='generic',w_gap=0.5, h_gap=0.5, type="obs",stl=''):
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    keylist = G.keys();
    idx = np.argsort([int(v) for v in keylist])
    keylist = [keylist[i] for i in idx]
    for v in keylist:
        print " ","& ".join(["\\node["+type+"]{"+v+"};" for i in range(1,s)])+"\\\\"
    matrix_end()

def matrix_edges(G,s,mname='generic'):
    nodes = G.keys()
    idx = np.argsort([int(v) for v in nodes])
    nodes = [nodes[i] for i in idx]    
    for v in nodes:
        idx1 = nodes.index(v)+1
        for u in G[v]:
            if G[v][u].intersection([(edge_type['directed'],1)]):
                idx2 = nodes.index(u)+1
                print '\\foreach \\x in{'+','.join(map(str,range(1,s-1)))+'}{'
                print '  \\pgfmathtruncatemacro{\\xn}{\\x+1}'
                print '  \\draw[pil] ('+mname+'-'+str(idx1)+'-\\x) -- ('\
                    +mname+'-'+str(idx2)+'-\\xn);'
                print '};'
def emacs_vars():
    print '%%% Local Variables:'
    print '%%% mode: latex'
    print '%%% TeX-master: "../master"'
    print '%%% End:'

def dbnprint(G,s,mname='generic',w_gap=0.5, h_gap=0.5,type="obs",stl=''):
    matrix_grid(G,s,mname=mname,w_gap=w_gap, h_gap=h_gap, type=type,stl=stl)
    matrix_edges(G,s,mname=mname)

from scipy import array,cos,sin,deg2rad,rad2deg
import math

def getangle(A,B):
    """
    When A and  B are two angles around the clock  returns an angle of
    the line that is connecting them.
    """
    x = array([cos(deg2rad(A)),sin(deg2rad(A))])
    y=array([cos(deg2rad(B)),sin(deg2rad(B))])
    d = y-x
    return rad2deg(math.atan2(d[1],d[0]))

def cdbnprint(G,mtype="obs",bend=5,curve=5,R=1):
    """
    Prints  out  a  compressed  dbn  repesentation  of  the  graph  in
    TikZ/Latex format
    """
    output = StringIO.StringIO()
    BE = set()
    n = len(G)
    nodes = G.keys()
    idx = np.argsort([int(v) for v in nodes])
    nodes = [nodes[i] for i in idx]

    g = dict2graph(ecj.cloneBfree(G))
    paintSCC(g,colors)


    for i in range(0,n):
        node = g.vs[i]['label']
        rc = g.vs[i]["color"]
        print >>output, "{ \\definecolor{mycolor}{RGB}{"\
            +str(rc[0])+","+str(rc[1])+","+str(rc[2])+"}"
        mcolor = "fill = {rgb: red,"+str(rc[0])+"; green,"+str(rc[1])+\
            "; blue,"+str(rc[2])+"}"
        print >>output, "\\node["+mtype+", fill=mycolor] ("+node+") at ("+\
            str(-i*360/n+180)+":"+str(R)+") {"+node+"};}"

#    print >>output,"\\foreach \\name/\\angle in {"+",".join(
#        [nodes[i]+"/"+str(-i*360/n+180) for i in range(0,n)])+"}"
#    print >>output,"\\node["+mtype+"] (\\name) at (\\angle:"\
#        +str(R)+") {\\name};"

    for i in range(0,n):
        v = nodes[i]
        ll=[v+'/'+u for u in G[v]]
        for l in ll:
            a,b = l.split('/')
            if G[a][b].intersection([(edge_type['bidirected'],0)]):
                if not(BE.intersection([(a,b)]) or BE.intersection([(b,a)])):
                    ang_a = -nodes.index(a)*360/n+180
                    ang_b = -nodes.index(b)*360/n+180
                    print >>output,'  \\draw[pilip, on layer=back] ('+a+') -- ('+b+');'
            if G[a][b].intersection([(edge_type['directed'],1)]):
                ang_a = -nodes.index(a)*360/n+180
                ang_b = -nodes.index(b)*360/n+180
                if a == b:
                    print >>output,"\\path[overlay,draw,pil] ("+a+")" +\
                        " .. controls +("+ "%.5f" % (bend+ang_a) +\
                        ":"+ fstr % (2*curve)+"mm) and +("+\
                        "%.5f" % (ang_a-bend)+\
                        ":"+"%.5f" % (2*curve)+"mm) .. ("+b+");"
                else:
                    print >>output,"\\path[overlay,draw,pil] ("+a+")" +\
                        " .. controls +("+"%.5f" % (bend+getangle(ang_a,ang_b)) +\
                        ":"+fstr % (curve)+"mm) and +("+\
                        fstr % (getangle(ang_b,ang_a)-bend)+\
                        ":"+fstr % (curve)+"mm) .. ("+b+");"
    return output

def gprint(G,mtype="obs",bend=5,curve=5,R=1,layout=None, scale=5):
    """
    Prints out an automatically layout compressed dbn repesentation of
    the graph in TikZ/Latex format
    """
    output = StringIO.StringIO()
    BE = set()
    n = len(G)
    if not layout:
        g = dict2graph(ecj.cloneBfree(G))
        layout = g.layout_fruchterman_reingold(maxiter=50000, coolexp=1.1)
        #layout = g.layout_graphopt(niter=50000, node_charge=0.08)
        layout.center([0,0])
        layout.scale(float(1/scipy.absolute(layout.coords).max()))
        layout.scale(R)
        cc = scipy.round_(array(layout.coords),decimals=4)
    else:
        g = dict2graph(ecj.cloneBfree(G))
        cc = array(layout.coords)
    paintSCC(g,colors)
    for i in range(0,n):
        node = g.vs[i]['label']
        rc = g.vs[i]["color"]
        print >>output, "{ \\definecolor{mycolor}{RGB}{"\
            +str(rc[0])+","+str(rc[1])+","+str(rc[2])+"}"
        mcolor = "fill = {rgb: red,"+str(rc[0])+"; green,"+str(rc[1])+\
            "; blue,"+str(rc[2])+"}"
        print >>output, "\\node["+mtype+", fill=mycolor] ("+node+") at ("+\
            str(cc[i][0])+","+str(cc[i][1])+") {"+node+"};}"

    for i in range(0,n):
        v = g.vs[i]['label']
        ll=[v+'/'+u for u in G[v]]
        for l in ll:
            a,b = l.split('/')
            if G[a][b].intersection([(edge_type['bidirected'],0)]):
                if not(BE.intersection([(a,b)]) or BE.intersection([(b,a)])):
                    print >>output,'  \\draw[pilip, on layer=back] ('+\
                        a+') -- ('+b+');'
            if G[a][b].intersection([(edge_type['directed'],1)]):
                if a == b:
                    dff = cc[g.vs['label'].index(a)] - scipy.mean(cc,0)
                    ang = scipy.arctan2(dff[1],dff[0])
                    ang_a = scipy.rad2deg(ang)
                    print >>output,"\\path[overlay,draw,pil] ("+a+")" +\
                        " .. controls +("+ "%.5f" % (bend+ang_a) +\
                        ":"+ fstr % (2*curve)+"mm) and +("+\
                        "%.5f" % (ang_a-bend)+\
                        ":"+"%.5f" % (2*curve)+"mm) .. ("+b+");"
                else:
                    dff = cc[g.vs['label'].index(b)] \
                        - cc[g.vs['label'].index(a)]
                    ang = scipy.arctan2(dff[1],dff[0])
                    ang_a = scipy.rad2deg(ang)
                    ang_b = ang_a+180
                    print >>output,"\\path[overlay,draw,pil] ("+a+")" +\
                        " .. controls +("+\
                        "%.5f" % (bend+ang_a) +\
                        ":"+fstr % (curve)+"mm) and +("+\
                        fstr % (ang_b-bend)+\
                        ":"+fstr % (curve)+"mm) .. ("+b+");"
    return output

def cdbnwrap(G,u,name='AAA',R=1,gap=0.5):
    output = StringIO.StringIO()
    print >>output,"\\node[right="+str(gap)+"cm of "+name+str(u-1)\
        +",scale=0.7]("+name+str(u)+"){"
    print >>output,"\\begin{tikzpicture}"
    s = cdbnprint(undersample(G,u),mtype='lahid',bend=25,curve=10,R=R)
    print >>output,s.getvalue()
    s.close()
    print >>output,"\\end{tikzpicture}"
    print >>output,"};"
    return output

def cdbnsingle(g,scale=0.7,R=1,gap=0.5,mtype="lahid"):
    output = StringIO.StringIO()
    print >>output,"\\node[scale="+str(scale)+"](){"
    print >>output,"\\begin{tikzpicture}"
    s = cdbnprint(g,mtype=mtype,bend=25,curve=10,R=R)
    print >>output,s.getvalue()
    s.close()
    print >>output,"\\end{tikzpicture}"
    print >>output,"};"
    return output

def cdbn_single(G,u,scale=0.7,R=1,gap=0.5,mtype="lahid"):
    return cdbnsingle(undersample(G,u),scale=scale,R=R,gap=gap,mtype=mtype)

def gsingle(g,scale=0.7,R=1,gap=0.5,mtype="lahid",layout=None):
    output = StringIO.StringIO()
    print >>output,"\\node[scale="+str(scale)+"](){"
    print >>output,"\\begin{tikzpicture}"
    s = gprint(g,mtype=mtype,bend=25,curve=6,R=R,layout=layout)
    print >>output,s.getvalue()
    s.close()
    print >>output,"\\end{tikzpicture}"
    print >>output,"};"
    return output

def g_single(G,u,scale=0.7,R=1,gap=0.5,mtype="lahid",layout=None):
    g = undersample(G,u)
    return gsingle(g,scale=scale,R=R,gap=gap,mtype=mtype,layout=layout)

def unfoldplot(G,steps=7,repeats=5,gap=0.5,R=1,hg=0.1,wgap=0.7,name='AAA',stl=''):
    u = 0
    dbnprint(undersample(G,u), repeats, w_gap=wgap,
             h_gap=hg, mname=name+str(u), type='hid', stl=stl)
    print "\\node[left="+str(gap)+"cm of "+name+str(u)+",scale=0.7] (C) {"
    print "\\begin{tikzpicture}"
    cdbnprint(G,mtype='hid',bend=15,curve=5,R=R)
    print "\\end{tikzpicture}"
    print "};"
    for u in range(1,steps):
        dbnprint(undersample(G,u),repeats,w_gap=wgap,h_gap=hg,mname=name+\
                     str(u),type='ahid',stl=', below=0.25cm of '+name+str(u-1))

        print "\\node[left="+str(gap)+"cm of "+name+str(u)+",scale=0.7] () {"
        print "\\begin{tikzpicture}"
        cdbnprint(undersample(G,u),mtype='lahid',bend=15,curve=5,R=R)
        print "\\end{tikzpicture}"
        print "};"
    emacs_vars()

def foldplot(G,steps=7,gap=0.5,R=1,hg=0.1,name='AAA',stl=''):
    u = 0
    print '\\node[scale=0.7'+stl+'] ('+name+str(u)+') {'
    print "\\begin{tikzpicture}"
    s=cdbnprint(G,mtype='hid',bend=25,curve=10,R=R)
    print s.getvalue()
    s.close()
    print "\\end{tikzpicture}"
    print "};"
    for u in range(1,steps):
        s=cdbnwrap(G,u,R=R,name=name)
        print s.getvalue()
        s.close()
    emacs_vars()

def matrix_fold(G,m,n,R=1,mname='g',w_gap=0.5, h_gap=0.5, stl=''):
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    for i in range(0,n):
        S = []
        for j in range(0,m):
            u = m*i+j
            if u == 0:
                s=cdbn_single(G,u,R=R,mtype="hid")
            else:
                s=cdbn_single(G,u,R=R)
            S.append(s.getvalue())
            s.close()
        print " ","& ".join(S)+"\\\\"
    matrix_end()

def gmatrix_fold(G,m,n,R=1,mname='g',w_gap=0.5, h_gap=0.5, stl='', shift=0):
    #g = dict2graph(ecj.cloneBfree(G))
    #layout = g.layout_graphopt(niter=50000, node_charge=0.02)
    #layout.center([0,0])
    #layout.scale(0.01)
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    for i in range(0,n):
        S = []
        for j in range(0,m):
            u = m*i+j+shift
            s=g_single(G,u,R=R,mtype="hid")
            S.append(s.getvalue())
            s.close()
        print " ","& ".join(S)+"\\\\"
    matrix_end()

def addselfedge(G,v):
    G[v].update({v:set([(edge_type['directed'],1)])})

def gmatrix_list(l,m,n,R=1,mname='g',w_gap=0.5,h_gap=0.5,stl='',shift=0):
    """
    Given a list of graphs prints them out as latex laying out each graph in a force-based layout
    """
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    for i in range(0,n):
        S = []
        for j in range(0,m):
            u = m*i+j+shift
            if u >len(l)-1: break
            s=gsingle(l[u],R=R,mtype="hid")
            S.append(s.getvalue())
            s.close()
        print " ","& ".join(S)+"\\\\"
    matrix_end()

def matrix_list(l,m,n,R=1,mname='g',w_gap=0.5,h_gap=0.5,stl='',shift=0):
    """
    Given a list of graphs prints them out as latex laying out each in a circle layout
    """
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    for i in range(0,n):
        S = []
        for j in range(0,m):
            u = m*i+j+shift
            if u >len(l)-1: break
            s=cdbnsingle(l[u],R=R,mtype="hid")
            S.append(s.getvalue())
            s.close()
        print " ","& ".join(S)+"\\\\"
    matrix_end()

class WritableObject:
    def __init__(self):
        self.content = []
    def write(self, string):
        self.content.append(string)
N = {
'a': set('bcd'),
'b': set('ce'),
'c': set('d'),
'd': set('e'),
'e': set('f'),
'f': set('cgh'),
'g': set('fh'),
'h': set('fg')
}

N2 = {
'a': set('bcf'),
'b': set('c'),
'c': set('de'),
'd': set('e'),
'e': set(),
'f': set(),
}

N4 = {
'a': set('bdc'),
'b': set('d'),
'c': set('d'),
'd': set('efg'),
'e': set('g'),
'f': set('g'),
'g': set(),
}

"""
Each  edge is a  tuple (t,w),  where t  signifies its  type and  w its
integer time length:
0 - an isochronal edge
1 - from slice t to t+1 (the regular case)
2 - from t to t+2
etc
"""
edge_type = {'directed':0,'undirected':1,'bidirected':2}
iedge_type = {0:'directed',1:'undirected',2:'bidirected'}

LG = { # loop graph
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'e':set([(0,1)])},
'e': {'a':set([(0,1)])}
}

# 2 cycles
C = {
    'a': {'b':set([(0,1)]),'d':set([(0,1)])},
    'b': {'c':set([(0,1)])},
    'c': {'d':set([(0,1)])},
    'd': {'a':set([(0,1)])}
}

# a cycle of two
N = {
'a': {'b':set([(0,1)])},
'b': {'a':set([(0,1)]),'d':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'c':set([(0,1)])},
}

N1 = {
'a': {'d':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'b':set([(0,1)])},
'd': {'b':set([(0,1)]),'a':set([(0,1)])},
}

# a cycle of three
A = {
'a': {'b':set([(0,1)]),'e':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'a':set([(0,1)])},
'd': {'a':set([(0,1)])},
'e': {'d':set([(0,1)])},
}

# no cycles
U = {
'a': {'d':set([(0,1)]),'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {},
'd': {'c':set([(0,1)])},
}

# DAG
D = {
'a': {'b':set([(0,1)]),'c':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {}
}

D2 = {
'a': {},
'b': {},
'c': {'d':set([(0,1)]),'h':set([(0,1)]),'a':set([(0,1)]),'b':set([(0,1)])},
'd': {},
'e': {'d':set([(0,1)]),'f':set([(0,1)])},
'f': {'g':set([(0,1)]),'h':set([(0,1)])},
'g': {'c':set([(0,1)])},
'h': {},
}

# 3 SCC DAG
D3 = {
'a': {'b':set([(0,1)]),'c':set([(0,1)])},
'b': {'e':set([(0,1)]),'i':set([(0,1)]),'d':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'a':set([(0,1)])},
'e': {'f':set([(0,1)])},
'f': {'g':set([(0,1)])},
'g': {'e':set([(0,1)]),'h':set([(0,1)])},
'h': {'i':set([(0,1)])},
'i': {'h':set([(0,1)])},
}

P = {
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'a':set([(0,1)]),'b':set([(0,1)])},
}

L = {
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'e':set([(0,1)])},
'e': {'f':set([(0,1)]),'a':set([(0,1)])},
'f': {'a':set([(0,1)])},
'g': {'h':set([0,1])},
'h': {'i':set([0,1])},
'i': {'g':set([0,1])},
}

UL = {
'1': {'3':set([(0,1)]),'9':set([(0,1)])},
'3': {'4':set([(0,1)])},
'4': {'5':set([(0,1)])},
'5': {'1':set([(0,1)]),'6':set([(0,1)])},
'6': {'5':set([(0,1)]),'7':set([(0,1)])},
'7': {'8':set([(0,1)])},
'8': {'6':set([(0,1)])},
'9': {'a':set([(0,1)])},
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'1':set([(0,1)])},
}

UL11 = {
'1': {'2':set([(0,1)]),'3':set([(0,1)])},
'2': {'4':set([(0,1)])},
'3': {'4':set([(0,1)]),'3':set([(0,1)])},
'4': {},
}

# tree
T2 = {
'a': {'b':set([(0,1)]),'o':set([(0,1)])},
'b': {'f':set([(0,1)]),'c':set([(0,1)])},
'c': {'d':set([(0,1)]),'e':set([(0,1)])},
'd': {},
'e': {},
'f': {'g':set([(0,1)]),'h':set([(0,1)])},
'g': {},
'h': {},
'i': {},
'j': {},
'k': {'i':set([(0,1)]),'j':set([(0,1)])},
'l': {},
'm': {},
'n': {'m':set([(0,1)]),'l':set([(0,1)])},
'o': {'n':set([(0,1)]),'k':set([(0,1)])},
}

# one loop
L1 = {
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'e':set([(0,1)])},
'e': {'a':set([(0,1)])},
}

L25 = {
'1': {'2':set([(0,1)]),'1':set([(0,1)])},
'2': {'3':set([(0,1)])},
'3': {'5':set([(0,1)]),'6':set([(0,1)])},
#'4': {'5':set([(0,1)])},
'5': {'6':set([(0,1)])},
'6': {'7':set([(0,1)]),'2':set([(0,1)])},
'7': {'1':set([(0,1)])},
}
#'8': {'9':set([(0,1)])},
#'9': {'01':set([(0,1)])},b
#'01': {'1':set([(0,1)])},
# '02': {'03':set([(0,1)]),'6':set([(0,1)])},
# '03': {'1':set([(0,1)])},
# }

L01 = {
'1': {'2':set([(0,1)])},
'2': {'3':set([(0,1)]),'7':set([(0,1)])},
'3': {'4':set([(0,1)]),'7':set([(0,1)])},
'4': {'5':set([(0,1)]),'6':set([(0,1)])},
'5': {'6':set([(0,1)])},
'6': {'3':set([(0,1)])},
'7': {'1':set([(0,1)])},
}
L02 = {
'1': {'2':set([(0,1)]),'1':set([(0,1)])},
'2': {'3':set([(0,1)])},
'3': {'4':set([(0,1)]),'6':set([(0,1)])},
'4': {'5':set([(0,1)])},
'5': {'6':set([(0,1)]),'4':set([(0,1)])},
'6': {'7':set([(0,1)])},
'7': {'1':set([(0,1)]),'7':set([(0,1)])},
}

L03 = {
'1': {'2':set([(0,1)]),'1':set([(0,1)])},
'2': {'3':set([(0,1)])},
'3': {'4':set([(0,1)])},
'4': {'5':set([(0,1)])},
'5': {'6':set([(0,1)])},
'6': {'7':set([(0,1)])},
'7': {'1':set([(0,1)])},
}


TT = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)]),'01':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'01':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T10 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'01':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T22 = {
'01': {'02':set([(0,1)]),'10':set([(0,1)])},
'02': {'01':set([(0,1)])},
'03': {'04':set([(0,1)]),'10':set([(0,1)])},
'04': {'03':set([(0,1)])},
'05': {'06':set([(0,1)]),'10':set([(0,1)])},
'06': {'05':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'07':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T63 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'01':set([(0,1)]),'10':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'07':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T631 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'01':set([(0,1)]),'10':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'07':set([(0,1)]),'10':set([(0,1)]),'09':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T33 = {
'01': {'02':set([(0,1)]),'10':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'01':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'04':set([(0,1)]),'10':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'07':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}


OO = {
'01': {'02':set([(0,1)]),'11':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'10':set([(0,1)])},
'10': {'01':set([(0,1)])},
'11': {'11':set([(0,1)])},
'12': {'05':set([(0,1)]),'06':set([(0,1)]),'12':set([(0,1)])},
}

OO5 = {
'01': {'03':set([(0,1)]),'11':set([(0,1)])},
'02': {'04':set([(0,1)])},
'03': {'05':set([(0,1)])},
'04': {'06':set([(0,1)])},
'05': {'07':set([(0,1)])},
'06': {'08':set([(0,1)])},
'07': {'09':set([(0,1)])},
'08': {'10':set([(0,1)])},
'09': {'01':set([(0,1)])},
'10': {'02':set([(0,1)]),'11':set([(0,1)])},
'11': {'11':set([(0,1)])},
'12': {'05':set([(0,1)]),'06':set([(0,1)]),'12':set([(0,1)])},
}

MV = {
'1': {'4':set([(0,1)])},
'2': {},
'3': {},
'4': {'6':set([(0,1)])},
'5': {'4':set([(0,1)]),'2':set([(0,1)])},
'6': {'3':set([(0,1)]), '5':set([(0,1)])},
}


MV_ = {
'1': {'2':set([(0,1)]),'3':set([(0,1)])},
'2': {'3':set([(2,0)])},
'3': {'2':set([(2,0)])},
}

G96 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'01':set([(0,1)]),'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)]),'08':set([(0,1)])},
'07': {'04':set([(0,1)])},
'08': {'09':set([(0,1)]),'24':set([(0,1)])},
'09': {'10':set([(0,1)])},
'10': {'11':set([(0,1)])},
'11': {'12':set([(0,1)])},
'12': {'13':set([(0,1)])},
'13': {'14':set([(0,1)])},
'14': {'06':set([(0,1)])},

# '15': {'16':set([(0,1)]),'15':set([(0,1)])},
# '16': {'17':set([(0,1)])},
# '17': {'18':set([(0,1)])},
# '18': {'19':set([(0,1)])},
# '19': {'20':set([(0,1)])},
# '20': {'21':set([(0,1)])},
# '21': {'21':set([(0,1)])},
'22': {'01':set([(0,1)]),'22':set([(0,1)])},



'24': {'24':set([(0,1)])},
#'25': {'26':set([(0,1)])},
#'26': {'26':set([(0,1)])},
}

G960 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)]),'01':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'01':set([(0,1)])},
}

DDD1 = {
    '01': {'02':set([(0,1)]),'01':set([(0,1)])},
    '02': {'03':set([(0,1)])},
    '03': {'04':set([(0,1)]), '06':set([(0,1)])},
    '04': {'05':set([(0,1)])},
    '05': {'06':set([(0,1)])},
    '06': {'01':set([(0,1)])}
}

DDD2 = {
    '01': {'02':set([(0,1)]),'01':set([(0,1)])},
    '02': {'03':set([(0,1)])},
    '03': {'06':set([(0,1)])},
    '04': {'03':set([(0,1)])},
    '05': {'04':set([(0,1)])},
    '06': {'01':set([(0,1)]), '05':set([(0,1)])}
}
"""Generic object pickler and compressor

This module saves and reloads compressed representations of generic Python
objects to and from the disk.
"""

__author__ = "Bill McNeill <billmcn@speakeasy.net>"
__version__ = "1.0"

import cPickle
import gzip

def save(object, filename, protocol = -1):
    """Save an object to a compressed disk file.
       Works well with huge objects.
    """
    with gzip.GzipFile(filename, 'wb') as file:
        cPickle.dump(object, file, protocol)

def load(filename):
    """Loads a compressed object from disk
    """
    with gzip.GzipFile(filename, 'rb') as file:
        object = cPickle.load(file)
    return object

if __name__ == "__main__":
	import sys
	import os.path
	
	class Object:
		x = 7
		y = "This is an object."
	
	filename = sys.argv[1]
	if os.path.isfile(filename):
		o = load(filename)
		print "Loaded %s" % o
	else:
		o = Object()
		save(o, filename)
		print "Saved %s" % o
import sys

sys.path.append('./tools/')
import graphkit as gk
import numpy as np
from ortools.constraint_solver import pywrapcp
import ipdb


def grelabel(g):
    """
    Relabel graph nodes to numerical contiguous labels 1 through len(g)
    :param g: input graph
    :return: relabeled graph and name mapping
    """
    gg = {}
    nodemap = {x[0]: x[1] + 1 for x in zip(g.keys(), range(len(g)))}
    for v in g:
        gg[nodemap[v]] = {}
        for w in g[v]:
            gg[nodemap[v]][nodemap[w]] = g[v][w]
    return gg, {v: k for k, v in nodemap.items()}


def crelabel(c, map):
    newcliq = set()
    for e in c:
        newcliq.add((map[e[0] + 1], map[e[1] + 1]))
    return newcliq


def printedges(g):
    l = gk.edgelist(g)
    for e in l:
        print e[0] - 1, '->', e[1] - 1


def edge_exists(i, j, g):
    return j + 1 in g[i + 1]


def numkids(i, edges, solver):
    N = len(edges)
    return solver.Sum([edges[i][k] for k in range(N)])


def numparents(i, edges, solver):
    N = len(edges)
    return solver.Sum([edges[k][i] for k in range(N)])


def clique_constrain(solver, parents, children, edges, g):
    N = len(parents)

    # declare constraints
    solver.Add(solver.Sum(parents) > 0)
    solver.Add(solver.Sum(children) > 0)

    for i in range(N):
        # makes sure that there exists at least one outgoing edge from node i if it is marked in parents for this b-clique
        solver.Add((parents[i] == 1) == (solver.Sum([edges[i][k] for k in range(N)]) >= 1))
        # this makes sure that there exists at least one incoming edge to node i if it is marked as a child
        solver.Add((children[i] == 1) == (solver.Sum([edges[k][i] for k in range(N)]) >= 1))

        solver.Add((children[i] == 1) == (solver.Sum(parents) <= numparents(i, edges, solver)))
        solver.Add((parents[i] == 1) == (solver.Sum(children) <= numkids(i, edges, solver)))

        for j in range(N):
            # edge existence constraints
            if not edge_exists(i, j, g):
                solver.Add(edges[i][j] == 0)
            else:
                solver.Add((edges[i][j] == 1) == (parents[i] * children[j] == 1))


def bcliques(g, verbose=False):
    solver = pywrapcp.Solver("b-clique")
    g, mp = grelabel(g)
    # declare variables
    edges = []
    N = len(g)
    for i in range(N):
        e = []
        for j in range(N):
            e.append(solver.IntVar(0, 1, "%i -> %i" % (i, j)))
        edges.append(e)

    parents = [solver.IntVar(0, 1, "%i" % (i)) for i in range(N)]
    children = [solver.IntVar(0, 1, "%i" % (i)) for i in range(N)]

    # declare constraints
    clique_constrain(solver, parents, children, edges, g)

    # run the solver
    solution = solver.Assignment()
    solution.Add([edges[i][j] for i in range(N) for j in range(N)])
    solution.Add(parents)
    solution.Add(children)

    collector = solver.AllSolutionCollector(solution)
    solver.Solve(solver.Phase([edges[i][j] for i in range(N) for j in range(N)] + children + parents,
                              solver.CHOOSE_FIRST_UNBOUND,
                              solver.ASSIGN_MAX_VALUE),
                 [collector])
    num_solutions = collector.SolutionCount()

    # output solutions
    if verbose:
        print "num_solutions:", num_solutions
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()

    cliques = []
    pts = set()
    for s in range(num_solutions):
        c = set()
        qval = [collector.Value(s, edges[i][j]) for i in range(N) for j in range(N)]
        pval = [collector.Value(s, parents[i]) for i in range(N)]
        cval = [collector.Value(s, children[i]) for i in range(N)]
        if verbose:
            print " ----------------- ", s
        if verbose:
            print np.asarray([pval, cval])
        for i in range(len(qval)):
            if qval[i]:
                e = np.unravel_index(i, [N, N])
                c.add((e[0], e[1]))
                if verbose:
                    print e[0], "->", e[1]
        if not (True in map(lambda x: c.issubset(x), cliques)):
            if not c.issubset(pts):
                pts = pts.union(c)
                cliques.append(c)
        if verbose:
            print "---"

    # check for b-cliques for which all edges are covered in other b-cliques
    cc = []
    for i in range(len(cliques)):
        bcl = cliques.pop()
        pts = set()
        for j in range(len(cliques)):
            pts = pts.union(cliques[j])
        for j in range(len(cc)):
            pts = pts.union(cc[j])

        if not bcl.issubset(pts):
            cc.append(bcl)
    return map(lambda x: crelabel(x, mp), cc)
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import ecj, bfutils
import bfutils as bfu
import graphkit as gk
import warnings
from statsmodels.tsa.api import VAR
from sympy.matrices import SparseMatrix
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed
from scipy import linalg, optimize
import numpy as np
import scipy
import ipdb

def symchol(M): # symbolic Cholesky
    B = SparseMatrix(M)
    t = B.row_structure_symbolic_cholesky()
    B = np.asarray(B)*0
    for i in range(B.shape[0]): B[i,t[i]] = 1
    return B

def G2SVAR(G):
    n = len(G)
    A,B = npG2SVAR(G)
    P,L,U = linalg.lu(B)
    A = linalg.inv(L).tolist()
    B = B.tolist()
    A = listplace(A, 0.0, 0.0)
    for i in range(0,n): A[i][i] = 1
    B = listplace(B, 0.0, 'e')
    for i in range(0,n): B[i][i] = 'e'
    return A,B,P

def G2AH(G):
    n = len(G)
    A,B = npG2SVAR(G)
    P,L,U = linalg.lu(B)
    A = linalg.inv(L).tolist()
    B = B.tolist()
    A = listplace(A, 0.0, 0.0)
    for i in range(0,n): A[i][i] = 1
    B = listplace(B, 0.0, 'e')
    for i in range(0,n): B[i][i] = 'e'
    return A,B,P

def bnf2CG(fname):
    d = eval(open(fname).read())
    G = {}
    for v in d:
        G[v] = {u: set((0,1)) for u in d[v]['pars']}
    G = ecj.tr(G)
    for v in G:
        ld = {u: set([(0,1)]) for u in G[v]}
        G[v] = ld
    return G

def npG2SVAR(G):
    n = len(G)
    A = bfu.graph2adj(G)
    B = np.tril(bfu.graph2badj(G))
    np.fill_diagonal(B,1)
    B = symchol(B)
    return A,B

def x2M(x, A, B, aidx, bidx):
    A[aidx] = x[:len(aidx[0])]
    B[bidx] = x[len(aidx[0]):]
    #B[(bidx[1],bidx[0])] = x[len(aidx[0]):]
    return A, B

def nllf(x, A, B, Y, aidx, bidx): # negative log likelihood
    A,B = x2M(x, A, B, aidx, bidx)
    T = Y.shape[1]
    X = Y[:,1:] - np.dot(A, Y[:,:-1])
    ldB = T*np.log(abs(1./linalg.det(B)))
    return ldB + 0.5*np.trace( np.dot(np.dot(B.T, B), np.dot(X,X.T)))

def nllf2(x, A, B, YY, XX, YX, T, aidx, bidx): # negative log likelihood
    A,B = x2M(x, A, B, aidx, bidx)
    AYX = np.dot(A, YX.T)
    S = YY - AYX - AYX.T + np.dot(np.dot(A,XX), A.T)
    ldB = T*np.log(abs(1./linalg.det(B)))
    return 0.5*np.dot(np.dot(B.T, B).T.flat,S.flat) + ldB
    #return ldB + 0.5*np.trace( np.dot(np.dot(B.T, B), S))

def VARbic(nllf, K, T):
    return 2*nllf + K*np.log(T)

def listplace(l, a, b):
    return [listplace(x,a,b) if not np.isscalar(x) else b if x != a  else x for x in l]

# -------------------------------------------------------------------
# data generation
# -------------------------------------------------------------------

def randweights(n, c=0.1, factor=9):
    rw = scipy.random.randn(n)
    idx = scipy.where(abs(rw) < factor*c)
    if idx:
        rw[idx] = rw[idx]+scipy.sign(rw[idx])*c*factor
    return rw

def transitionMatrix(cg, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = scipy.where(A==1)
    A[edges] = randweights(edges[0].shape[0], c=minstrength)
    l = linalg.eig(A)[0]
    c = 0
    pbar = ProgressBar(widgets=['Searching for weights: ', Percentage(), ' '], maxval=10000).start()
    while max(l*scipy.conj(l)) > 1:
        A[edges] = randweights(edges[0].shape[0], c=c)
        c += 1
        l = linalg.eig(A)[0]
        pbar.update(c)
    pbar.finish()
    return A

def sampleWeights(n, minstrength=0.1):
    r = scipy.randn(n)
    s = minstrength/np.min(np.abs(r))
    r = s*r
    return r

def transitionMatrix2(cg, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = scipy.where(A==1)
    A[edges] = sampleWeights(edges[0].shape[0], minstrength=minstrength)
    l = linalg.eig(A)[0]
    c = 0
    pbar = ProgressBar(widgets=['Searching for weights: ', Percentage(), ' '], maxval=10000).start()
    while max(l*scipy.conj(l)) > 1:
        A[edges] = sampleWeights(edges[0].shape[0], minstrength=minstrength)
        c += 1
        l = linalg.eig(A)[0]
        if c>pbar.maxval:
            raise ValueError
        pbar.update(c)
    pbar.finish()
    return A

def transitionMatrix3(cg, x0=None, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = scipy.where(A==1)

    try:
        s = x0.shape
        x = x0
    except AttributeError:
        A = initRandomMatrix(A, edges)
        x = A[edges]

    def objective(x):
        A[edges] = np.real(x)
        l = linalg.eig(A)[0]
        m = np.max(np.real(l*scipy.conj(l)))-0.99
        n = np.min(np.min(np.abs(x)),minstrength)-minstrength
        return m*m + 0.1*n*n

    o = np.zeros(len(edges))
    while np.min(np.abs(o[0])) < 0.8*minstrength:
        rpt = True
        while rpt:
            try:
                try:
                    o = optimize.fmin_bfgs(objective, x,
                                           gtol=1e-10, maxiter=100,
                                           disp=False, full_output=True)
                    A[edges]=np.real(o[0])
                    l = linalg.eig(A)[0]
                    if np.max(np.real(l*scipy.conj(l))) < 1:
                        rpt = False

                except:
                    rpt = True
            except Warning:
                x = scipy.randn(len(edges[0]))
                rpt = True
    A[edges]=np.real(o[0])
    return A

def initRandomMatrix(A, edges, maxtries=100, distribution='beta', stable=True):
    '''
    possible distributions:
    flat
    flatsigned
    beta
    normal
    uniform
    '''
    s = 2.0

    def init():
        if distribution=='flat':
            x = np.ones(len(edges[0]))
        elif distribution=='flatsigned':
            x = np.sign(scipy.randn(len(edges[0])))*scipy.ones(len(edges[0]))
        elif distribution=='beta':
            x = np.random.beta(0.5,0.5,len(edges[0]))*3-1.5
        elif distribution=='normal':
            x = scipy.randn(len(edges[0]))
        elif distribution=='uniform':
            x = np.sign(scipy.randn(len(edges[0])))*scipy.rand(len(edges[0]))
        else:
             raise ValueError('Wrong option!')
        return x

    def eigenvalue(A):
        l = linalg.eig(A)[0]
        s = np.max(np.real(l*scipy.conj(l)))
        return s

    x = init()
    A[edges] = x
    s = eigenvalue(A)
    alpha = np.random.rand()*(0.99-0.8)+0.8
    A = A/(alpha*s)
    s = eigenvalue(A)

    return A

def transitionMatrix4(g, minstrength=0.1, distribution='normal', maxtries=1000):
    A = gk.CG2adj(g)
    edges = np.where(A==1)
    s = 2.0
    c = 0
    pbar = ProgressBar(widgets=['Searching for weights: ',
                                Percentage(), ' '],
                       maxval=maxtries).start()
    while s > 1.0:
        minstrength -= 0.001
        A = initRandomMatrix(A, edges, distribution=distribution)
        x = A[edges]
        delta = minstrength/np.min(np.abs(x))
        A[edges] = delta*x
        l = linalg.eig(A)[0]
        s = np.max(np.real(l*scipy.conj(l)))
        c += 1
        if c > maxtries:
            return None
        pbar.update(c)
    pbar.finish()

    return A

def drawsamplesLG(A, nstd=0.1, samples=100):
    n = A.shape[0]
    data = scipy.zeros([n, samples])
    data[:,0] = nstd*scipy.random.randn(A.shape[0])
    for i in range(1,samples):
        data[:,i] = scipy.dot(A,data[:,i-1]) \
                    + nstd*scipy.random.randn(A.shape[0])
    return data


def getAgraph(n, mp=2, st=0.5, verbose=True):
    keeptrying = True
    while keeptrying:
        G = gk.rnd_cg(n, maxindegree=mp, force_connected=True)
        try:
            A = transitionMarix2(G, minstrength=st)
            keeptrying = False
        except ValueError as e:
            if verbose:
                print "!!! Unable to find strong links for a stable matrix !!!"
                print "*** trying a different graph"
    return {'graph':      G,
            'transition': A,
            'converges':  len(bfutils.call_undersamples(G))}

def getAring(n, density=0.1, st=0.5, verbose=True, dist='flatsigned', permute=False):
    keeptrying = True
    plusedges = bfutils.dens2edgenum(density,n)
    while keeptrying:
        G = gk.ringmore(n, plusedges, permute=permute)
        try:
            A = transitionMatrix4(G, minstrength=st, distribution=dist)
            try:
                s = A.shape
                keeptrying = False
            except AttributeError:
                keeptrying = True
        except ValueError:
            if verbose:
                print "!!! Unable to find strong links for a stable matrix !!!"
                print "*** trying a different graph"
    return {'graph':      G,
            'transition': A,
            'converges':  len(bfutils.call_undersamples(G))}



# -------------------------------------------------------------------
# estimation
# -------------------------------------------------------------------

def scoreAGraph(G, data, x0 = None):
    A,B = npG2SVAR(G)
    n = len(G)
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    K = scipy.sum(len(a_idx[0])+len(b_idx[0])/2)

    if x0:
        o = optimize.fmin_bfgs(nllf, x0, args=(A, B, data, a_idx, b_idx),
                               disp=False, full_output=True)
    else:
        o = optimize.fmin_bfgs(nllf, scipy.randn(len(a_idx[0])+len(b_idx[0])),
                               args=(np.double(A), np.double(B),
                                     data, a_idx, b_idx),
                               disp=False, full_output=True)
    ipdb.set_trace()
    return 2*o[1] + K*np.log(data.shape[1]) #VARbic(o[1],K,data.shape[1])

def scoreAGraph2(G, data, x0 = None):
    A,B = npG2SVAR(G)
    n = len(G)
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)

    T = data.shape[1]
    YY = np.dot(data[:,1:],data[:,1:].T)
    XX = np.dot(data[:,:-1],data[:,:-1].T)
    YX = np.dot(data[:,1:],data[:,:-1].T)

    K = scipy.sum(len(a_idx[0])+len(b_idx[0])/2)

    if x0:
        o = optimize.fmin_bfgs(nllf2, x0,
                               args=(np.double(A), np.double(B),
                                     YY,XX,YX,T,a_idx, b_idx),
                                     gtol=1e-12, maxiter=500,
                                     disp=False, full_output=True)
    else:
        o = optimize.fmin_bfgs(nllf2, scipy.randn(len(a_idx[0])+len(b_idx[0])),
                               args=(np.double(A), np.double(B),
                                     YY,XX,YX,T,a_idx, b_idx),
                                gtol=1e-12, maxiter=500,
                                disp=False, full_output=True)

    return 2*o[1] + K*np.log(data.shape[1]) #VARbic(o[1],K,data.shape[1])

def estimateG(G,YY,XX,YX,T,x0=None):
    A,B = npG2SVAR(G)
    K = scipy.sum(abs(A)+abs(B))
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    try:
        s = x0.shape
        x = x0
    except AttributeError:
        x = scipy.randn(K)
    o = optimize.fmin_bfgs(nllf2, x,
                           args=(np.double(A), np.double(B),
                                 YY,XX,YX,T,a_idx, b_idx),
                           disp=False, full_output=True)
    A,B = x2M(o[0], np.double(A), np.double(B), a_idx, b_idx)
    return  A,B


def data2AB(data,x0=None):
    n = data.shape[0]
    T = data.shape[1]
    YY = np.dot(data[:,1:],data[:,1:].T)
    XX = np.dot(data[:,:-1],data[:,:-1].T)
    YX = np.dot(data[:,1:],data[:,:-1].T)

    model = VAR(data.T)
    r = model.fit(1)
    A = r.coefs[0,:,:]

    #A = np.ones((n,n))
    B = np.ones((n,n))
    np.fill_diagonal(B,0)
    B[np.triu_indices(n)] = 0
    K = np.int(scipy.sum(abs(B)))#abs(A)+abs(B)))

    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    np.fill_diagonal(B,1)

    try:
        s = x0.shape
        x = x0
    except AttributeError:
        x = np.r_[A.flatten(),0.1*scipy.randn(K)]
    o = optimize.fmin_bfgs(nllf2, x,
                           args=(np.double(A), np.double(B),
                                 YY,XX,YX,T,a_idx, b_idx),
                           gtol=1e-12, maxiter=500,
                           disp=False, full_output=True)
    ipdb.set_trace()
    A,B = x2M(o[0], np.double(A), np.double(B), a_idx, b_idx)
    B = B+B.T
    return  A,B

def amap(f, a):
     v = np.vectorize(f)
     return v(a)

def AB2intAB(A,B, th=0.09):
    A[amap(lambda x: abs(x) > th, A)] = 1
    A[amap(lambda x: abs(x) < 1, A)] = 0
    B[amap(lambda x: abs(x) > th, B)] = 1
    B[amap(lambda x: np.abs(x) < 1, B)] = 0
    return A,B

def intAB2graph(A,B):
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if A[j,i]: g[str(i+1)][str(j+1)] = set([(0,1)])

    for i in range(n):
        for j in range(n):
            if B[j,i] and j!=i:
                if str(j+1) in g[str(i+1)]:
                    g[str(i+1)][str(j+1)].add((2,0))
                else:
                    g[str(i+1)][str(j+1)] = set([(2,0)])
    return g

def data2graph(data,x0=None, th=0.09):
    A,B = data2AB(data,x0=x0)
    Ab,Bb = AB2intAB(A,B,th=th)
    return intAB2graph(Ab,Bb)

def data2VARgraph(data, pval=0.05):
    model = VAR(data.T)
    r = model.fit(1)
    A = r.coefs[0,:,:]
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if np.abs(A[j,i]) > pval: g[str(i+1)][str(j+1)] = set([(0,1)])
    return g

def data2VARgraph_model(data, pval=0.05):
    model = VAR(data.T)
    r = model.fit(1)
    A = r.coefs[0,:,:]
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if np.abs(A[j,i]) > pval: g[str(i+1)][str(j+1)] = set([(0,1)])
    return g, r
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))
import traversal as trv
import bfutils as bfu
import graphkit as gk
import zickle as zkl
import numpy as np
import itertools as iter
import statsmodels.api as sm
import linear_model as lm
import ipdb
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore, norm

def independent(y,X,pval=0.05):
    X  = sm.add_constant(X)
    est  = sm.OLS(y,X).fit()
    return est.pvalues[1] > pval

def kernel(z):
    if np.abs(z) > 1.: return 0.
    return .5

def residuals_(x,y,z):
    N = len(x)
    pd = map(kernel,pdist(z.T))
    PD = squareform(pd)
    sumsx = np.dot(PD,x) + 0.5*x
    sumsy = np.dot(PD,y) + 0.5*y
    weights = np.sum(PD, axis=1) + 0.5
    residualsx = x - sumsx/weights
    residualsy = y - sumsy/weights

    return residualsx, residualsy

def moment22(x,y): return np.dot(x*x,y*y)/len(x)

def fdr(alpha, pvalues):
    m = len(pvalues)
    c = np.cumsum(1./(np.arange(m)+1))
    pcompare = (pvalues <= alpha*(np.arange(1,m+1)/(c*(m+1))))
    idx = np.where(pcompare == True)[0]
    if len(idx) == 0: return -1
    return idx[-1]

def fdrQ(alpha,pvalues):
    pvalues = np.sort(pvalues)
    pvalues = pvalues[~np.isnan(pvalues)]
    min  = np.nan if len(pvalues) == 0 else pvalues[0]
    high = 1.0
    low  = 0.
    q    = alpha
    while (high - low) > 1e-5:
        midpoint = (high + low)/2.0
        q = midpoint
        cutoff = pvalues[fdr(q, pvalues)]

        if cutoff < min:
            low = midpoint
        elif cutoff > min:
            high = midpoint
        else:
            low  = midpoint
            high = midpoint
    return q

def fdrCutoff(alpha, pvalues):
    pvalues = np.sort(pvalues)
    k = fdr(alpha,pvalues)
    if k < 0: return 0
    return pvalues[k]

def np_fisherZ(x,y,r):
    z = 0.5 * (np.log(1.0 + r) - np.log(1.0 - r))
    w = np.sqrt(len(x)) * z
    x_ = zscore(x)
    y_ = zscore(y)
    t2 = moment22(x_,y_)
    t = np.sqrt(t2)
    p = 2. * (1. - norm.cdf(np.abs(w), 0.0, t))
    return p

def independent_(x,y, alpha = 0.05):
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]

    # For PC, should not remove the edge for this reason.
    if len(x) < 10: return False

    ps = []
    for i in range(15):
        for j in range(15):
            x_ = x**i
            y_ = y**j
            r = np.corrcoef(x_,y_)[0,1]
            #r = max(min(r,1),-1) # Tetrad had this
            p = np_fisherZ(x_,y_,r)
            #if not np.isnan(p):
            ps.append(p)

    if not ps: return True
    return fdrCutoff(alpha, ps) > alpha


def addallb(g):
    n = len(g)
    for i in range(n):
        for j in range(n):
            if str(j+1) in g[str(i+1)]:
                g[str(i+1)][str(j+1)].add((2,0))
            else:
                g[str(i+1)][str(j+1)] = set([(2,0)])
    return g

def dpc(data, pval=0.05):
    n = data.shape[0]
    # stack the data: first n rows is t-1 slice, the next n are slice t
    data = np.asarray(np.r_[data[:,:-1],data[:,1:]])

    def tetrad_cind_(y,x,condset=[], alpha=0.01, shift=0):
        y = data[n+int(y)-1,:]
        x = data[shift+int(x)-1,:]
        if condset:
            X  = data[condset,:]
            ry, rx = residuals_(y,x,X)
        else:
            ry, rx = [y,x]
        return independent_(ry, rx, alpha = alpha)

    def cind_(y,x, condset=[], pval=pval, shift=0):
        yd = data[n+int(y)-1,:].T
        X  = data[[shift+int(x)-1]+condset,:].T
        return independent(yd, X, pval=pval)

    def cindependent(y, x, counter, parents=[], pval=pval):
        for S in [j for j in iter.combinations(parents,counter)]:
            if cind_(y, x, condset=list(S), pval=pval): return True
            #if tetrad_cind_(x, y, condset=list(S), alpha=pval): return True                
        return False

    def bindependent(y, x, parents=[], pval=pval):
        return cind_(y, x, condset=parents, pval=pval, shift=n)
        #return tetrad_cind_(y, x, condset=parents, alpha=pval, shift=n)

    def prune(elist, mask, g):
        for e in mask:
            g[e[0]][e[1]].remove((0,1))
            elist.remove(e)
        gk.clean_leaf_nodes(g)

    g  = gk.superclique(n)
    gtr= bfu.gtranspose(g)

    el = gk.edgelist(g)
    for counter in range(n):
        to_remove = []
        for e in el:
            ppp = [int(k)-1 for k in gtr[e[1]] if k != e[0]]
            if counter <= len(ppp):
                if cindependent(e[1], e[0], counter, parents=ppp, pval=pval):
                    to_remove.append(e)
                    gtr[e[1]].pop(e[0],None)
        prune(el, to_remove, g)

    bel = [map(lambda k: str(k+1), x) for x in iter.combinations(range(n),2)]
    for e in bel:
        ppp = list(set(gtr[e[0]].keys()) | set(gtr[e[1]].keys()))
        ppp = map(lambda x: int(x)-1, ppp)
        if bindependent(e[0], e[1], parents=ppp, pval=pval):
            g[e[0]][e[1]].remove((2,0))
            g[e[1]][e[0]].remove((2,0))
    gk.clean_leaf_nodes(g)

    return g


# Local Variables:
# mode: python
# python-indent-offset: 4
# End:
import networkx
from bfutils import g2num, num2CG
from comparison import graph2nx
from ecj import undersample
from numpy import argsort

def simple_loops(g, u):
    """
    iterator over the list of simple loops of graph g at the undersample rate u
    """
    gx = graph2nx(num2CG(g2num(undersample(g,u)), len(g)))
    for l in networkx.simple_cycles(gx):
        yield l

def print_loops(g, u):
    l = [x for x in simple_loops(g,u)]
    lens = map(len, l)    
    for i in argsort(lens):
        print l[i]

def ul(l):
    """
    returns list elements that are present only once
    """
    u, r = set(), set()
    for e in l:
        if e not in u:
            u.add(e)
        else:
            r.add(e)
    return u.difference(r)
import sys

sys.path.append('./tools/')

import bclique as bq
import latents as lt
import pathtreetools as ptt
from pathtree import PathTree


def can_add_loop(pt, num, elements):
    r = False
    if ptt.isptelement(PathTree(num - pt.preset, pre=pt.preset), num):
        r = True
    for e in pt.loopset:
        if type(e) is int:
            can_add_loop(PathTree(pt.preset))

def learn_path_tree(pt):
    elements = ptt.pt2seq(pt, 1)
    newpt = PathTree(set(), pre={elements[0]})

    def rpath(elements, npt):
        if not ptt.isptelement(npt, element[0]):


    return newpt

def bpts(bc):
    """
    Given a b-clique returns a set of PathTrees (one for each edge in bclique) such that loops are identified across edges in bclique
    :param bc: b-clique
    :return:
    """

def getagenerator(g):
    bclqs = bq.bcliques(g)
    for clq in bclqs:
        bpts(clq)
    return bclqs
import ecj
import scipy
import numpy
import operator
import networkx as nx
#from progressbar import ProgressBar, Percentage
numpy.random.RandomState()
import bfutils as bfu
import numpy as np
import gmpy as gmp

def num2CG(num,n):
    """num2CG - converts a number  whose binary representaion encodes edge
    presence/absence into a compressed graph representaion

    """
    n2 = n*n
    G = {'%i'%(i+1):{} for i in xrange(n)}
    if num == 0: return G
    bl = gmp.bit_length(num)
    idx = [n2-i-1 for i in xrange(bl) if num & (1<<i)]
    idx = np.unravel_index(idx,(n,n))
    x = idx[0]+1
    y = idx[1]+1
    for i in xrange(len(x)):
        G['%i' % x[i]]['%i' % y[i]] = set([(0,1)])
    return G

def hasSelfLoops(G):
    for u in G:
        if G[u].has_key(u):
            return True
    return False

def randSCC(n):
    G = num2CG(scipy.random.randint(2**(n**2)),n)
    while (len(ecj.scc(G)) > 1) or gcd4scc(G)>1:
        G = num2CG(scipy.random.randint(2**(n**2)),n)
    return G

def SM_fixed(Gstar,G, iter=5):
    compat = []
    for j in range(0,iter):
        if Gstar == ecj.undersample(G,j):
            compat.append(j)
    return compat
def SM_converging(Gstar,G):
    """Gstar is the undersampled reference graph, while G is the starting
    graph. The  code searches  over all undersampled  version of  G to
    find all matches with Gstar
    """
    compat = []
    GG = G
    Gprev = G
    if G == Gstar: return [0]
    j = 1
    G = ecj.undersample(GG,j)
    while not (G == Gprev):
        if Gstar == G: compat.append(j)
        j += 1
        Gprev = G
        G = ecj.undersample(GG,j)
    return compat
def searchMatch(Gstar,G, iter=5):
    if gcd4scc(G) >1: return SM_fixed(Gstar, G, iter=iter)
    return SM_converging(Gstar, G)
def hasSink(G):
    return not reduce(operator.and_, [bool(G[n]) for n in G], True)
def hasRoot(G): return hasSink(ecj.tr(G))
def isSclique(G):
    n = len(G)
    for v in G:
        if sum([(0,1) in G[v][w] for w in G[v]]) < n: return False
        if sum([(2,0) in G[v][w] for w in G[v]]) < n-1: return False
    return True

def graph2nx(G):
    g = nx.DiGraph()
    for v in G:
        g.add_edges_from([(v,x) for x in G[v] if (0,1) in G[v][x]])
    return g
def nx2graph(G):
    g = {str(n+1):{} for n in G}
    for n in G:
        g['%i' % (n+1)] = {'%i' % (x+1):set([(0,1)]) for x in G[n]}
    return g

def gcd4scc(SCC):
    g = graph2nx(SCC)
    return ecj.listgcd(map(lambda x: len(x)-1, nx.simple_cycles(g)))

def compatibleAtU(uGstar):
    compat = []
    n = len(uGstar)
    numG = 2**(n**2)
    #pbar = Percentage()
    for i in range(1,numG):
        G = num2CG(i,n)
        #pbar.update(i+1)
        if len(ecj.scc(G)) > 1: continue
        l = searchMatch(uGstar,G, iter = 5)
        if l: compat.append((l,G))
    #pbar.finish()
    return compat
from copy import deepcopy
all_acceptable = []

def find_next_graph(graph_g):
    next_graph = {}
    for i in graph_g:next_graph[i] = {}
    for i in graph_g:
        for j in graph_g[i]:
            for m in graph_g[j]:
                if not m in next_graph[i]:next_graph[i][m] = set([(0,1)])
                elif next_graph[i][m] == set([(0,2)]) or next_graph[i][m] == set([(0,3)]):next_graph[i][m] = set([(0,3)])
            for n in graph_g[i]:
                if j != n:
                    if not n in next_graph[j]:next_graph[j][n] = set([(0,2)])
                    elif next_graph[j][n] == set([(0,1)]) or next_graph[j][n] == set([(0,3)]):next_graph[j][n] = set([(0,3)])
    return next_graph

def compare(g_two_star,graph_g):
    for i in g_two_star:
        for j in g_two_star[i]:
            if j not in graph_g[i]:
                return False
            elif g_two_star[i][j] != graph_g[i][j] and  graph_g[i][j] != set([(0,3)]):
                return False
    return True

def try_next_double_edge(g_de, g_one_star_de, stack_de):
    global all_acceptable
    if stack_de != []:
        j = stack_de.pop()
        i = stack_de.pop()
        for k in g_one_star_de:
            boolean_first_edge_has_value = False
            boolean_second_edge_has_value = False
            if i in g_one_star_de[k]:
                boolean_first_edge_has_value = True
            if j in g_one_star_de[k]:
                boolean_second_edge_has_value = True
            g_one_star_de[k][i] = 1
            g_one_star_de[k][j] = 1
            g_two_star_de = find_next_graph(g_one_star_de)
            if compare(g_two_star_de, g_de):
                try_next_double_edge(g_de, g_one_star_de, stack_de)
            if not boolean_first_edge_has_value:
                del g_one_star_de[k][i]
            if not boolean_second_edge_has_value:
                del g_one_star_de[k][j]
        stack_de.append(i)
        stack_de.append(j)
    else:
        in_all_accept = False
        for xxxx in all_acceptable:
            if(xxxx == g_one_star_de):
                in_all_accept = True
        if (not in_all_accept):
            all_acceptable.append(deepcopy(g_one_star_de))
            print "Added 1."
#print g_one_star_de
#if g_one_star_de not in all_acceptable:
#all_acceptable.append(g_one_star_de)
        #print "G2 --------------------------"
        #print find_next_graph(g_one_star_de) == g_de
        #print "End-------------------------"


def sample_g_double_edge(g_de, g_one_star_de):
    stack_de = []
    for i in g_de:
        for j in g_de[i]:
            if g_de[i][j] != set([(0, 1)]):
                if int(i) < int(j):
                    stack_de.append(i)
                    stack_de.append(j)
    try_next_double_edge(g_de, g_one_star_de, stack_de)

def try_next_single_edge(g_se, g_one_star_se, stack_se):
    if stack_se != []:
        j = stack_se.pop()
        i = stack_se.pop()
        for k in g_one_star_se:
            boolean_first_edge_has_value = False
            boolean_second_edge_has_value = False
            if k in g_one_star_se[i]:
                boolean_first_edge_has_value = True
            if j in g_one_star_se[k]:
                boolean_second_edge_has_value = True
            g_one_star_se[i][k] = 1
            g_one_star_se[k][j] = 1
            g_two_star_se = find_next_graph(g_one_star_se)
            if compare(g_two_star_se, g_se):
                try_next_single_edge(g_se, g_one_star_se, stack_se)
            if not boolean_first_edge_has_value:
                del g_one_star_se[i][k]
            if not boolean_second_edge_has_value and j in g_one_star_se[k]:
                del g_one_star_se[k][j]
        stack_se.append(i)
        stack_se.append(j)
    else:
        sample_g_double_edge(g_se,g_one_star_se)

def sample_g_single_edge(g_se, g_one_star_se):
    stack_se = []
    for i in g_se:
        for j in g_se[i]:
            if g_se[i][j] != set([(0, 2)]):
                stack_se.append(i)
                stack_se.append(j)
    try_next_single_edge(g_se, g_one_star_se, stack_se)



def main(gStart):
    global all_acceptable
    all_acceptable = []
    g_one_star = {}
    for i in gStart:
        g_one_star[i] = {};
    sample_g_single_edge(gStart, g_one_star)
    return all_acceptableimport scipy
import itertools
# from progressbar import ProgressBar, Percentage
from multiprocessing import Pool, Array, Process, Manager
from numpy.random import randint
import numpy as np
# import ipdb
import networkx as nx
# local
import ecj
import zickle as zkl
import graphkit as gk
from comparison import num2CG, nx2graph, isSclique
import itertools
#import rpy2


def pure_directed_inc(G, D):
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in D[v]:
            if G[w]:
                for e in G[w]: G_un[v][e] = set([(0, 1)])
    return G_un


def directed_inc(G, D):
    G_un = {}
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in D[v]:
            if G[w] and (0, 1) in D[v][w]:
                for e in G[w]: G_un[v][e] = set([(0, 1)])
    return G_un


def bidirected_inc(G, D):
    # bidirected edges
    for w in G:
        # transfer old bidirected edges
        for l in D[w]:
            if (2, 0) in D[w][l]:
                G[w].setdefault(l, set()).add((2, 0))
        # new bidirected edges
        l = [e for e in D[w] if (0, 1) in D[w][e]]
        for pair in itertools.permutations(l, 2):
            G[pair[0]].setdefault(pair[1], set()).add((2, 0))
    return G


def increment(g):
    """
    Undersample the input graph g by 2
    Only works for g1 to g2 directed

    Args:
        g: an input graph in the dictionary format

    Returns:
        r: a graph in the dictionary format
    """
    r = {n: {} for n in g}

    for n in g:
        for h in g[n]:
            for e in g[h]:
                if not e in r[n]:
                    r[n][e] = set([(0, 1)])

    for n in g:
        for pair in itertools.combinations(g[n], 2):

            if pair[1] in r[pair[0]]:
                r[pair[0]][pair[1]].add((2, 0))
            else:
                r[pair[0]][pair[1]] = set([(2, 0)])

            if pair[0] in r[pair[1]]:
                r[pair[1]][pair[0]].add((2, 0))
            else:
                r[pair[1]][pair[0]] = set([(2, 0)])

    return r


def dincrement_u(G_star, G_u):
    # directed edges
    G_un = pure_directed_inc(G_star, G_u)
    return G_un


def increment_u(G_star, G_u):
    # directed edges
    G_un = directed_inc(G_star, G_u)
    # bidirected edges
    G_un = bidirected_inc(G_un, G_u)
    return G_un


def undersample(G, u):
    Gu = G
    for i in range(u):
        Gu = increment_u(G, Gu)
    return Gu


def all_undersamples(G_star):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if ecj.isSclique(g): return glist  # superclique convergence
        # this will (may be) capture DAGs and oscillations
        if g in glist: return glist
        glist.append(g)
    return glist


def graph2adj(G):
    n = len(G)
    A = scipy.zeros((n, n), dtype=np.int8)
    for v in G:
        A[int(v) - 1, [int(w) - 1 for w in G[v] if (0, 1) in G[v][w]]] = 1
    return A


def graph2badj(G):
    n = len(G)
    A = scipy.zeros((n, n), dtype=np.int8)
    for v in G:
        A[int(v) - 1, [int(w) - 1 for w in G[v] if (2, 0) in G[v][w]]] = 1
    return A


def adj2graph(A):
    G = {str(i): {} for i in range(1, A.shape[0] + 1)}
    idx = np.where(A == 1)
    for i in range(len(idx[0])):
        G['%i' % (idx[0][i] + 1)]['%i' % (idx[1][i] + 1)] = set([(0, 1)])
    return G


def adjs2graph(A, B):
    names = [str(i) for i in range(1, A.shape[0] + 1)]
    G = {}
    for name in names:
        G[name] = {}
    for i in range(A.shape[0]):
        for name in map(str, np.where(A[i, :] == 1)[0] + 1):
            G[str(i + 1)][name] = set([(0, 1)])

    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            if B[i, j]:
                if str(j + 1) in G[str(i + 1)]:
                    G[str(i + 1)][str(j + 1)].add((2, 0))
                else:
                    G[str(i + 1)][str(j + 1)] = set([(2, 0)])
    return G


def g2vec(g):
    A = graph2adj(g)
    B = graph2badj(g)
    return np.r_[A.flatten(), B[np.triu_indices(B.shape[0])]]


def vec2adj(v, n):
    A = np.zeros((n, n))
    B = np.zeros((n, n))
    A[:] = v[:n ** 2].reshape(n, n)
    B[np.triu_indices(n)] = v[n ** 2:]
    B = B + B.T
    return A, B


def vec2g(v, n):
    A, B = vec2adj(v, n)
    return adjs2graph(A, B)


# tried mutable ctypes buffer - not faster :(
def graph2str(G):
    n = len(G)
    d = {((0, 1),): '1', ((2, 0),): '0', ((2, 0), (0, 1),): '1', ((0, 1), (2, 0),): '1'}
    A = ['0'] * (n * n)
    for v in G:
        for w in G[v]:
            A[n * (int(v, 10) - 1) + int(w, 10) - 1] = d[tuple(G[v][w])]
    return ''.join(A)


def graph2bstr(G):
    n = len(G)
    d = {((0, 1),): '0', ((2, 0),): '1', ((2, 0), (0, 1),): '1', ((0, 1), (2, 0),): '1'}
    A = ['0'] * (n * n)
    for v in G:
        for w in G[v]:
            A[n * (int(v, 10) - 1) + int(w, 10) - 1] = d[tuple(G[v][w])]
    return ''.join(A)


def adj2num(A):
    s = reduce(lambda y, x: y + str(x),
               A.flatten().tolist(), '')
    return int(s, 2)


# def g2num(G): return int(graph2str(G),2) #adj2num(graph2adj(G))
def g2num(g):
    n = len(g)
    n2 = n ** 2 + n
    num = 0
    for v in range(1, n + 1):
        for w in g[str(v)]:
            num |= (1 << (n2 - v * n - int(w, 10)))
    return num


def ug2num(g):
    n = len(g)
    n2 = n ** 2 + n
    num = 0
    mask = 0
    num2 = 0
    for v in g:
        for w in g[v]:
            if (0, 1) in g[v][w]:
                mask = (1 << (n2 - int(v, 10) * n - int(w, 10)))
                num |= mask
            if (2, 0) in g[v][w]: num2 |= mask
    return num, num2


def bg2num(g):
    n = len(g)
    n2 = n ** 2 + n
    num = 0
    for v in g:
        for w in g[v]:
            if (2, 0) in g[v][w]:
                num = num | (1 << (n2 - int(v) * n - int(w)))
    return num


# def bg2num(G): return int(graph2bstr(G),2)#adj2num(graph2badj(G))
# def ug2num(G): return (g2num(G),bg2num(G))#(adj2num(graph2adj(G)),adj2num(graph2badj(G)))

def num2adj(num, n):
    l = list(bin(num)[2:])
    l = ['0' for i in range(0, n ** 2 - len(l))] + l
    return scipy.reshape(map(int, l), [n, n])


def add_bd_by_adj(G, adj):
    c = 0
    for e in adj:
        for v in range(len(e)):
            if e[v] == 1:
                try:
                    G[str(c + 1)][str(v + 1)].add((2, 0))
                except KeyError:
                    G[str(c + 1)][str(v + 1)] = set([(2, 0)])
        c += 1
    return G


def tuple2graph(t, n):
    g = num2CG(t[0], n)
    return add_bd_by_adj(g, num2adj(t[1], n))


def call_undersamples(G_star):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if g in glist: return glist
        glist.append(g)
    return glist


def overshoot(G_star, H):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if isSclique(g): return False
        if gk.isedgesubset(H, g): return True
        if g in glist: return False
        glist.append(g)
    return False


def forms_loop(G_star, loop):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if (g2num(gk.digonly(g)) & loop) == loop:
            return True
        if g in glist: return False
        glist.append(g)
    return False


def call_u_conflicts_d(G_star, H, checkrate=0):
    glist = [G_star]
    while True:
        g = dincrement_u(G_star, glist[-1])
        if gk.isedgesubset(g, H): return False
        if g in glist: return True
        glist.append(g)
    return True


def call_u_conflicts(G_star, H, checkrate=0):
    glist = [G_star]
    while True:
        # g = increment_u(G_star, glist[-1])
        g = directed_inc(G_star, glist[-1])
        if gk.isdedgesubset(g, H): return False
        g = bidirected_inc(g, glist[-1])
        if gk.isedgesubset(g, H): return False
        if g in glist: return True
        glist.append(g)
    return True


def call_u_conflicts2(G_star, H):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if gk.isedgesubset(g, H): return False, glist
        if g in glist: return True, glist
        glist.append(g)
    return True, glist


def call_u_equals2(G_star, glist, H):
    while True:
        g = increment_u(G_star, glist[-1])
        if g == H: return True
        if g in glist: return False
        glist.append(g)
    return False


def call_u_equals(G_star, H):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if g == H: return True
        if g in glist: return False
        glist.append(g)
    return False


def compact_call_undersamples(G_star, steps=None):
    glist = [ug2num(G_star)]
    lastgraph = G_star
    while True:
        g = increment_u(G_star, lastgraph)
        if ug2num(g) in glist: return glist
        glist.append(ug2num(g))
        lastgraph = g
    return glist


def cc_undersamples(G_star, steps=1):
    glist = [ug2num(G_star)]
    lastgraph = G_star
    for i in xrange(steps):
        g = increment_u(G_star, lastgraph)
        n = ug2num(g)
        if n in glist: return []
        glist.append(n)
        lastgraph = g
    return glist[-1]


def compatible(d1, d2):
    idx = scipy.where(scipy.array([[r == l for l in d2] for r in d1]))
    return idx


def compat(G):
    n = len(G)
    num = g2num(G)
    # sample all the graph for gStar
    star_l = all_undersamples(G)
    hits = {}
    # brute force all graphs
    for i in range(0, 2 ** (n ** 2)):
        tmpG = num2CG(i, n)
        tmp_l = all_undersamples(tmpG)
        c = compatible(tmp_l, star_l)
        if len(sum(c)) > 0: hits[i] = c
    return hits


def icompat(i, nodes):
    print i
    g = num2CG(i, nodes)
    return compat(g)


def ilength(i, nodes):
    print i
    g = num2CG(i, nodes)
    return len(call_undersamples(g))


def iall(i, nodes):
    print i
    g = num2CG(i, nodes)
    return compact_call_undersamples(g)


def cc_all(i, nodes, steps):
    # print i
    g = num2CG(i, nodes)
    return i, cc_undersamples(g, steps=steps)


def make_rect(l):
    max_seq = max(map(len, l))
    nl = []
    for e in l:
        e += [e[-1]] * (max_seq - len(e))
        nl.append(e)
    return nl


def uniqseq(l):
    s = []
    ltr = map(lambda *a: list(a), *l)
    for i in range(len(ltr)):
        s.append(len(np.unique(ltr[i])))


def loadgraphs(fname):
    g = zkl.load(fname)
    return g


def savegraphs(l, fname):
    zkl.save(l, fname)


def jason2graph(g):
    r = {}
    d = {1: set([(0, 1)]),
         2: set([(2, 0)]),
         3: set([(0, 1), (2, 0)])}
    for head in g:
        r[head] = {}
        for tail in g[head]:
            r[head][tail] = d[g[head][tail]]
    return r


def graph2jason(g):
    r = {}
    for head in g:
        r[head] = {}
        for tail in g[head]:
            if g[head][tail] == set([(0, 1)]):
                r[head][tail] = 1
            elif g[head][tail] == set([(2, 0)]):
                r[head][tail] = 2
            elif g[head][tail] == set([(0, 1), (2, 0)]):
                r[head][tail] = 3
    return r


def ring(n):
    g = {}
    for i in range(1, n):
        g[str(i)] = {str(i + 1): set([(0, 1)])}
    g[str(n)] = {'1': set([(0, 1)])}
    return g


def addAring(g):
    for i in range(1, len(g)):
        if str(i + 1) in g[str(i)]:
            g[str(i)][str(i + 1)].add((0, 1))
        else:
            g[str(i)][str(i + 1)] = set([(0, 1)])
    if '1' in g[str(len(g))]:
        g[str(len(g))]['1'].add((0, 1))
    else:
        g[str(len(g))]['1'] = set([(0, 1)])


def upairs(n, k):
    '''
    n unique nonsequential pairs
    '''
    s = set()
    for p in randint(n, size=(3 * k, 2)):
        if p[1] - p[0] == 1: continue
        s.add(tuple(p))
    l = [e for e in s]
    return l[:k]


def ringarcs(g, n):
    for edge in upairs(len(g), n):
        g[str(edge[0] + 1)][str(edge[1] + 1)] = set([(0, 1)])
    return g


def ringmore(n, m):
    return ringarcs(ring(n), m)


# talking about extra edges on top of the ring
def dens2edgenum(d, n=10): return int(d * n ** 2) - n


def edgenum2dens(e, n=10): return np.double(e + n) / n ** 2


def gtranspose(G):  # Transpose (rev. edges of) G
    GT = {u: {} for u in G}
    for u in G:
        for v in G[u]:
            GT[v][u] = set([(0, 1)])  # Add all reverse edges
    return GT


def scale_free(n, alpha=0.7, beta=0.25,
               delta_in=0.2, delta_out=0.2):
    g = nx.scale_free_graph(n, alpha=alpha,
                            beta=beta,
                            delta_in=delta_in, delta_out=delta_out)
    g = nx2graph(g)
    g = gtranspose(g)
    addAring(g)
    return g


def randH(n, d1, d2):
    g = bfu.ringmore(n, d1)
    pairs = [x for x in itertools.combinations(g.keys(), 2)]
    for p in np.random.permutation(pairs)[:d2]:
        g[p[0]].setdefault(p[1], set()).add((2, 0))
        g[p[1]].setdefault(p[0], set()).add((2, 0))
    return g
import functools
import sys
from collections import Counter
sys.path.append('tools')

# local package
from bfutils import *

nodes = 5 # number of nodes in the graph
PNUM = 75 # number of processes to use

template = "{0:9}{1:9}{2:9}{3:9}{4:9}{5:10}"

def wrapper_c(fold):
    return icompat(fold, nodes)
def wrapper_l(fold):
    return ilength(fold, nodes)
def wrapper_list(fold):
    return iall(fold, nodes)

def wrapper_u(fold, steps):
    return cc_all(fold, nodes, steps)

def make_rect(l):
    max_seq = max(map(len,l))
    for e in l:
        e += [e[-1]] * (max_seq - len(e))
    return l

def wrapper_unique(fold):
    counter = 0
    for i in results[number]:
        if i not in unique_appeared_graph[counter]:
            unique_appeared_graph[counter].append(i)
        counter += 1
    return 1

def wrapper_non_eliminate(fold):
	return tuple(results[fold][counter])
resultsmp=True
results=[]
pool=Pool(processes=PNUM)

print template.format(*('        u', ' all_uniq', '   unique', '  seen_Gu', ' converge',  ' uconverge'))
cumset = set()
clen = 0
for s in xrange(100):

    results=pool.map(functools.partial(wrapper_u, steps=s),
                     range(2**(nodes**2)))
    
    converged = len([e for e in results if e[1]==[]])
    notconverged = len(results) - converged

    if notconverged < 100 and notconverged > 0:
        survivors = [e for e in results if e[1]]
    if notconverged == 0: break
    results = filter(lambda (x): x[1] != [], results)

    r = set(results)
    d = r.difference(cumset)
    cumset = cumset.union(r)

    cl = 2**(nodes**2) - len(results) - clen
    clen += cl

    print template.format(*(s, len(r), len(d), len(cumset),
                            converged, notconverged))

pool.close()
pool.join()
import sys, os
sys.path.append('./tools/')
import traversal, bfutils, graphkit
import unknownrate as ur
from multiprocessing import Pool,Process, Queue, cpu_count, current_process, active_children
import functools
import zickle as zkl
import time, socket
import scipy

KEY='rasl_il_u2'
UMAX = 2
INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 1000 # stop traversing after growing equivalence class tothis size
REPEATS = 100
if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=60
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'mars':
    PNUM=20
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'hooke':
    PNUM=21
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'saturn':
    PNUM=12
    PNUM=max((1,PNUM/INPNUM))
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count())
    PNUM=max((1,PNUM/INPNUM))
print 'processes: ',PNUM, INPNUM

def multiprocess(argslist, ncpu):
    total = len(argslist)
    done = 0
    result_queue = Queue()
    jobs = []

    def ra_wrapper_(fold, n=10, k=10):
        scipy.random.seed()
        l = {}
        while True:
            try:
                g = bfutils.ringmore(n,k) # random ring of given density
                gs= bfutils.call_undersamples(g)
                for u in range(1,min([len(gs),UMAX])):
                    g2 = bfutils.undersample(g,u)
                    print fold,': ',traversal.density(g),':',
                    startTime = int(round(time.time() * 1000))
                    s = ur.liteqclass(g2, verbose=False, capsize=CAPSIZE)
                    endTime = int(round(time.time() * 1000))
                    print len(s), u
                    l[u] = {'eq':s,'ms':endTime-startTime}
            except MemoryError:
                print 'memory error... retrying'
                continue
            break
        result_queue.put( {'gt':g,'solutions':l} )

    while argslist != [] and done<10 :
        if len(active_children()) < ncpu:
            p = Process(target=ra_wrapper_,args=(argslist.pop(),))
            jobs.append(p)
            p.start()
            done+=1
            print "\r",float(done)/total,"%",
    #get results here
    res = [result_queue.get() for p in jobs]
    print res

def ra_wrapper(fold, n=10, k=10):
    scipy.random.seed()
    l = {}
    while True:
        try:
            g = bfutils.ringmore(n,k) # random ring of given density
            gs= bfutils.call_undersamples(g)
            for u in range(1,min([len(gs),UMAX])):
                g2 = bfutils.undersample(g,u)
                print fold,': ',traversal.density(g),':', traversal.density(g2),':',
                startTime = int(round(time.time() * 1000))
                #s = ur.liteqclass(g2, verbose=False, capsize=CAPSIZE)
                s = ur.eqclass(g2)
                endTime = int(round(time.time() * 1000))
                print len(s), u
                l[u] = {'eq':s,'ms':endTime-startTime}
        except MemoryError:
            print 'memory error... retrying'
            continue
        break
    return {'gt':g,'solutions':l}

def ra_wrapper_preset(fold, glist=[]):
    scipy.random.seed()
    l = {}
    while True:
        try:
            g = glist[fold]
            gs= bfutils.call_undersamples(g)
            for u in range(1,min([len(gs),UMAX])):
                g2 = bfutils.undersample(g,u)
                print fold,': ',traversal.density(g),':',
                startTime = int(round(time.time() * 1000))
                s = ur.liteqclass(g2, verbose=False, capsize=CAPSIZE)
                endTime = int(round(time.time() * 1000))
                print len(s), u
                l[u] = {'eq':s,'ms':endTime-startTime}
        except MemoryError:
            print 'memory error... retrying'
            continue
        break
    return {'gt':g,'solutions':l}

def killall(l):
    for e in l:
        e.join(timeout=0.001)
        if not e.is_alive():
            #print 'first result'
            for p in l:
                if p != e:
                    #print 'terminating ', p.name
                    p.terminate()
                    p.join()
                else:
                    p.join()
            return True
    return False

def fan_wrapper(fold,n=10,k=10):
    scipy.random.seed()
    curr_proc=current_process()
    curr_proc.daemon=False
    output = Queue()
    while True:
        try:
            g = bfutils.ringmore(n,k)
            gdens = traversal.density(g)
            g2 = bfutils.increment_u(g,g)
            #g2 = bfutils.undersample(g,2)
            def inside_wrapper():
                scipy.random.seed()
                try:
                    startTime = int(round(time.time() * 1000))
                    #s = traversal.v2g22g1(g2, capsize=CAPSIZE)
                    s = traversal.backtrack_more2(g2, rate=2, capsize=CAPSIZE)
                    endTime = int(round(time.time() * 1000))
                    print "{:2}: {:8} : {:4}  {:10} seconds".\
                        format(fold, round(gdens,3), len(s),
                               round((endTime-startTime)/1000.,3))
                    output.put({'gt':g,'eq':s,'ms':endTime-startTime})
                except MemoryError:
                    print 'memory error...'
		    raise
            pl = [Process(target=inside_wrapper) for x in range(INPNUM)]
            for e in pl: e.start()
            while True:
                if killall(pl): break
            r = output.get()
        except MemoryError:
            print 'memory error... retrying'
            for p in pl:
                p.terminate()
                p.join()
            continue
        break
    for p in pl: p.join()
    return r

densities = {5: [0.2],
             6: [0.2, .25, .3],
             7: [0.2, .25, .3],
             8: [0.15, 0.2, 0.25, 0.3],
             9: [.2],
             10:[.2],
             15:[0.2],
             20:[0.1],# 0.15, 0.2, 0.25, 0.3],
             25:[0.1],
             30:[0.1],
             35:[0.1],
             40:[0.1],
             50:[0.05, 0.1],
             60:[0.05, 0.1]}

for nodes in [5]:
    z = {}
    #pool=Pool(processes=PNUM)
    for dens in densities[nodes]:
        print "{:2}: {:8} : {:10} : {:10}  {:10}".format('id', 'densityi(G)', 'density(H)', 'eq class', 'time')
        e = bfutils.dens2edgenum(dens, n=nodes)

        eqclasses = []
        for x in pool.imap(functools.partial(ra_wrapper, n=nodes, k=e),
                           range(REPEATS)):
            eqclasses.append(x)
            z[dens] = eqclasses
            zkl.save(eqclasses,
                     socket.gethostname().split('.')[0]+\
                     '_nodes_'+str(nodes)+'_density_'+\
                     str(dens)+'_'+KEY+'_.zkl')

        print ''
        print '----'
        print ''
        #pool.close()
        #pool.join()
        #zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_'+KEY+'_.zkl')
import seaborn as sb
import zickle as zkl
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append('./tools/')
densities = [.15, 0.2, 0.25, 0.3, 0.35]

def gettimes(d):
    t = [x['ms'] for x in d]
    time  = map(lambda x: x/1000./60., t)
    return time

l = ['leibnitz_nodes_15_density_0.1_newp_.zkl',
     'leibnitz_nodes_20_density_0.1_newp_.zkl',
     'leibnitz_nodes_25_density_0.1_newp_.zkl',
     'leibnitz_nodes_30_density_0.1_newp_.zkl',
     'leibnitz_nodes_35_density_0.1_newp_.zkl']

alltimes_new = []
for fname in l:
    d = zkl.load(fname)
    alltimes_new.append(gettimes(d))

shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

plt.figure(figsize=[10,2])


g = sb.boxplot(alltimes_new,names=map(lambda x: str(int(x*100))+"",
                                      densities),
               widths=wds, color="Reds",fliersize=fliersz,
               linewidth=lwd,
               **{'positions':np.arange(len(densities))+shift,
                  'label':'MSL'})

# plt.plot(np.arange(len(densities))-shift,
#         map(np.median,alltimes_old), 'ro-', lw=0.5, mec='k')
# plt.plot(np.arange(len(densities))+shift,
#         map(np.median,alltimes_new), 'bo-', lw=0.5, mec='k')
g.figure.get_axes()[0].set_yscale('log')
plt.xlabel('number of nodes in a graph')
plt.ylabel('computation time (minutes)')
#plt.title('100 6 node graphs per density\n$G_2 \\rightarrow G_1$',
#          multialignment='center')
plt.subplots_adjust(right=0.99, left=0.2)
plt.legend(loc=0)
plt.show()
import pandas as pd
import pylab as pl
import seaborn as sb
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np

SBDIR = '~/soft/src/dev/tools/stackedBarGraph/'
GFDIR = '/na/home/splis/soft/src/dev/craft/gunfolds/tools/'
import sys, os

sys.path.append(os.path.expanduser(SBDIR))
sys.path.append(os.path.expanduser(GFDIR))

import zickle as zkl
from stackedBarGraph import StackedBarGrapher
SBG = StackedBarGrapher()

def gettimes(d):
    t = [x['ms'] for x in d]
    time  = map(lambda x: x/1000./60., t)
    return time

l = [(0.15, 'leibnitz_nodes_15_density_0.1_newp_.zkl'),
     (0.20, 'leibnitz_nodes_20_density_0.1_newp_.zkl'),
     (0.25, 'leibnitz_nodes_25_density_0.1_newp_.zkl'),
     (0.30, 'leibnitz_nodes_30_density_0.1_newp_.zkl'),
     (0.35, 'leibnitz_nodes_35_density_0.1_newp_.zkl')]

fig = pl.figure(figsize=[10,3])
#Read in data & create total column

d = zkl.load("hooke_nodes_6_g32g1_.zkl")#hooke_nodes_35_newp_.zkl")
densities = [.15, .20 ,.25, .30, .35]
d = {}
for fname in l:
    d[fname[0]] = zkl.load(fname[1])

def get_counts(d):
    eqc = [len(x['eq']) for x in d]
    keys = np.sort(np.unique(eqc))
    c = {}
    for k in keys:
        c[k] = len(np.where(eqc == k)[0])
    return c

# unique size
usz = set()
dc = {}
for u in densities:
    dc[u] = get_counts(d[u])
    for v in dc[u]:
        usz.add(v)

for u in densities:
    for c in usz:
        if not c in dc[u]:
            dc[u][c] = 0

A = []
for u in densities:
    A.append([dc[u][x] for x in np.sort(dc[u].keys())])
    
#print A
#A = np.array(A)
pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.color_palette("Paired",len(usz)))
#pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.dark_palette("#5178C7",len(usz)))
#pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.blend_palette(["mediumseagreen", "ghostwhite", "#4168B7"],len(usz)))
scalarMap = mpl.cm.ScalarMappable(norm = lambda x: x/np.double(len(usz)),
                                  cmap=pp)

d_widths = [.5]*len(densities)
d_labels = map(lambda x: str(int(x*100))+"%",densities)
#u = np.sort(list(usz))
d_colors = [scalarMap.to_rgba(i) for i in range(len(A[0]))]
#d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', '#fee090']

#ax = fig.add_subplot(211)
ax = plt.subplot2grid((3,1), (0, 0), rowspan=2)

SBG.stackedBarPlot(ax,
                   A,
                   d_colors,
                   xLabels=d_labels,
                   yTicks=3,
                   widths=d_widths,
                   gap = 0.005,
                   scale=False
)

for i in range(len(A)):
    Ai = [x for x in A[i] if x>0]
    y = [x/2.0 for x in Ai]
    for j in range(len(Ai)):
        if j>0:
            yy = y[j]+np.sum(Ai[0:j])
        else:
            yy = y[j]        
        pl.text(0.5*i-0.02,yy-1.2,str(Ai[j]),fontsize=12,zorder=10)
    
# #Set general plot properties
# sns.set_style("white")
# sns.set_context({"figure.figsize": (24, 10)})

# for i in np.sort(list(usz))[::-1]:
#     y = [100-dc[u][i] for u in np.sort(dc.keys())]
#     bottom_plot=sns.barplot(x=np.asarray(densities)*100, y=y)
#    #                         color=scalarMap.to_rgba(i))
#     #y = (sbd[i+1]-sbd[i])/2.+sbd[i]scala
#     #for j in range(len(sbd.Density)):
#     #    pl.text(j-0.1,y[j],'1',fontsize=16,zorder=i)

# #Optional code - Make plot look nicer
sns.despine(left=True)
# #Set fonts to consistent 16pt size
ax.set(xticklabels="")
for item in ([ax.xaxis.label, ax.yaxis.label] +
             #ax.get_xticklabels() +
             ax.get_yticklabels()):
    item.set_fontsize(12)



alltimes_new = []
for fname in l:
    dp = zkl.load(fname[1])
    alltimes_new.append(gettimes(dp))

shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

ax = plt.subplot2grid((3,1), (2, 0))
g = sb.boxplot(alltimes_new,names=map(lambda x: str(int(x*100))+"",
                                      densities),
               widths=wds, color="Reds",fliersize=fliersz,
               linewidth=lwd,
               **{'positions':np.arange(len(densities))+shift,
                  'label':'MSL'})

# plt.plot(np.arange(len(densities))-shift,
#         map(np.median,alltimes_old), 'ro-', lw=0.5, mec='k')
# plt.plot(np.arange(len(densities))+shift,
#         map(np.median,alltimes_new), 'bo-', lw=0.5, mec='k')
g.figure.get_axes()[1].set_yscale('log')
plt.xlabel('number of nodes in a graph')
plt.ylabel('computation time\n(minutes)')
#plt.title('100 6 node graphs per density\n$G_2 \\rightarrow G_1$',
#          multialignment='center')
#plt.subplots_adjust(right=0.99, left=0.2)
plt.legend(loc=0)
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() +
             ax.get_yticklabels()):
    item.set_fontsize(12)

pl.subplots_adjust(bottom=0.1,hspace=0.01,top=0.98)
# plt.show()

pl.show()
import seaborn as sb
import zickle as zkl
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.append('./tools/')
densities = [.15, 0.2, 0.25, 0.3, 0.35]

def gettimes(d):
    t = [x['ms'] for x in d]
    time  = map(lambda x: x/1000./60., t)
    return time

def getalltimes(data, densities=densities):
    alltimes = []
    for dens in densities:        
        alltimes.append(gettimes(data[dens]))
    return alltimes

def timesfromfile(fname):
    d = zkl.load(fname)
    return getalltimes(d)

#alltimes_old = timesfromfile("hooke_nodes_8_old_nomem.zkl")
#alltimes_new = timesfromfile("hooke_nodes_8_newp_.zkl")

#alltimes_old = timesfromfile("hooke_nodes_10_newp_.zkl")
#alltimes_old = timesfromfile("hooke_nodes_6_g32g1_.zkl")
#alltimes_new = timesfromfile("leibnitz_nodes_35_newp_.zkl")


shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

plt.figure(figsize=[10,2])

# g = sb.boxplot(alltimes_old,names=map(lambda x: str(int(x*100))+"%",
#                                      densities),
#               widths=wds, color="Reds", fliersize=fliersz, linewidth=lwd,
#               **{'positions':np.arange(len(densities))-shift,
#                  'label':'naive approach'})

g = sb.boxplot(alltimes_new,names=map(lambda x: str(int(x*100))+"",
                                      densities),
               widths=wds, color="Blues",fliersize=fliersz,
               linewidth=lwd,
               **{'positions':np.arange(len(densities))+shift,
                  'label':'MSL'})

# plt.plot(np.arange(len(densities))-shift,
#         map(np.median,alltimes_old), 'ro-', lw=0.5, mec='k')
# plt.plot(np.arange(len(densities))+shift,
#         map(np.median,alltimes_new), 'bo-', lw=0.5, mec='k')
g.figure.get_axes()[0].set_yscale('log')
plt.xlabel('density (% of 36 total possible edges)')
plt.ylabel('computation time (minutes)')
plt.title('100 6 node graphs per density\n$G_2 \\rightarrow G_1$',
          multialignment='center')
plt.subplots_adjust(right=0.99, left=0.2)
plt.legend(loc=0)
plt.show()
import pandas as pd
import pylab as pl
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np

import sys
sys.path.append('/na/home/splis/soft/src/dev/craft/gunfolds/tools/')
sys.path.append('/na/homes/splis/soft/src/dev/tools/stackedBarGraph/')
import zickle as zkl
from stackedBarGraph import StackedBarGrapher
SBG = StackedBarGrapher()

fig = pl.figure(figsize=[10,1.3])
#Read in data & create total column

d = zkl.load("hooke_nodes_6_g32g1_.zkl")#hooke_nodes_35_newp_.zkl")
densities = np.sort(d.keys())

def get_counts(d):
    eqc = [len(x['eq']) for x in d]
    keys = np.sort(np.unique(eqc))
    c = {}
    for k in keys:
        c[k] = len(np.where(eqc == k)[0])
    return c

# unique size
usz = set()
dc = {}
for u in densities:
    dc[u] = get_counts(d[u])
    for v in dc[u]:
        usz.add(v)

for u in densities:
    for c in usz:
        if not c in dc[u]:
            dc[u][c] = 0

A = []
for u in densities:
    A.append([dc[u][x] for x in np.sort(dc[u].keys())])
    
#print A
#A = np.array(A)
pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.color_palette("Paired",len(usz)))
#pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.dark_palette("#5178C7",len(usz)))
#pp = mpl.colors.LinearSegmentedColormap.from_list("t",sns.blend_palette(["mediumseagreen", "ghostwhite", "#4168B7"],len(usz)))
scalarMap = mpl.cm.ScalarMappable(norm = lambda x: x/np.double(len(usz)),
                                  cmap=pp)

d_widths = [.5]*len(densities)
d_labels = map(lambda x: str(int(x*100))+"%",densities)
#u = np.sort(list(usz))
d_colors = [scalarMap.to_rgba(i) for i in range(len(A[0]))]
#d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', '#fee090']

ax = fig.add_subplot(111)
SBG.stackedBarPlot(ax,
                   A,
                   d_colors,
                   xLabels=d_labels,
                   yTicks=3,
                   widths=d_widths,
                   gap = 0.005,
                   scale=False
)

for i in range(len(A)):
    Ai = [x for x in A[i] if x>0]
    y = [x/2.0 for x in Ai]
    for j in range(len(Ai)):
        if j>0:
            yy = y[j]+np.sum(Ai[0:j])
        else:
            yy = y[j]        
        pl.text(0.5*i-0.02,yy-1.2,str(Ai[j]),fontsize=12,zorder=10)
    
# #Set general plot properties
# sns.set_style("white")
# sns.set_context({"figure.figsize": (24, 10)})

# for i in np.sort(list(usz))[::-1]:
#     y = [100-dc[u][i] for u in np.sort(dc.keys())]
#     bottom_plot=sns.barplot(x=np.asarray(densities)*100, y=y)
#    #                         color=scalarMap.to_rgba(i))
#     #y = (sbd[i+1]-sbd[i])/2.+sbd[i]scala
#     #for j in range(len(sbd.Density)):
#     #    pl.text(j-0.1,y[j],'1',fontsize=16,zorder=i)

# #Optional code - Make plot look nicer
sns.despine(left=True)
# #Set fonts to consistent 16pt size
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(12)
pl.subplots_adjust(bottom=0.2)
# plt.show()
pl.show()
import sys, os
sys.path.append('./tools/')
import traversal, bfutils, graphkit
import unknownrate as ur
from multiprocessing import Pool,Process, Queue, cpu_count, current_process
import functools
import zickle as zkl
import time, socket
import scipy
import clingo as cg

UMAX = 1
INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 1000 # stop traversing after growing equivalence class tothis size
REPEATS = 100
if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=60
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'mars':
    PNUM=12
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'hooke':
    PNUM=21
    PNUM=max((1,PNUM/INPNUM))
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count())
    PNUM=max((1,PNUM/INPNUM))
print 'processes: ',PNUM, INPNUM

def wrapper_rate_agnostic(fold, n=10, k=10):
    scipy.random.seed()
    l = {}
    while True:
        try:
            g = bfutils.ringmore(n,k) # random ring of given density
            gs= bfutils.call_undersamples(g)
            for u in range(1,min([len(gs),UMAX])):
                g2 = bfutils.undersample(g,u)
                print fold,': ',traversal.density(g),':',
                startTime = int(round(time.time() * 1000))
                s = ur.iteqclass(g2, verbose=False)
                endTime = int(round(time.time() * 1000))
                print len(s)
                l[u] = {'eq':s,'ms':endTime-startTime}
        except MemoryError:
            print 'memory error... retrying'
            continue
        break
    return {'gt':g,'solutions':l}

def killall(l):
    for e in l:
        e.join(timeout=0.001)
        if not e.is_alive():
            #print 'first result'
            for p in l:
                if p != e:
                    #print 'terminating ', p.name
                    p.terminate()
                    p.join()
                else:
                    p.join()
            return True
    return False

def fan_wrapper(fold,n=10,k=10):
    scipy.random.seed()
    curr_proc=current_process()
    curr_proc.daemon=False
    output = Queue()
    while True:
        try:
            g = bfutils.ringmore(n,k)
            gdens = traversal.density(g)
            g2 = bfutils.increment_u(g,g)
            #g2 = bfutils.undersample(g,2)
            def inside_wrapper():
                scipy.random.seed()
                try:
                    startTime = int(round(time.time() * 1000))
                    s = traversal.v2g22g1(g2, capsize=CAPSIZE)
                    #s = traversal.backtrack_more2(g2, rate=2, capsize=CAPSIZE)
                    endTime = int(round(time.time() * 1000))
                    print "{:2}: {:8} : {:4}  {:10} seconds".\
                        format(fold, round(gdens,3), len(s),
                               round((endTime-startTime)/1000.,3))
                    output.put({'gt':g,'eq':s,'ms':endTime-startTime})
                except MemoryError:
                    print 'memory error...'
		    raise
            pl = [Process(target=inside_wrapper) for x in range(INPNUM)]
            for e in pl: e.start()
            while True:
                if killall(pl): break
            r = output.get()
        except MemoryError:
            print 'memory error... retrying'
            for p in pl:
                p.terminate()
                p.join()
            continue
        break
    for p in pl: p.join()
    return r

densities = {6: [0.2, 0.25, 0.3, 0.35],
             8: [0.3],
             10:[0.1],# 0.15, 0.2, 0.25, 0.3],
             15:[0.25, 0.3],
             20:[0.1],# 0.15, 0.2, 0.25, 0.3],
             25:[0.1],
             30:[0.1],
             35:[0.1],
             40:[0.1],
             50:[0.05, 0.1],
             60:[0.05, 0.1]}

for nodes in [15]:
    z = {}
    pool=Pool(processes=PNUM)
    for dens in densities[nodes]:
        print "{:2}: {:8} : {:10}  {:10}".format('id', 'density', 'eq class', 'time')
        e = bfutils.dens2edgenum(dens, n=nodes)
        eqclasses = pool.map(functools.partial(fan_wrapper, n=nodes, k=e),
                             range(REPEATS))
        z[dens] = eqclasses
        zkl.save(z[dens],
                 socket.gethostname().split('.')[0]+\
                     '_nodes_'+str(nodes)+'_density_'+str(dens)+'_newp_.zkl')
        print ''
        print '----'
        print ''
    pool.close()
    pool.join()
    zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_newp_.zkl')
# This is a use-case of the tools in the tools directory. The example defines a graph and shows how to generate a figure that shows the graph at different undersampling rates. Running the file in python (python dbnplot.py) generates a figure in figures folder: shipfig.pdf


# system packages
import os, sys
import numpy as np
from random import random
sys.path.append('tools')
import zickle as zkl
# local packages
import dbn2latex as d2l
import bfutils as bfu
from bfutils import jason2graph

def ring(n):
    g = {}
    g['1'] = {'1': {(0,1)}, '2': {(0,1)}}
    for i in range(2,n):
        g[str(i)] = {str(i+1): {(0,1)}}
    g[str(n)] = {'1': {(0,1)}}
    return g

def listplot(fname, mname='JJ', stl='', width=5, R=2):
    l = zkl.load(fname)
    l = l[:17*10]
    y = min(width,len(l))
    x = np.int(np.ceil(len(l)/float(y)))
    d2l.matrix_list(l,y,x, R=R, w_gap=1, h_gap=2, mname=mname, stl=stl)


g = {'1': {'2': set([(0, 1)]), '7': set([(0, 1)])},
     '10': {'1': set([(0, 1)]), '5': set([(0, 1)]), '9': set([(0, 1)])},
     '2': {'3': set([(0, 1)]),
           '4': set([(0, 1)]),
           '6': set([(0, 1)]),
           '7': set([(0, 1)])},
     '3': {'4': set([(0, 1)])},
     '4': {'1': set([(0, 1)]), '4': set([(0, 1)]), '5': set([(0, 1)])},
     '5': {'10': set([(0, 1)]),
           '5': set([(0, 1)]),
           '6': set([(0, 1)]),
           '8': set([(0, 1)]),
           '9': set([(0, 1)])},
     '6': {'2': set([(0, 1)]), '7': set([(0, 1)])},
     '7': {'8': set([(0, 1)])},
     '8': {'4': set([(0, 1)]),
           '7': set([(0, 1)]),
           '8': set([(0, 1)]),
           '9': set([(0, 1)])},
     '9': {'1': set([(0, 1)]),
           '10': set([(0, 1)]),
           '2': set([(0, 1)]),
           '6': set([(0, 1)]),
           '7': set([(0, 1)])}}



# generation of the output
g = {'1': {'2':set([(0,1)])},
     '2': {'3':set([(0,1)]), '1':set([(0,1)])},
     '3': {'4':set([(0,1)])},
     '4': {'1':set([(0,1)])},
}


#d2l.matrix_unfold(l[0],2,1,R=5, w_gap=1, h_gap=2, mname='TT1')

n = 20
dens = 0.07

for i in range(10):
    print i
    g = bfu.ringmore(n,bfu.dens2edgenum(dens,n))
    gl = bfu.all_undersamples(g)
    zkl.save(gl, 'list1.zkl')

    # output file
    foo = open('figures/shipfig_figure.tex', 'wb')
    sys.stdout = foo
    listplot('list1.zkl', width=17, R=6)
    sys.stdout = sys.__stdout__              # remember to reset sys.stdout!
    foo.flush()
    foo.close()
    PPP = os.getcwd()
    os.chdir('figures')
    os.system('pdflatex --shell-escape shipfig.tex 2>&1 > /dev/null')
    os.system('mv shipfig.pdf /tmp/all_graphs_'+str(n)+'nodes_density_'+str(dens)+'_'+str(i)+'.pdf')
    os.chdir(PPP)
import signal
import pprint
import time, socket
import numpy as np
import scipy
import functools, itertools
import progressbar as pb
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))
import copy

from multiprocessing import Pool,Process, Queue, cpu_count, current_process
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed
import linear_model as lm
import traversal as trv
import bfutils as bfu
import graphkit as gk
import zickle as zkl
import pc
import pylab as plt
import unknownrate as ur
NOISE_STD = '0.1'
DEPTH=2
URATE=3
DIST='beta'
BURNIN=100
SAMPLESIZE=2000
PARALLEL=True
POSTFIX='_rasl_u'+str(URATE)
EST = 'svar'
INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 100 # stop traversing after growing equivalence class tothis size
REPEATS = 20
if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=30
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'mars':
    PNUM=21
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'saturn':
    PNUM=12
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'hooke':
    PNUM=22
    PNUM=max((1,PNUM/INPNUM))
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count())
    PNUM=max((1,PNUM/INPNUM))
print 'processes: ',PNUM, INPNUM

def timeout(func, args=(), kwargs={}, timeout_duration=1, default=None):
    import signal

    class TimeoutError(Exception):
        pass

    def handler(signum, frame):
        raise TimeoutError()

    # set the timeout handler
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout_duration)
    try:
        result = func(*args, **kwargs)
    except TimeoutError as exc:
        result = default
    finally:
        signal.alarm(0)

    return result

def hamming_neighbors(v, step):
    l = []
    for e in itertools.combinations(range(len(v)),step):
        b = copy.copy(v)
        for i in e: b[i] = int(not b[i])
        l.append(b)
    return l

def find_nearest_reachable(g2, max_depth=4):
	s = ur.liteqclass(g2, capsize=CAPSIZE, verbose=False)	
	if s: return s
	step = 1
	n = len(g2)
	v = bfu.g2vec(g2)
	while True:
		l = hamming_neighbors(v,step)
		pbar = ProgressBar(widgets=['neighbors checked @ step '+str(step)+': ', Percentage(), ' '], maxval=len(l)).start()
		c = 0
		for e in l:
			g = bfu.vec2g(e,n)
			if not gk.scc_unreachable(g):
				s = ur.liteqclass(g, capsize=CAPSIZE, verbose=False)
			else:
				s = set()
			if s: return s
			pbar.update(c)
			c += 1
		pbar.finish()
		if step > max_depth:
			return set()
		step += 1


def wrapper(fold,n=10,dens=0.1, urate=URATE):
    scipy.random.seed()
    rate = urate

    r = None
    s = set()
    counter = 0
    while not s:
        scipy.random.seed()
        sst = 0.9
        r = None
        while not r:
            r = lm.getAring(n, dens, sst, False, dist=DIST)
            print sst,
            sys.stdout.flush()
            if sst < 0.03:
                sst -= 0.001
            else:
                sst -= 0.01
            if sst < 0: sst = 0.02
        #pprint.pprint(r['transition'].round(2),width=200)
        #d = zkl.load('leibnitz_nodes_'+str(n)+'_OCE_model_.zkl')
        #r = d[dens][fold]
        g = r['graph']
        true_g2 = bfu.undersample(g, rate-1)
        data = lm.drawsamplesLG(r['transition'], samples=BURNIN+SAMPLESIZE*2,
                                nstd=np.double(NOISE_STD))
        data = data[:,BURNIN:]
        if np.max(data) > 1000.:
            pprint.pprint(r['transition'].round(2),width=200)
            #raise ValueError
        startTime = int(round(time.time() * 1000))
        if EST=='pc':
            g2 = pc.dpc(data[:,::rate], pval=0.0001)
        elif EST=='svar':
            g2 = lm.data2graph(data[:,::rate])
        if trv.density(g2) < 0.7:
            print gk.OCE(g2,true_g2)
            #s = examine_bidirected_flips(g2, depth=DEPTH)
            s = find_nearest_reachable(g2, max_depth=1)
            #s = trv.v2g22g1(g2, capsize=CAPSIZE, verbose=False)
            #s = trv.edge_backtrack2g1_directed(g2, capsize=CAPSIZE)
            #s = timeout(trv.v2g22g1,
            #s = timeout(trv.edge_backtrack2g1_directed,
            #            args=(g2,CAPSIZE),
            #            timeout_duration=1000, default=set())
            print 'o',
            sys.stdout.flush()
            if -1 in s: s=set()
        endTime = int(round(time.time() * 1000))
        #if counter > 3:
        #    print 'not found'
        #    return None
        counter += 1
    print ''
    oce = [gk.OCE(bfu.num2CG(x,n),g) for x in s]
    cum_oce = [sum(x['directed'])+sum(x['bidirected']) for x in oce]
    idx = np.argmin(cum_oce)
    print "{:2}: {:8} : {:4}  {:10} seconds".\
          format(fold, round(dens,3), cum_oce[idx],
                 round((endTime-startTime)/1000.,3))
    #np.set_printoptions(formatter={'float': lambda x: format(x, '6.3f')+", "})
    #pprint.pprint(r['transition'].round(2))
    #np.set_printoptions()

    return {'gt':r,
            'eq':s,
            'OCE':oce[idx],
            'tries_till_found': counter,
            'estimate': g2,
            'graphs_tried': counter,
            'strength':sst+0.01,
            'ms':endTime-startTime}


def wrapgen(fold,n=10,dens=0.1):
    scipy.random.seed()
    rate = 2

    s = set()
    sst = 0.06
    r = None
    while not r:
        r = timeout(lm.getAring, args=(n, dens, sst, False),
                    timeout_duration=3)
        print sst,
        if sst < 0.03:
            sst -= 0.002
        else:
            sst -= 0.01
        if sst < 0: break
    print 'model '+str(fold)+' found \n'+str(r['transition'].round(2))
    sys.stdout.flush()
    return r

densities = {5: [0.25, 0.3, 0.35],
			 6: [0.25, 0.3, 0.35],
             8: [.15, .2, 0.25, 0.3],
             10:[0.3],
             15:[0.1],
             20:[0.1],
             25:[0.1],
             30:[0.1],
             35:[0.1]}

wrp = wrapper

for nodes in [6]:
    z = {}
    pool=Pool(processes=PNUM)
    for dens in densities[nodes]:
        print "{:2}: {:8} : {:10}  {:10}".format('id', 'density', 'OCE', 'time')

        if PARALLEL:
            errors = pool.map(functools.partial(wrp, n=nodes,
                                                dens=dens),
                              range(REPEATS))
            print 'done'
        else:
            errors = []
            for i in range(REPEATS):
                errors.append(wrp(i,n=nodes,dens=dens))
        print 'computed'
        z[dens] = errors
        zkl.save(z[dens],
                 socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_samples_'+str(SAMPLESIZE)+'_density_'+str(dens)+'_noise_'+NOISE_STD+'_OCE_b_'+EST+'_'+DIST+POSTFIX+'.zkl')
        print ''
        print '----'
        print ''
    pool.close()
    pool.join()
    zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_samples_'+str(SAMPLESIZE)+'_noise_'+NOISE_STD+'_OCE_b_'+EST+'_'+DIST+POSTFIX+'.zkl')
import sys
sys.path.append('./tools/')
import seaborn as sb
import zickle as zkl
import graphkit as gk
import bfutils as bfu
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


#

uplimit = 0.25

#

#d05 = zkl.load('leibnitz_nodes_5_samples_2000_noise_0.1_OCE_b_svar_beta_rasl_more.zkl')
#d05 = zkl.load('oranos_nodes_5_samples_2000_noise_0.1_OCE_b_svar_beta_rasl_more.zkl')
#d05 = zkl.load('leibnitz_nodes_6_samples_2000_noise_0.1_OCE_b_svar_beta_rasl.zkl')
d05 = zkl.load('leibnitz_nodes_6_samples_2000_noise_0.1_OCE_b_svar_beta_rasl.zkl')

def estOE(d):
    gt= d['gt']['graph']
    gt=bfu.undersample(gt,1)
    e = gk.OCE(d['estimate'],gt)
    N = np.double(len(gk.edgelist(gt))) +\
        np.double(len(gk.bedgelist(gt)))
    return (e['directed'][0]+e['bidirected'][0])/N

def estCOE(d):
    gt= d['gt']['graph']
    gt=bfu.undersample(gt,1)
    e = gk.OCE(d['estimate'],gt)
    n = len(gt)
    N = np.double(n**2+(n-1)**2/2.0\
                  -len(gk.edgelist(gt))
                  -len(gk.bedgelist(gt)))
    return (e['directed'][1]+e['bidirected'][1])/N

d = d05#d08_01
density = np.sort(d.keys())
n = len(d[density[0]][0]['gt']['graph'])
OE = [[gk.oerror(x) for x in d[dd]] for dd in density]
COE = [[gk.cerror(x) for x in d[dd]] for dd in density]

eOE = [[estOE(x) for x in d[dd]] for dd in density]
eCOE = [[estCOE(x) for x in d[dd]] for dd in density]

samplesnum = 20
denscol = [0.25]*samplesnum+[0.30]*samplesnum+[0.35]*samplesnum
OE = OE[0]+OE[1]+OE[2]
eOE = eOE[0]+eOE[1]+eOE[2]
COE = COE[0]+COE[1]+COE[2]
eCOE = eCOE[0]+eCOE[1]+eCOE[2]


OE = pd.DataFrame(np.asarray([denscol+denscol, OE+eOE, pd.Categorical(['RASL']*samplesnum*3+['SVAR']*samplesnum*3)]).T, columns=['density', 'time', 'OE'])

COE = pd.DataFrame(np.asarray([denscol+denscol, COE+eCOE, pd.Categorical(['RASL']*samplesnum*3+['SVAR']*samplesnum*3)]).T, columns=['density', 'time', 'COE'])


shift = 0.15
wds = 0.3
fliersz = 2
lwd = 1

plt.figure(figsize=[14,6])

plt.subplot(121)
ax = sb.boxplot(x="density", y="time", hue="OE",
                data=OE,
                palette="Set2",
                linewidth=lwd,
                width=wds,
                fliersize=fliersz)
sb.stripplot(x="density", y="time", hue="OE",
                data=OE,
                palette='Set2',
                size=4, jitter=True, edgecolor="gray")
#ax.figure.get_axes()[0].set_yscale('log')
plt.xticks([0,1,2],('25%','30%','35%'))

plt.ylim([-0.02,uplimit])
plt.xlabel('density (% of '+str(n**2)+' total possible edges)')
plt.ylabel('Edge omission error')
plt.title(str(samplesnum)+' '+str(n)+'-node graphs per density',
          multialignment='center')
plt.legend(loc=0)


plt.subplot(122)
ax = sb.boxplot(x="density", y="time", hue="COE",
                data=COE,
                palette="Set2",
                linewidth=lwd,
                width=wds,
                fliersize=fliersz)
sb.stripplot(x="density", y="time", hue="COE",
                data=COE,
                palette='Set2',
                linewidth=lwd,
                size=4, jitter=True, edgecolor="gray")

#ax.figure.get_axes()[0].set_yscale('log')
plt.xticks([0,1,2],('25%','30%','35%'))

plt.ylim([-0.02,uplimit])
plt.xlabel('density (% of '+str(n**2)+' total possible edges)')
plt.ylabel('Edge comission error')
plt.title(str(samplesnum)+' '+str(n)+'-node graphs per density',
          multialignment='center')
plt.legend(loc=2)

sb.set_context('poster')
plt.savefig('/tmp/RASL_simulation.svgz',transparent=False, bbox_inches='tight', pad_inches=0)
plt.show()
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
import sys

sys.path.append('./tools/')
from pathtree import PathTree
from ortools.constraint_solver import pywrapcp
from matplotlib.cbook import flatten
from functools import wraps
import numpy as np
import bisect
from sortedcontainers import SortedDict
import ipdb


class SolutionNotFoundInTime(Exception):
    pass


def ptloopnum(pt):
    """
    Given a PathTree object returns the number of loops in it
    :param pt: PathTree object
    :return: number of loops (n)
    """

    def ptn(pt, n=0):
        for e in pt.loopset:
            if type(e) is int:
                n += 1
                continue
            n += ptn(e, n=1)
        return n

    return ptn(pt)


def ptnodenum(pt):
    """
    Given a PathTree object returns the number of latents that comprise it
    :param pt: PathTree object
    :return: number of nodes (n)
    """
    n = pt.preset - 1

    def ptn(pt, n=0):
        for e in pt.loopset:
            if type(e) is int:
                n += e - 1
                continue
            n += ptn(e, n=1)
        return n

    return n + ptn(pt)


def ptelement(pt, w):
    """
    An element generated by a PathTree with a given weight setting
    :param pt: PathTree
    :param w: a list of weights
    :return: an integer
    """
    n = pt.preset

    def sumloops(pt, w):
        n = 0
        ls = list(pt.loopset)
        for i in range(len(ls)):
            if type(ls[i]) is int:
                n += w[i] * ls[i]
                continue
            n += w[i][0] * ls[i].preset \
                 + min(1, w[i][0]) * sumloops(ls[i], w[i][1])
        return n

    return n + sumloops(pt, w)


def weights_pt(pt, weights):
    c = [0]

    def crawl(pt, w, c):
        wl = []
        for e in pt.loopset:
            if type(e) is int:
                wl.append(w[c[0]])
                c[0] += 1
                continue
            ww = w[c[0]]
            c[0] += 1
            wl.append([ww, crawl(e, w, c)])
        return wl

    return crawl(pt, weights, c)


def extraloops_pt(pt, loops):  # loops are tuples (loop, weight)
    c = [0]

    def crawl(pt, l, c):
        first = [l[c[0]]]
        wl = []
        for e in pt.loopset:
            c[0] += 1
            if type(e) is int:
                wl.append(l[c[0]])
                continue
            wl.append(crawl(e, l, c))
        return first + [wl]

    return crawl(pt, loops, c)


def ptelement_extraloop(pt, w, eloops):
    """
    An element generated by a PathTree with a given weight setting and extra loops on each level
    :param pt: PathTree
    :param w: a list of list of weights
    :param eloops: a list of tuples with lengths of extra loops and their weights
    :return: an integer
    """
    n = pt.preset + eloops[0][0] * eloops[0][1]

    def sumloops(pt, w, lps):
        ls = list(pt.loopset)
        n = 0
        for i in range(len(ls)):
            if type(ls[i]) is int:
                n += w[i] * ls[i] + min(1, w[i]) * lps[i][0] * lps[i][1]
                continue
            n += w[i][0] * ls[i].preset \
                 + min(1, w[i][0]) * (lps[i][0][0] * lps[i][0][1] + sumloops(ls[i], w[i][1], lps[i][1]))
        return n

    return n + sumloops(pt, w, eloops[1])


def isptelement_el(el, pt, w, eloops):
    return el == ptelement_extraloop(pt, w, eloops)


def isptsubset_el(elist, pt, w, eloops):
    for i in range(elist[-1]):
        if isptelement_el(i, pt, w, eloops):
            if not i in elist:
                return False
    return True


def isrightpt(el, elist, pt, w, eloops):
    for i in range(elist[-1]):
        if isptelement_el(i, pt, w, eloops):
            if not i in elist:
                return False
        if i == el and not isptelement_el(i, pt, w, eloops):
            return False
    return True


def ptelements(pt, seqlen=100, verbose=False, maxloop=100):
    """
    Generate first `seqlen` elements from a pathtree
    :param pt: a path tree object from pathtree.py
    :param seqlen: number of elements to generate in ascending order
    :param verbose: whether to print debugging information
    :return: a list of elements
    """
    solver = pywrapcp.Solver("pt-elements")

    # declare variables
    weights = []
    N = ptloopnum(pt)
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))

    # declare constraints
    # solver.Add()

    # run the solver
    solution = solver.Assignment()
    solution.Add(weights)
    db = solver.Phase(weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    num_solutions = 0
    els = set()
    while solver.NextSolution():
        w = [x.Value() for x in weights]
        num_solutions += 1
        els.add(ptelement(pt, w))
        if len(els) == seqlen:
            break
    solver.EndSearch()

    # output solutions
    if verbose:
        print "num_solutions:", num_solutions
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()

    return list(els)


def isptelement(pt, element, verbose=False, maxloop=100):
    """
    Check if an integer element is in the weight set represented by the path tree
    :param pt: a path tree object from pathtree.py
    :param element: an integer to check for presence in the weight
    :param verbose: whether to print debugging information
    :return: True or False
    """
    solver = pywrapcp.Solver("isptelement")

    # declare variables
    weights = []
    N = ptloopnum(pt)
    if not N:
        return element == pt.preset
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))

    wpt = weights_pt(pt, weights)

    # declare constraints
    solver.Add(element == ptelement(pt, wpt))

    # run the solver
    solution = solver.Assignment()
    solution.Add(weights)
    db = solver.Phase(weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    solution_exists = False
    while solver.NextSolution():
        solution_exists = True
        break
    solver.EndSearch()

    # output solutions
    if verbose:
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()

    return solution_exists


def loops_and_weights(solver, loops, weights):
    """
    Add constraints to solver that make sure loops are not generated if subtree is not active due to a zero weight upstream
    :param solver:
    :param loops:
    :param weights:
    :return:
    """

    def recurse(s, l, w):
        for ww, ll in zip(w, l):
            if type(ww) is list:
                for e in flatten(ll):
                    s.Add((ww[0] == 0) <= (e == 0))
                recurse(s, ll[1:], ww[1:])
            else:
                for e in flatten(ll):
                    s.Add((ww == 0) <= (e == 0))

    recurse(solver, loops[1], weights)


def eloops_simplify(eloops):
    l = []
    for e in eloops:
        if type(e) is list:
            l.append(eloops_simplify(e))
        else:
            l.append(int(e[0].Value()))
    return l


def ptaugmented(pt, eloops):
    def augment(pt, ls):
        pre = pt.preset
        loop = pt.loopset
        s = set()
        if ls[0]:
            s.add(ls[0])
        for l, el in zip(loop, ls[1]):
            if type(l) is int:
                if not el:
                    s.add(l)
                else:
                    s.add(PathTree({el}, pre=l))
                continue
            s.add(augment(l, el))

        return PathTree(s, pre=pre)

    t = augment(pt, eloops)

    return t


def ptsubset(pt, elist):
    for i in range(elist[-1]):
        if isptelement(pt, i) and not i in elist:
            return False
    return True


def smallest_pt(ptlist):
    if ptlist:
        idx = np.argsort(map(ptnodenum, ptlist))
        sol = ptlist[idx[0]]
    else:
        sol = None
    return sol


def pairprint(pt1, pt2, k=40):
    print np.c_[pt2seq(pt1, k), pt2seq(pt2, k)]


def etesteq(pt1, pt2, k=100):
    a1 = np.asarray(pt2seq(pt1, k))
    a2 = np.asarray(pt2seq(pt2, k))
    return np.sum(a1 - a2) == 0


def keeptreegrow(pt, e, seq, cutoff=10, cap=1000):
    t = None
    while t is None:
        t = growtree(pt, e, seq, cutoff=cutoff)
        cutoff += 10
        if cutoff > cap:
            raise SolutionNotFoundInTime("Cannot keep the tree growing")
    return t


def add_element(d, pt):
    """
    Add a PathTree to dictionary d such that it is either appended to the list or added anew
    Args:
        d: a dictionary
        pt: a PathTree

    Returns:

    """
    key = ptnodenum(pt)
    if key in d:
        d[key].append(pt)
    else:
        d[key] = pt


def del_element(d, pt, key=None):
    """
    Delete a PathTree from dictionary d such that it is either removed from the list or the list that only contains one element is removed
    Args:
        d: a dictionary
        pt: a PathTree

    Returns:

    """
    if key is None:
        key = ptnodenum(pt)
    if len(d[key]) == 1:
        del d[key]
    else:
        d[key].remove(pt)


def swap_elements(d, pt1, pt2, key=None):
    del_element(d, pt1, key=key)
    add_element(d, pt2)


def seq2pt(seq, verbose=False, cutoff=100):
    if not seq:
        return None

    pt = PathTree({}, pre=seq[0])
    pts = SortedDict()  # PathTrees
    pts[ptnodenum(pt)] = [pt]

    for e in seq[1:]:
        e_is_in = False
        for key in pts:
            for pt in pts[key]:
                if verbose:
                    print e
                try:
                    newpt = keeptreegrow(pt, e, seq, cutoff=cutoff)
                    swap_elements(pts, pt, newpt, key=key)
                    e_is_in = True
                    break
                except SolutionNotFoundInTime:
                    continue
        if not e_is_in:
            newpt = PathTree({}, pre=e)
            add_element(d, newpt)

    return pt


def growtree(pt, element, ref_elements, verbose=False, maxloop=100, cutoff=100):
    """
    Add a loop with the minimal length to a path tree to enable it to generate a given element and still be a subset of a given list
    :param pt: a path tree object from pathtree.py
    :param element: an integer to check for presence in the weight
    :param ref_elements: a (finite) list that should be a superset of numbers generated by the new path tree, for numbers smaller than tosubset[-1]
    :param verbose: whether to print debugging information
    :return: a PathTree augmented with a new loop
    """
    solver = pywrapcp.Solver("loop_an_element")

    # PathTree already can generate that number. Just to foolproof
    if isptelement(pt, element):
        return pt

    # declare variables
    weights = []  # weights denoting how many times a loop is active (marginalized)
    loops = []  # extra loops that can be potentially added
    lweights = []  # weights for the extra loops (marginalized out in the end)
    ltuples = []  # tuple list to hold loops and weights together

    N = ptloopnum(pt)  # number of loops in the PathTree
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))

    for i in range(N + 1):
        w = solver.IntVar(0, maxloop, "lw[%04i]" % i)
        l = solver.IntVar(0, maxloop, "l[%04i]" % i)
        lweights.append(w)  # loop related weight
        loops.append(l)
        ltuples.append((l, w))

    eloops = extraloops_pt(pt, ltuples)
    ws = weights_pt(pt, weights)

    # declare constraints
    solver.Add(solver.MemberCt(ptelement_extraloop(pt, ws, eloops), ref_elements))
    solver.Add(element == ptelement_extraloop(pt, ws, eloops))  # make sure the element can be generated
    solver.Add(solver.Count(loops, 0, len(loops) - 1))  # only one loop is on
    solver.Add(solver.Count(lweights, 0, len(lweights) - 1))  # only one loop is weighted
    for i in range(len(lweights)):
        solver.Add((lweights[i] == 0) <= (loops[i] == 0))  # if a loop has weight zero then it can't be active
        # solver.Add(lweights[i] >= loops[i])
    loops_and_weights(solver, eloops, ws)  # if a subtree is off (weight zero) no need to add loops

    # run the solver
    solution = solver.Assignment()
    solution.Add(loops)
    db = solver.Phase(loops + lweights + weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    numsol = 0
    pts = []
    while solver.NextSolution():
        # print numsol,
        new_pt = ptaugmented(pt, eloops_simplify(eloops))
        if verbose:
            print "trying PathTree: ", new_pt
        if ptsubset(new_pt, ref_elements):
            pts.append(new_pt)
            if verbose:
                print "OK PathTree: ", pts[-1]
        numsol += 1
        if numsol >= cutoff:
            break
    solver.EndSearch()

    # output solutions
    if verbose:
        print "solutions:", numsol
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()
        print "for ", element, "solutions found ", numsol

    return smallest_pt(pts)


def pt2seq(pt, num):
    if not pt.loopset:
        return [pt.preset]
    i = 0
    s = set()
    while len(s) < num:
        if isptelement(pt, i, maxloop=10 * num):
            s.add(i)
        i += 1
    l = list(s)
    l.sort()
    return l


def s2spt(s):  # convert edge set to pt
    ss = set()
    for e in s:
        if type(e) is int:
            ss.add(PathTree({0}, pre={e}))
            continue
        ss.add(e)
    return ss


def spt_elements(spt, num):
    """
    Generate numbers from a set of PathTrees
    :param spt: set of PathTrees
    :param num: number of elements (from the first) to generate
    :return: list of num numbers
    """
    i = 0
    s = set()
    while len(s) < num:
        if issptelement(spt, i):
            s.add(i)
        i += 1
    return list(s)


def issptelement(spt, element):
    a = False
    for pt in s2spt(spt):
        a = a or isptelement(pt, element)
    return a
import sys

sys.path.append('./tools/')


def osumnum(s, num):
    return set(num + x for x in s)


def osumset(s1, s2):
    s = set()
    for e in s1:
        s.update(osumnum(s2, e))
    return s


class PathTree:
    def __init__(self, lset, pre=0):
        self.loopset = lset
        self.preset = pre

    def __add__(self, other):
        if type(other) is int:
            return PathTree(self.loopset, pre=self.preset + other)
        if type(other) is set:
            return PathTree(self.loopset, pre=osumnum(other, self.preset))
        return PathTree(self.loopset.union(other.loopset), pre=self.preset + other.preset)

    def __radd__(self, other):
        return self.__add__(other)

    def __str__(self):
        s = '('
        comma = False
        if not self.preset == 0:
            s += str(self.preset) + ' '
            comma = True
        if not self.loopset == {0}:
            if comma:
                s += ', '
            s += '<' + ', '.join(map(str, self.loopset)) + '>'
        s += ')'
        return s
import sys

sys.path.append('./tools/')
import traversal, bfutils
import numpy as np
from ortools.constraint_solver import pywrapcp

U = 2 # undersampling rate: 1 means no undersampling
N = 10 # number of nodes
k = 25 # number of extra edges

solver = pywrapcp.Solver("MSL")

# generate a random graph and undersample
g = bfutils.ringmore(N,k)
gdens = traversal.density(g)
g2 = bfutils.undersample(g,U-1)

# undersampled edges
dedgeu = {}
bedgeu = {}
for i in range(N):
    for j in range(N):
        dedgeu[(i,j)] = 0
        bedgeu[(i,j)] = 0
        v = str(i+1)
        w = str(j+1)
        if w in g2[v]:
            if (0,1) in g2[v][w]: dedgeu[(i,j)] = 1
            if (2,0) in g2[v][w]: bedgeu[(i,j)] = 1
        

# declare variables
edges = []
for i in range(N):
    e = []
    for j in range(N):
        e.append(solver.IntVar(0,1,"%i -> %i" % (i,j)))
    edges.append(e)

# path constraint
def apath(i,j,u, e=edges):
    n = len(e)
    if u <= 1: return [e[i][k]*e[k][j] for k in range(n)]
    l = []
    for k in range(n):
        for z in range(n):
            l.extend(map(lambda x: e[i][k]*x*e[z][j], apath(k,z,u-1)))
    return l                

# directed path constraint
def dcons(i, j, u, e=edges, s=solver):
    return s.Sum(apath(i,j,u,e=e))
# bidirected edge constraint (a balanced trek constraint)
def bcons(i, j, u, e=edges, s=solver):
    n = len(e)
    l = [e[k][i]*e[k][j] for k in range(n)]
    for ui in range(1,u):
        for k in range(n):
            l.extend([x*y for x in apath(k,i,ui,e=e) for y in apath(k,j,ui,e=e)])
    return s.Sum(l)
    
# declare constraints
for i in range(N):
    for j in range(N):
        # directed edge constraints
        de = dcons(i,j,U-1) 
        if dedgeu[(i,j)] == 1:
            solver.Add( 0 < de )
        else:
            solver.Add( 0 == de)

# bidirected edge constraints
for i in range(N):
    for j in range(i,N):
        if j == i: continue        
        be = bcons(i,j,U-1) #solver.Sum([edges[k][i]*edges[k][j] for k in range(N)])
        if bedgeu[(i,j)] == 1:
            solver.Add( 0 < be )
        else:
            solver.Add( 0 == be)


# run the solver
solution = solver.Assignment()
solution.Add([edges[i][j] for i in range(N) for j in range(N)])
collector = solver.AllSolutionCollector(solution)
solver.Solve(solver.Phase([edges[i][j] for i in range(N) for j in range(N)],
                          solver.CHOOSE_FIRST_UNBOUND,
                          solver.ASSIGN_MIN_VALUE),
                          [collector])
num_solutions = collector.SolutionCount()

# output solutions
print "num_solutions:", num_solutions
print "failures:", solver.Failures()
print "branches:", solver.Branches()
print "WallTime:", solver.WallTime()

if num_solutions > 0 and num_solutions < 5:
    for s in range(num_solutions):
        qval = [collector.Value(s, edges[i][j]) for i in range(N) for j in range(N)]
        for i in range(len(qval)):
            if qval[i]:
                e = np.unravel_index(i,[N,N])
                print e[0],"->",e[1]
        print
        print
import sys, os
import sys, os

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import graphkit as gk
import numpy as np
from random import shuffle
import ecj
from copy import deepcopy
from pathtree import PathTree, osumset
import ipdb


def g2lg(g):
    """
    Convert a data structure encoding the MSL-type graph into a structure encoding latents graph
    :return: a graph with integer indices and sets as weights

    Args:
        g (MSL-graph): 
    """
    edge_type = {(0, 1): 1, (2, 0): 2}
    edge_weight = {(0, 1): 1, (2, 0): 0}
    lg = {int(e): {int(c): {edge_type[w]: {edge_weight[w]}
                            for w in g[e][c]}
                   for c in g[e]}
          for e in g}
    return fix_selfloops(lg)


def fix_selfloops(g):
    for v in g:
        if v in g[v]:
            g[v][v] = {1: {PathTree({1})}}
    return g


def gtranspose(G):  # Transpose (rev. edges of) G
    GT = {u: {} for u in G}
    for u in G:
        for v in G[u]:
            if 1 in G[u][v]:
                GT[v][u] = {1: {1}}  # Add all reverse edges
    return GT


def parents(g, N):  # an inefficient way to do it
    t = gtranspose(g)  # Transpose the graph
    return {n: g[n][N][1]
            for n in t[N] if n != N}  # Return all children


def children(g, N):
    return {n: g[N][n][1]
            for n in g[N] if n != N}


def remove_node(g, N):
    del g[N]
    for V in g:
        if N in g[V]:
            del g[V][N]


def iterate_ws(ws):
    starts = []
    for e in ws:
        if type(e) is int:
            starts.append(e)
            continue
        if type(e) is set:
            starts.extend(list(e))
            continue
        starts.append(e.preset)
        starts.append([e.preset, iterate_ws(e.loopset)])
    return starts


def iterate_pt(pt):  # iterate over a path tree
    starts = [pt.preset]
    for e in pt.loopset:
        if type(e) is int:
            starts.append(e)
            continue
        if type(e) is set:
            starts.extend(list(e))
            continue
        starts.append(e.preset)
        starts.append([e.preset, iterate_ws(e.loopset)])
    return starts


def merge_weightsets(ab, ah, hb, hh):
    ws = osumset(ah, hb)
    if hh:
        ws = osumset(ws, hh)
    if not type(ws) is set:
        ws = {ws}
    return ws.union(ab)


def hide_node(g, H):
    """
    Removes a node from a graph taking care of the edges and weights as if the node was not observed
    :param g: input graph
    :param H: node to hide
    :return: the graph
    """

    gg = deepcopy(g)

    if not H in g:
        raise KeyError
    ch = children(g, H)
    pa = parents(g, H)
    if H in g[H]:
        sl = g[H][H][1]
    else:
        sl = set()
    remove_node(gg, H)

    for p in pa:
        for c in ch:
            if c in gg[p]:
                ab = gg[p][c][1]  # self loop
            else:
                ab = set()
            w = merge_weightsets(ab, pa[p], ch[c], sl)  # new weight set
            if c == p:  # a new loop is forming
                w = {PathTree(w)}
            gg[p][c] = {1: w}

    return gg


def degrees(nodes, g):
    return [len(parents(g, v)) + len(children(g, v)) for v in nodes]


def sortbydegree(nodes, g):
    idx = np.argsort(degrees(nodes, g))
    return list(np.asarray(nodes)[idx])


def hide_nodes(g, nodelist, dosort=True):
    nodeset = set()  # make sure not to delete a node twice
    if dosort: nodelist = sortbydegree(nodelist, g)
    gg = deepcopy(g)
    for n in nodelist:
        if n in nodeset: continue
        gg = hide_node(gg, n)
        nodeset.add(n)
    return gg


def hide_random(g, ratio):
    """
    Hire random modes in the `ratio` proportion from graph g
    :param g: input graph
    :param ratio: what percentage of nodes to hide
    :return: the graph with hidden variables
    """
    nodes = g.keys()
    shuffle(nodes)
    return hide_nodes(g, nodes[:int(len(g) * ratio)])


def print_ws(ws):
    print '{',
    for e in ws:
        print e, ', ',
    print '}'


def test_osumnum():
    assert osumnum(set(range(5)), 1) == set(range(1, 5 + 1))


def testcase(n):
    g1 = {1: {2: {1: {1}}, 4: {1: {1}}},
          2: {3: {1: {1}}, 7: {1: {1}}},
          3: {},
          4: {5: {1: {1}}},
          5: {3: {1: {1}}, 6: {1: {1}}},
          6: {5: {1: {1}}},
          7: {8: {1: {1}}},
          8: {2: {1: {1}}}}

    g2 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 6: {1: {1}}},
          3: {4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {},
          6: {7: {1: {1}}},
          7: {3: {1: {1}}}}

    g3 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 6: {1: {1}}},
          3: {4: {1: {1}}, 8: {1: {1}}},
          4: {5: {1: {1}}},
          5: {},
          6: {7: {1: {1}}},
          7: {2: {1: {1}}},
          8: {9: {1: {1}}},
          9: {10: {1: {1}}},
          10: {11: {1: {1}}},
          11: {3: {1: {1}}}}

    g4 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {}}

    g5 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    g6 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}, 6: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    g7 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}, 6: {1: {1}}},
          7: {8: {1: {1}}, 7: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    g8 = {1: {2: {1: {1}}, 5: {1: {1}}},
          2: {3: {1: {1}}, 2: {1: {1}}},
          3: {4: {1: {1}}},
          4: {8: {1: {1}}},
          5: {6: {1: {1}}},
          6: {7: {1: {1}}},
          7: {4: {1: {1}}},
          8: {9: {1: {1}}},
          9: {9: {1: {1}}, 10: {1: {1}}},
          10: {}}

    cases = [g1, g2, g3, g4, g5, g6, g7, g8]

    return fix_selfloops(cases[n])
import graph_tool as gt
from graph_tool import draw as gtd
import numpy as np

def lg2gt(g):
    gr = gt.Graph()
    vlabel = gr.new_vertex_property("string")
    verts = {}
    edges = {}
    for v in g:
        verts[v] = gr.add_vertex()
        vlabel[verts[v]] = str(v)
    gr.vertex_properties["label"] = vlabel
    for v in g:
        for w in g[v]:
            edges[(v,w)] = gr.add_edge(verts[v], verts[w])
    return gr

def plotg(g, layout='sfdp', pos=True):
    gg = lg2gt(g)
    if not pos:
        if layout=='fr':
            pos = gtd.fruchterman_reingold_layout(gg)
        else:
            pos = gtd.sfdp_layout(gg)
    else:
        pos = gg.new_vertex_property("vector<double>")
        n = gg.num_vertices()
        s = 2.0*np.pi/n
        for v in range(gg.num_vertices()):
            idx = int(gg.vertex_properties['label'][gg.vertex(v)]) - 1
            pos[gg.vertex(v)] = (n * np.cos(s * idx),
                                 n * np.sin(s * idx))

    gtd.graph_draw(gg, pos,
               vertex_text=gg.vertex_properties['label'],
               vertex_font_size=32,
               edge_pen_width=1,
               edge_marker_size=15,
               vertex_pen_width=1,
               vertex_fill_color=[0.62109375,
                                  0.875     ,
                                  0.23828125,
                                  1])
# tools to construct (random) graphs
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import ecj
import bfutils as bfu
import traversal as trv

import random as std_random
import numpy as np
import scipy
import networkx as nx
import jgraph as igr
from numpy.random import randint
from comparison import nx2graph

def edgelist(g): # directed
    '''
    return a list of tuples for edges of g
    '''
    l = []
    for n in g:
        l.extend([(n,e) for e in g[n] if (0,1) in g[n][e]])
    return l

def edgenumber(g):
    return sum([sum([len(g[y][x]) for x in g[y]]) for y in g])

def iedgelist(g): # directed iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    for v in g:
        for w in g[v]:
            if (0,1) in g[v][w]: yield (v,w)
def inedgelist(g): # missing directed iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    n = len(g)
    for v in g:
        for i in xrange(1,n+1):
            w = str(i)
            if not w in g[v]:
                yield (v,w)
            elif not (0,1) in g[v][w]:
                yield (v,w)
def ibedgelist(g): # directed iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    for v in g:
       for w in g[v]:
           if (2,0) in g[v][w]: yield (v,w)
def inbedgelist(g): # missing bidirected iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    for v in g:
       for w in g:
           if v!=w:
               if not w in g[v]:
                   yield (v,w)
               elif not (2,0) in g[v][w]:
                   yield (v,w)

def bedgelist(g): # bidirected edge list with flips
    l = []
    for n in g:
        l.extend([tuple(sorted((n,e))) for e in g[n] if (2,0) in g[n][e]])
    l = list(set(l))
    l = l + map(lambda x: (x[1],x[0]), l)
    return l

def rnd_edges(n):
    """ generate a random uniformly distributed mask
    """
    rnum = std_random.getrandbits(n**2)
    l = list(bin(rnum)[2:])
    l = ['0' for i in range(0,n**2 - len(l))] + l
    return l

def list2dbn(l):
    """ convert list of edge presences/absences (0,1) to a DBN graph
    """
    n = scipy.sqrt(len(l))
    l = scipy.reshape(map(int, l), [n, n])
    G = ecj.adj2DBN(l)
    return G

def list2CG(l):
    """ convert list of edge presences/absences (0,1) to a compressed
    graph (CG) representation
    """
    n = scipy.sqrt(len(l))
    l = scipy.reshape(map(int, l), [n, n])
    G = bfu.adj2graph(l)
    return G

def rnd_dbn(n): return list2dbn(rnd_edges(n))
def rnd_cg(n):  return list2CG(rnd_edges(n))

def rnd_adj(n, maxindegree=5):
    l = scipy.zeros([n,n])
    for u in range(0,n):
        cap = scipy.random.randint(min([n,maxindegree+1]))
        idx = scipy.random.randint(n,size=cap)
        l[u, idx] = 1
    return l

def sp_rnd_edges(n, maxindegree=5):
    '''
    a sparse set of random edges
    '''
    l = rnd_adj(n, maxindegree=maxdegree)
    return scipy.reshape(l, n**2)

def sp_rnd_dbn(n, maxindegree=3):
    '''
    a sparse random DBN graph
    '''
    l = sp_rnd_edges(n)
    return list2dbn(l)

def emptyG(n):
    A = [[0 for j in range(n)] for i in range(n)]
    return bfu.adj2graph(np.asarray(A))

def fullG(n):
    A = [[1 for j in range(n)] for i in range(n)]
    return bfu.adj2graph(np.asarray(A))


def CG2uCG(cg):
    """
    convert to an undirected graph
    """
    G = {}
    for u in cg:
        G[u] = cg[u].copy()
    for u in cg:
        for v in cg[u]:
            G[v][u] = cg[u][v]
    return G

def connected(cg):
    n = len(cg)
    return sum(1 for _ in ecj.traverse(CG2uCG(cg),'1')) == n

def sp_rnd_CG(n, maxindegree=3, force_connected=False):
    l = sp_rnd_edges(n, maxindegree=maxindegree)
    cg = list2CG(l)
    if force_connected:
        while not connected(cg):
            cg = list2CG(sp_rnd_edges(n, maxindegree=maxindegree))
    return cg

def CG2adj(G):
    n = len(G)
    A = [[0 for i in range(0,n)] for j in range(0,n)]
    for v in G:
        if G[v]:
            directed = [w for w in G[v] if (0,1) in G[v][w]]
            for w in directed:
                A[int(w)-1][int(v)-1] = 1
    A = np.double(np.asarray(A))
    return A

def g2ig(g):
    """
    Converts our graph represenataion to an igraph for plotting
    """
    t = scipy.where(CG2adj(g)==1)
    l = zip(t[0],t[1])
    ig = igr.Graph(l,directed=True)
    ig.vs["name"] = scipy.sort([u for u in g])
    ig.vs["label"] = ig.vs["name"]
    return ig

def superclique(n):
    g = {}
    for i in range(n):
        g[str(i+1)] = {str(j+1):set([(0,1),(2,0)])
                       for j in range(n) if j!=i}
        g[str(i+1)][str(i+1)] = set([(0,1)])
    return g

def complement(g):
    n = len(g)
    sq = superclique(n)
    for v in g:
        for w in g[v]:
            sq[v][w].difference_update(g[v][w])
            if not sq[v][w]: sq[v].pop(w)
    return sq

def gtranspose(G):                      # Transpose (rev. edges of) G
    GT = {u:{} for u in G}
    for u in G:
        for v in G[u]:
            if (0,1) in G[u][v]:
                GT[v][u] = set([(0,1)])        # Add all reverse edges
    return GT

def scale_free(n, alpha=0.7, beta=0.25,
               delta_in=0.2, delta_out=0.2):
    g = nx.scale_free_graph(n, alpha=alpha,
                            beta=beta,
                            delta_in=delta_in, delta_out=delta_out)
    g = nx2graph(g)
    g = gtranspose(g)
    addAring(g)
    return g

def ring(n, permute=False):
    g = {}
    names = [str(x+1) for x in range(n)]
    if permute: names = np.random.permutation(names) 
    for i in range(n-1):
        g[names[i]] = {names[i+1]: set([(0,1)])}
    g[names[n-1]] = {names[0]: set([(0,1)])}
    return g

def addAring(g):
    for i in range(1,len(g)):
        if str(i+1) in g[str(i)]:
            g[str(i)][str(i+1)].add((0,1))
        else:
            g[str(i)][str(i+1)] = set([(0,1)])
    if '1' in g[str(len(g))]:
        g[str(len(g))]['1'].add((0,1))
    else:
        g[str(len(g))]['1'] = set([(0,1)])

def upairs(n,k):
    '''
    n unique nonsequential pairs
    '''
    s = set()
    for p in randint(n, size=(3*k, 2)):
        if p[1]-p[0] == 1: continue
        s.add(tuple(p))
    l = [e for e in s]
    return l[:k]

def ringarcs(g,n):
    for edge in upairs(len(g),n):
        g[str(edge[0]+1)][str(edge[1]+1)] = set([(0,1)])
    return g
def ringmore(n,m, permute=False):
    return ringarcs(ring(n,permute=permute),m)

def digonly(H):
    """returns a subgraph of H contatining all directed edges of H

    Arguments:
    - `H`: undersampled graph
    """
    g = {n:{} for n in H}
    for v in g:
        g[v] = {w:set([(0,1)]) for w in H[v] if not H[v][w] == set([(2,0)])}
    return g

# Justin's ternary representation: 1 = directed edge; 2 = bidirected; 3 = both
def justin2graph(g):
    r = {}
    d = {1: set([(0,1)]),
         2: set([(2,0)]),
         3: set([(0,1),(2,0)]) }
    for head in g:
        r[head] = {}
        for tail in g[head]:
            r[head][tail] = d[g[head][tail]]
    return r

def graph2justin(g):
    r = {}
    for head in g:
        r[head] = {}
        for tail in g[head]:
            if g[head][tail] == set([(0,1)]):
                r[head][tail] = 1
            elif g[head][tail] == set([(2,0)]):
                r[head][tail] = 2
            elif g[head][tail] == set([(0,1),(2,0)]):
                r[head][tail] = 3
    return r

def OCE(g1,g2):
    '''
    omission/commision error of g1 referenced to g2
    '''
    s1 = set(edgelist(g1))
    s2 = set(edgelist(g2))
    omitted = len(s2 - s1)
    comitted = len(s1 - s2)

    s1 = set(bedgelist(g1))
    s2 = set(bedgelist(g2))
    bomitted = len(s2 - s1)
    bcomitted = len(s1 - s2)

    return {'directed': (omitted, comitted),
            'bidirected': (bomitted, bcomitted)}

def clean_leaf_nodes(g):
    for v in g: g[v] = {w:g[v][w] for w in g[v] if g[v][w]}

def cerror(d):
    return d['OCE']['directed'][1]/np.double(len(d['gt']['graph'])**2-len(edgelist(d['gt']['graph'])))

def oerror(d):
    return d['OCE']['directed'][0]/np.double(len(edgelist(d['gt']['graph'])))

def bidirected_no_fork(g):
    be = bedgelist(g)
    T = gtranspose(g)
    for e in be:
        if not set(T[e[0]].keys())&set(T[e[1]].keys()):
            return True
    return False

def fork_mismatch(g):
    be = bedgelist(g)
    benum = len(be)/2
    forknum = 0
    for v in g:
        fn = len([w for w in g[v] if (0,1) in g[v][w]])
        forknum += fn*(fn-1)/2.
    if benum < len(g)*(len(g)-1)/2.:
        return (forknum-benum) > benum
    else:
        return False

def no_parents(g):
    T = gtranspose(g)
    for n in T:
        if not T[n]: return True
    return False

def no_children(g):
    for n in g:
        if not g[n]: return True
    return False

def scc_unreachable(g):
    if bidirected_no_fork(g): return True
    if no_parents(g): return True
    if no_children(g): return True
    return False

# unlike functions from traversal package these do no checking
def addanedge(g,e): g[e[0]][e[1]] =  set([(0,1)])
def delanedge(g,e): g[e[0]].pop(e[1], None)
def addedges(g,es):
    for e in es: addanedge(g,e)
def deledges(g,es):
    for e in es: delanedge(g,e)

def checkequality(H,G_test, au = None):
    if not au:
        allundersamples = bfu.call_undersamples(G_test)
    else:
        allundersamples = au
    for graph in allundersamples:
        if graph == H: return True
    return False

def isdedgesubset(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                if not (0,1) in g2[n][h]:
                    return False
            else:
                    return False
    return True

def isedgesubset(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                #if not (0,1) in g2[n][h]:
                if not g2star[n][h].issubset(g2[n][h]):
                    return False
            else:
                    return False
    return True

def isedgesubset_(g,H):
    '''
    check if g edges are a subset of those of H
    '''
    for e in inbedgelist(H):
        if e[1] in g[e[0]] and (2,0) in g[e[0]][e[1]]: return False
    for e in inedgelist(H):
        if e[1] in g[e[0]] and (0,1) in g[e[0]][e[1]]: return False
    return True

def checkconflict(H,G_test, au = None):
    if not au:
        allundersamples = bfu.call_undersamples(G_test)
    else:
        allundersamples = au
    for graph in allundersamples:
        if isedgesubset(graph,H): return False
    return True

def checkconflict_(Hnum, G_test, au = None):
    if not au:
        allundersamples = bfu.call_undersamples(G_test)
    else:
        allundersamples = au
    #Hnum = bfu.ug2num(H)
    for graph in allundersamples:
        gnum = bfu.ug2num(graph)
        if gnum[0]&Hnum[0] == gnum[0] and gnum[1]&Hnum[1] == gnum[1]:
            return False
    return True
import zickle as zkl
alloops = zkl.load('allloops.zkl')
import scipy
from numpy.lib.arraysetops import intersect1d
from itertools import combinations
from testgraphs import *

def walk(G, s, S=set()):
    P, Q = dict(), set()
    P[s] = None
    Q.add(s)
    while Q:
        u = Q.pop()
        for v in G[u].difference(P,S):
            Q.add(v)
            P[v] = u
    return P

def traverse(G, s, qtype=set):
    S, Q = set(), qtype()
    Q.add(s)
    while Q:
        u = Q.pop()
        if u in S: continue
        S.add(u)
        for v in G[u]:
            Q.add(v)
        yield u

def dfs_topsort(G):
    S, res = set(), []
    def recurse(u):
        if u in S: return
        S.add(u)
        for v in G[u]:
            recurse(v)
        res.append(u)
    for u in G:
        recurse(u)
    res.reverse()
    return res

def tr(G):                      # Transpose (rev. edges of) G
    GT = {}
    for u in G: GT[u] = set()   # Get all the nodes in there
    for u in G:
        for v in G[u]:
            GT[v].add(u)        # Add all reverse edges
    return GT

def scc(G):                   # Kosaraju's algorithm
    GT = tr(G)                # Get the transposed graph
    sccs, seen = [], set()
    for u in dfs_topsort(G):   # DFS starting points
        if u in seen: continue # Ignore covered nodes
        C = walk(GT, u, seen)  # Don't go "backward" (seen)
        seen.update(C)         # We've now seen C
        sccs.append(C)         # Another SCC found
    return sccs


TT = {
    'A': ['B', 'C'],
    'B': ['D','E'],
    'C': [],
    'D': [],
    'E': []
}
def bfs_print_tree(tree,r):
    """
    A modified single list solution
    """
    Q = []
    idx = 1
    print str(idx)+': '+r
    Q.extend(tree[r])
    while Q:
        idx += 1
        print str(idx)+':',
        for u in range(0,len(Q)):
            e = Q.pop(0)
            print e,
            Q.extend(tree[e])
        print ''

def bfs_dict(tree,r):
    """
    Bob's suggested dictionary based solution
    """
    D = {}
    idx = 1
    D[idx] = [r]
    while D[idx]:
        idx += 1
        D[idx] = []
        for u in D[idx-1]: D[idx].extend(tree[u])
    D.pop(idx) # the last dictionary element is empty - must go
    for idx in D: print str(idx)+': '+' '.join(D[idx])


def cloneBfree(G):
    D = {}
    for v in G:
        D[v] = {}
        for u in G[v]:            
            if not (len(G[v][u].intersection(set([(0,1)]))) == 0):
                D[v][u] = set([(0,1)])
    return D

def clrbi(G):    
    for v in G:
        d = []
        for u in G[v]:
            try:
                G[v][u].remove((edge_type['bidirected'],0))
                if len(G[v][u]) == 0:
                    d.append(u)
            except KeyError:
                pass
        for e in d:
            G[v].pop(e)

def ecj(G,s,sccs=set()):           # elementary circuit by Johnson
    blocked = {v:False for v in G} # unblock all
    B = {v:[] for v in G}
    stack = []
    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            B[u].remove(w)
            if blocked[w]: unblock(w)
    def circuit(v, stack):
        f = False
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                f = True
                #print stack
                sccs.add(len(stack))
            elif not blocked[u]:
                if circuit(u, stack):
                    f = True
        if f:
            unblock(v)
        else:
            for w in G[v]:
                if v not in B[w]:
                    B[w].append(v)
        stack.pop()
        return f
    circuit(s,stack)
    return sccs

def ecj_loops(G,s,sl=set()):           # elementary circuit by Johnson
    blocked = {v:False for v in G} # unblock all
    B = {v:[] for v in G}
    stack = []
    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            B[u].remove(w)
            if blocked[w]: unblock(w)
    def circuit(v, stack):
        f = False
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                f = True
                #print scipy.sort(stack)
                sl.add(tuple(scipy.sort(stack)))
            elif not blocked[u]:
                if circuit(u, stack):
                    f = True
        if f:
            unblock(v)
        else:
            for w in G[v]:
                if v not in B[w]:
                    B[w].append(v)
        stack.pop()
        return f
    circuit(s,stack)
    return sl

def gcd(a, b):
    while b != 0: a, b = b, a%b
    return a

def listgcd(l):
    if len(l)>0:
        return gcd(l[0],listgcd(l[1:]))
    else:
        return 0
def lcm(a, b): return a*b/gcd(a,b)
def chmatch(n,m,delta):
    m,n = scipy.sort([n,m])
    sq = scipy.mod(range(n,lcm(n,m)+1,n),m)
    return scipy.mod(delta,m) in sq

def reachable(s, G, g):
    S, Q = set(), []
    Q.append(s)
    while Q:
        u = Q.pop()
        if u in S: continue
        if g in G[u]: return True
        S.add(u)
        Q.extend(G[u])
    return False

def allpaths(G, s, g, S=[]):
    if S is None: S = []
    S.append(s)
    if s == g:
        print S
    else:
        for u in G[s]:
            if u in S: continue
            allpaths(G,u,g,S)
    S.remove(s)

def lz_ecj(G,s,sccs=set()):           # elementary circuit by Johnson
    blocked = {v:False for v in G}    # unblock all
    B = {v:set() for v in G}
    stack = []
    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            if blocked[w]: unblock(w)
        B[u].clear()
    def circuit(v, stack):
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                print 'bottom'
                unblock(v)
                yield len(stack)
            elif not blocked[u]:
                print 'recurse'
                for x in circuit(u, stack):
                    unblock(v)
                    yield x
            else:
                print 'unmet'
                for w in G[v]: B[w].add(v)
        stack.pop()
#    circuit(s,stack)
    for v in circuit(s,stack): yield v

def iterate_allpaths(G, s, g, d=0, S=[], c=True):
    if S is None: S = []
    S.append(s)
    d += 1
    if s == g:
        if c:
            yield d-1
        else:
            yield list(S)
    else:
        for u in G[s]:
            if u in S: continue
            for v in iterate_allpaths(G,u,g,d,S,c):
                yield v
    S.remove(s)

def iddfs(G, s, g): # iterative depening DFS paths
    yielded = set()
    def recurse(G, s, g, d, S=None):
        if s not in yielded:
            yielded.add(s)
        if d == 0: return
        if S is None: S = []
        S.append(s)
        if s == g:
            yield list(S)
        else:
            for u in G[s]:
                if u in S: continue
                for v in recurse(G, u, g, d-1, S):
                    yield v
        S.remove(s)
    n = len(G)
    for d in range(n):
        #if len(yielded) == n: break
        for u in recurse(G, s, g, d):
            yield u

def reached_at_step(G, s, d):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    yielded = set()
    def recurse(G, s, d, B=None):
        if d == 0:
            if s not in yielded: # this avoids yielding duplicates
                yielded.add(s)
                yield s
            return
        if B is None: B = [] # black - backed out of this path
        for u in G[s]:
            #if u in B: continue
            if G[s][u] == (edge_type['bidirected'],0): continue
            for v in recurse(G, u, d-1, B):
                yield v
        B.append(s)
    for u in recurse(G, s, d):
            yield u

def d_trek(h,G,a,b,d):
    """
    Does there exist a trek with head h connecting a and b in d steps.
    """
    return set([a,b]).issubset(reached_at_step(G,h,d))

def d_biegde(G,a,b,d):
    """
    Do  a and  b  become connected  by  a bidirectional  edge after  d
    undersamples
    """
    for i in range(1,d+1):
        for u in G:
            if d_trek(u,G,a,b,i):
                return True
    return False

def undersample(G,d,bid=True):
    """
    
    """
    N = {}
    for u in G: 
        N.update({u:{v:set([(0,1)]) for v in reached_at_step(G,u,d+1)}})
    if bid:
        items = G.keys()
        for i in range(len(items)):
            for j in range(i+1, len(items)):
                u,v = items[i],items[j]
                if d_biegde(G,u,v,d):
                    try:
                        N[u][v].add((edge_type['bidirected'],0))
                    except KeyError:
                        N[u].update({v:set([(edge_type['bidirected'],0)])})
                    try:
                        N[v][u].add((edge_type['bidirected'],0))
                    except KeyError:
                        N[v].update({u:set([(edge_type['bidirected'],0)])})
    return N

import itertools as itt
def exist_equal_paths(h,G,a,b):
    Sa, Sb, Dff = set(), set(), set()
    ag=iterate_allpaths(G,h,a,0,[],True)
    bg=iterate_allpaths(G,h,b,0,[],True)
    for v in itt.izip_longest(ag,bg):
        print v
        Sa.add(v[0])
        Sb.add(v[1])
        if v[0] in Sb or v[1] in Sa: return True
    return False

# checks if there  exist exact length paths from the  head node to the
# nodes at question  by iterative deepining to avoid  oing through all
# paths
def iexist_equal_paths(h,G,a,b):
    Sa, Sb = set(), set()
    Pa, Pb = [], []
    ag=iddfs(G,h,a)
    bg=iddfs(G,h,b)
    for v in itt.izip(ag,bg):
        print v
        Sa.add(len(v[0])); Pa.append(v[0])
        Sb.add(len(v[1])); Pb.append(v[1])
        if len(v[0]) in Sb or len(v[1]) in Sa: return True
    return False

# check  if two  unequal  length  paths can  be  compensated by  their
# elementary cycles

def has_unit_cycle(G,path):
    for v in path:
        if v in G[v]: return True
    return False
def ecj_compat(G,p1,p2):
    n = len(p1)
    m = len(p2)
    p2, p1 = [[p1,p2][i] for i in scipy.argsort([n,m])]
    m,n = scipy.sort([n,m])
    delta = n - m
    if not delta: return True # redundant check for 0
    if has_unit_cycle(G,p2): return True # equivalent
    # if the shorter path does not have cycles they are not compatible
    # if the  longer path does not  have cycles: check  if the shorter
    #                                            path    has    cycles
    #                                            divisible by delta

    # otherwise start checking
    print p1, p2, n, m, delta


def wc(n):
    n = n*3
    a={str(v):set([str(v+1),str(v+2),str(v+3)]) for v in range(1,n,3)}
    b={str(v):set([str(v+2)]) for v in range(2,n,3)}
    c={str(v):set([str(v+1)]) for v in range(3,n,3)}
    a.update(b)
    a.update(c)
    a.update({str(n):set()})
    return a

# Frobenius number from here: http://cgi.gladman.plus.com/wp/?page_id=563 
def residue_table(a):
  n = [0] + [None] * (a[0] - 1)
  for i in range(1, len(a)):
    d = gcd(a[0], a[i])
    for r in range(d):
      try:
        nn = min(n[q] for q in range(r, a[0], d) if n[q] != None)
      except:
        continue
      if nn != None: 
        for c in range(a[0] // d):
          nn += a[i]
          p = nn % a[0]
          nn = min(nn, n[p]) if n[p] != None else nn
          n[p] = nn
  return n
 
def frobenius_number(a):
  return max(residue_table(sorted(a))) - min(a)
 
def isSclique(G):
    n = len(G)
    for v in G:
        if sum([(0,1) in G[v][w] for w in G[v]]) < n: return False
        if sum([(2,0) in G[v][w] for w in G[v]]) < n-1: return False
    return True

# Jianyu does not use bidirected edges
def isJclique(G):
    return (sum([len(G[w].keys()) for w in G]) == len(G)**2)

def directed_inc(G,D):
    G_un = {}
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in [el for el in D[v] if (0,1) in D[v][el]]:
            for e in G[w]:
                G_un[v][e] = set([(0,1)])
    return G_un
def bidirected_inc(G,D):
    # bidirected edges
    for w in D:
        # transfer old bidirected edges
        l = [e for e in D[w] if (2,0) in D[w][e]]
        for p in l:
            try: 
                 G[w][p].add((2,0))
            except KeyError: 
                 G[w][p] = set([(2,0)])
        # new bidirected dges
        l = [e for e in D[w] if (0,1) in D[w][e]]
        for p in list(combinations(l, 2)):
            try: 
                 G[p[0]][p[1]].add((2,0))
            except KeyError: 
                 G[p[0]][p[1]] = set([(2,0)])
            try: 
                G[p[1]][p[0]].add((2,0))
            except KeyError: 
                G[p[1]][p[0]] = set([(2,0)])
    return G
def increment_u(G_star, G_u):
    # directed edges
    G_un = directed_inc(G_star,G_u)
    # bidirected edges
    G_un = bidirected_inc(G_un,G_u)
    return G_un

def sample_graph(graph_g,steps=5):
    graph_g_list = [graph_g]
    for i in range(0,steps):
        g = increment_u(graph_g,graph_g_list[-1])
        graph_g_list.append(g)
    return graph_g_list
#BFS implementation of subgraph and supergraph
#Gu to G1 algorithm
from itertools import combinations, permutations
from functools import wraps
import copy
import time
import sys,os
import numpy as np
import ipdb
import operator
from scipy.misc import comb
import math
import gmpy as gmp
import gmpy as gmp
from scipy.misc import comb
import zickle as zkl
import simpleloops as sls
import math
import load_loops
from matplotlib.cbook import flatten
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed

TOOLSPATH='./tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))


circp = zkl.load('circular_p.zkl')
alloops = load_loops.alloops

import pprint
import bfutils as bfu
import traversal as trv
import graphkit as gk
import comparison as cmp
import simpleloops as sl

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = trv.gsig(args[0])         # Signature: just the g
        #s = tool.signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def prune_conflicts(H, g, elist):
    """checks if adding an edge from the list to graph g causes a
    conflict with respect to H and if it does removes the edge
    from the list

    Arguments:
    - `H`: the undersampled graph
    - `g`: a graph under construction
    - `elist`: list of edges to check
    """
    l  = []
    for e in elist:
        gk.addanedge(g,e)
        if not bfu.call_u_conflicts(g, H): l.append(e)
        gk.delanedge(g,e)
    return l

def eqclass(H):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    g = {n:{} for n in H}
    s = set()

    @memo
    def addedges(g,H,edges):
        if edges:
            nedges = prune_conflicts(H, g, edges)
            n = len(nedges)

            if n == 0: return None

            for i in range(n):
                gk.addanedge(g,nedges[i])
                if bfu.call_u_equals(g, H): s.add(bfu.g2num(g))
                addedges(g,H,nedges[:i]+nedges[i+1:])
                gk.delanedge(g,nedges[i])
    edges = gk.edgelist(gk.complement(g))
    addedges(g,H,edges)

    return s-set([None])

# these two functions come fromt his answer:
# http://stackoverflow.com/a/12174125
def set_bit(value, bit): return value | (1<<bit)
def clear_bit(value, bit): return value & ~(1<<bit)
def e2num(e,n): return (1<<(n*n +n - int(e[0],10)*n - int(e[1],10)))
def le2num(elist,n):
    num = 0
    for e in elist:
        num |= e2num(e,n)
    return num
def ekey2e(ekey,n):
    idx = np.unravel_index(n*n - ekey .bit_length() - 1 + 1,(n,n))
    idx = tuple([x+1 for x in idx])
    return ('%i %i'%idx).split(' ')

def cacheconflicts(num, cache):
    """Given a number representation of a graph and an iterable of
    conflicting subgraphs return True if the graph conflicts with any
    of them and false otherwise

    Arguments:
    - `num`: the number representation of a graph
    - `cache`: an iterable of number representations of conflicting
      graphs
    """
    conflict = False
    for c in cache:
        if num & c == c:
            return True
    return False

class nobar:
    def update(self,c): return None
    def finish(self): return None

def start_progress_bar(iter, n, verbose = True):
    if verbose:
        pbar = ProgressBar(widgets=['%3s' % str(iter) +
                                '%10s' % str(n)+' ',
                                Bar('-'), ' '],
                        maxval=n).start()
    else:
        pbar = nobar()
    return pbar

def add2set_loop(ds, H, cp, ccf, iter=1, verbose=True,
                 capsize=100, currsize=0):
    n = len(H)
    n2 = n*n +n
    dsr = {}
    s = set()
    ss = set()

    pbar = start_progress_bar(iter, len(ds), verbose = verbose)

    c = 0

    for gnum in ds:
        c += 1
        pbar.update(c)
        gset = set()
        eset = set()
        for sloop in ds[gnum]:
            if sloop & gnum == sloop: continue
            num = sloop | gnum
            if sloop in ccf and skip_conflictors(num,ccf[sloop]):continue
            if not num in s:
                g = bfu.num2CG(num, n)
                if not bfu.call_u_conflicts(g, H):
                    s.add(num)
                    gset.add((num,sloop))
                    eset.add(sloop)
                    if bfu.call_u_equals(g, H):
                        ss.add(num)
                        if capsize <= len(ss)+currsize: return dsr, ss

        #for gn,e in gset:
        #   if e in cp:
        #       dsr[gn] = eset - cp[e] - set([e])
        #   else:
        #       dsr[gn] = eset - set([e])
        for gn in gset: dsr[gn[0]] = eset - set([gn[1]])
    pbar.finish()
    return dsr, ss

def add2set_(ds, H, cp, ccf, iter=1, verbose=True, capsize=100):
    n = len(H)
    n2 = n*n +n
    dsr = {}
    s = set()
    ss = set()

    pbar = start_progress_bar(iter, len(ds), verbose = verbose)

    c = 0

    for gnum in ds:
        g = bfu.num2CG(gnum, n)
        c += 1
        pbar.update(c)
        glist = []
        elist = []
        eset = set()
        for e in ds[gnum]:
            if not e[1] in g[e[0]]:
                gk.addanedge(g,e)
                num = bfu.g2num(g)
                ekey = (1<<(n2 - int(e[0],10)*n - int(e[1],10)))
                if ekey in ccf and skip_conflictors(num,ccf[ekey]):
                    gk.delanedge(g,e)
                    continue
                if not num in s:
                    s.add(num)
                    if not bfu.call_u_conflicts(g, H):
                        #cf, gl2 = bfu.call_u_conflicts2(g, H)
                        #if not cf:
                        glist.append((num,ekey))
                        elist.append(e)
                        eset.add(ekey)
                        if bfu.call_u_equals(g, H):
                            ss.add(num)
                            #if bfu.call_u_equals2(g, gl2, H): ss.add(num)
                        if capsize <= len(ss): break
                gk.delanedge(g,e)

        for gn,e in glist:
            if e in cp:
                dsr[gn] = [ekey2e(k,n) for k in eset - cp[e]]
            else:
                dsr[gn] = elist
        if capsize <= len(ss): return dsr, ss

    pbar.finish()
    return dsr, ss

def skip_conflictors(gnum, ccf):
    pss = False
    for xx in ccf:
        if xx&gnum == xx:
            pss = True
            break
    return pss

def bconflictor(e,H):
    n = len(H)
    s = set()
    for v in H:
        s.add(e2num((v,e[0]),n)|e2num((v,e[1]),n))
    return s

def conflictor(e,H):
    n = len(H)
    def pairs(e):
        ekey = e2num(e,n)
        return [ekey|e2num((e[0],e[0]),n),
                ekey|e2num((e[1],e[1]),n)]

    def trios(e,H):
        s = set()
        for v in H:
            if not v in e:
                s.add(e2num((e[0],v),n)|
                      e2num((v,e[1]),n)|
                      e2num((e[1],e[1]),n))
                s.add(e2num((e[0],e[0]),n)|
                      e2num((e[0],v),n)|
                      e2num((v,e[1]),n))
                s.add(e2num((e[0],v),n)|
                      e2num((v,v),n)|
                      e2num((v,e[1]),n))
        return s

    return trios(e,H).union(pairs(e))

def conflictor_set(H):
    s = set()
    for x in gk.inedgelist(H):  s = s | conflictor(x,H)
    for x in gk.inbedgelist(H): s = s | bconflictor(x,H)
    return s

def conflictors(H):
    s = conflictor_set(H)
    ds = {}
    num = reduce(operator.or_,s)
    for i in xrange(gmp.bit_length(num)):
        if num & 1<<i:
            ds[1<<i] = [x for x in s if x&(1<<i)]
    return ds

def may_be_true_selfloop(n,H):
    for v in H[n]:
        if v == n: continue
        if (0,1) in H[n][v] and not ((2,0) in H[n][v]): return False
    return True

def issingleloop(num):
    bl = gmp.bit_length(num)
    idx = [1 for i in xrange(bl) if num & (1<<i)]
    return len(idx) == 1

def nonbarren(H):
    for v in H:
        if H[v]: return v
    return False

def prune_loops(loops, H):
    l = []
    n = len(H)
    for loop in loops:
        g = bfu.num2CG(loop, n)
        x = [k for k in g if g[k]]
        if len(x) == 1:
            s = reduce(lambda x, s: s.union(x),
                       [H[x[0]][w] for w in H[x[0]]])
            if not (2,0) in s: continue
        if not bfu.call_u_conflicts_d(g, H): l.append(loop)        
    return l

def lconflictors(H, sloops=None):
    if not sloops: sloops = prune_loops(allsloops(len(H)),H)
    s = conflictor_set(H)
    ds = {}
    num = reduce(operator.or_,s)
    for i in xrange(gmp.bit_length(num)):
        if num & 1<<i:
            cset = [x for x in s if x&(1<<i)]
            for sloop in sloops:
                if sloop & 1<<i:
                    ds.setdefault(sloop,[]).extend(cset)
    return ds

def confpairs(H):
    n = len(H)
    g = {n:{} for n in H}
    d = {}

    edges = gk.edgelist(gk.complement(g))
    edges = prune_conflicts(H, g, edges)

    for p in combinations(edges,2):
        gk.addedges(g,p)
        if bfu.call_u_conflicts(g, H):
            n1 = e2num(p[0],n)
            n2 = e2num(p[1],n)
            d.setdefault(n1,set()).add(n2)
            d.setdefault(n2,set()).add(n1)
        gk.deledges(g,p)

    return d

def lconfpairs(H, cap=10, sloops=None):
    n = len(H)
    d = {}
    if not sloops: sloops = prune_loops(allsloops(len(H)),H)
    c = 0
    for p in combinations(sloops,2):
        g = bfu.num2CG(p[0]|p[1], n)
        if bfu.call_u_conflicts(g, H):
            d.setdefault(p[0],set()).add(p[1])
            d.setdefault(p[1],set()).add(p[0])
        if c >= cap: break
        c +=1
    return d


def iteqclass(H, verbose=True, capsize=100):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    if cmp.isSclique(H):
        print 'not running on superclique'
        return None
    g = {n:{} for n in H}
    s = set()
    Hnum = bfu.ug2num(H)
    if Hnum[1]==0: s.add(Hnum[0])

    cp = confpairs(H)
    ccf = conflictors(H)

    edges = gk.edgelist(gk.complement(g))
    ds = {bfu.g2num(g): edges}

    if verbose: print '%3s'%'i'+'%10s'%' graphs'
    for i in range(len(H)**2):
        ds, ss = add2set_(ds, H, cp, ccf, iter=i,
                            verbose=verbose,
                            capsize=capsize)
        s = s | ss
        if capsize <= len(ss): break
        if not ds: break

    return s

def liteqclass(H, verbose=True, capsize=100, asl=None):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set([-1])
    g = {n:{} for n in H}
    s = set()



    if asl:
        sloops = asl
    else:
        sloops = prune_loops(allsloops(len(H)),H)

    cp  = []#lconfpairs(H, sloops=sloops)
    ccf = lconflictors(H, sloops=sloops)
    ds = {0: sloops}

    if verbose: print '%3s'%'i'+'%10s'%' graphs'
    i=0
    while ds:
        ds, ss = add2set_loop(ds, H, cp, ccf, iter=i,
                              verbose=verbose,
                              capsize=capsize,
                              currsize=len(s))
        s = s | ss
        i += 1
        if capsize <= len(s): break

    return s

def edgemask(gl,H, cds):
    """given a list of encoded graphs and observed undersampled graph
    H returns a matrix with -1 on diagonal, 0 at the conflicting graph
    combination and encoded graph at non-conflicted
    positions. Furthermore, returns a set of graphs that are in the
    equivalence class of H

    Arguments:
    - `gl`: list of integer encoded graphs
    - `H`: the observed undersampled graph
    """
    n = len(H)
    nl= len(gl)
    s = set()
    mask = np.zeros((nl,nl),'int')
    np.fill_diagonal(mask,-1)

    for i in xrange(nl):
        for j in xrange(i+1,nl):

            if gl[i] & gl[j]: continue
            if skip_conflict(gl[i], gl[j], cds): continue

            gnum = gl[i] | gl[j]
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                if bfu.call_u_equals(g, H): s.add(gnum)
                mask[i,j] = gnum
                mask[j,i] = gnum
    return mask, s

def ledgemask(gl,H, cds):
    """given a list of encoded graphs and observed undersampled graph
    H returns a matrix with -1 on diagonal, 0 at the conflicting graph
    combination and encoded graph at non-conflicted
    positions. Furthermore, returns a set of graphs that are in the
    equivalence class of H

    Arguments:
    - `gl`: list of integer encoded graphs
    - `H`: the observed undersampled graph
    """
    n = len(H)
    nl= len(gl)
    s = set()
    mask = np.zeros((nl,nl),'int')
    np.fill_diagonal(mask,-1)

    for i in xrange(nl):
        for j in xrange(i+1,nl):

            if gl[i] & gl[j]: continue
            gnum = gl[i] | gl[j]
            if skip_conflictors(gnum, cds): continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                if bfu.call_u_equals(g, H): s.add(gnum)
                mask[i,j] = gnum
                mask[j,i] = gnum
    return mask, s

def edgeds(mask):
    """construct an edge dictionary from the mask matrix

    Arguments:
    - `mask`:
    """
    ds = {}
    nl = mask.shape[0]
    idx = np.triu_indices(nl,1)
    for i,j in zip(idx[0], idx[1]):
        if mask[i,j]:
            ds[(i,j)] = set()
            conf = set([i,j])
            conf = conf.union(np.where(mask[i,:]==0)[0])
            conf = conf.union(np.where(mask[j,:]==0)[0])
            for k,m in zip(idx[0], idx[1]):
                if not mask[k,m]: continue
                if k in conf: continue
                if m in conf: continue
                if not (k,m) in ds: ds[(i,j)].add(mask[k,m])
            if not ds[(i,j)]: ds.pop((i,j))
    return ds

def edgedsg(mask):
    """construct an edge dictionary from the mask matrix

    Arguments:
    - `mask`:
    """
    ds = {}
    nl = mask.shape[0]
    idx = np.triu_indices(nl,1)
    for i,j in zip(idx[0], idx[1]):
        if mask[i,j]:
            ds[mask[i,j]] = set()
            conf = set([i,j])
            conf = conf.union(np.where(mask[i,:]==0)[0])
            conf = conf.union(np.where(mask[j,:]==0)[0])
            for k,m in zip(idx[0], idx[1]):
                if not mask[k,m]: continue
                if k in conf: continue
                if m in conf: continue
                if not (k,m) in ds: ds[mask[i,j]].add(mask[k,m])
            if not ds[mask[i,j]]: ds.pop(mask[i,j])
    return ds


def quadlister(glist, H, cds):
    n = len(H)
    s = set()
    cache = {}

    def edgemask(gl, H, cds):
        nl= len(gl)
        ss = set()
        mask = np.zeros((nl,nl),'int')
        np.fill_diagonal(mask,-1)
        idx = np.triu_indices(nl,1)
        for i,j in zip(idx[0], idx[1]):
            if gl[i] & gl[j]:
                mask[i,j] = -1
                mask[j,i] = -1
                continue
            if skip_conflict(gl[i], gl[j], cds):
                gnum = gl[i] | gl[j]
                cache[gnum] = False
                continue

            gnum = gl[i] | gl[j]
            if gnum in cache:
                if cache[gnum]:
                    mask[i,j] = gnum
                    mask[j,i] = gnum
            else:
                cache[gnum] = False
                g = bfu.num2CG(gnum,n)
                if not bfu.call_u_conflicts(g, H):
                    if bfu.call_u_equals(g, H): ss.add(gnum)
                    mask[i,j] = gnum
                    mask[j,i] = gnum
                    cache[gnum] = True
        return mask, ss


    def quadmerger(gl, H, cds):
        mask, ss = edgemask(gl, H, cds)
        ds = edgeds(mask)
        #ipdb.set_trace()
        return [[mask[x]]+list(ds[x]) for x in ds], ss

    l = []
    for gl in glist:
        ll, ss = quadmerger(gl, H, cds)
        l.extend(ll)
        s = s|ss

    return l, s


def dceqc(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cds = confpairs(H)

    glist =  [2**np.arange(n**2)]
    i = 1
    #for i in range(int(np.log2(n**2))):
    while glist != []:
        print i, np.max(map(len,glist)), len(glist)
        glist, ss = quadlister(glist, H, cds)
        s = s|ss
        i += 1
    return s

def quadmerge(gl, H, cds):
    n = len(H)
    l = set()
    s = set()
    mask, ss = edgemask(gl, H, cds)
    s = s | ss
    ds = edgeds(mask)

    #pp = pprint.PrettyPrinter(indent=1)
    #pp.pprint(ds)

    for idx in ds:
        for gn in ds[idx]:
            if mask[idx]&gn: continue
            if skip_conflict(mask[idx], gn, cds): continue
            gnum = mask[idx] | gn
            if gnum in l or gnum in ss: continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                l.add(gnum)
                if bfu.call_u_equals(g, H): s.add(gnum)

    return list(l), s

def skip_conflict(g1, g2, ds):
    pss = False
    for ekey in ds:
        if (g1 & ekey) == ekey:
            if ekey in ds and cacheconflicts(g2,ds[ekey]):
                pss = True
                break
    return pss

def edgemask2(gl,H, cds):
    n = len(H)
    nl= len(gl)
    s = set()
    o = set()
    mask = np.zeros((nl,nl),'int')
    np.fill_diagonal(mask,-1)
    for i in xrange(nl):
        for j in xrange(i+1,nl):
            if gl[i] & gl[j]: continue
            gnum = gl[i] | gl[j]
            if skip_conflictors(gnum, cds): continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                if bfu.call_u_equals(g, H): s.add(gnum)
                mask[i,j] = gnum
                mask[j,i] = gnum
            elif bfu.overshoot(g, H): o.add(gnum)
    return mask, s, o # mask, found eqc members, overshoots

def patchmerge(ds, H, cds):
    n = len(H)
    l = set()
    s = set()
    o = set()
    for gkey in ds:
        for num in ds[gkey]:
            if gkey & num: continue
            gnum = gkey | num
            if gnum is s: continue
            if skip_conflictors(gnum, cds): continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                l.add(gnum)
                if bfu.call_u_equals(g, H): s.add(gnum)
            elif not gnum in o and bfu.overshoot(g, H): o.add(gnum)
    return l, s, o

def quadmerge2(gl, H, cds):
    n = len(H)

    mask, s, o = edgemask2(gl, H, cds)
    #ipdb.set_trace()
    ds = edgedsg(mask)
    l, ss, oo = patchmerge(ds, H, cds)

    o = o | oo
    s = s | ss

    print 'overshoots: ', len(o)

    return list(l), s

def quadmerge21(gl, H, cds):
    n = len(H)
    l = set()

    mask, ss, o = edgemask2(gl, H, cds)
    idx = np.triu_indices(mask.shape[0], 1)
    print len(o)
    for i in range(len(idx[0])):
        if mask[idx[0][i],idx[1][i]]: l.add(mask[idx[0][i],idx[1][i]])

    return list(l), ss

def dceqclass2(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cp = confpairs(H)
    confs = conflictor_set(H)
    ccf = conflictors(H)

    def prune_loops(gl, H):
        l = []
        for e in gl:
            if e[0] == e[1] and not (e[1] in H[e[0]] and (1,0) in H[e[0]][e[1]]): continue
            l.append(e)
        return l
    edges = gk.edgelist(gk.complement(bfu.num2CG(0,n)))
    edges = prune_loops(edges, H)
    glist = map(lambda x: e2num(x,n),edges)

    #glist =  list(2**np.arange(n**2))
    i = 0
    while glist != []:
        print 2**i, len(glist)
        glist_prev = glist
        glist, ss = quadmerge21(glist, H, confs)
        s = s|ss
        i += 1


    ds = {x: edges for x in glist_prev}

    for j in range(i, len(H)**2):
        ds, ss = add2set_(ds, H, cp, ccf, iter=j, verbose=True)
        s = s | ss
        if not ds: break

    return s


def dceqclass(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cds = confpairs(H)

    glist =  [0]+list(2**np.arange(n**2))
    i = 1
    while glist != []:
        print i, len(glist)
        glist, ss = quadmerge(glist, H, cds)
        s = s|ss
        i += 1
    return s

def ldceqclass(H,asl=None):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cds = lconfpairs(H)
    if asl:
        sloops = asl
    else:
        sloops = prune_loops(allsloops(len(H)),H)

    glist =  sloops
    i = 1
    while glist != []:
        print i, len(glist)
        glist, ss = lquadmerge(glist, H, cds)
        s = s|ss
        i += 1
    return s

def lquadmerge(gl, H, cds):
    n = len(H)
    l = set()
    s = set()
    mask, ss = ledgemask(gl, H, cds)
    s = s | ss
    ds = edgeds(mask)

    #pp = pprint.PrettyPrinter(indent=1)
    #pp.pprint(ds)

    for idx in ds:
        for gn in ds[idx]:
            if mask[idx]&gn: continue
            if skip_conflictors(mask[idx], gn, cds): continue
            gnum = mask[idx] | gn
            if gnum in l or gnum in ss: continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                l.add(gnum)
                if bfu.call_u_equals(g, H): s.add(gnum)

    return list(l), s

def quadmerge_(glist, H, ds):
    n = len(H)
    gl = set()
    ss = set()
    conflicts = set()
    for gi in combinations(glist, 2):
        if gi[0] & gi[1]: continue
        #if skip_conflict(gi[0], gi[1], ds): continue
        gnum = gi[0] | gi[1]
        if gnum in conflicts: continue
        if skip_conflictors(gnum, ds):
            conflicts.add(gnum)
            continue
        if gnum in gl: continue
        g = bfu.num2CG(gnum,n)
        if not bfu.call_u_conflicts(g, H):
            gl.add(gnum)
            if bfu.call_u_equals(g, H): ss.add(gnum)
        else:
            conflicts.add(gnum)
    return gl, ss

def ecmerge(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return None
    n = len(H)
    s = set()
    ds = confpairs(H)
    ccf = conflictors(H)
    cset = set()
    for e in ccf:
        cset = cset.union(ccf[e])

    glist =  np.r_[[0],2**np.arange(n**2)]
    #glist =  2**np.arange(n**2)

    #glist, ss = quadmerge(glist,H)

    for i in range(int(2*np.log2(n))):
        print i, len(glist)
        glist, ss = quadmerge_(glist,H, cset)
        s = s | ss
    return s

def getrates(g,H):
    n = len(H)
    au = bfu.call_undersamples(g)
    return list(np.where(map(lambda x: x == H, au))[0])

def withrates(s,H):
    n = len(H)
    d = {g:set() for g in s}
    for g in s:
        d[g] = getrates(bfu.num2CG(g,n),H)
    return d

def add2set(gset, elist, H):
    n = len(H)

    s = set()
    ss = set()

    eremove = {e: True for e in elist}

    for gnum in gset:
        g = bfu.num2CG(gnum, n)
        for e in elist:
            if not e[1] in g[e[0]]:
                gk.addanedge(g,e)
                num = bfu.g2num(g)
                if not num in s:
                    au = bfu.call_undersamples(g)
                    if not gk.checkconflict(H, g, au=au):
                        eremove[e] = False
                        s.add(num)
                        if gk.checkequality(H, g, au=au): ss.add(num)
                gk.delanedge(g,e)

    for e in eremove:
        if eremove[e]: elist.remove(e)

    return s, ss, elist

def eqclass_list(H):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    g = {n:{} for n in H}
    s = set()

    edges = gk.edgelist(gk.complement(g))
    #edges = prune_conflicts(H, g, edges)

    gset = set([bfu.g2num(g)])
    for i in range(len(H)**2):
        print i
        gset, ss, edges = add2set(gset, edges, H)
        s = s | ss
        if not edges: break

    return s

def loop2graph(l,n):
    g = {str(i):{} for i in range(1,n+1)}
    for i in range(len(l)-1):
        g[l[i]][l[i+1]] = set([(0,1)])
    g[l[-1]][l[0]] = set([(0,1)])
    return g

def set_loop(loop, graph):
    for i in range(0,len(loop)-1):
        graph[loop[i]][loop[i+1]] = set([(0,1)])
    graph[loop[-1]][loop[0]] = set([(0,1)])

def rotate(l,n): return l[n:] + l[:n]
def get_perm(loop1, loop2, n=None):
    if not n: n = len(loop1)
    basel = [str(i) for i in xrange(1,n+1)]
    diff1 = set(basel) - set(loop1)
    diff2 = set(basel) - set(loop2)
    if loop1[0] in loop2:
        l2 = rotate(loop2, loop2.index(loop1[0]))
    else:
        l2 = loop2
    mp = {}
    for x,y in zip(loop1+list(diff1),l2+list(diff2)):
        mp[x] = y
    return mp

def permute(g, perm):
    gn = {x:{} for x in g}
    for e in g:
        gn[perm[e]] = {perm[x]:g[e][x] for x in g[e]}
    return gn

def permuteAset(s, perm):
    n = len(perm)
    ns = set()
    for e in s:
        ns.add(bfu.g2num(permute(bfu.num2CG(e,n),perm)))
    return ns

def noverlap_loops(loops):
    d = {}
    for l in loops:
        el = []
        for k in loops:
            if not set(l) & set(k):
                el.append(tuple(k))
                #d.setdefault(tuple(l),set()).add(tuple(k))
        d[tuple(l)] = noverlap_loops(el)
    return d

def loop_combinations(loops):
    s = set()
    d = noverlap_loops(loops)
    def dfs_traverse(d, gs):
        if d:
            for e in d:
                dfs_traverse(d[e], gs|set([e]))
        else:
            s.add(frozenset(gs))
    for e in d:
        dfs_traverse(d[e],set([e]))
    return list(s)

def sorted_loops(g):
    l = [x for x in sl.simple_loops(g,0)]
    s = {}
    for e in l:
        s.setdefault(len(e),[]).append(e)
    return s

def loopgroups(g, n=None):
    d = sorted_loops(g)
    if n:
        return loop_combinations(d[n])
    else:
        l=[]
        for key in d:
            l.append(loop_combinations(d[key]))
        return l

def count_loops(n):
    s = 0
    for i in range(1,n+1):
        s += comb(n,i) * math.factorial(i-1)
    return s

def perm_cyclic(l): return [tuple(l[i:]+l[:i]) for i in range(len(l))]
def hashloop(l):
    t = [int(x) for x in l]
    idx = np.argmin(t)
    return tuple(l[idx:]+l[:idx])
def perm_circular_slow2(l):
    s = [tuple(l)]
    c = {}
    c[hashloop(l)] = True
    for e in permutations(l):
        if not hashloop(e) in c:
            s.append(e)
            c[hashloop(e)] = True
    return s

def perm_circular_slow(l):
    s = [tuple(l)]
    c = set(perm_cyclic(l))
    for e in permutations(l):
        if not e in c:
            s.append(e)
            c = c | set(perm_cyclic(e))
    return s
def perm_circular(l, cp=circp):
    r = []
    n = len(l)
    for e in cp[n]:
        r.append([l[i] for i in e])
    return r

def gen_loops(n):
    l = [str(i) for i in range(1,n+1)]
    s = []
    for i in range(1,n+1):
        for e in combinations(l,i):
            s.extend(perm_circular(e))
    return s

def allsloops(n, asl = alloops):
    if asl: return asl[n]
    s = []
    l = gen_loops(n)
    for e in l:
        s.append(bfu.g2num(loop2graph(e,n)))
    return s

def reverse(H, verbose=True, capsize=1000):
    n = len(H)
    s = set()

    g      = gk.superclique(n)
    sloops = set(allsloops(n))

    ds = {bfu.g2num(g): sloops}

    if verbose: print '%3s'%'i'+'%10s'%' graphs'
    i=0
    while ds:
        ds, ss = del_loop(ds, H, iter=i,
                          verbose=verbose,
                          capsize=capsize)
        s = s | ss
        i += 1
        if capsize <= len(s): break

    return s

# ----------------------

def build_loop_step(ds, loop, n, iter=1):
    n2 = n*n +n
    dsr = {}
    s = set()
    ss = set()

    pbar = start_progress_bar(iter, len(ds))

    c = 0

    for gnum in ds:
        c += 1
        pbar.update(c)
        gset = set()
        eset = set()
        for sloop in ds[gnum]:
            num = sloop | gnum
            if not num in s:
                g = bfu.num2CG(num, n)
                s.add(num)
                if bfu.forms_loop(g, loop):
                    ss.add(num)
                else:
                    gset.add((num,sloop))
                    eset.add(sloop)

        for gn,e in gset:
            dsr[gn] = eset - set([e])

    pbar.finish()
    return dsr, ss

def forward_loop_match(loop, n):
    """start with an empty graph and keep adding simple loops until
    the loop is generated at some undersampling rate

    Arguments:
    - `loop`: binary encoding of the loop
    - `n`: number of nodes in the graph
    """
    s = set()
    sloops = allsloops(n)
    ds = {0: sloops}

    i=0
    while ds:
        ds, ss = build_loop_step(ds, loop, n, iter=i)
        s = s | ss
        i += 1

    return s

def delAloop(g, loop):
    n = len(g)
    l = []
    l = [bfu.g2num(ur.loop2graph(s,n)) for s in sls.simple_loops(g,0)]
    l = [num for num in l if not num == loop ]
    print loop, ': ',  l
    return bfu.num2CG(reduce(operator.or_, l),n)

def reverse_loop_match(g,loop):
    """start with a graph and keep removing loops while the loop is still matched

    Arguments:
    - `g`: graph that generates the loop
    - `loop`: the reference loop
    """
    s = set()
    n = len(g)

    def prune(g):
        numh = bfu.g2num(g)
        cannotprune = True
        for l in sls.simple_loops(gk.digonly(g),0):
            gg = delAloop(g,bfu.g2num(loop2graph(l,n)))
            if bfu.forms_loop(gg, loop):
                cannotprune = False
                prune(gg)
        if cannotprune:
            print 'one'
            s.add(g)

    prune(g)
    return s

def reverse_edge_match(g,loop):
    """start with a graph and keep removing loops while the loop is still matched

    Arguments:
    - `g`: graph that generates the loop
    - `loop`: the reference loop
    """
    s = set()
    n = len(g)

    def prune(g):
        numh = bfu.g2num(g)
        cannotprune = True
        for l in gk.edgelist(gk.digonly(g)):
            gk.delanedge(g,l)
            if bfu.forms_loop(g, loop):
                cannotprune = False
                prune(g)
            gk.addanedge(g,l)
        if cannotprune: s.add(bfu.g2num(g))

    prune(g)
    return s

def matchAloop(loop, n):
    """returns a set of minimal graphs that generate this loop

    Arguments:
    - `loop`: binary encoding of the loop
    - `n`: number of nodes in the graph
    """
    s = set()
    l = forward_loop_match(loop,n)
    print len(l)
    for g in l:
        s = s | reverse_edge_match(bfu.num2CG(g,n),loop)

    return s

# ----------------------

def del_loop(ds, H, iter=0, verbose=True, capsize=1000):
    n = len(H)

    dsr = {}
    s = set()
    ss = set()
    print iter,
    for gnum in ds:
        gset = []
        s = set()
        for sloop in ds[gnum]:
            rset = ds[gnum] - set([sloop])
            num = reduce(operator.or_, rset)
            if not num in s:
                g = bfu.num2CG(num, n)
                if bfu.overshoot(g, H):
                    s.add(num)
                    gset.append((num,rset))

        if gset == []:
            print '.',
            ss.add(gnum)

        for gn in gset: dsr[gn[0]] = gn[1]
    print ''
    return dsr, ss

def main():
    g = bfu.ringmore(6,1);
    H = bfu.undersample(g,1);
    ss = liteqclass(H)
    print ss

if __name__ == "__main__":
    main()
from networkx import strongly_connected_components
from functools import wraps
import scipy
import numpy as np
import itertools, copy, time
import random
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

# local packages
import bfutils as bfu
import graphkit as gk
import comparison
import ecj


def isedgesubsetD(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if (0,1) in g2star[n][h]:
                if h in g2[n]:
                    if not (0,1) in g2[n][h]:
                        return False
                else:
                    return False
    return True

def isedgesubset(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                #if not (0,1) in g2[n][h]:
                if not g2star[n][h].issubset(g2[n][h]):
                    return False
            else:
                    return False
    return True

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def purgepath(path, l):
    for i in range(1,len(path)-1):
        l.remove((path[i],path[i+1]))

def next_or_none(it):
    try:
        n = it.next()
    except StopIteration:
        return None
    return n

def try_till_d_path(g,d,gt,order=None):
    k = []
    i = 1
    while not k:
        if order:
            k = [x for x in length_d_paths(g,order(i),d)]
        else:
            k = [x for x in length_d_paths(g,str(i),d)]
        i += 1
        if i > len(g): return []

    ld = []
    for i in range(min(10, len(k))):
        ld.append(len(checkApath(['2']+k[i], gt)))
    idx = np.argmin(ld)
    return k[0]

def try_till_path(g, gt):
    d = len(g)-1
    gx = comparison.graph2nx(g)
    sccl = [x for x in strongly_connected_components(gx)]
    # take the largest
    ln = [len(x) for x in sccl]
    idx = np.argsort(ln)
    d = len(sccl[idx[-1]])-1
    sccl = [sccl[i] for i in idx[::-1]]
    order = [item for sublist in sccl for item in sublist]
    k = []
    while not k:
        if d < 5: return []
        k = try_till_d_path(g,d,gt)
        d -= 1
    return k

def gpurgepath(g,path):
    for i in range(1,len(path)-1):
        del g[path[i]][path[i+1]]

def forks(n,c,el,bl,doempty=lambda x,y: True):
    '''
    INPUT:
    n - single node string
    c - mutable list of children of node n (will be changed as a side effect)
    el - list of edges available for construction (side effect change)
    bl - list of bidirected edges

    OUTPUT:
    l - list of forks
    '''
    l = []
    r = set()
    for p in [x for x in itertools.combinations(c,2)]:
        if doempty(p, bl) and (n,p[0]) in el and (n,p[1]) in el:
            l.append((n,)+p)
            el.remove((n,p[0]))
            el.remove((n,p[1]))
            r.add(p[0])
            r.add(p[1])
    for e in r: c.remove(e)
    return l

def childrenedges(n,c,el,bl):
    '''
    INPUT:
    n - single node string
    c - mutable list of children of node n (will be changed as a side effect)
    el - list of edges available for construction (side effect change)
    bl - list of bidirected edges

    OUTPUT:
    l - list of
    '''
    l = []
    r = set()
    for ch in c:
        if (n,ch) in el:
            l.append((n,ch))
            el.remove((n,ch))
            r.add(ch)
    for e in r: c.remove(e)
    return l

# an empty fork is a for without the bidirected edge
def make_emptyforks(n,c,el,bl):
    return forks(n,c,el,bl,doempty=lambda x,y: not x in y)
def make_fullforks(n,c,el,bl):
    return forks(n,c,el,bl,doempty=lambda x,y: x in y)

def make_longpaths(g, el):
    '''
    el - list of edges that is modified as a side effect
    '''
    l = []
    gc = copy.deepcopy(g)
    for i in range(16):
        k = try_till_path(gc,g)
        if len(k) < 5: break
        if k:
            l.append(('2',)+tuple(k))
            purgepath(l[-1],el)
            gpurgepath(gc,l[-1])
        else:
            break
    return l

def make_allforks_and_rest(g,el,bl, dofullforks=True):
    '''
    el - list of edges that is modified as a side effect
    '''
    l = []
    r = []
    nodes = [n for n in g]
    random.shuffle(nodes)
    for n in nodes:

        c = [e for e in g[n] if (0,1) in g[n][e]]# all children
        if len(c) == 1:
            if (n,c[0]) in el:
                r.append((n,c[0]))
                el.remove((n,c[0]))
        elif len(c) > 1:
            l.extend(make_emptyforks(n,c,el,bl))
            if dofullforks: l.extend(make_fullforks(n,c,el,bl))
            r.extend(childrenedges(n,c,el,bl))
    return l, r

def vedgelist(g, pathtoo=False):
    """ Return a list of tuples for edges of g and forks
    a superugly organically grown function that badly needs refactoring
    """
    l = []
    el = gk.edgelist(g)
    bl = gk.bedgelist(g)

    if pathtoo: l.extend(make_longpaths(g,el))
    l2,r = make_allforks_and_rest(g,el,bl,dofullforks=True)
    l.extend(l2)

    A, singles = makechains(r)

    if singles:
        B, singles = makesinks(singles)
    else:
        B, singles = [], []

    l = longpaths_pick(l)+threedges_pick(l) + A + B + singles
    return l

def twoedges_pick(l):  return [e for e in l if len(e)==2]
def threedges_pick(l): return [e for e in l if len(e)==3]
def longpaths_pick(l): return [e for e in l if len(e)>3 and e[0]=='2']
def makechains(l):
    """ Greedily construct 2 edge chains from edge list
    """
    ends = {e[1]:e for e in l}
    starts = {e[0]:e for e in l}
    r = []
    singles = []
    while l:

        e = l.pop()
        if e[1] in starts and e[0] != e[1] and starts[e[1]] in l:
            r.append(('0', e[0],)+starts[e[1]])
            l.remove(starts[e[1]])
        elif e[0] in ends and e[0] != e[1] and ends[e[0]] in l:
            r.append(('0',)+ends[e[0]]+(e[1],))
            l.remove(ends[e[0]])
        else:
            singles.append(e)
    return r, singles
def makesink(es): return ('1', es[0][0],) + es[1]
def makesinks(l):
    """ Greedily construct 2 edge sinks ( b->a<-c ) from edge list
    """
    sinks = {}
    for e in l:
        if e[1] in sinks:
            sinks[e[1]].append(e)
        else:
            sinks[e[1]] = [e]
    r = []
    singles = []
    for e in sinks:
        if len(sinks[e])>1:
            for es in chunks(sinks[e],2):
                if len(es)==2:
                    r.append(makesink(es))
                else:
                    singles.append(es[0])
        else:
            singles.append(sinks[e][0])
    return r, singles

def selfloop(n,g):
    return n in g[n]

def selfloops(l,g):
    return reduce(lambda x,y: x and y, map(lambda x: selfloop(x,g), l))

def checkbedges(v,bel,g2):
    r = []
    for e in bel:
        if e == tuple(v[1:]) and not selfloops(e, g2):
            r.append(e)
        if e == (v[2],v[1]) and not selfloops(e, g2):
            r.append(e)
    for e in r: bel.remove(e)
    return bel

def checkedge(e, g2):
    if e[0] == e[1]:
        l = [n for n in g2 if n in g2[n]]
        s = set()
        for v in g2[e[0]]: s=s.union(g2[e[0]][v])
        if not (2,0) in s:
            l.remove(e[0])
        return l
    else:
        return g2.keys()

def single_nodes(v,g2):
    """ Returns a list of singleton nodes allowed for merging with v
    """
    l = [(n,n) for n in g2 if not n in v and len(g2[n])>1]
    return l

def checkvedge(v, g2):
    """ Nodes to check to merge the virtual nodes of v ( b<-a->c )
    """
    l = gk.bedgelist(g2)
    if (v[1],v[2]) in l:
        l = single_nodes(v,g2) + checkbedges(v,l,g2)
        for n in v:
            if n in g2[n]: l.append((n,n))
    else:
        l = checkbedges(v,l,g2)
    return list(set(l))

def checkAedge(v, g2):
    """ Nodes to check to merge the virtual nodes of A ( b->a<-c )
    """
    l = []
    # try all pairs but the sources
    for pair in itertools.combinations(g2,2):
        #if pair == (v[1],v[2]): continue
        #if pair == (v[2],v[1]): continue
        l.append(pair)
        l.append(pair[::-1])
    for n in g2:
        l.append((n,n))
    return l

def checkcedge(c, g2):
    """ Nodes to check to merge the virtual nodes of c ( a->b->c )
    """
    l = gk.edgelist(g2)
    return list(set(l))

def checkApath(p, g2):
    sl = [x for x in g2 if selfloop(x,g2)]
    d = len(p) - 2
    l = []
    for n in g2:
        l.extend([tuple(x) for x in length_d_loopy_paths(g2, n, d, p[1:])])
    #k = prunepaths_1D(g2, p, l)
    return l


def isedge(v):  return len(v) == 2 # a->b
def isvedge(v): return len(v) == 3 # b<-a->c
def isCedge(v): return len(v) == 4 and v[0] == '0' # a->b->c
def isAedge(v): return len(v) == 4 and v[0] == '1'# a->c<-b
def isApath(v):  return len(v) >= 4 and v[0] == '2'# a->b->...->z

def checker(n,ee):
    g = bfu.ringmore(n,ee)
    g2 = bfu.increment(g)
    d = checkable(g2)
    t = [len(d[x]) for x in d]
    r = []
    n = len(g2)
    ee= len(gk.edgelist(g2))
    for i in range(1,len(t)):
        r.append(sum(np.log10(t[:i])) - ee*np.log10(n))
    return r

def checkerDS(n,ee):
    g = bfu.ringmore(n,ee)
    g2 = bfu.increment(g)
    gg = checkable(g2)
    d,p,idx = conformanceDS(g2,gg,gg.keys())
    t = [len(x) for x in p]
    r = []
    n = len(g2)
    ee= len(gk.edgelist(g2))
    for i in range(1,len(t)):
        r.append(sum(np.log10(t[:i])) - ee*np.log10(n))
    return r

def fordens(n,denslist, repeats=100):
    rl={}
    for d in denslist:
        ee = bfu.dens2edgenum(d,n)
        l=[checker(n,ee)[-1] for i in range(repeats)]
        rl[d] = (round(scipy.mean(l),3),round(scipy.std(l),3))
    return rl

def fordensDS(n,denslist, repeats=100):
    rl={}
    for d in denslist:
        print d
        ee = bfu.dens2edgenum(d,n)
        l=[checkerDS(n,ee)[-1] for i in range(repeats)]
        rl[d] = (round(scipy.mean(l),3),round(scipy.std(l),3))
    return rl

def checkable(g2):
    d = {}
    g = cloneempty(g2)
    vlist = vedgelist(g2,pathtoo=False)
    for v in vlist:
        if isvedge(v):
            d[v] = checkvedge(v,g2)
        elif isCedge(v):
            d[v] = checkcedge(v,g2)
        elif isAedge(v):
            d[v] = checkAedge(v,g2)
        elif isApath(v):
            d[v] = checkApath(v,g2)
        else:
            d[v] = checkedge(v,g2)

    # # check if some of the otherwise permissible nodes still fail
    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge),
         (addaAedge,delaAedge),
         (addapath,delapath)]
    c = [ok2add2edges,
         ok2addavedge,
         ok2addacedge,
         ok2addaAedge,
         ok2addapath]

    for e in d:
        adder, remover = f[edge_function_idx(e)]
        checks_ok = c[edge_function_idx(e)]
        for n in d[e]:
            if not checks_ok(e,n,g,g2):
                d[e].remove(n)

    return d

def inorder_check2(e1, e2, j1, j2, g2, f=[], c=[]):
    g = cloneempty(g2) # the graph to be used for checking

    if f==[]:
        f = [(add2edges,del2edges,mask2edges),
             (addavedge,delavedge,maskavedge),
             (addacedge,delacedge,maskaCedge),
             (addaAedge,delaAedge,maskaAedge),
             (addapath,delapath,maskapath)]

    if c==[]:
        c = [ok2add2edges,
             ok2addavedge,
             ok2addacedge,
             ok2addaAedge,
             ok2addapath]

    adder, remover, masker = f[edge_function_idx(e1)]
    checks_ok = c[edge_function_idx(e2)]

    d = {}
    s1 = set()
    s2 = set()
    for c1 in j1: # for each connector
        mask = adder(g,e1,c1)
        d[c1] = set()
        for c2 in j2:
            if checks_ok(e2,c2,g,g2):
               d[c1].add(c2)
               s2.add(c2)
        remover(g,e1,c1,mask)
        if d[c1]: s1.add(c1)
    return d,s1,s2

def check3(e1, e2, e3, j1, j2, j3, g2, f=[], c=[]):
    g = cloneempty(g2) # the graph to be used for checking
    if f==[]:
        f = [(add2edges,del2edges,mask2edges),
             (addavedge,delavedge,maskavedge),
             (addacedge,delacedge,maskaCedge),
             (addaAedge,delaAedge,maskaAedge),
             (addapath,delapath,maskapath)]
    if c==[]:
        c = [ok2add2edges,
             ok2addavedge,
             ok2addacedge,
             ok2addaAedge,
             ok2addapath]

    adder1, remover1, masker1 = f[edge_function_idx(e1)]
    adder2, remover2, masker2 = f[edge_function_idx(e2)]

    checks_ok2 = c[edge_function_idx(e2)]
    checks_ok3 = c[edge_function_idx(e3)]

    d = {}
    s1 = set()
    s2 = set()
    s3 = set()

    for c1 in j1: # for each connector
        mask1 = adder1(g,e1,c1)
        append_set1 = False
        for c2 in j2:
            append_set2 = False
            if checks_ok2(e2,c2,g,g2):
                mask2 = adder2(g,e2,c2)
                for c3 in j3:
                    if checks_ok3(e3,c3,g,g2):
                        append_set1 = append_set2 = True
                        s3.add(c3)
                remover2(g,e2,c2,mask2)
            if append_set2: s2.add(c2)
        if append_set1: s1.add(c1)
        remover1(g,e1,c1,mask1)
    return s1, s2, s3

def del_empty(d):
    l = [e for e in d]
    for e in l:
        if d[e]==set(): del d[e]
    return d
def inorder_checks(g2, gg):
    #idx = np.argsort([len(gg[x]) for x in gg.keys()])
    #ee = [gg.keys()[i] for i in idx] # to preserve the order
    ee = [e for e in gg] # to preserve the order
    #cds = conformanceDS(g2, ee)
    #oo = new_order(g2, ee, repeats=100, cds=None)
    #ee = oo[0]
    random.shuffle(ee)
    d = {} # new datastructure
    d[ee[0]] = {('0'):gg[ee[0]]}
    for i in range(len(ee)-1):
        d[ee[i+1]] = del_empty(inorder_check2(ee[i], ee[i+1],
                                              gg[ee[i]], gg[ee[i+1]], g2)[0])
    return ee, d

def cloneempty(g): return {n:{} for n in g} # return a graph with no edges

def ok2addanedge1(s, e, g, g2,rate=1):
    """
    s - start,
    e - end
    """
    # directed edges
    # self-loop
    if s == e and not e in g2[s]: return False
    for u in g: # Pa(s) -> e
        if s in g[u] and not (e in g2[u] and (0,1) in g2[u][e]):
            return False
    for u in g[e]: # s -> Ch(e)
        if not (u in g2[s] and (0,1) in g2[s][u]):return False
    # bidirected edges
    for u in g[s]: # e <-> Ch(s)
        if u!=e and not (u in g2[e] and (2,0) in g2[e][u]):return False
    return True

def ok2addanedge2(s, e, g, g2, rate=1):
    mask = addanedge(g,(s,e))
    value = bfu.undersample(g,rate) == g2
    delanedge(g,(s,e),mask)
    return value

def ok2addanedge_sub(s, e, g, g2, rate=1):
    mask = addanedge(g,(s,e))
    value = isedgesubset(bfu.undersample(g,rate),g2)
    delanedge(g,(s,e),mask)
    return value

def ok2addanedge(s, e, g, g2, rate=1):
    f = [ok2addanedge1, ok2addanedge2]
    return f[min([1,rate-1])](s,e,g,g2,rate=rate)

def ok2addanedge_(s, e, g, g2, rate=1):
    f = [ok2addanedge1, ok2addanedge_sub]
    return f[min([1,rate-1])](s,e,g,g2,rate=rate)


def ok2add2edges(e,p,g,g2): return edge_increment_ok(e[0],p,e[1],g,g2)

def maskanedge(g,e): return [e[1] in g[e[0]]]
def mask2edges(g,e,p): return [p in g[e[0]], e[1] in g[p]]
def maskavedge(g,e,p):
    return [p[0] in g[e[0]], p[1] in g[e[0]],
            e[1] in g[p[0]], e[2] in g[p[1]]]
def maskaAedge(g,e,p):
    return [p[0] in g[e[1]], p[1] in g[e[2]],
            e[3] in g[p[0]], e[3] in g[p[1]]]
def maskaCedge(g,e,p):
    return [p[0] in g[e[1]], e[2] in g[p[0]],
            p[1] in g[e[2]], e[3] in g[p[1]]]
def maskapath(g,e,p):
    mask = []
    for i in range(len(p)):
        mask.append(p[i] in g[e[i+1]])
        mask.append(e[i+2] in g[p[i]])
    return mask
def maskaVpath(g,e,p):
    mask = []
    mask.extend([p[0] in g[e[0]], e[1] in g[p[-1]]])
    for i in range(1,len(p)):
        mask.append(p[i] in g[p[i-1]])
    return mask

def addanedge(g,e):
    '''
    add edge e[0] -> e[1] to g
    '''
    mask = maskanedge(g,e)
    g[e[0]][e[1]] =  set([(0,1)])
    return mask
def delanedge(g,e,mask):
    '''
    delete edge e[0] -> e[1] from g if it was not there before
    '''
    if not mask[0]: g[e[0]].pop(e[1], None)


def add2edges(g,e,p):
    '''
    break edge e[0] -> e[1] into two pieces
    e[0] -> p and p -> e[1]
    and add them to g
    '''
    mask = mask2edges(g,e,p)
    g[e[0]][p] = g[p][e[1]] = set([(0,1)])
    return mask

def del2edges(g,e,p,mask):
    '''
    restore the graph as it was before adding e[0]->p and p->e[1]
    '''
    if not mask[0]: g[e[0]].pop(p, None)
    if not mask[1]: g[p].pop(e[1], None)

def ok2addavedge(e,p,g,g2):
    if p[1] == e[0]:
        if p[0] != p[1] and p[0] != e[2] and not (e[2] in g2[p[0]] and (2,0) in g2[p[0]][e[2]]):
                return False
        if p[0] == p[1] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False
        if p[0] == e[1] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False

    if p[0] == e[0]:
        if p[0] != p[1] and p[1] != e[1] and not (e[1] in g2[p[1]] and (2,0) in g2[p[1]][e[1]]):
                return False
        if p[0] == p[1] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False
        if p[1] == e[2] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False

    if p[0] == e[1] and p[1] == e[2] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False
    if p[0] == e[2] and not (e[1] in g2[p[1]] and (0,1) in g2[p[1]][e[1]]):
            return False
    if p[1] == e[1] and not (e[2] in g2[p[0]] and (0,1) in g2[p[0]][e[2]]):
            return False
    if p[0] == p[1] == e[0] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False

    if not edge_increment_ok(e[0],p[0],e[1],g,g2): return False
    if not edge_increment_ok(e[0],p[1],e[2],g,g2): return False

    return  True

def addavedge(g,v,b):
    mask = maskavedge(g,v,b)
    g[v[0]][b[0]] = g[v[0]][b[1]] = g[b[0]][v[1]] = g[b[1]][v[2]] = set([(0,1)])
    return mask

def delavedge(g,v,b,mask):
    if not mask[0]: g[v[0]].pop(b[0], None)
    if not mask[1]: g[v[0]].pop(b[1], None)
    if not mask[2]: g[b[0]].pop(v[1], None)
    if not mask[3]: g[b[1]].pop(v[2], None)

def ok2addaAedge(e,p,g,g2):
    if p[1] == e[1] and not (p[0] in g2[e[2]] and (0,1) in g2[e[2]][p[0]]): return False
    if p[0] == e[2] and not (p[1] in g2[e[1]] and (0,1) in g2[e[1]][p[1]]): return False

    if not edge_increment_ok(e[1],p[0],e[3],g,g2): return False
    if not edge_increment_ok(e[2],p[1],e[3],g,g2): return False

    return True

def addaAedge(g,v,b):
    mask = maskaAedge(g,v,b)
    g[v[1]][b[0]] = g[v[2]][b[1]] = g[b[0]][v[3]] = g[b[1]][v[3]] = set([(0,1)])
    return mask

def delaAedge(g,v,b,mask):
    if not mask[0]: g[v[1]].pop(b[0], None)
    if not mask[1]: g[v[2]].pop(b[1], None)
    if not mask[2]: g[b[0]].pop(v[3], None)
    if not mask[3]: g[b[1]].pop(v[3], None)

def cleanedges(e,p,g, mask):
    i = 0
    for m in mask:
        if not m[0]: g[e[i+1]].pop(p[i], None)
        if not m[1]: g[p[i]].pop(e[i+2], None)
        i += 1
def cleanVedges(g, e,p, mask):

    if mask:
        if not mask[0]: g[e[0]].pop(p[0], None)
        if not mask[1]: g[p[-1]].pop(e[1], None)

        i = 0
        for m in mask[2:]:
            if not m: g[p[i]].pop(p[i+1], None)
            i += 1

def ok2addapath(e,p,g,g2):
    mask = []
    for i in range(len(p)):
        if not edge_increment_ok(e[i+1],p[i],e[i+2],g,g2):
            cleanedges(e,p,g,mask)
            return False
        mask.append(add2edges(g,(e[i+1],e[i+2]),p[i]))
    cleanedges(e,p,g,mask)
    return True

def ok2addaVpath(e,p,g,g2,rate=2):
    mask = addaVpath(g,e,p)
    if not isedgesubset(bfu.undersample(g,rate), g2):
        cleanVedges(g,e,p,mask)
        return False
    cleanVedges(g,e,p,mask)
    return True

def ok2addapath1(e,p,g,g2):
    for i in range(len(p)):
        if not edge_increment_ok(e[i+1],p[i],e[i+2],g,g2):
            return False
    return True

def addapath(g,v,b):

    mask = maskapath(g,v,b)
    s = set([(0,1)])
    for i in range(len(b)):
        g[v[i+1]][b[i]] = g[b[i]][v[i+2]] = s

    return mask

def addaVpath(g,v,b):
    mask = maskaVpath(g,v,b)
    s = set([(0,1)])
    l = [v[0]] + list(b) + [v[1]]
    for i in range(len(l)-1):
        g[l[i]][l[i+1]] = s
    return mask

def delaVpath(g, v, b, mask):
    cleanVedges(g, v, b, mask)

def delapath(g, v, b, mask):
    for i in range(len(b)):
        if not mask[2*i]: g[v[i+1]].pop(b[i], None)
        if not mask[2*i+1]:g[b[i]].pop(v[i+2], None)

def prunepaths_1D(g2, path, conn):
    c = []
    g = cloneempty(g2)
    for p in conn:
        mask = addapath(g,path,p)
        if isedgesubset(bfu.increment(g), g2): c.append(tuple(p))
        delapath(g,path,p,mask)
    return c

def ok2addacedge(e,p,g,g2):

    if p[0] == p[1]:
        if not e[2] in g2[e[2]]: return False
        if not p[0] in g2[p[0]]: return False
        if not (e[3] in g2[e[1]] and (0,1) in g2[e[1]][e[3]]): return False

    if not edge_increment_ok(e[1],p[0],e[2],g,g2): return False
    if not edge_increment_ok(e[2],p[1],e[3],g,g2): return False

    return True

def addacedge(g,v,b): # chain
    mask = maskaCedge(g,v,b)
    g[v[1]][b[0]] = g[v[2]][b[1]] = g[b[0]][v[2]] = g[b[1]][v[3]] = set([(0,1)])
    return mask

def delacedge(g,v,b,mask):
    if not mask[0]: g[v[1]].pop(b[0], None)
    if not mask[1]: g[b[0]].pop(v[2], None)
    if not mask[2]: g[v[2]].pop(b[1], None)
    if not mask[3]: g[b[1]].pop(v[3], None)

def rotate(l): return l[1:] + l[:1] # rotate a list
def density(g): return len(gk.edgelist(g))/np.double(len(g)**2)
def udensity(g): return (len(gk.edgelist(g))+len(gk.bedgelist(g))/2.)/np.double(len(g)**2 + len(g)*(len(g)-1)/2.)

def esig(l,n):
    '''
    turns edge list into a hash string
    '''
    z = len(str(n))
    n = map(lambda x: ''.join(map(lambda y: y.zfill(z),x)), l)
    n.sort()
    n = ''.join(n[::-1])
    return int('1'+n)

def gsig(g):
    '''
    turns input graph g into a hash string using edges
    '''
    return bfu.g2num(g)

def signature(g, edges): return (gsig(g),esig(edges,len(g)))

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def memo1(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = gsig(args[0])             # Signature: g
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def memo2(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return set()#cache[s]               # Return the cached solution
    return wrap

def eqsearch(g2, rate=1):
    '''Find  all  g  that are also in  the equivalence
    class with respect to g2 and the rate.
    '''

    s = set()
    noop = set()

    @memo1
    def addnodes(g,g2,edges):
        if edges:
            masks  = []
            for e in edges:
                if ok2addanedge_(e[0],e[1],g,g2,rate=rate):
                    masks.append(True)
                else:
                    masks.append(False)
            nedges = [edges[i] for i in range(len(edges)) if masks[i]]
            n = len(nedges)
            if n:
                for i in range(n):
                    mask = addanedge(g,nedges[i])
                    if bfu.undersample(g,rate) == g2: s.add(bfu.g2num(g))
                    addnodes(g,g2,nedges[:i]+nedges[i+1:])
                    delanedge(g,nedges[i],mask)
                return s
            else:
                return noop
        else:
            return noop


    g = cloneempty(g2)
    edges = gk.edgelist(gk.complement(g))
    addnodes(g,g2,edges)
    return s


def supergraphs_in_eq(g, g2, rate=1):
    '''Find  all supergraphs of g  that are also in  the same equivalence
    class with respect to g2 and the rate.
    Currently works only for bfu.undersample by 1
    '''
    if bfu.undersample(g,rate) != g2:
        raise ValueError('g is not in equivalence class of g2')

    s = set()

    def addnodes(g,g2,edges):
        if edges:
            masks  = []
            for e in edges:
                if ok2addanedge(e[0],e[1],g,g2,rate=rate):
                    masks.append(True)
                else:
                    masks.append(False)
            nedges = [edges[i] for i in range(len(edges)) if masks[i]]
            n = len(nedges)
            if n:
                for i in range(n):
                    mask = addanedge(g,nedges[i])
                    s.add(bfu.g2num(g))
                    addnodes(g,g2,nedges[:i]+nedges[i+1:])
                    delanedge(g,nedges[i],mask)

    edges = gk.edgelist(gk.complement(g))
    addnodes(g,g2,edges)
    return s

def edge_backtrack2g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    single_cache = {}

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            e = edges.pop()
            ln = [n for n in g2]
            for n in ln:
                if (n,e) in single_cache: continue
                mask = add2edges(g,e,n)
                if isedgesubset(bfu.increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and bfu.increment(r)==g2:
                        s.add(bfu.g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements in eqclass')
                del2edges(g,e,n,mask)
            edges.append(e)
        else:
            return g
    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:
            mask = add2edges(g,e,n)
            if not isedgesubset(bfu.increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def edge_backtrack2g1_directed(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    single_cache = {}

    def edgeset(g):
        return set(gk.edgelist(g))
    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            e = edges.pop()
            ln = [n for n in g2]
            for n in ln:
                if (n,e) in single_cache: continue
                mask = add2edges(g,e,n)
                if isedgesubsetD(bfu.increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and edgeset(bfu.increment(r))==edgeset(g2):
                        s.add(bfu.g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements in eqclass')
                del2edges(g,e,n,mask)
            edges.append(e)
        else:
            return g
    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:
            mask = add2edges(g,e,n)
            if not isedgesubsetD(bfu.increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def g22g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    single_cache = {}

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            if bfu.increment(g) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                return g
            e = edges[0]
            for n in g2:

                if (n,e) in single_cache: continue
                if not edge_increment_ok(e[0],n,e[1],g,g2): continue

                mask = add2edges(g,e,n)
                r = nodesearch(g,g2,edges[1:],s)
                del2edges(g,e,n,mask)

        elif bfu.increment(g)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements in eqclass')
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:

            mask = add2edges(g,e,n)
            if not isedgesubset(bfu.increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def backtrack_more(g2, rate=1, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    single_cache = {}
    if rate == 1:
        ln = [n for n in g2]
    else:
        ln = []
        for x in itertools.combinations_with_replacement(g2.keys(),rate):
            ln.extend(itertools.permutations(x,rate))
        ln = set(ln)

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            if bfu.undersample(g,rate) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                return g
            e = edges[0]
            for n in ln:

                if (n,e) in single_cache: continue
                if not ok2addaVpath(e, n, g, g2, rate=rate): continue

                mask = addaVpath(g,e,n)
                r = nodesearch(g,g2,edges[1:],s)
                delaVpath(g,e,n,mask)

        elif bfu.undersample(g,rate)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements in eqclass')
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in ln:

            mask = addaVpath(g,e,n)
            if not isedgesubset(bfu.undersample(g,rate), g2):
                single_cache[(n,e)] = False
            delaVpath(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def backtrackup2u(H,umax=2):
    s = set()
    for i in xrange(1,umax+1):
        s = s | backtrack_more(H,rate=i)
    return s

def vg22g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge),
         (addaAedge,delaAedge),
         (addapath,delapath)]
    c = [ok2add2edges,
         ok2addavedge,
         ok2addacedge,
         ok2addaAedge,
         ok2addapath]
    @memo2 # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            #key, checklist = edges.popitem()
            key = random.choice(edges.keys())
            checklist = edges.pop(key)
            adder, remover = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]
            for n in checklist:
                mask = adder(g,key,n)
                if isedgesubset(bfu.increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and bfu.increment(r)==g2:
                        s.add(bfu.g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements')
                remover(g,key,n,mask)
            edges[key] = checklist
        else:
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    chlist = checkable(g2)
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g,g2,chlist,s)
    except ValueError:
        s.add(0)
    return s

def edge_function_idx(edge):
    return min(4,len(edge))-2+min(max(3,len(edge))-3,1)*int(edge[0])

def v2g22g1(g2, capsize=None, verbose=True):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    f = [(add2edges,del2edges,mask2edges),
         (addavedge,delavedge,maskavedge),
         (addacedge,delacedge,maskaCedge),
         (addaAedge,delaAedge,maskaAedge),
         (addapath,delapath,maskapath)]
    c = [ok2add2edges,
         ok2addavedge,
         ok2addacedge,
         ok2addaAedge,
         ok2addapath]

    def predictive_check(g,g2,pool,checks_ok, key):
        s = set()
        for u in pool:
            if not checks_ok(key,u,g,g2): continue
            s.add(u)
        return s

    @memo2 # memoize the search
    def nodesearch(g, g2, order, inlist, s, cds, pool, pc):
        if order:
            if bfu.increment(g) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                s.update(supergraphs_in_eq(g, g2))
                return g

            key = order[0]
            if pc:
                tocheck = [x for x in pc if x in cds[len(inlist)-1][inlist[0]]]
            else:
                tocheck = cds[len(inlist)-1][inlist[0]]

            if len(order) > 1:
                kk = order[1]
                pc = predictive_check(g,g2,pool[len(inlist)],
                                      c[edge_function_idx(kk)],kk)
            else:
                pc = set()

            adder, remover, masker = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]

            for n in tocheck:
                if not checks_ok(key,n,g,g2): continue
                masked = np.prod(masker(g,key,n))
                if masked:
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                else:
                    mask = adder(g,key,n)
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                    remover(g,key,n,mask)

        elif bfu.increment(g)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements')
            return g

    @memo2 # memoize the search
    def nodesearch0(g, g2, order, inlist, s, cds):

        if order:
            key = order.pop(0)
            tocheck = cds[len(inlist)-1][inlist[0]]

            adder, remover, masker = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]

            if len(tocheck) > 1:
                for n in tocheck:
                    if not checks_ok(key,n,g,g2): continue
                    mask = masker(g,key,n)
                    if not np.prod(mask):
                        mask = adder(g,key,n)
                        r = nodesearch0(g,g2,order, [n]+inlist, s, cds)
                        if r and bfu.increment(r)==g2:
                            s.add(bfu.g2num(r))
                            if capsize and len(s)>capsize:
                                raise ValueError('Too many elements')
                        remover(g,key,n,mask)
                    else:
                        r = nodesearch0(g,g2,order, [n]+inlist, s, cds)
                        if r and bfu.increment(r)==g2:
                            s.add(bfu.g2num(r))
                            if capsize and len(s)>capsize:
                                raise ValueError('Too many elements')
            elif tocheck:
                (n,) = tocheck
                mask = adder(g,key,n)
                r = nodesearch0(g,g2, order, [n]+inlist, s, cds)
                if r and bfu.increment(r) == g2:
                    s.add(bfu.g2num(r))
                    if capsize and len(s)>capsize:
                        raise ValueError('Too many elements')
                remover(g,key,n,mask)

            order.insert(0,key)

        else:
            return g

    # find all directed g1's not conflicting with g2

    startTime = int(round(time.time() * 1000))
    gg = checkable(g2)

    idx = np.argsort([len(gg[x]) for x in gg.keys()])
    keys = [gg.keys()[i] for i in idx]

    cds, order, idx = conformanceDS(g2, gg, keys)
    endTime = int(round(time.time() * 1000))
    if verbose:
        print "precomputed in {:10} seconds".format(round((endTime-startTime)/1000.,3))
    if 0 in [len(x) for x in order]:
        return set()
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g, g2, [keys[i] for i in idx], ['0'], s, cds, order, set())
        #nodesearch0(g, g2, [gg.keys()[i] for i in idx], ['0'], s, cds)
    except ValueError, e:
        print e
        s.add(0)
    return s

def backtrack_more2(g2, rate=2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    f = [(addaVpath,delaVpath,maskaVpath)]
    c = [ok2addaVpath]

    def predictive_check(g,g2,pool,checks_ok, key):
        s = set()
        for u in pool:
            if not checks_ok(key,u,g,g2,rate=rate): continue
            s.add(u)
        return s

    @memo2 # memoize the search
    def nodesearch(g, g2, order, inlist, s, cds, pool, pc):
        if order:
            if bfu.undersample(g,rate) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                s.update(supergraphs_in_eq(g, g2, rate=rate))
                return g

            key = order[0]
            if pc:
                tocheck = [x for x in pc if x in cds[len(inlist)-1][inlist[0]]]
            else:
                tocheck = cds[len(inlist)-1][inlist[0]]

            if len(order) > 1:
                kk = order[1]
                pc = predictive_check(g,g2,pool[len(inlist)],
                                      c[edge_function_idx(kk)],kk)
            else:
                pc = set()

            adder, remover, masker = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]

            for n in tocheck:
                if not checks_ok(key,n,g,g2,rate=rate): continue
                masked = np.prod(masker(g,key,n))
                if masked:
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                else:
                    mask = adder(g,key,n)
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                    remover(g,key,n,mask)

        elif bfu.undersample(g,rate)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements')
            return g

    # find all directed g1's not conflicting with g2

    startTime = int(round(time.time() * 1000))
    ln = [x for x in itertools.permutations(g2.keys(),rate)] + \
         [(n,n) for n in g2]
    gg = {x:ln for x in gk.edgelist(g2)}
    keys = gg.keys()
    cds, order, idx = conformanceDS(g2, gg, gg.keys(), f=f, c=c)
    endTime = int(round(time.time() * 1000))
    print "precomputed in {:10} seconds".format(round((endTime-startTime)/1000.,3))
    if 0 in [len(x) for x in order]:
        return set()
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g, g2, [keys[i] for i in idx], ['0'], s, cds, order, set())
    except ValueError, e:
        print e
        s.add(0)
    return s


def unionpool(idx, cds):
    s = set()
    for u in cds[idx]:
        for v in cds[idx][u]:
            s = s.union(cds[idx][u][v])
    return s

def conformanceDS(g2, gg, order, f=[], c=[]):
    CDS = {}
    pool = {}

    CDS[0] = set(gg[order[0]])
    pool = [set(gg[order[i]]) for i in range(len(order))]

    for x in itertools.combinations(range(len(order)),2):

        d, s_i1, s_i2 = inorder_check2(order[x[0]], order[x[1]],
                                       pool[x[0]], pool[x[1]],
                                       g2, f=f, c=c)

        pool[x[0]] = pool[x[0]].intersection(s_i1)
        pool[x[1]] = pool[x[1]].intersection(s_i2)

        d = del_empty(d)
        if not x[1] in CDS:
            CDS[x[1]] = {}
            CDS[x[1]][x[0]] = d
        else:
            CDS[x[1]][x[0]] = d
    if density(g2) > 0.35:
        itr3 = [x for x in itertools.combinations(range(len(order)),3)]
        for x in random.sample(itr3, min(10,np.int(scipy.misc.comb(len(order),3)))):
            s1, s2, s3 = check3(order[x[0]], order[x[1]], order[x[2]],
                                pool[x[0]], pool[x[1]], pool[x[2]],
                                g2, f=f, c=c)

            pool[x[0]] = pool[x[0]].intersection(s1)
            pool[x[1]] = pool[x[1]].intersection(s2)
            pool[x[2]] = pool[x[2]].intersection(s3)

    return prune_sort_CDS(CDS, pool)

def prune_modify_CDS(cds, pool):
    ds = {}
    ds[0]={}
    ds[0]['0'] = pool[0]
    for i in range(1,len(pool)):
        ds[i] = {}
        for j in cds[i].keys():
            for e in pool[i-1].intersection(cds[i][j].keys()):
                ds[i][e] = pool[i].intersection(cds[i][j][e])
    return ds, pool, range(len(pool))

def prune_sort_CDS(cds, pool):

    idx = np.argsort([len(x) for x in pool])
    p = [pool[i] for i in idx]

    ds = {}
    ds[0]={}
    ds[0]['0'] = pool[idx[0]]

    for i in range(1,len(idx)):
        ds[i] = {}
        for j in range(i):
            if idx[j] > idx[i]:
                dd = invertCDSelement(cds[idx[j]][idx[i]])
            else:
                dd = cds[idx[i]][idx[j]]
            for e in pool[idx[i-1]].intersection(dd.keys()):
                ds[i][e] = pool[idx[i]].intersection(dd[e])

    return ds, p, idx

def invertCDSelement(d_i):
    d = {}
    for e in d_i:
        for v in d_i[e]:
            if v in d:
                d[v].add(e)
            else:
                d[v]=set([e])
    return d

def conformant(cds, inlist):

    if inlist[len(inlist)-2] in cds[len(inlist)-1][0]:
        s = cds[len(inlist)-1][0][inlist[len(inlist)-2]]
    else:
        return set()
    for i in range(1,len(inlist)-1):
        if inlist[len(inlist)-i-2] in cds[len(inlist)-1][i]:
            s = s.intersection(cds[len(inlist)-1][i][inlist[len(inlist)-i-2]])
        else:
            return set()
    return s

def length_d_paths(G, s, d):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    def recurse(G, s, d, path=[]):
        if d == 0:
            yield path
            return

        for u in G[s]:
            if G[s][u] == set([(2,0)]) or u in path: continue
            for v in recurse(G, u, d-1, path+[u]):
                yield v

    for u in recurse(G, s, d, [s]):
            yield u

def edge_increment_ok(s,m,e,g,g2):
    """
    s - start,
    m - middle,
    e - end
    """
    # bidirected edges
    for u in g[s]:
        if u!=m and not (m in g2[u] and (2,0) in g2[u][m]):return False

    # directed edges
    if s == e:
        if not (m in g2[m] and (0,1) in g2[m][m]):return False
        if not (s in g2[s] and (0,1) in g2[s][s]):return False
    for u in g[m]:
        if not (u in g2[s] and (0,1) in g2[s][u]):return False
        # bidirected edges
        if u!=e and not (e in g2[u] and (2,0) in g2[u][e]):return False
    for u in g[e]:
        if not (u in g2[m] and (0,1) in g2[m][u]):return False

    for u in g:
        if s in g[u] and not (m in g2[u] and (0,1) in g2[u][m]):
            return False
        if m in g[u] and not (e in g2[u] and (0,1) in g2[u][e]):
            return False

    return True


def length_d_loopy_paths(G, s, dt, p):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    g1 = cloneempty(G)

    def recurse(g, g2, s, d, path=[]):

        if edge_increment_ok(p[-d-2],s,p[-d-1],g,g2):

            if d == 0:
                yield path
                return

            mask = add2edges(g,(p[-d-2],p[-d-1]),s)
            for u in g2[s]:
                if g2[s][u] == set([(2,0)]): continue
                for v in recurse(g, g2, u, d-1, path+[u]):
                    yield v
            del2edges(g,(p[-d-2],p[-d-1]),s, mask)

    for u in recurse(g1, G, s, dt-1, [s]):
            yield u
# system packages
import StringIO
import scipy
import sys
import os
import igraph
import numpy as np

# local packages
import zickle
import ecj
from ecj import undersample
import testgraphs
reload(testgraphs)
from testgraphs import *

colors = zickle.load('data/colors.zkl') # load the colors

def graph2dict(g):
    D = {}
    for v in range(0,len(g.vs)):
        D[g.vs[v]["label"]] = {}
        for u in g.neighbors(v,mode="OUT"):
            D[g.vs[v]["label"]][g.vs[u]["label"]] = set([(0,1)])
    return D

def paintSCC(g, cm):
    nameidx = {}
    D = graph2dict(g)
    scc = ecj.scc(D)
    for i in range(0,len(scc)):
        for v in scc[i]:
            g.vs[g.vs['label'].index(v)]["color"] = cm[i]
fstr = "%.5f" # precision to use
def dict2graph(D):
    A = scipy.zeros([len(D),len(D)])
    nodes = []
    indices = {}
    c = 0
    for v in D:
        nodes.append(v)
    idx = np.argsort([int(v) for v in nodes])
    nodes = [nodes[i] for i in idx]
    for v in nodes:
        indices[v] = c
        c +=1
    for v in nodes:
        for u in D[v]:
            if  not ((2,0) in D[v][u]):
                A[indices[v],indices[u]] = 1
    g = igraph.Graph.Adjacency(A.tolist())
    g.vs["label"] = nodes
    return g

def unroll(G,steps):
    N = {}
    for i in range(0,steps):
        N.update({v+str(i):set([u+str(i+1) for u in G[v]]) for v in G})
    N.update({v+str(steps):set() for v in G})
    return N

def unroll_undersample(G,steps):
    # does not provide isochronal bidirectional edges
    N = {}
    steps += 2
    U = unroll(G,steps)
    nodes = G.keys()
    for v in G:
        N.update({v:set([nodes[k] for k in scipy.where([ecj.reachable(v+'0',U,u+str(steps-1)) for u in G])[0]])})
    return N

def matrix_start(mname='generic',w_gap=0.45,h_gap=0.5,stl=''):
    print "\matrix ("+mname+") [matrix of nodes, row sep="\
        +str(h_gap)+"cm,column sep="+str(w_gap)+"cm"+stl+"]"
    print "{"
def matrix_end():
    print "};"

def matrix_grid(G,s,mname='generic',w_gap=0.5, h_gap=0.5, type="obs",stl=''):
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    keylist = G.keys();
    idx = np.argsort([int(v) for v in keylist])
    keylist = [keylist[i] for i in idx]
    for v in keylist:
        print " ","& ".join(["\\node["+type+"]{"+v+"};" for i in range(1,s)])+"\\\\"
    matrix_end()

def matrix_edges(G,s,mname='generic'):
    nodes = G.keys()
    idx = np.argsort([int(v) for v in nodes])
    nodes = [nodes[i] for i in idx]    
    for v in nodes:
        idx1 = nodes.index(v)+1
        for u in G[v]:
            if G[v][u].intersection([(edge_type['directed'],1)]):
                idx2 = nodes.index(u)+1
                print '\\foreach \\x in{'+','.join(map(str,range(1,s-1)))+'}{'
                print '  \\pgfmathtruncatemacro{\\xn}{\\x+1}'
                print '  \\draw[pil] ('+mname+'-'+str(idx1)+'-\\x) -- ('\
                    +mname+'-'+str(idx2)+'-\\xn);'
                print '};'
def emacs_vars():
    print '%%% Local Variables:'
    print '%%% mode: latex'
    print '%%% TeX-master: "../master"'
    print '%%% End:'

def dbnprint(G,s,mname='generic',w_gap=0.5, h_gap=0.5,type="obs",stl=''):
    matrix_grid(G,s,mname=mname,w_gap=w_gap, h_gap=h_gap, type=type,stl=stl)
    matrix_edges(G,s,mname=mname)

from scipy import array,cos,sin,deg2rad,rad2deg
import math

def getangle(A,B):
    """
    When A and  B are two angles around the clock  returns an angle of
    the line that is connecting them.
    """
    x = array([cos(deg2rad(A)),sin(deg2rad(A))])
    y=array([cos(deg2rad(B)),sin(deg2rad(B))])
    d = y-x
    return rad2deg(math.atan2(d[1],d[0]))

def cdbnprint(G,mtype="obs",bend=5,curve=5,R=1):
    """
    Prints  out  a  compressed  dbn  repesentation  of  the  graph  in
    TikZ/Latex format
    """
    output = StringIO.StringIO()
    BE = set()
    n = len(G)
    nodes = G.keys()
    idx = np.argsort([int(v) for v in nodes])
    nodes = [nodes[i] for i in idx]

    g = dict2graph(ecj.cloneBfree(G))
    paintSCC(g,colors)


    for i in range(0,n):
        node = g.vs[i]['label']
        rc = g.vs[i]["color"]
        print >>output, "{ \\definecolor{mycolor}{RGB}{"\
            +str(rc[0])+","+str(rc[1])+","+str(rc[2])+"}"
        mcolor = "fill = {rgb: red,"+str(rc[0])+"; green,"+str(rc[1])+\
            "; blue,"+str(rc[2])+"}"
        print >>output, "\\node["+mtype+", fill=mycolor] ("+node+") at ("+\
            str(-i*360/n+180)+":"+str(R)+") {"+node+"};}"

#    print >>output,"\\foreach \\name/\\angle in {"+",".join(
#        [nodes[i]+"/"+str(-i*360/n+180) for i in range(0,n)])+"}"
#    print >>output,"\\node["+mtype+"] (\\name) at (\\angle:"\
#        +str(R)+") {\\name};"

    for i in range(0,n):
        v = nodes[i]
        ll=[v+'/'+u for u in G[v]]
        for l in ll:
            a,b = l.split('/')
            if G[a][b].intersection([(edge_type['bidirected'],0)]):
                if not(BE.intersection([(a,b)]) or BE.intersection([(b,a)])):
                    ang_a = -nodes.index(a)*360/n+180
                    ang_b = -nodes.index(b)*360/n+180
                    print >>output,'  \\draw[pilip, on layer=back] ('+a+') -- ('+b+');'
            if G[a][b].intersection([(edge_type['directed'],1)]):
                ang_a = -nodes.index(a)*360/n+180
                ang_b = -nodes.index(b)*360/n+180
                if a == b:
                    print >>output,"\\path[overlay,draw,pil] ("+a+")" +\
                        " .. controls +("+ "%.5f" % (bend+ang_a) +\
                        ":"+ fstr % (2*curve)+"mm) and +("+\
                        "%.5f" % (ang_a-bend)+\
                        ":"+"%.5f" % (2*curve)+"mm) .. ("+b+");"
                else:
                    print >>output,"\\path[overlay,draw,pil] ("+a+")" +\
                        " .. controls +("+"%.5f" % (bend+getangle(ang_a,ang_b)) +\
                        ":"+fstr % (curve)+"mm) and +("+\
                        fstr % (getangle(ang_b,ang_a)-bend)+\
                        ":"+fstr % (curve)+"mm) .. ("+b+");"
    return output

def gprint(G,mtype="obs",bend=5,curve=5,R=1,layout=None, scale=5):
    """
    Prints out an automatically layout compressed dbn repesentation of
    the graph in TikZ/Latex format
    """
    output = StringIO.StringIO()
    BE = set()
    n = len(G)
    if not layout:
        g = dict2graph(ecj.cloneBfree(G))
        layout = g.layout_fruchterman_reingold(maxiter=50000, coolexp=1.1)
        #layout = g.layout_graphopt(niter=50000, node_charge=0.08)
        layout.center([0,0])
        layout.scale(float(1/scipy.absolute(layout.coords).max()))
        layout.scale(R)
        cc = scipy.round_(array(layout.coords),decimals=4)
    else:
        g = dict2graph(ecj.cloneBfree(G))
        cc = array(layout.coords)
    paintSCC(g,colors)
    for i in range(0,n):
        node = g.vs[i]['label']
        rc = g.vs[i]["color"]
        print >>output, "{ \\definecolor{mycolor}{RGB}{"\
            +str(rc[0])+","+str(rc[1])+","+str(rc[2])+"}"
        mcolor = "fill = {rgb: red,"+str(rc[0])+"; green,"+str(rc[1])+\
            "; blue,"+str(rc[2])+"}"
        print >>output, "\\node["+mtype+", fill=mycolor] ("+node+") at ("+\
            str(cc[i][0])+","+str(cc[i][1])+") {"+node+"};}"

    for i in range(0,n):
        v = g.vs[i]['label']
        ll=[v+'/'+u for u in G[v]]
        for l in ll:
            a,b = l.split('/')
            if G[a][b].intersection([(edge_type['bidirected'],0)]):
                if not(BE.intersection([(a,b)]) or BE.intersection([(b,a)])):
                    print >>output,'  \\draw[pilip, on layer=back] ('+\
                        a+') -- ('+b+');'
            if G[a][b].intersection([(edge_type['directed'],1)]):
                if a == b:
                    dff = cc[g.vs['label'].index(a)] - scipy.mean(cc,0)
                    ang = scipy.arctan2(dff[1],dff[0])
                    ang_a = scipy.rad2deg(ang)
                    print >>output,"\\path[overlay,draw,pil] ("+a+")" +\
                        " .. controls +("+ "%.5f" % (bend+ang_a) +\
                        ":"+ fstr % (2*curve)+"mm) and +("+\
                        "%.5f" % (ang_a-bend)+\
                        ":"+"%.5f" % (2*curve)+"mm) .. ("+b+");"
                else:
                    dff = cc[g.vs['label'].index(b)] \
                        - cc[g.vs['label'].index(a)]
                    ang = scipy.arctan2(dff[1],dff[0])
                    ang_a = scipy.rad2deg(ang)
                    ang_b = ang_a+180
                    print >>output,"\\path[overlay,draw,pil] ("+a+")" +\
                        " .. controls +("+\
                        "%.5f" % (bend+ang_a) +\
                        ":"+fstr % (curve)+"mm) and +("+\
                        fstr % (ang_b-bend)+\
                        ":"+fstr % (curve)+"mm) .. ("+b+");"
    return output

def cdbnwrap(G,u,name='AAA',R=1,gap=0.5):
    output = StringIO.StringIO()
    print >>output,"\\node[right="+str(gap)+"cm of "+name+str(u-1)\
        +",scale=0.7]("+name+str(u)+"){"
    print >>output,"\\begin{tikzpicture}"
    s = cdbnprint(undersample(G,u),mtype='lahid',bend=25,curve=10,R=R)
    print >>output,s.getvalue()
    s.close()
    print >>output,"\\end{tikzpicture}"
    print >>output,"};"
    return output

def cdbnsingle(g,scale=0.7,R=1,gap=0.5,mtype="lahid"):
    output = StringIO.StringIO()
    print >>output,"\\node[scale="+str(scale)+"](){"
    print >>output,"\\begin{tikzpicture}"
    s = cdbnprint(g,mtype=mtype,bend=25,curve=10,R=R)
    print >>output,s.getvalue()
    s.close()
    print >>output,"\\end{tikzpicture}"
    print >>output,"};"
    return output

def cdbn_single(G,u,scale=0.7,R=1,gap=0.5,mtype="lahid"):
    return cdbnsingle(undersample(G,u),scale=scale,R=R,gap=gap,mtype=mtype)

def gsingle(g,scale=0.7,R=1,gap=0.5,mtype="lahid",layout=None):
    output = StringIO.StringIO()
    print >>output,"\\node[scale="+str(scale)+"](){"
    print >>output,"\\begin{tikzpicture}"
    s = gprint(g,mtype=mtype,bend=25,curve=6,R=R,layout=layout)
    print >>output,s.getvalue()
    s.close()
    print >>output,"\\end{tikzpicture}"
    print >>output,"};"
    return output

def g_single(G,u,scale=0.7,R=1,gap=0.5,mtype="lahid",layout=None):
    g = undersample(G,u)
    return gsingle(g,scale=scale,R=R,gap=gap,mtype=mtype,layout=layout)

def unfoldplot(G,steps=7,repeats=5,gap=0.5,R=1,hg=0.1,wgap=0.7,name='AAA',stl=''):
    u = 0
    dbnprint(undersample(G,u), repeats, w_gap=wgap,
             h_gap=hg, mname=name+str(u), type='hid', stl=stl)
    print "\\node[left="+str(gap)+"cm of "+name+str(u)+",scale=0.7] (C) {"
    print "\\begin{tikzpicture}"
    cdbnprint(G,mtype='hid',bend=15,curve=5,R=R)
    print "\\end{tikzpicture}"
    print "};"
    for u in range(1,steps):
        dbnprint(undersample(G,u),repeats,w_gap=wgap,h_gap=hg,mname=name+\
                     str(u),type='ahid',stl=', below=0.25cm of '+name+str(u-1))

        print "\\node[left="+str(gap)+"cm of "+name+str(u)+",scale=0.7] () {"
        print "\\begin{tikzpicture}"
        cdbnprint(undersample(G,u),mtype='lahid',bend=15,curve=5,R=R)
        print "\\end{tikzpicture}"
        print "};"
    emacs_vars()

def foldplot(G,steps=7,gap=0.5,R=1,hg=0.1,name='AAA',stl=''):
    u = 0
    print '\\node[scale=0.7'+stl+'] ('+name+str(u)+') {'
    print "\\begin{tikzpicture}"
    s=cdbnprint(G,mtype='hid',bend=25,curve=10,R=R)
    print s.getvalue()
    s.close()
    print "\\end{tikzpicture}"
    print "};"
    for u in range(1,steps):
        s=cdbnwrap(G,u,R=R,name=name)
        print s.getvalue()
        s.close()
    emacs_vars()

def matrix_fold(G,m,n,R=1,mname='g',w_gap=0.5, h_gap=0.5, stl=''):
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    for i in range(0,n):
        S = []
        for j in range(0,m):
            u = m*i+j
            if u == 0:
                s=cdbn_single(G,u,R=R,mtype="hid")
            else:
                s=cdbn_single(G,u,R=R)
            S.append(s.getvalue())
            s.close()
        print " ","& ".join(S)+"\\\\"
    matrix_end()

def gmatrix_fold(G,m,n,R=1,mname='g',w_gap=0.5, h_gap=0.5, stl='', shift=0):
    #g = dict2graph(ecj.cloneBfree(G))
    #layout = g.layout_graphopt(niter=50000, node_charge=0.02)
    #layout.center([0,0])
    #layout.scale(0.01)
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    for i in range(0,n):
        S = []
        for j in range(0,m):
            u = m*i+j+shift
            s=g_single(G,u,R=R,mtype="hid")
            S.append(s.getvalue())
            s.close()
        print " ","& ".join(S)+"\\\\"
    matrix_end()

def addselfedge(G,v):
    G[v].update({v:set([(edge_type['directed'],1)])})

def gmatrix_list(l,m,n,R=1,mname='g',w_gap=0.5,h_gap=0.5,stl='',shift=0):
    """
    Given a list of graphs prints them out as latex laying out each graph in a force-based layout
    """
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    for i in range(0,n):
        S = []
        for j in range(0,m):
            u = m*i+j+shift
            if u >len(l)-1: break
            s=gsingle(l[u],R=R,mtype="hid")
            S.append(s.getvalue())
            s.close()
        print " ","& ".join(S)+"\\\\"
    matrix_end()

def matrix_list(l,m,n,R=1,mname='g',w_gap=0.5,h_gap=0.5,stl='',shift=0):
    """
    Given a list of graphs prints them out as latex laying out each in a circle layout
    """
    matrix_start(mname=mname,w_gap=w_gap,h_gap=h_gap,stl=stl)
    for i in range(0,n):
        S = []
        for j in range(0,m):
            u = m*i+j+shift
            if u >len(l)-1: break
            s=cdbnsingle(l[u],R=R,mtype="hid")
            S.append(s.getvalue())
            s.close()
        print " ","& ".join(S)+"\\\\"
    matrix_end()

class WritableObject:
    def __init__(self):
        self.content = []
    def write(self, string):
        self.content.append(string)
N = {
'a': set('bcd'),
'b': set('ce'),
'c': set('d'),
'd': set('e'),
'e': set('f'),
'f': set('cgh'),
'g': set('fh'),
'h': set('fg')
}

N2 = {
'a': set('bcf'),
'b': set('c'),
'c': set('de'),
'd': set('e'),
'e': set(),
'f': set(),
}

N4 = {
'a': set('bdc'),
'b': set('d'),
'c': set('d'),
'd': set('efg'),
'e': set('g'),
'f': set('g'),
'g': set(),
}

"""
Each  edge is a  tuple (t,w),  where t  signifies its  type and  w its
integer time length:
0 - an isochronal edge
1 - from slice t to t+1 (the regular case)
2 - from t to t+2
etc
"""
edge_type = {'directed':0,'undirected':1,'bidirected':2}
iedge_type = {0:'directed',1:'undirected',2:'bidirected'}

LG = { # loop graph
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'e':set([(0,1)])},
'e': {'a':set([(0,1)])}
}

# 2 cycles
C = {
    'a': {'b':set([(0,1)]),'d':set([(0,1)])},
    'b': {'c':set([(0,1)])},
    'c': {'d':set([(0,1)])},
    'd': {'a':set([(0,1)])}
}

# a cycle of two
N = {
'a': {'b':set([(0,1)])},
'b': {'a':set([(0,1)]),'d':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'c':set([(0,1)])},
}

N1 = {
'a': {'d':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'b':set([(0,1)])},
'd': {'b':set([(0,1)]),'a':set([(0,1)])},
}

# a cycle of three
A = {
'a': {'b':set([(0,1)]),'e':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'a':set([(0,1)])},
'd': {'a':set([(0,1)])},
'e': {'d':set([(0,1)])},
}

# no cycles
U = {
'a': {'d':set([(0,1)]),'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {},
'd': {'c':set([(0,1)])},
}

# DAG
D = {
'a': {'b':set([(0,1)]),'c':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {}
}

D2 = {
'a': {},
'b': {},
'c': {'d':set([(0,1)]),'h':set([(0,1)]),'a':set([(0,1)]),'b':set([(0,1)])},
'd': {},
'e': {'d':set([(0,1)]),'f':set([(0,1)])},
'f': {'g':set([(0,1)]),'h':set([(0,1)])},
'g': {'c':set([(0,1)])},
'h': {},
}

# 3 SCC DAG
D3 = {
'a': {'b':set([(0,1)]),'c':set([(0,1)])},
'b': {'e':set([(0,1)]),'i':set([(0,1)]),'d':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'a':set([(0,1)])},
'e': {'f':set([(0,1)])},
'f': {'g':set([(0,1)])},
'g': {'e':set([(0,1)]),'h':set([(0,1)])},
'h': {'i':set([(0,1)])},
'i': {'h':set([(0,1)])},
}

P = {
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'a':set([(0,1)]),'b':set([(0,1)])},
}

L = {
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'e':set([(0,1)])},
'e': {'f':set([(0,1)]),'a':set([(0,1)])},
'f': {'a':set([(0,1)])},
'g': {'h':set([0,1])},
'h': {'i':set([0,1])},
'i': {'g':set([0,1])},
}

UL = {
'1': {'3':set([(0,1)]),'9':set([(0,1)])},
'3': {'4':set([(0,1)])},
'4': {'5':set([(0,1)])},
'5': {'1':set([(0,1)]),'6':set([(0,1)])},
'6': {'5':set([(0,1)]),'7':set([(0,1)])},
'7': {'8':set([(0,1)])},
'8': {'6':set([(0,1)])},
'9': {'a':set([(0,1)])},
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'1':set([(0,1)])},
}

UL11 = {
'1': {'2':set([(0,1)]),'3':set([(0,1)])},
'2': {'4':set([(0,1)])},
'3': {'4':set([(0,1)]),'3':set([(0,1)])},
'4': {},
}

# tree
T2 = {
'a': {'b':set([(0,1)]),'o':set([(0,1)])},
'b': {'f':set([(0,1)]),'c':set([(0,1)])},
'c': {'d':set([(0,1)]),'e':set([(0,1)])},
'd': {},
'e': {},
'f': {'g':set([(0,1)]),'h':set([(0,1)])},
'g': {},
'h': {},
'i': {},
'j': {},
'k': {'i':set([(0,1)]),'j':set([(0,1)])},
'l': {},
'm': {},
'n': {'m':set([(0,1)]),'l':set([(0,1)])},
'o': {'n':set([(0,1)]),'k':set([(0,1)])},
}

# one loop
L1 = {
'a': {'b':set([(0,1)])},
'b': {'c':set([(0,1)])},
'c': {'d':set([(0,1)])},
'd': {'e':set([(0,1)])},
'e': {'a':set([(0,1)])},
}

L25 = {
'1': {'2':set([(0,1)]),'1':set([(0,1)])},
'2': {'3':set([(0,1)])},
'3': {'5':set([(0,1)]),'6':set([(0,1)])},
#'4': {'5':set([(0,1)])},
'5': {'6':set([(0,1)])},
'6': {'7':set([(0,1)]),'2':set([(0,1)])},
'7': {'1':set([(0,1)])},
}
#'8': {'9':set([(0,1)])},
#'9': {'01':set([(0,1)])},b
#'01': {'1':set([(0,1)])},
# '02': {'03':set([(0,1)]),'6':set([(0,1)])},
# '03': {'1':set([(0,1)])},
# }

L01 = {
'1': {'2':set([(0,1)])},
'2': {'3':set([(0,1)]),'7':set([(0,1)])},
'3': {'4':set([(0,1)]),'7':set([(0,1)])},
'4': {'5':set([(0,1)]),'6':set([(0,1)])},
'5': {'6':set([(0,1)])},
'6': {'3':set([(0,1)])},
'7': {'1':set([(0,1)])},
}
L02 = {
'1': {'2':set([(0,1)]),'1':set([(0,1)])},
'2': {'3':set([(0,1)])},
'3': {'4':set([(0,1)]),'6':set([(0,1)])},
'4': {'5':set([(0,1)])},
'5': {'6':set([(0,1)]),'4':set([(0,1)])},
'6': {'7':set([(0,1)])},
'7': {'1':set([(0,1)]),'7':set([(0,1)])},
}

L03 = {
'1': {'2':set([(0,1)]),'1':set([(0,1)])},
'2': {'3':set([(0,1)])},
'3': {'4':set([(0,1)])},
'4': {'5':set([(0,1)])},
'5': {'6':set([(0,1)])},
'6': {'7':set([(0,1)])},
'7': {'1':set([(0,1)])},
}


TT = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)]),'01':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'01':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T10 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'01':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T22 = {
'01': {'02':set([(0,1)]),'10':set([(0,1)])},
'02': {'01':set([(0,1)])},
'03': {'04':set([(0,1)]),'10':set([(0,1)])},
'04': {'03':set([(0,1)])},
'05': {'06':set([(0,1)]),'10':set([(0,1)])},
'06': {'05':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'07':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T63 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'01':set([(0,1)]),'10':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'07':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T631 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'01':set([(0,1)]),'10':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'07':set([(0,1)]),'10':set([(0,1)]),'09':set([(0,1)])},
'10': {'10':set([(0,1)])},
}

T33 = {
'01': {'02':set([(0,1)]),'10':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'01':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'04':set([(0,1)]),'10':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'07':set([(0,1)]),'10':set([(0,1)])},
'10': {'10':set([(0,1)])},
}


OO = {
'01': {'02':set([(0,1)]),'11':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'10':set([(0,1)])},
'10': {'01':set([(0,1)])},
'11': {'11':set([(0,1)])},
'12': {'05':set([(0,1)]),'06':set([(0,1)]),'12':set([(0,1)])},
}

OO5 = {
'01': {'03':set([(0,1)]),'11':set([(0,1)])},
'02': {'04':set([(0,1)])},
'03': {'05':set([(0,1)])},
'04': {'06':set([(0,1)])},
'05': {'07':set([(0,1)])},
'06': {'08':set([(0,1)])},
'07': {'09':set([(0,1)])},
'08': {'10':set([(0,1)])},
'09': {'01':set([(0,1)])},
'10': {'02':set([(0,1)]),'11':set([(0,1)])},
'11': {'11':set([(0,1)])},
'12': {'05':set([(0,1)]),'06':set([(0,1)]),'12':set([(0,1)])},
}

MV = {
'1': {'4':set([(0,1)])},
'2': {},
'3': {},
'4': {'6':set([(0,1)])},
'5': {'4':set([(0,1)]),'2':set([(0,1)])},
'6': {'3':set([(0,1)]), '5':set([(0,1)])},
}


MV_ = {
'1': {'2':set([(0,1)]),'3':set([(0,1)])},
'2': {'3':set([(2,0)])},
'3': {'2':set([(2,0)])},
}

G96 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'01':set([(0,1)]),'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)]),'08':set([(0,1)])},
'07': {'04':set([(0,1)])},
'08': {'09':set([(0,1)]),'24':set([(0,1)])},
'09': {'10':set([(0,1)])},
'10': {'11':set([(0,1)])},
'11': {'12':set([(0,1)])},
'12': {'13':set([(0,1)])},
'13': {'14':set([(0,1)])},
'14': {'06':set([(0,1)])},

# '15': {'16':set([(0,1)]),'15':set([(0,1)])},
# '16': {'17':set([(0,1)])},
# '17': {'18':set([(0,1)])},
# '18': {'19':set([(0,1)])},
# '19': {'20':set([(0,1)])},
# '20': {'21':set([(0,1)])},
# '21': {'21':set([(0,1)])},
'22': {'01':set([(0,1)]),'22':set([(0,1)])},



'24': {'24':set([(0,1)])},
#'25': {'26':set([(0,1)])},
#'26': {'26':set([(0,1)])},
}

G960 = {
'01': {'02':set([(0,1)])},
'02': {'03':set([(0,1)])},
'03': {'04':set([(0,1)])},
'04': {'05':set([(0,1)])},
'05': {'06':set([(0,1)])},
'06': {'07':set([(0,1)]),'01':set([(0,1)])},
'07': {'08':set([(0,1)])},
'08': {'09':set([(0,1)])},
'09': {'01':set([(0,1)])},
}

DDD1 = {
    '01': {'02':set([(0,1)]),'01':set([(0,1)])},
    '02': {'03':set([(0,1)])},
    '03': {'04':set([(0,1)]), '06':set([(0,1)])},
    '04': {'05':set([(0,1)])},
    '05': {'06':set([(0,1)])},
    '06': {'01':set([(0,1)])}
}

DDD2 = {
    '01': {'02':set([(0,1)]),'01':set([(0,1)])},
    '02': {'03':set([(0,1)])},
    '03': {'06':set([(0,1)])},
    '04': {'03':set([(0,1)])},
    '05': {'04':set([(0,1)])},
    '06': {'01':set([(0,1)]), '05':set([(0,1)])}
}
"""Generic object pickler and compressor

This module saves and reloads compressed representations of generic Python
objects to and from the disk.
"""

__author__ = "Bill McNeill <billmcn@speakeasy.net>"
__version__ = "1.0"

import cPickle
import gzip

def save(object, filename, protocol = -1):
    """Save an object to a compressed disk file.
       Works well with huge objects.
    """
    with gzip.GzipFile(filename, 'wb') as file:
        cPickle.dump(object, file, protocol)

def load(filename):
    """Loads a compressed object from disk
    """
    with gzip.GzipFile(filename, 'rb') as file:
        object = cPickle.load(file)
    return object

if __name__ == "__main__":
	import sys
	import os.path
	
	class Object:
		x = 7
		y = "This is an object."
	
	filename = sys.argv[1]
	if os.path.isfile(filename):
		o = load(filename)
		print "Loaded %s" % o
	else:
		o = Object()
		save(o, filename)
		print "Saved %s" % o
import sys

sys.path.append('./tools/')
import graphkit as gk
import numpy as np
from ortools.constraint_solver import pywrapcp
import ipdb


def grelabel(g):
    """
    Relabel graph nodes to numerical contiguous labels 1 through len(g)
    :param g: input graph
    :return: relabeled graph and name mapping
    """
    gg = {}
    nodemap = {x[0]: x[1] + 1 for x in zip(g.keys(), range(len(g)))}
    for v in g:
        gg[nodemap[v]] = {}
        for w in g[v]:
            gg[nodemap[v]][nodemap[w]] = g[v][w]
    return gg, {v: k for k, v in nodemap.items()}


def crelabel(c, map):
    newcliq = set()
    for e in c:
        newcliq.add((map[e[0] + 1], map[e[1] + 1]))
    return newcliq


def printedges(g):
    l = gk.edgelist(g)
    for e in l:
        print e[0] - 1, '->', e[1] - 1


def edge_exists(i, j, g):
    return j + 1 in g[i + 1]


def numkids(i, edges, solver):
    N = len(edges)
    return solver.Sum([edges[i][k] for k in range(N)])


def numparents(i, edges, solver):
    N = len(edges)
    return solver.Sum([edges[k][i] for k in range(N)])


def clique_constrain(solver, parents, children, edges, g):
    N = len(parents)

    # declare constraints
    solver.Add(solver.Sum(parents) > 0)
    solver.Add(solver.Sum(children) > 0)

    for i in range(N):
        # makes sure that there exists at least one outgoing edge from node i if it is marked in parents for this b-clique
        solver.Add((parents[i] == 1) == (solver.Sum([edges[i][k] for k in range(N)]) >= 1))
        # this makes sure that there exists at least one incoming edge to node i if it is marked as a child
        solver.Add((children[i] == 1) == (solver.Sum([edges[k][i] for k in range(N)]) >= 1))

        solver.Add((children[i] == 1) == (solver.Sum(parents) <= numparents(i, edges, solver)))
        solver.Add((parents[i] == 1) == (solver.Sum(children) <= numkids(i, edges, solver)))

        for j in range(N):
            # edge existence constraints
            if not edge_exists(i, j, g):
                solver.Add(edges[i][j] == 0)
            else:
                solver.Add((edges[i][j] == 1) == (parents[i] * children[j] == 1))


def bcliques(g, verbose=False):
    solver = pywrapcp.Solver("b-clique")
    g, mp = grelabel(g)
    # declare variables
    edges = []
    N = len(g)
    for i in range(N):
        e = []
        for j in range(N):
            e.append(solver.IntVar(0, 1, "%i -> %i" % (i, j)))
        edges.append(e)

    parents = [solver.IntVar(0, 1, "%i" % (i)) for i in range(N)]
    children = [solver.IntVar(0, 1, "%i" % (i)) for i in range(N)]

    # declare constraints
    clique_constrain(solver, parents, children, edges, g)

    # run the solver
    solution = solver.Assignment()
    solution.Add([edges[i][j] for i in range(N) for j in range(N)])
    solution.Add(parents)
    solution.Add(children)

    collector = solver.AllSolutionCollector(solution)
    solver.Solve(solver.Phase([edges[i][j] for i in range(N) for j in range(N)] + children + parents,
                              solver.CHOOSE_FIRST_UNBOUND,
                              solver.ASSIGN_MAX_VALUE),
                 [collector])
    num_solutions = collector.SolutionCount()

    # output solutions
    if verbose:
        print "num_solutions:", num_solutions
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()

    cliques = []
    pts = set()
    for s in range(num_solutions):
        c = set()
        qval = [collector.Value(s, edges[i][j]) for i in range(N) for j in range(N)]
        pval = [collector.Value(s, parents[i]) for i in range(N)]
        cval = [collector.Value(s, children[i]) for i in range(N)]
        if verbose:
            print " ----------------- ", s
        if verbose:
            print np.asarray([pval, cval])
        for i in range(len(qval)):
            if qval[i]:
                e = np.unravel_index(i, [N, N])
                c.add((e[0], e[1]))
                if verbose:
                    print e[0], "->", e[1]
        if not (True in map(lambda x: c.issubset(x), cliques)):
            if not c.issubset(pts):
                pts = pts.union(c)
                cliques.append(c)
        if verbose:
            print "---"

    # check for b-cliques for which all edges are covered in other b-cliques
    cc = []
    for i in range(len(cliques)):
        bcl = cliques.pop()
        pts = set()
        for j in range(len(cliques)):
            pts = pts.union(cliques[j])
        for j in range(len(cc)):
            pts = pts.union(cc[j])

        if not bcl.issubset(pts):
            cc.append(bcl)
    return map(lambda x: crelabel(x, mp), cc)
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import ecj, bfutils
import bfutils as bfu
import graphkit as gk
import warnings
from statsmodels.tsa.api import VAR
from sympy.matrices import SparseMatrix
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed
from scipy import linalg, optimize
import numpy as np
import scipy
import ipdb

def symchol(M): # symbolic Cholesky
    B = SparseMatrix(M)
    t = B.row_structure_symbolic_cholesky()
    B = np.asarray(B)*0
    for i in range(B.shape[0]): B[i,t[i]] = 1
    return B

def G2SVAR(G):
    n = len(G)
    A,B = npG2SVAR(G)
    P,L,U = linalg.lu(B)
    A = linalg.inv(L).tolist()
    B = B.tolist()
    A = listplace(A, 0.0, 0.0)
    for i in range(0,n): A[i][i] = 1
    B = listplace(B, 0.0, 'e')
    for i in range(0,n): B[i][i] = 'e'
    return A,B,P

def G2AH(G):
    n = len(G)
    A,B = npG2SVAR(G)
    P,L,U = linalg.lu(B)
    A = linalg.inv(L).tolist()
    B = B.tolist()
    A = listplace(A, 0.0, 0.0)
    for i in range(0,n): A[i][i] = 1
    B = listplace(B, 0.0, 'e')
    for i in range(0,n): B[i][i] = 'e'
    return A,B,P

def bnf2CG(fname):
    d = eval(open(fname).read())
    G = {}
    for v in d:
        G[v] = {u: set((0,1)) for u in d[v]['pars']}
    G = ecj.tr(G)
    for v in G:
        ld = {u: set([(0,1)]) for u in G[v]}
        G[v] = ld
    return G

def npG2SVAR(G):
    n = len(G)
    A = bfu.graph2adj(G)
    B = np.tril(bfu.graph2badj(G))
    np.fill_diagonal(B,1)
    B = symchol(B)
    return A,B

def x2M(x, A, B, aidx, bidx):
    A[aidx] = x[:len(aidx[0])]
    B[bidx] = x[len(aidx[0]):]
    #B[(bidx[1],bidx[0])] = x[len(aidx[0]):]
    return A, B

def nllf(x, A, B, Y, aidx, bidx): # negative log likelihood
    A,B = x2M(x, A, B, aidx, bidx)
    T = Y.shape[1]
    X = Y[:,1:] - np.dot(A, Y[:,:-1])
    ldB = T*np.log(abs(1./linalg.det(B)))
    return ldB + 0.5*np.trace( np.dot(np.dot(B.T, B), np.dot(X,X.T)))

def nllf2(x, A, B, YY, XX, YX, T, aidx, bidx): # negative log likelihood
    A,B = x2M(x, A, B, aidx, bidx)
    AYX = np.dot(A, YX.T)
    S = YY - AYX - AYX.T + np.dot(np.dot(A,XX), A.T)
    ldB = T*np.log(abs(1./linalg.det(B)))
    return 0.5*np.dot(np.dot(B.T, B).T.flat,S.flat) + ldB
    #return ldB + 0.5*np.trace( np.dot(np.dot(B.T, B), S))

def VARbic(nllf, K, T):
    return 2*nllf + K*np.log(T)

def listplace(l, a, b):
    return [listplace(x,a,b) if not np.isscalar(x) else b if x != a  else x for x in l]

# -------------------------------------------------------------------
# data generation
# -------------------------------------------------------------------

def randweights(n, c=0.1, factor=9):
    rw = scipy.random.randn(n)
    idx = scipy.where(abs(rw) < factor*c)
    if idx:
        rw[idx] = rw[idx]+scipy.sign(rw[idx])*c*factor
    return rw

def transitionMatrix(cg, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = scipy.where(A==1)
    A[edges] = randweights(edges[0].shape[0], c=minstrength)
    l = linalg.eig(A)[0]
    c = 0
    pbar = ProgressBar(widgets=['Searching for weights: ', Percentage(), ' '], maxval=10000).start()
    while max(l*scipy.conj(l)) > 1:
        A[edges] = randweights(edges[0].shape[0], c=c)
        c += 1
        l = linalg.eig(A)[0]
        pbar.update(c)
    pbar.finish()
    return A

def sampleWeights(n, minstrength=0.1):
    r = scipy.randn(n)
    s = minstrength/np.min(np.abs(r))
    r = s*r
    return r

def transitionMatrix2(cg, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = scipy.where(A==1)
    A[edges] = sampleWeights(edges[0].shape[0], minstrength=minstrength)
    l = linalg.eig(A)[0]
    c = 0
    pbar = ProgressBar(widgets=['Searching for weights: ', Percentage(), ' '], maxval=10000).start()
    while max(l*scipy.conj(l)) > 1:
        A[edges] = sampleWeights(edges[0].shape[0], minstrength=minstrength)
        c += 1
        l = linalg.eig(A)[0]
        if c>pbar.maxval:
            raise ValueError
        pbar.update(c)
    pbar.finish()
    return A

def transitionMatrix3(cg, x0=None, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = scipy.where(A==1)

    try:
        s = x0.shape
        x = x0
    except AttributeError:
        A = initRandomMatrix(A, edges)
        x = A[edges]

    def objective(x):
        A[edges] = np.real(x)
        l = linalg.eig(A)[0]
        m = np.max(np.real(l*scipy.conj(l)))-0.99
        n = np.min(np.min(np.abs(x)),minstrength)-minstrength
        return m*m + 0.1*n*n

    o = np.zeros(len(edges))
    while np.min(np.abs(o[0])) < 0.8*minstrength:
        rpt = True
        while rpt:
            try:
                try:
                    o = optimize.fmin_bfgs(objective, x,
                                           gtol=1e-10, maxiter=100,
                                           disp=False, full_output=True)
                    A[edges]=np.real(o[0])
                    l = linalg.eig(A)[0]
                    if np.max(np.real(l*scipy.conj(l))) < 1:
                        rpt = False

                except:
                    rpt = True
            except Warning:
                x = scipy.randn(len(edges[0]))
                rpt = True
    A[edges]=np.real(o[0])
    return A

def initRandomMatrix(A, edges, maxtries=100, distribution='beta', stable=True):
    '''
    possible distributions:
    flat
    flatsigned
    beta
    normal
    uniform
    '''
    s = 2.0

    def init():
        if distribution=='flat':
            x = np.ones(len(edges[0]))
        elif distribution=='flatsigned':
            x = np.sign(scipy.randn(len(edges[0])))*scipy.ones(len(edges[0]))
        elif distribution=='beta':
            x = np.random.beta(0.5,0.5,len(edges[0]))*3-1.5
        elif distribution=='normal':
            x = scipy.randn(len(edges[0]))
        elif distribution=='uniform':
            x = np.sign(scipy.randn(len(edges[0])))*scipy.rand(len(edges[0]))
        else:
             raise ValueError('Wrong option!')
        return x

    def eigenvalue(A):
        l = linalg.eig(A)[0]
        s = np.max(np.real(l*scipy.conj(l)))
        return s

    x = init()
    A[edges] = x
    s = eigenvalue(A)
    alpha = np.random.rand()*(0.99-0.8)+0.8
    A = A/(alpha*s)
    s = eigenvalue(A)

    return A

def transitionMatrix4(g, minstrength=0.1, distribution='normal', maxtries=1000):
    A = gk.CG2adj(g)
    edges = np.where(A==1)
    s = 2.0
    c = 0
    pbar = ProgressBar(widgets=['Searching for weights: ',
                                Percentage(), ' '],
                       maxval=maxtries).start()
    while s > 1.0:
        minstrength -= 0.001
        A = initRandomMatrix(A, edges, distribution=distribution)
        x = A[edges]
        delta = minstrength/np.min(np.abs(x))
        A[edges] = delta*x
        l = linalg.eig(A)[0]
        s = np.max(np.real(l*scipy.conj(l)))
        c += 1
        if c > maxtries:
            return None
        pbar.update(c)
    pbar.finish()

    return A

def drawsamplesLG(A, nstd=0.1, samples=100):
    n = A.shape[0]
    data = scipy.zeros([n, samples])
    data[:,0] = nstd*scipy.random.randn(A.shape[0])
    for i in range(1,samples):
        data[:,i] = scipy.dot(A,data[:,i-1]) \
                    + nstd*scipy.random.randn(A.shape[0])
    return data


def getAgraph(n, mp=2, st=0.5, verbose=True):
    keeptrying = True
    while keeptrying:
        G = gk.rnd_cg(n, maxindegree=mp, force_connected=True)
        try:
            A = transitionMarix2(G, minstrength=st)
            keeptrying = False
        except ValueError as e:
            if verbose:
                print "!!! Unable to find strong links for a stable matrix !!!"
                print "*** trying a different graph"
    return {'graph':      G,
            'transition': A,
            'converges':  len(bfutils.call_undersamples(G))}

def getAring(n, density=0.1, st=0.5, verbose=True, dist='flatsigned', permute=False):
    keeptrying = True
    plusedges = bfutils.dens2edgenum(density,n)
    while keeptrying:
        G = gk.ringmore(n, plusedges, permute=permute)
        try:
            A = transitionMatrix4(G, minstrength=st, distribution=dist)
            try:
                s = A.shape
                keeptrying = False
            except AttributeError:
                keeptrying = True
        except ValueError:
            if verbose:
                print "!!! Unable to find strong links for a stable matrix !!!"
                print "*** trying a different graph"
    return {'graph':      G,
            'transition': A,
            'converges':  len(bfutils.call_undersamples(G))}



# -------------------------------------------------------------------
# estimation
# -------------------------------------------------------------------

def scoreAGraph(G, data, x0 = None):
    A,B = npG2SVAR(G)
    n = len(G)
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    K = scipy.sum(len(a_idx[0])+len(b_idx[0])/2)

    if x0:
        o = optimize.fmin_bfgs(nllf, x0, args=(A, B, data, a_idx, b_idx),
                               disp=False, full_output=True)
    else:
        o = optimize.fmin_bfgs(nllf, scipy.randn(len(a_idx[0])+len(b_idx[0])),
                               args=(np.double(A), np.double(B),
                                     data, a_idx, b_idx),
                               disp=False, full_output=True)
    ipdb.set_trace()
    return 2*o[1] + K*np.log(data.shape[1]) #VARbic(o[1],K,data.shape[1])

def scoreAGraph2(G, data, x0 = None):
    A,B = npG2SVAR(G)
    n = len(G)
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)

    T = data.shape[1]
    YY = np.dot(data[:,1:],data[:,1:].T)
    XX = np.dot(data[:,:-1],data[:,:-1].T)
    YX = np.dot(data[:,1:],data[:,:-1].T)

    K = scipy.sum(len(a_idx[0])+len(b_idx[0])/2)

    if x0:
        o = optimize.fmin_bfgs(nllf2, x0,
                               args=(np.double(A), np.double(B),
                                     YY,XX,YX,T,a_idx, b_idx),
                                     gtol=1e-12, maxiter=500,
                                     disp=False, full_output=True)
    else:
        o = optimize.fmin_bfgs(nllf2, scipy.randn(len(a_idx[0])+len(b_idx[0])),
                               args=(np.double(A), np.double(B),
                                     YY,XX,YX,T,a_idx, b_idx),
                                gtol=1e-12, maxiter=500,
                                disp=False, full_output=True)

    return 2*o[1] + K*np.log(data.shape[1]) #VARbic(o[1],K,data.shape[1])

def estimateG(G,YY,XX,YX,T,x0=None):
    A,B = npG2SVAR(G)
    K = scipy.sum(abs(A)+abs(B))
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    try:
        s = x0.shape
        x = x0
    except AttributeError:
        x = scipy.randn(K)
    o = optimize.fmin_bfgs(nllf2, x,
                           args=(np.double(A), np.double(B),
                                 YY,XX,YX,T,a_idx, b_idx),
                           disp=False, full_output=True)
    A,B = x2M(o[0], np.double(A), np.double(B), a_idx, b_idx)
    return  A,B


def data2AB(data,x0=None):
    n = data.shape[0]
    T = data.shape[1]
    YY = np.dot(data[:,1:],data[:,1:].T)
    XX = np.dot(data[:,:-1],data[:,:-1].T)
    YX = np.dot(data[:,1:],data[:,:-1].T)

    model = VAR(data.T)
    r = model.fit(1)
    A = r.coefs[0,:,:]

    #A = np.ones((n,n))
    B = np.ones((n,n))
    np.fill_diagonal(B,0)
    B[np.triu_indices(n)] = 0
    K = np.int(scipy.sum(abs(B)))#abs(A)+abs(B)))

    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    np.fill_diagonal(B,1)

    try:
        s = x0.shape
        x = x0
    except AttributeError:
        x = np.r_[A.flatten(),0.1*scipy.randn(K)]
    o = optimize.fmin_bfgs(nllf2, x,
                           args=(np.double(A), np.double(B),
                                 YY,XX,YX,T,a_idx, b_idx),
                           gtol=1e-12, maxiter=500,
                           disp=False, full_output=True)
    ipdb.set_trace()
    A,B = x2M(o[0], np.double(A), np.double(B), a_idx, b_idx)
    B = B+B.T
    return  A,B

def amap(f, a):
     v = np.vectorize(f)
     return v(a)

def AB2intAB(A,B, th=0.09):
    A[amap(lambda x: abs(x) > th, A)] = 1
    A[amap(lambda x: abs(x) < 1, A)] = 0
    B[amap(lambda x: abs(x) > th, B)] = 1
    B[amap(lambda x: np.abs(x) < 1, B)] = 0
    return A,B

def intAB2graph(A,B):
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if A[j,i]: g[str(i+1)][str(j+1)] = set([(0,1)])

    for i in range(n):
        for j in range(n):
            if B[j,i] and j!=i:
                if str(j+1) in g[str(i+1)]:
                    g[str(i+1)][str(j+1)].add((2,0))
                else:
                    g[str(i+1)][str(j+1)] = set([(2,0)])
    return g

def data2graph(data,x0=None, th=0.09):
    A,B = data2AB(data,x0=x0)
    Ab,Bb = AB2intAB(A,B,th=th)
    return intAB2graph(Ab,Bb)

def data2VARgraph(data, pval=0.05):
    model = VAR(data.T)
    r = model.fit(1)
    A = r.coefs[0,:,:]
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if np.abs(A[j,i]) > pval: g[str(i+1)][str(j+1)] = set([(0,1)])
    return g

def data2VARgraph_model(data, pval=0.05):
    model = VAR(data.T)
    r = model.fit(1)
    A = r.coefs[0,:,:]
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if np.abs(A[j,i]) > pval: g[str(i+1)][str(j+1)] = set([(0,1)])
    return g, r
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))
import traversal as trv
import bfutils as bfu
import graphkit as gk
import zickle as zkl
import numpy as np
import itertools as iter
import statsmodels.api as sm
import linear_model as lm
import ipdb
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore, norm

def independent(y,X,pval=0.05):
    X  = sm.add_constant(X)
    est  = sm.OLS(y,X).fit()
    return est.pvalues[1] > pval

def kernel(z):
    if np.abs(z) > 1.: return 0.
    return .5

def residuals_(x,y,z):
    N = len(x)
    pd = map(kernel,pdist(z.T))
    PD = squareform(pd)
    sumsx = np.dot(PD,x) + 0.5*x
    sumsy = np.dot(PD,y) + 0.5*y
    weights = np.sum(PD, axis=1) + 0.5
    residualsx = x - sumsx/weights
    residualsy = y - sumsy/weights

    return residualsx, residualsy

def moment22(x,y): return np.dot(x*x,y*y)/len(x)

def fdr(alpha, pvalues):
    m = len(pvalues)
    c = np.cumsum(1./(np.arange(m)+1))
    pcompare = (pvalues <= alpha*(np.arange(1,m+1)/(c*(m+1))))
    idx = np.where(pcompare == True)[0]
    if len(idx) == 0: return -1
    return idx[-1]

def fdrQ(alpha,pvalues):
    pvalues = np.sort(pvalues)
    pvalues = pvalues[~np.isnan(pvalues)]
    min  = np.nan if len(pvalues) == 0 else pvalues[0]
    high = 1.0
    low  = 0.
    q    = alpha
    while (high - low) > 1e-5:
        midpoint = (high + low)/2.0
        q = midpoint
        cutoff = pvalues[fdr(q, pvalues)]

        if cutoff < min:
            low = midpoint
        elif cutoff > min:
            high = midpoint
        else:
            low  = midpoint
            high = midpoint
    return q

def fdrCutoff(alpha, pvalues):
    pvalues = np.sort(pvalues)
    k = fdr(alpha,pvalues)
    if k < 0: return 0
    return pvalues[k]

def np_fisherZ(x,y,r):
    z = 0.5 * (np.log(1.0 + r) - np.log(1.0 - r))
    w = np.sqrt(len(x)) * z
    x_ = zscore(x)
    y_ = zscore(y)
    t2 = moment22(x_,y_)
    t = np.sqrt(t2)
    p = 2. * (1. - norm.cdf(np.abs(w), 0.0, t))
    return p

def independent_(x,y, alpha = 0.05):
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]

    # For PC, should not remove the edge for this reason.
    if len(x) < 10: return False

    ps = []
    for i in range(15):
        for j in range(15):
            x_ = x**i
            y_ = y**j
            r = np.corrcoef(x_,y_)[0,1]
            #r = max(min(r,1),-1) # Tetrad had this
            p = np_fisherZ(x_,y_,r)
            #if not np.isnan(p):
            ps.append(p)

    if not ps: return True
    return fdrCutoff(alpha, ps) > alpha


def addallb(g):
    n = len(g)
    for i in range(n):
        for j in range(n):
            if str(j+1) in g[str(i+1)]:
                g[str(i+1)][str(j+1)].add((2,0))
            else:
                g[str(i+1)][str(j+1)] = set([(2,0)])
    return g

def dpc(data, pval=0.05):
    n = data.shape[0]
    # stack the data: first n rows is t-1 slice, the next n are slice t
    data = np.asarray(np.r_[data[:,:-1],data[:,1:]])

    def tetrad_cind_(y,x,condset=[], alpha=0.01, shift=0):
        y = data[n+int(y)-1,:]
        x = data[shift+int(x)-1,:]
        if condset:
            X  = data[condset,:]
            ry, rx = residuals_(y,x,X)
        else:
            ry, rx = [y,x]
        return independent_(ry, rx, alpha = alpha)

    def cind_(y,x, condset=[], pval=pval, shift=0):
        yd = data[n+int(y)-1,:].T
        X  = data[[shift+int(x)-1]+condset,:].T
        return independent(yd, X, pval=pval)

    def cindependent(y, x, counter, parents=[], pval=pval):
        for S in [j for j in iter.combinations(parents,counter)]:
            if cind_(y, x, condset=list(S), pval=pval): return True
            #if tetrad_cind_(x, y, condset=list(S), alpha=pval): return True                
        return False

    def bindependent(y, x, parents=[], pval=pval):
        return cind_(y, x, condset=parents, pval=pval, shift=n)
        #return tetrad_cind_(y, x, condset=parents, alpha=pval, shift=n)

    def prune(elist, mask, g):
        for e in mask:
            g[e[0]][e[1]].remove((0,1))
            elist.remove(e)
        gk.clean_leaf_nodes(g)

    g  = gk.superclique(n)
    gtr= bfu.gtranspose(g)

    el = gk.edgelist(g)
    for counter in range(n):
        to_remove = []
        for e in el:
            ppp = [int(k)-1 for k in gtr[e[1]] if k != e[0]]
            if counter <= len(ppp):
                if cindependent(e[1], e[0], counter, parents=ppp, pval=pval):
                    to_remove.append(e)
                    gtr[e[1]].pop(e[0],None)
        prune(el, to_remove, g)

    bel = [map(lambda k: str(k+1), x) for x in iter.combinations(range(n),2)]
    for e in bel:
        ppp = list(set(gtr[e[0]].keys()) | set(gtr[e[1]].keys()))
        ppp = map(lambda x: int(x)-1, ppp)
        if bindependent(e[0], e[1], parents=ppp, pval=pval):
            g[e[0]][e[1]].remove((2,0))
            g[e[1]][e[0]].remove((2,0))
    gk.clean_leaf_nodes(g)

    return g


# Local Variables:
# mode: python
# python-indent-offset: 4
# End:
import networkx
from bfutils import g2num, num2CG
from comparison import graph2nx
from ecj import undersample
from numpy import argsort

def simple_loops(g, u):
    """
    iterator over the list of simple loops of graph g at the undersample rate u
    """
    gx = graph2nx(num2CG(g2num(undersample(g,u)), len(g)))
    for l in networkx.simple_cycles(gx):
        yield l

def print_loops(g, u):
    l = [x for x in simple_loops(g,u)]
    lens = map(len, l)    
    for i in argsort(lens):
        print l[i]

def ul(l):
    """
    returns list elements that are present only once
    """
    u, r = set(), set()
    for e in l:
        if e not in u:
            u.add(e)
        else:
            r.add(e)
    return u.difference(r)
import sys

sys.path.append('./tools/')

import bclique as bq
import latents as lt
import pathtreetools as ptt
from pathtree import PathTree


def can_add_loop(pt, num, elements):
    r = False
    if ptt.isptelement(PathTree(num - pt.preset, pre=pt.preset), num):
        r = True
    for e in pt.loopset:
        if type(e) is int:
            can_add_loop(PathTree(pt.preset))

def learn_path_tree(pt):
    elements = ptt.pt2seq(pt, 1)
    newpt = PathTree(set(), pre={elements[0]})

    def rpath(elements, npt):
        if not ptt.isptelement(npt, element[0]):


    return newpt

def bpts(bc):
    """
    Given a b-clique returns a set of PathTrees (one for each edge in bclique) such that loops are identified across edges in bclique
    :param bc: b-clique
    :return:
    """

def getagenerator(g):
    bclqs = bq.bcliques(g)
    for clq in bclqs:
        bpts(clq)
    return bclqs
import ecj
import scipy
import numpy
import operator
import networkx as nx
#from progressbar import ProgressBar, Percentage
numpy.random.RandomState()
import bfutils as bfu
import numpy as np
import gmpy as gmp

def num2CG(num,n):
    """num2CG - converts a number  whose binary representaion encodes edge
    presence/absence into a compressed graph representaion

    """
    n2 = n*n
    G = {'%i'%(i+1):{} for i in xrange(n)}
    if num == 0: return G
    bl = gmp.bit_length(num)
    idx = [n2-i-1 for i in xrange(bl) if num & (1<<i)]
    idx = np.unravel_index(idx,(n,n))
    x = idx[0]+1
    y = idx[1]+1
    for i in xrange(len(x)):
        G['%i' % x[i]]['%i' % y[i]] = set([(0,1)])
    return G

def hasSelfLoops(G):
    for u in G:
        if G[u].has_key(u):
            return True
    return False

def randSCC(n):
    G = num2CG(scipy.random.randint(2**(n**2)),n)
    while (len(ecj.scc(G)) > 1) or gcd4scc(G)>1:
        G = num2CG(scipy.random.randint(2**(n**2)),n)
    return G

def SM_fixed(Gstar,G, iter=5):
    compat = []
    for j in range(0,iter):
        if Gstar == ecj.undersample(G,j):
            compat.append(j)
    return compat
def SM_converging(Gstar,G):
    """Gstar is the undersampled reference graph, while G is the starting
    graph. The  code searches  over all undersampled  version of  G to
    find all matches with Gstar
    """
    compat = []
    GG = G
    Gprev = G
    if G == Gstar: return [0]
    j = 1
    G = ecj.undersample(GG,j)
    while not (G == Gprev):
        if Gstar == G: compat.append(j)
        j += 1
        Gprev = G
        G = ecj.undersample(GG,j)
    return compat
def searchMatch(Gstar,G, iter=5):
    if gcd4scc(G) >1: return SM_fixed(Gstar, G, iter=iter)
    return SM_converging(Gstar, G)
def hasSink(G):
    return not reduce(operator.and_, [bool(G[n]) for n in G], True)
def hasRoot(G): return hasSink(ecj.tr(G))
def isSclique(G):
    n = len(G)
    for v in G:
        if sum([(0,1) in G[v][w] for w in G[v]]) < n: return False
        if sum([(2,0) in G[v][w] for w in G[v]]) < n-1: return False
    return True

def graph2nx(G):
    g = nx.DiGraph()
    for v in G:
        g.add_edges_from([(v,x) for x in G[v] if (0,1) in G[v][x]])
    return g
def nx2graph(G):
    g = {str(n+1):{} for n in G}
    for n in G:
        g['%i' % (n+1)] = {'%i' % (x+1):set([(0,1)]) for x in G[n]}
    return g

def gcd4scc(SCC):
    g = graph2nx(SCC)
    return ecj.listgcd(map(lambda x: len(x)-1, nx.simple_cycles(g)))

def compatibleAtU(uGstar):
    compat = []
    n = len(uGstar)
    numG = 2**(n**2)
    #pbar = Percentage()
    for i in range(1,numG):
        G = num2CG(i,n)
        #pbar.update(i+1)
        if len(ecj.scc(G)) > 1: continue
        l = searchMatch(uGstar,G, iter = 5)
        if l: compat.append((l,G))
    #pbar.finish()
    return compat
from copy import deepcopy
all_acceptable = []

def find_next_graph(graph_g):
    next_graph = {}
    for i in graph_g:next_graph[i] = {}
    for i in graph_g:
        for j in graph_g[i]:
            for m in graph_g[j]:
                if not m in next_graph[i]:next_graph[i][m] = set([(0,1)])
                elif next_graph[i][m] == set([(0,2)]) or next_graph[i][m] == set([(0,3)]):next_graph[i][m] = set([(0,3)])
            for n in graph_g[i]:
                if j != n:
                    if not n in next_graph[j]:next_graph[j][n] = set([(0,2)])
                    elif next_graph[j][n] == set([(0,1)]) or next_graph[j][n] == set([(0,3)]):next_graph[j][n] = set([(0,3)])
    return next_graph

def compare(g_two_star,graph_g):
    for i in g_two_star:
        for j in g_two_star[i]:
            if j not in graph_g[i]:
                return False
            elif g_two_star[i][j] != graph_g[i][j] and  graph_g[i][j] != set([(0,3)]):
                return False
    return True

def try_next_double_edge(g_de, g_one_star_de, stack_de):
    global all_acceptable
    if stack_de != []:
        j = stack_de.pop()
        i = stack_de.pop()
        for k in g_one_star_de:
            boolean_first_edge_has_value = False
            boolean_second_edge_has_value = False
            if i in g_one_star_de[k]:
                boolean_first_edge_has_value = True
            if j in g_one_star_de[k]:
                boolean_second_edge_has_value = True
            g_one_star_de[k][i] = 1
            g_one_star_de[k][j] = 1
            g_two_star_de = find_next_graph(g_one_star_de)
            if compare(g_two_star_de, g_de):
                try_next_double_edge(g_de, g_one_star_de, stack_de)
            if not boolean_first_edge_has_value:
                del g_one_star_de[k][i]
            if not boolean_second_edge_has_value:
                del g_one_star_de[k][j]
        stack_de.append(i)
        stack_de.append(j)
    else:
        in_all_accept = False
        for xxxx in all_acceptable:
            if(xxxx == g_one_star_de):
                in_all_accept = True
        if (not in_all_accept):
            all_acceptable.append(deepcopy(g_one_star_de))
            print "Added 1."
#print g_one_star_de
#if g_one_star_de not in all_acceptable:
#all_acceptable.append(g_one_star_de)
        #print "G2 --------------------------"
        #print find_next_graph(g_one_star_de) == g_de
        #print "End-------------------------"


def sample_g_double_edge(g_de, g_one_star_de):
    stack_de = []
    for i in g_de:
        for j in g_de[i]:
            if g_de[i][j] != set([(0, 1)]):
                if int(i) < int(j):
                    stack_de.append(i)
                    stack_de.append(j)
    try_next_double_edge(g_de, g_one_star_de, stack_de)

def try_next_single_edge(g_se, g_one_star_se, stack_se):
    if stack_se != []:
        j = stack_se.pop()
        i = stack_se.pop()
        for k in g_one_star_se:
            boolean_first_edge_has_value = False
            boolean_second_edge_has_value = False
            if k in g_one_star_se[i]:
                boolean_first_edge_has_value = True
            if j in g_one_star_se[k]:
                boolean_second_edge_has_value = True
            g_one_star_se[i][k] = 1
            g_one_star_se[k][j] = 1
            g_two_star_se = find_next_graph(g_one_star_se)
            if compare(g_two_star_se, g_se):
                try_next_single_edge(g_se, g_one_star_se, stack_se)
            if not boolean_first_edge_has_value:
                del g_one_star_se[i][k]
            if not boolean_second_edge_has_value and j in g_one_star_se[k]:
                del g_one_star_se[k][j]
        stack_se.append(i)
        stack_se.append(j)
    else:
        sample_g_double_edge(g_se,g_one_star_se)

def sample_g_single_edge(g_se, g_one_star_se):
    stack_se = []
    for i in g_se:
        for j in g_se[i]:
            if g_se[i][j] != set([(0, 2)]):
                stack_se.append(i)
                stack_se.append(j)
    try_next_single_edge(g_se, g_one_star_se, stack_se)



def main(gStart):
    global all_acceptable
    all_acceptable = []
    g_one_star = {}
    for i in gStart:
        g_one_star[i] = {};
    sample_g_single_edge(gStart, g_one_star)
    return all_acceptableimport scipy
import itertools
# from progressbar import ProgressBar, Percentage
from multiprocessing import Pool, Array, Process, Manager
from numpy.random import randint
import numpy as np
# import ipdb
import networkx as nx
# local
import ecj
import zickle as zkl
import graphkit as gk
from comparison import num2CG, nx2graph, isSclique
import itertools
#import rpy2


def pure_directed_inc(G, D):
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in D[v]:
            if G[w]:
                for e in G[w]: G_un[v][e] = set([(0, 1)])
    return G_un


def directed_inc(G, D):
    G_un = {}
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in D[v]:
            if G[w] and (0, 1) in D[v][w]:
                for e in G[w]: G_un[v][e] = set([(0, 1)])
    return G_un


def bidirected_inc(G, D):
    # bidirected edges
    for w in G:
        # transfer old bidirected edges
        for l in D[w]:
            if (2, 0) in D[w][l]:
                G[w].setdefault(l, set()).add((2, 0))
        # new bidirected edges
        l = [e for e in D[w] if (0, 1) in D[w][e]]
        for pair in itertools.permutations(l, 2):
            G[pair[0]].setdefault(pair[1], set()).add((2, 0))
    return G


def increment(g):
    """
    Undersample the input graph g by 2
    Only works for g1 to g2 directed

    Args:
        g: an input graph in the dictionary format

    Returns:
        r: a graph in the dictionary format
    """
    r = {n: {} for n in g}

    for n in g:
        for h in g[n]:
            for e in g[h]:
                if not e in r[n]:
                    r[n][e] = set([(0, 1)])

    for n in g:
        for pair in itertools.combinations(g[n], 2):

            if pair[1] in r[pair[0]]:
                r[pair[0]][pair[1]].add((2, 0))
            else:
                r[pair[0]][pair[1]] = set([(2, 0)])

            if pair[0] in r[pair[1]]:
                r[pair[1]][pair[0]].add((2, 0))
            else:
                r[pair[1]][pair[0]] = set([(2, 0)])

    return r


def dincrement_u(G_star, G_u):
    # directed edges
    G_un = pure_directed_inc(G_star, G_u)
    return G_un


def increment_u(G_star, G_u):
    # directed edges
    G_un = directed_inc(G_star, G_u)
    # bidirected edges
    G_un = bidirected_inc(G_un, G_u)
    return G_un


def undersample(G, u):
    Gu = G
    for i in range(u):
        Gu = increment_u(G, Gu)
    return Gu


def all_undersamples(G_star):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if ecj.isSclique(g): return glist  # superclique convergence
        # this will (may be) capture DAGs and oscillations
        if g in glist: return glist
        glist.append(g)
    return glist


def graph2adj(G):
    n = len(G)
    A = scipy.zeros((n, n), dtype=np.int8)
    for v in G:
        A[int(v) - 1, [int(w) - 1 for w in G[v] if (0, 1) in G[v][w]]] = 1
    return A


def graph2badj(G):
    n = len(G)
    A = scipy.zeros((n, n), dtype=np.int8)
    for v in G:
        A[int(v) - 1, [int(w) - 1 for w in G[v] if (2, 0) in G[v][w]]] = 1
    return A


def adj2graph(A):
    G = {str(i): {} for i in range(1, A.shape[0] + 1)}
    idx = np.where(A == 1)
    for i in range(len(idx[0])):
        G['%i' % (idx[0][i] + 1)]['%i' % (idx[1][i] + 1)] = set([(0, 1)])
    return G


def adjs2graph(A, B):
    names = [str(i) for i in range(1, A.shape[0] + 1)]
    G = {}
    for name in names:
        G[name] = {}
    for i in range(A.shape[0]):
        for name in map(str, np.where(A[i, :] == 1)[0] + 1):
            G[str(i + 1)][name] = set([(0, 1)])

    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            if B[i, j]:
                if str(j + 1) in G[str(i + 1)]:
                    G[str(i + 1)][str(j + 1)].add((2, 0))
                else:
                    G[str(i + 1)][str(j + 1)] = set([(2, 0)])
    return G


def g2vec(g):
    A = graph2adj(g)
    B = graph2badj(g)
    return np.r_[A.flatten(), B[np.triu_indices(B.shape[0])]]


def vec2adj(v, n):
    A = np.zeros((n, n))
    B = np.zeros((n, n))
    A[:] = v[:n ** 2].reshape(n, n)
    B[np.triu_indices(n)] = v[n ** 2:]
    B = B + B.T
    return A, B


def vec2g(v, n):
    A, B = vec2adj(v, n)
    return adjs2graph(A, B)


# tried mutable ctypes buffer - not faster :(
def graph2str(G):
    n = len(G)
    d = {((0, 1),): '1', ((2, 0),): '0', ((2, 0), (0, 1),): '1', ((0, 1), (2, 0),): '1'}
    A = ['0'] * (n * n)
    for v in G:
        for w in G[v]:
            A[n * (int(v, 10) - 1) + int(w, 10) - 1] = d[tuple(G[v][w])]
    return ''.join(A)


def graph2bstr(G):
    n = len(G)
    d = {((0, 1),): '0', ((2, 0),): '1', ((2, 0), (0, 1),): '1', ((0, 1), (2, 0),): '1'}
    A = ['0'] * (n * n)
    for v in G:
        for w in G[v]:
            A[n * (int(v, 10) - 1) + int(w, 10) - 1] = d[tuple(G[v][w])]
    return ''.join(A)


def adj2num(A):
    s = reduce(lambda y, x: y + str(x),
               A.flatten().tolist(), '')
    return int(s, 2)


# def g2num(G): return int(graph2str(G),2) #adj2num(graph2adj(G))
def g2num(g):
    n = len(g)
    n2 = n ** 2 + n
    num = 0
    for v in range(1, n + 1):
        for w in g[str(v)]:
            num |= (1 << (n2 - v * n - int(w, 10)))
    return num


def ug2num(g):
    n = len(g)
    n2 = n ** 2 + n
    num = 0
    mask = 0
    num2 = 0
    for v in g:
        for w in g[v]:
            if (0, 1) in g[v][w]:
                mask = (1 << (n2 - int(v, 10) * n - int(w, 10)))
                num |= mask
            if (2, 0) in g[v][w]: num2 |= mask
    return num, num2


def bg2num(g):
    n = len(g)
    n2 = n ** 2 + n
    num = 0
    for v in g:
        for w in g[v]:
            if (2, 0) in g[v][w]:
                num = num | (1 << (n2 - int(v) * n - int(w)))
    return num


# def bg2num(G): return int(graph2bstr(G),2)#adj2num(graph2badj(G))
# def ug2num(G): return (g2num(G),bg2num(G))#(adj2num(graph2adj(G)),adj2num(graph2badj(G)))

def num2adj(num, n):
    l = list(bin(num)[2:])
    l = ['0' for i in range(0, n ** 2 - len(l))] + l
    return scipy.reshape(map(int, l), [n, n])


def add_bd_by_adj(G, adj):
    c = 0
    for e in adj:
        for v in range(len(e)):
            if e[v] == 1:
                try:
                    G[str(c + 1)][str(v + 1)].add((2, 0))
                except KeyError:
                    G[str(c + 1)][str(v + 1)] = set([(2, 0)])
        c += 1
    return G


def tuple2graph(t, n):
    g = num2CG(t[0], n)
    return add_bd_by_adj(g, num2adj(t[1], n))


def call_undersamples(G_star):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if g in glist: return glist
        glist.append(g)
    return glist


def overshoot(G_star, H):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if isSclique(g): return False
        if gk.isedgesubset(H, g): return True
        if g in glist: return False
        glist.append(g)
    return False


def forms_loop(G_star, loop):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if (g2num(gk.digonly(g)) & loop) == loop:
            return True
        if g in glist: return False
        glist.append(g)
    return False


def call_u_conflicts_d(G_star, H, checkrate=0):
    glist = [G_star]
    while True:
        g = dincrement_u(G_star, glist[-1])
        if gk.isedgesubset(g, H): return False
        if g in glist: return True
        glist.append(g)
    return True


def call_u_conflicts(G_star, H, checkrate=0):
    glist = [G_star]
    while True:
        # g = increment_u(G_star, glist[-1])
        g = directed_inc(G_star, glist[-1])
        if gk.isdedgesubset(g, H): return False
        g = bidirected_inc(g, glist[-1])
        if gk.isedgesubset(g, H): return False
        if g in glist: return True
        glist.append(g)
    return True


def call_u_conflicts2(G_star, H):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if gk.isedgesubset(g, H): return False, glist
        if g in glist: return True, glist
        glist.append(g)
    return True, glist


def call_u_equals2(G_star, glist, H):
    while True:
        g = increment_u(G_star, glist[-1])
        if g == H: return True
        if g in glist: return False
        glist.append(g)
    return False


def call_u_equals(G_star, H):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if g == H: return True
        if g in glist: return False
        glist.append(g)
    return False


def compact_call_undersamples(G_star, steps=None):
    glist = [ug2num(G_star)]
    lastgraph = G_star
    while True:
        g = increment_u(G_star, lastgraph)
        if ug2num(g) in glist: return glist
        glist.append(ug2num(g))
        lastgraph = g
    return glist


def cc_undersamples(G_star, steps=1):
    glist = [ug2num(G_star)]
    lastgraph = G_star
    for i in xrange(steps):
        g = increment_u(G_star, lastgraph)
        n = ug2num(g)
        if n in glist: return []
        glist.append(n)
        lastgraph = g
    return glist[-1]


def compatible(d1, d2):
    idx = scipy.where(scipy.array([[r == l for l in d2] for r in d1]))
    return idx


def compat(G):
    n = len(G)
    num = g2num(G)
    # sample all the graph for gStar
    star_l = all_undersamples(G)
    hits = {}
    # brute force all graphs
    for i in range(0, 2 ** (n ** 2)):
        tmpG = num2CG(i, n)
        tmp_l = all_undersamples(tmpG)
        c = compatible(tmp_l, star_l)
        if len(sum(c)) > 0: hits[i] = c
    return hits


def icompat(i, nodes):
    print i
    g = num2CG(i, nodes)
    return compat(g)


def ilength(i, nodes):
    print i
    g = num2CG(i, nodes)
    return len(call_undersamples(g))


def iall(i, nodes):
    print i
    g = num2CG(i, nodes)
    return compact_call_undersamples(g)


def cc_all(i, nodes, steps):
    # print i
    g = num2CG(i, nodes)
    return i, cc_undersamples(g, steps=steps)


def make_rect(l):
    max_seq = max(map(len, l))
    nl = []
    for e in l:
        e += [e[-1]] * (max_seq - len(e))
        nl.append(e)
    return nl


def uniqseq(l):
    s = []
    ltr = map(lambda *a: list(a), *l)
    for i in range(len(ltr)):
        s.append(len(np.unique(ltr[i])))


def loadgraphs(fname):
    g = zkl.load(fname)
    return g


def savegraphs(l, fname):
    zkl.save(l, fname)


def jason2graph(g):
    r = {}
    d = {1: set([(0, 1)]),
         2: set([(2, 0)]),
         3: set([(0, 1), (2, 0)])}
    for head in g:
        r[head] = {}
        for tail in g[head]:
            r[head][tail] = d[g[head][tail]]
    return r


def graph2jason(g):
    r = {}
    for head in g:
        r[head] = {}
        for tail in g[head]:
            if g[head][tail] == set([(0, 1)]):
                r[head][tail] = 1
            elif g[head][tail] == set([(2, 0)]):
                r[head][tail] = 2
            elif g[head][tail] == set([(0, 1), (2, 0)]):
                r[head][tail] = 3
    return r


def ring(n):
    g = {}
    for i in range(1, n):
        g[str(i)] = {str(i + 1): set([(0, 1)])}
    g[str(n)] = {'1': set([(0, 1)])}
    return g


def addAring(g):
    for i in range(1, len(g)):
        if str(i + 1) in g[str(i)]:
            g[str(i)][str(i + 1)].add((0, 1))
        else:
            g[str(i)][str(i + 1)] = set([(0, 1)])
    if '1' in g[str(len(g))]:
        g[str(len(g))]['1'].add((0, 1))
    else:
        g[str(len(g))]['1'] = set([(0, 1)])


def upairs(n, k):
    '''
    n unique nonsequential pairs
    '''
    s = set()
    for p in randint(n, size=(3 * k, 2)):
        if p[1] - p[0] == 1: continue
        s.add(tuple(p))
    l = [e for e in s]
    return l[:k]


def ringarcs(g, n):
    for edge in upairs(len(g), n):
        g[str(edge[0] + 1)][str(edge[1] + 1)] = set([(0, 1)])
    return g


def ringmore(n, m):
    return ringarcs(ring(n), m)


# talking about extra edges on top of the ring
def dens2edgenum(d, n=10): return int(d * n ** 2) - n


def edgenum2dens(e, n=10): return np.double(e + n) / n ** 2


def gtranspose(G):  # Transpose (rev. edges of) G
    GT = {u: {} for u in G}
    for u in G:
        for v in G[u]:
            GT[v][u] = set([(0, 1)])  # Add all reverse edges
    return GT


def scale_free(n, alpha=0.7, beta=0.25,
               delta_in=0.2, delta_out=0.2):
    g = nx.scale_free_graph(n, alpha=alpha,
                            beta=beta,
                            delta_in=delta_in, delta_out=delta_out)
    g = nx2graph(g)
    g = gtranspose(g)
    addAring(g)
    return g


def randH(n, d1, d2):
    g = bfu.ringmore(n, d1)
    pairs = [x for x in itertools.combinations(g.keys(), 2)]
    for p in np.random.permutation(pairs)[:d2]:
        g[p[0]].setdefault(p[1], set()).add((2, 0))
        g[p[1]].setdefault(p[0], set()).add((2, 0))
    return g
