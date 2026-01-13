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
