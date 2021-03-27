import json
import numpy as np
import time
import multiprocessing as mp
import collections
import sys
import os

os.system("taskset -p 0xff %d" % os.getpid())

threads = int(sys.argv[2])
splits = 60

filename = sys.argv[1]

clustlist = json.load(open(filename+'_clusters.json','r'))
clustsizes = json.load(open(filename+'_clustersizes.json','r'))

clusters = {}

def findhead_nochanges(r):
    stop = 0
    rcurrent = r
    while stop == 0:
        if isinstance(clustlist[rcurrent], int):
            clust = clustlist[rcurrent]
            stop = 1
        else:
            rcurrent = clustlist[rcurrent]
    return (rcurrent,clust)

nreads = len(clustlist.keys())

starttime = time.time()
count = 0
p_nreads = nreads/100
percent = 1
for key in clustlist.keys():
    rcurrent, clust = findhead_nochanges(key)
    clustlist[key] = clust
    count += 1
    if count > percent*p_nreads:
#        print(str(percent)+' percent done in '+str(time.time()-starttime)+' seconds.')
        percent += 1

json.dump(clustlist,open(filename+'_clusters_unchained.json','w'))

#clustlist = json.load(open('CAMI1_clusters_final445_unchained.json','r'))

#clusters_large = [k for k,v in clustersizes.items() if v>=20]
starttime = time.time()
clustlist_large = dict((k, v) for k, v in clustlist.items() if clustsizes[str(v)] >= 4)

#print('There are '+str(len(clustlist_large))+' clusters with at least 20 reads, out of a total of '+str(len(clustlist))+' clusters.')
#print(time.time()-starttime)

clustlist = clustlist_large
clustlist_large = None

starttime = time.time()
dictsize = int(len(clustlist.items())/threads)
dcl = collections.defaultdict(dict)
for i in range(threads):
    start_idx = i*dictsize
    end_idx = min((i+1)*dictsize,len(clustlist.items()))
    dcl[i] = dict(list(clustlist.items())[start_idx:end_idx])
#print('Splitting this huge dict into '+str(splits)+' takes '+str(time.time()-starttime)+' seconds.')

#sys.exit()

'''
testdict = {'a': 1, 'b':2, 'c':3, 'd':4, 'e':5, 'f':1,  'g':1, 'h':2, 'i':2, 'j':3, 'k':3, 'l':4, 'm':4, 'n':5, 'o':5, 'p':1, 'q':1, 'r':1, 's':2, 't':2, 'u':2, 'v':3, 'w':3, 'x':4, 'y':4, 'z':5}
d = collections.defaultdict(dict)
d[0] = dict(list(testdict.items())[int(len(testdict)/2):])
d[1] = dict(list(testdict.items())[:int(len(testdict)/2)])
'''

def clustering(d):
    clusters_sub = {}
    for val in np.unique(list(d.values())):
        clusters_sub[str(val)] = [k for k,v in d.items() if v==val]
    return clusters_sub

if __name__ == "__main__":
#    print('Start the clustering procedure')
    starttime = time.time()
    maplist = []
#    print(type(dcl))
#    for i in range(threads):
#        maplist.extend(dcl[i])
#    print(maplist)
    po = mp.Pool(threads)
    result = po.map(clustering,[dcl[0],dcl[1],dcl[2],dcl[3],dcl[4],dcl[5],dcl[6],dcl[7],dcl[8],dcl[9],dcl[10],dcl[11],dcl[12],dcl[13],dcl[14],dcl[15],dcl[16],dcl[17],dcl[18],dcl[19],dcl[20],dcl[21],dcl[22],dcl[23],dcl[24],dcl[25],dcl[26],dcl[27],dcl[28],dcl[29],dcl[30],dcl[31],dcl[32],dcl[33],dcl[34],dcl[35],dcl[36],dcl[37],dcl[38],dcl[39],dcl[40],dcl[41],dcl[42],dcl[43],dcl[44],dcl[45],dcl[46],dcl[47],dcl[48],dcl[49],dcl[50],dcl[51],dcl[52],dcl[53],dcl[54],dcl[55],dcl[56],dcl[57],dcl[58],dcl[59]])
    po.close()
    po.join()
    d = result[0]
    for i in range(1,len(result)):
        for k,v in result[i].items():
            if k in d.keys():
                d[k].extend(v)
            else:
                d[k] = v
#    print('Clustering took '+str(time.time()-starttime)+' seconds on '+str(threads)+' threads.')
    json.dump(d,open(filename+'_clusters_grouped.json','w'))

'''
# Do a sanity check by checking whether you have the same distribution of cluster sizes!!!!


starttime = time.time()
count = 0
nclusts = np.unique(list(clustlist.values()))
p_nclusts = nclusts/100
percent = 1
for val in np.unique(list(clustlist.values())):
    clusters[val] = [k for k,v in clustlist.items() if v==val]
    count += 1
    if count > percent*p_nreads:
#        print(str(percent)+' percent of clusters done in '+str(time.time()-starttime)+' seconds.')
        percent += 1


json.dump(clusters,open('clusters_final_grouped.json','w'))
'''
