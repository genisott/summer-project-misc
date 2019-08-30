import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time

import pysambar as sm
import condor as cd



def cut_histogram(net,typec):
    #Finds minimum on the histogram and prints threshold
    weights = list(net[typec])
    weights = [i for i in weights if i >0]
    u,b = np.histogram(weights,bins=50000)
    step  = (b[-1]-b[0])/len(b)
    st = u[int((1-b[0])/step):int((3.5-b[0])/step)]
    thr = b[list(u).index(min(st))]
    print(thr)
    #Reformats the node names so the network is truly bipartite and cuts it by the threshold.
    net["reg"] = "r"+net["reg"]
    cut = net[["tar","reg",typec]]
    cut = cut[cut[typec]>thr] #Cut-off
    return cut

def to_sign(filename):
    condor_output = pd.read_csv(filename,sep=",",index_col=0)
    num_com = max(condor_output["com"])
    print(num_com+1)

    base=dict()
    for i in range(0,num_com+1):
            base[i]=[]
    for _,r in condor_output.iterrows():
                base[r["com"]].append(r["tar"])
    ls = list()
    for i in range(0,num_com+1):
            fullStr = "Com"+str(i)+"\t"+'\t'.join(base[i])
            ls.append(fullStr)
    f = open("signPF.txt","w")
    for line in ls:
        f.write(line)
        f.write("\n")
    f.close()

def main(net_filename,mut_filename):
    t = time.time()
    print("Start")
    net = pd.read_csv(net_filename,index_col=0)
    print("Network read")
    net = cut_histogram(net,net.columns[2])
    print("Cut-off done")
    co = cd.condor_object(net)
    co = cd.initial_community(co)
    co = cd.brim(co)
    co["tar_memb"].to_csv("tar_memb.txt")
    print("Condor done")
    to_sign("tar_memb.txt")
    print("Sign file done")
    pt,cl = sm.sambar(mut_file=mut_filename,gmtfile="signPF.txt",gmtMSigDB=False,kmax=25)
    print("Sambar done")
    print("Everything done in ",time.time()-t)
    return pt,cl