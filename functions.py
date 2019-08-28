import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pysambar as sm



def silhouette(Dm,clusters):
    """Silhouette scores for each sample for a given clustering. """
    sil = list()
    for cluster in clusters:
        #Compute s(i) for every element of the cluster.
        for element in cluster:
            if len(cluster)==1: #If the cluster has only one element si=0
                sil.append([element,0])
                continue
            #a(i): 
            ai = (Dm[cluster].loc[element]).sum()/(len(cluster)-1)
            #b(i):
            im = np.array([Dm[element].loc[other].sum()/(len(other)) for other in clusters if other is not cluster])
            bi = min(im)
            si = (bi-ai)/max(ai,bi)
            sil.append([element,si])
    return sil

def silhouette_score(Dm,clusters):
	"""Silhouette width score for a clustering. 
	Average of individual silhouette scores."""
	sl = silhouette(Dm,clusters)
	return sum([u[1] for u in sl])/len(sl)

def silhouette_plot(Dm,cl):
    cl_l = cl.sort_values(by=[cl.index[0]],axis=1).columns
    nclust = cl.loc[cl.index[0]].max()+1
    gr = {j:[r for r in cl.columns if cl.loc[cl.index[0],r]==j-1] for j in range(1,nclust+1)}
    sl = silhouette(Dm,[gr[j] for j in range(1,nclust+1)])
    
    df=pd.DataFrame(sl)
    df = df.set_index(0)
    df["s"] = df.index
    df = df.loc[cl_l]
    
    df["g"] = cl.iloc[0]
    df.columns = ["sil","sam","gr"]
    
    df.sort_values(["gr","sil"],ascending=[False,True],inplace=True)
    ax = plt.figure(figsize=(20,20))
    sns.barplot(x="sil",y="sam",hue="gr",data=df,dodge=False)
    plt.title("Silhouette plot, k = "+cl.index[0]+", score ="+str(sum([u[1] for u in sl])/len(sl)))
    plt.yticks([])
    
    plt.ylabel("Samples")
    plt.xlabel("Silhouette score")
    
    
    return ax

def cut_histogram(net,typec):
	"""Function to cut network histogram in the minimum.
	It doesn't work as well as one would like but it was a temporary solution."""
    
    weights = list(net[typec])
    weights = [i for i in weights if i >0]
    u,b = np.histogram(weights,bins=50000)
    step  = (b[-1]-b[0])/len(b)
    st = u[int((1-b[0])/step):int((3.5-b[0])/step)]
    thr = b[list(u).index(min(st))]
    cut = net[["tar","reg",typec]]
    cut = cut[cut[typec]>thr] #Cut-off
    return cut


def to_sign(filename):
	"""Transform output file of condor to a sign_file accepted by pysambar."""

	condor_output = pd.read_csv(filename,sep=",",index_col=0)
    num_com = max(condor_output["com"])
    print(typec, num_com+1)

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