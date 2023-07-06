from matrix import *
import re
import matplotlib as ptl
import pandas as pd
import math



def read_GT(path: str) -> dict:

    '''
    Read ground truth from file

    Parameters
    ----------
    path: str
        Path to folder containing ground truth files

    Returns
    -------
    gt: dict
        Ground truth dict for all networks in the folder. The key is the network name, 
            the value is a list of ground truth communities
    '''

    gt = {}

    files = os.listdir(path)
    
    for file in files:

        with open(path + '/' + file, 'r') as f:
            lines = f.readlines()
            networtk_gt = []
            for line in lines:
                communities = [int(node) for node in line.split(' ')]               
                networtk_gt.append(communities)
            gt[file] = networtk_gt
    
    return gt

def read_rc_output(path: str) -> dict:

    '''
    Read communities from file

    Parameters
    ----------
    path: str
        Path to folder containing community files

    Returns
    -------
    communities: dict
        Communities dict for all networks in the folder. The key is the network name, 
            the value is a list of communities
    '''

    communities = {}

    files = os.listdir(path)
    
    for file in files:

        with open(path + '/' + file, 'r') as f:
            lines = f.readlines()
            clean_lines = []
            for line in lines:
                clean_lines.append(line.strip('\n').rstrip())
            lines = clean_lines
            networtk_communities = []
            for line in lines:
                communitie = [int(node) for node in line.split(' ')]               
                networtk_communities.append(communitie)
            communities[file] = networtk_communities
    
    return communities

def overlaping_detection(dict: dict) -> dict:
    
        '''
        Detect overlaping communities
    
        Parameters
        ----------
        dict: dict
            Communities dict for all networks in the folder. The key is the network name, 
                the value is a list of communities
    
        Returns
        -------
        overlaping: dict
            Overlaping communities dict for all networks. The key is the network name, 
                the value is a set of overlaping nodes.
        '''
    
        overlaping = {}
    
        for key, value in dict.items():
            nodes = set()
            overlaping[key] = set()
            for community in value:
                intersection = nodes.intersection(set(community))
                if len(intersection) > 0:
                    overlaping[key].update(intersection)
                    nodes.update(set(community))
                else:
                    nodes.update(set(community))
        
        return overlaping

def generate_pkl(path: str) -> None:

    '''
    Generate pkl files from txt files

    Parameters
    ----------
    path: str
        Path to folder containing txt files

    Returns
    -------
    None
    '''

    files = os.listdir('dataset/' + path)
    files.remove('GT')
    files.remove('README.txt')

    for file in files:
        G = nx.Graph()
        with open('dataset/' + path + '/' + file + '/' + file + '.dat', 'r') as f:
            lines = f.readlines()
                
        for line in lines:
            a, b = line.split('\t')
            G.add_edge(int(a) - 1, int(b) - 1)

        #     #nx.nx_pylab.draw(G, with_labels=True)
        #     #plt.show()
        pickle.dump(G, open('dataset/' + path + '/' + file + '/' + file + '.pkl', 'wb'))
        G.clear()

def runAlgorithmSimple(m, folder_version = 'NetsType_1.3'):

    for j in range(1, 12):

        m.G = pickle.load(open('dataset/' + folder_version + '/network'+ str(j) + '/network'+ str(j) + '.pkl', 'rb'))

        n = 0
        top = 1

        exportpath_Simple = folder_version

        for i in range(n, top):
            result = nx.algorithms.community.label_propagation.asyn_lpa_communities(m.G, seed=random.randint(0, 10000))
            communities = [list(x) for x in result]
            m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Lpa.txt', communities)
        
        for i in range(n, int(top)):
            result = nx.algorithms.community.greedy_modularity_communities(m.G, resolution= 1)        
            communities = [list(x) for x in result]
            m.export_Simple(exportpath_Simple, '/network'+ str(j) +'_Greedy.txt', communities)
        
        for i in range(n, top):
            result = nx.algorithms.community.louvain.louvain_communities(m.G, seed=random.randint(0, 10000))
            communities = [list(x) for x in result]
            m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Louvain.txt', communities)

    print('done')


def runAlgorithmSimpleTunning(m, init, step, top, folder_version = 'NetsType_1.1_Tunning'):

    exportpath_Simple = folder_version
    folder_version = folder_version.split('_Tunning')[0]

    for j in range(1, 12):

        m.G = pickle.load(open('dataset/' + folder_version + '/network'+ str(j) + '/network'+ str(j) + '.pkl', 'rb'))

        size = (top - init)/step
        resolution = [init+(step*i) for i in range(int(size)+1)]


        # for i in range(n, top):
        #     seed_i = random.randint(init, end)
        #     result = nx.algorithms.community.label_propagation.asyn_lpa_communities(m.G, seed=seed_i)
        #     communities = [list(x) for x in result]
        #     m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Lpa_' + str(i)  + '.txt', communities)

        
        for rs_i in resolution:
            result = nx.algorithms.community.greedy_modularity_communities(m.G, resolution = rs_i)        
            communities = [list(x) for x in result]
            m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Greedy_' + str(rs_i)  + '.txt', communities)
        
        
        # for rs_i in resolution:
        #     result = nx.algorithms.community.louvain.louvain_communities(m.G, resolution = rs_i, seed=random.randint(0, 10000))
        #     communities = [list(x) for x in result]
        #     m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Louvain_' + str(rs_i)  + '.txt', communities)

    print('done')


def runAlgorithmSimpleTunning(m, init, step, top, folder_version = 'NetsType_1.1_Tunning'):

    exportpath_Simple = folder_version
    folder_version = folder_version.split('_Tunning')[0]

    for j in range(1, 12):

        m.G = pickle.load(open('dataset/' + folder_version + '/network'+ str(j) + '/network'+ str(j) + '.pkl', 'rb'))

        size = (top - init)/step
        resolution = [init+(step*i) for i in range(int(size)+1)]


        # for i in range(n, top):
        #     seed_i = random.randint(init, end)
        #     result = nx.algorithms.community.label_propagation.asyn_lpa_communities(m.G, seed=seed_i)
        #     communities = [list(x) for x in result]
        #     m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Lpa_' + str(i)  + '.txt', communities)

        
        for rs_i in resolution:
            result = nx.algorithms.community.greedy_modularity_communities(m.G, resolution = rs_i)        
            communities = [list(x) for x in result]
            m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Greedy_' + str(rs_i)  + '.txt', communities)
        
        
        # for rs_i in resolution:
        #     result = nx.algorithms.community.louvain.louvain_communities(m.G, resolution = rs_i, seed=random.randint(0, 10000))
        #     communities = [list(x) for x in result]
        #     m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Louvain_' + str(rs_i)  + '.txt', communities)

    print('done')

def drawResultAlgorithm(folderpath, nameFile):

    print('Begin!!!!!!!!!!!')

    dictResult = pickle.load(open('output/' + folderpath + '/' + nameFile, 'rb'))

    df = pd.DataFrame()

    for _ , v in dictResult.items():
        df = df.append(v, ignore_index = True)
    
    columnsSorted = sorted(df.columns)
    nameColumns = columnsSorted[0]
    columnsSorted.remove(nameColumns)
    columnsSorted = sorted(columnsSorted, key=lambda x: int("".join([i for i in x if i.isdigit()])))
    columnsSorted.insert(0, nameColumns)
    
    df = df.reindex(columnsSorted, axis=1)
    df = df.set_index(nameColumns)
    
    print('created df done')
    print(df.columns)
    print(df)
    dfT = df.transpose()
    print(dfT)
    
    
    print('df transpose done')
    for item in dfT.columns:
        dfT[item].plot()
    
    plt.title('Run Algorithms and NMI accuracy' + ' Network: ' + folderpath)
    plt.xlabel('Nets')
    plt.ylabel('NMI Accuracy')
    plt.legend()
    plt.show()

   
if __name__ == "__main__":

    # create Data Structure
    m = Matrix([], {},[])
    # run Algorithm simple
    #runAlgorithmSimpleTunning(m, 5.5, 0.5, 10.0, 'NetsType_1.1_Tunning')
    # draw result
    drawResultAlgorithm('NetsType_1.6', 'NetsType_1.6_result.pkl')

    #G = pickle.load(open('dataset/NetsType_1.6/network10/network10.pkl', 'rb'))

    # result = nx.algorithms.community.louvain.louvain_communities(G, seed=random.randint(0, 10000), resolution=random.uniform(2,3.5)) # type: ignore

    # communities = [list(x) for x in result] # type: ignore

    # pickle.dump(communities, open('dataset/NetsType_1.6/network10/network10_Louvain.pkl', 'wb'))

    # from cdlib import evaluation, NodeClustering

    # communities = pickle.load(open('dataset/NetsType_1.6/network10/network10_Louvain.pkl', 'rb'))



    # nodes= []

    # with open('dataset/' + 'NetsType_1.6' + '/GT/community' + '10' + '_GT.dat', 'r') as f:
    #             lines = f.readlines()        
    #             for line in lines:
    #                 data = line.split(' ')
    #                 inter_data = [int(x) for x in data]
    #                 nodes.append(inter_data)

    # nodeClustA = NodeClustering(communities=nodes, graph=G, method_name='GT', method_parameters={}, overlap=True)
    # nodeClustB = NodeClustering(communities=communities, graph=G, method_name='Louvain', method_parameters={}, overlap=True)

    # match_resoult = evaluation.overlapping_normalized_mutual_information_MGH(nodeClustA, nodeClustB)

    # print(match_resoult)