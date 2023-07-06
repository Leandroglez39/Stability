import pandas as pd
import pickle
import multiprocessing
import time
from multiprocessing import Array
import sys
import concurrent.futures
from timeit import default_timer as timer


def temp():
    df = pd.read_csv('./data/matrix/0 file.csv')

    nodes = df.columns.to_list()    

    ady_list = [[] for _ in range(len(nodes))]

    for i in range(df.shape[0]):
        row = df.loc[i].to_list()
        for j in range(len(row)):
            if row[j] != 0:
                ady_list[i].append((j, row[j]))


    with open('./data/adym_0.pkl', 'wb') as f:
        pickle.dump(ady_list, f)

    a = pickle.load(open('./data/adym_0.pkl', 'rb'))

    print(a)

def cpu_bound(number):
    return sum(i * i for i in range(number))


def find_sums(numbers):
    with multiprocessing.Pool(8) as pool:
        pool.map(cpu_bound, numbers)

def mdouble(x):
    return x * 2

def lovain_concurrent(G, weight = 'weight', resolution = 1, threshold = 1e-07, seed = 1 , n = 10):

    '''
    This functiosn is for execute louvain algorithm in parallel.

    Parameters
    ----------
    G : NetworkX graph
    weight : string or None, optional (default="weight")
        The name of an edge attribute that holds the numerical value
        used as a weight. If None then each edge has weight 1.
    resolution : float, optional (default=1)
        If resolution is less than 1, the algorithm favors larger communities.
        Greater than 1 favors smaller communities
    threshold : float, optional (default=0.0000001)
        Modularity gain threshold for each level. If the gain of modularity
        between 2 levels of the algorithm is less than the given threshold
        then the algorithm stops and returns the resulting communities.
    seed : integer, random_state, or None (default)
        Indicator of random number generation state.
        See :ref:`Randomness<randomness>`.
    n :int, optional (default=10)
        Number of times to execute the algorithm.

    Returns
    -------
    list 
        A list of sets (partition of `G`). Each set represents one community and contains
        all the nodes that constitute it.
    '''

    import networkx.algorithms.community as nx_comm

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        communities = pool.starmap(nx_comm.louvain_communities, [(G, weight, resolution, threshold, seed) for _ in range(n)])
    return communities

def lpa_wrapper(G, weight = 'weight', seed = 1):

    import networkx.algorithms.community as nx_comm
    return list(nx_comm.asyn_lpa_communities(G, weight, seed))

def asyn_lpa_concurrent(G, weight = 'weight', seed = 1 , n = 10):
    

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        communities = pool.starmap(lpa_wrapper, [(G, weight, seed) for _ in range(n)])
    

    return [com for com in communities]

def greedy_modularity_concurrent(G, weight=None, resolution=1, cutoff=1, best_n=None, n = 10):
    
    '''
    This functiosn is for execute greedy modularity algorithm in parallel.

    Parameters
    ----------
    G : NetworkX graph

    weight : string or None, optional (default=None)
        The name of an edge attribute that holds the numerical value used
        as a weight.  If None, then each edge has weight 1.
        The degree is the sum of the edge weights adjacent to the node.

    resolution : float, optional (default=1)
        If resolution is less than 1, modularity favors larger communities.
        Greater than 1 favors smaller communities.

    cutoff : int, optional (default=1)
        A minimum number of communities below which the merging process stops.
        The process stops at this number of communities even if modularity
        is not maximized. The goal is to let the user stop the process early.
        The process stops before the cutoff if it finds a maximum of modularity.

    best_n : int or None, optional (default=None)
        A maximum number of communities above which the merging process will
        not stop. This forces community merging to continue after modularity
        starts to decrease until `best_n` communities remain.
        If ``None``, don't force it to continue beyond a maximum.

    n :int, optional (default=10) 
        Number of times to execute the algorithm. 
    
    Returns:
        list (frozenset): A list of sets (partition of G). Each set represents one community and contains all the nodes that constitute it.
    
    '''
    
    
    import networkx.algorithms.community as nx_comm

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        communities = pool.starmap(nx_comm.greedy_modularity_communities, [(G, weight ,resolution, cutoff,best_n) for _ in range(n)])
    return communities

def infomap_concurrent(G, n = 10):

    '''
    This functiosn is for execute infomap algorithm in parallel.

    Args:
        G (networkx.Graph): Graph to be clustered.
        n (int, optional): Number of times to execute the algorithm. Defaults to 10.
    Returns:
        list (cdlib.classes.node_clustering.NodeClustering): List of communities.
            
    NodeClustering type Properties:

        communities: List of communities
        graph: A networkx/igraph object
        method_name: Community discovery algorithm name
        method_parameters: Configuration for the community discovery algorithm used
        overlap: Boolean, whether the partition is overlapping or not

    '''
    
    from cdlib import algorithms

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        communities = pool.map(algorithms.infomap, [G for _ in range(n)])
    return communities

data_a = Array('i', 198034851, lock=False)

def match_count_parallel(tuple):
    pos = (19902 * tuple[0]) - sum([i for i in range(tuple[0])])
    pos = pos + tuple[1] - tuple[0] 
    data_a[pos] = data_a[pos] + 1

def match_count_parallel_return(tuple):
    pos = (19902 * tuple[0]) - sum([i for i in range(tuple[0])])
    pos = pos + tuple[1] - tuple[0] 
    return pos

if __name__ == '__main__':
    
    n = 11000
    #params = [(i,j) for i in range(n) for j in range( i + 1, n)]

    print('Params ready')
    

    import networkx as nx

    # Create a graph
    G = nx.Graph()
    G.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4)])

    # Create two subgraphs
    subgraph1 = G.subgraph([1, 2, 3])
    subgraph2 = G.subgraph([2, 3, 4])

    print(subgraph1.edges)
    print(subgraph2.edges)

    # Calculate the number of edges between the subgraphs
    
    a = nx.algorithms.community.quality.inter_community_edges(G, subgraph2.nodes)
    # Print the result
    print("Number of edges between subgraphs:", a)


    # concurrent
    # comment out to only run sequential
    start = timer()
    result = []

    # with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
    #     result = pool.map(match_count_parallel, params)
    
    # for i in params:
    #     match_count_parallel(i)
    
    # with concurrent.futures.ProcessPoolExecutor(max_workers= 20) as executor:
    #     futures = [executor.submit(match_count_parallel_return, i) for i in params]

    #     for i, future in enumerate(concurrent.futures.as_completed(futures)):
    #         if future.result():
    #             result.append(future.result())
    
    


    #print('Result 2:', result)
    print('Took: %.2f seconds.' % (timer() - start))

    print(len(result))
   

    print('Done')

    print(data_a)

    '''Region Concurrent Communities Algorithms'''
    
    #communities = lovain_concurrent(G)
    
    #communities = asyn_lpa_concurrent(G) 

    #communities = greedy_modularity_concurrent(G)

    #print(communities)

    

    #communities = infomap_concurrent(G)

    #print(communities[0].communities )

    ''' End Region Concurrent Communities Algorithms '''

    # from infomap import Infomap 

    # from cdlib import algorithms

    # nx_communities = algorithms.infomap(G)

    # print(nx_communities.communities)

    # numbers = [5000000 + x for x in range(20)]

    # start_time = time.time()
    # find_sums(numbers)
    # duration = time.time() - start_time
    # print(f"Duration {duration} seconds")

    
    
  