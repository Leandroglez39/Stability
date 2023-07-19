import time
import pandas as pd
import pickle
from dataclasses import dataclass
import networkx as nx
import networkx.algorithms.community as nx_comm
import matplotlib.pyplot as plt
import multiprocessing
import os
import datetime
import math
import numpy as np
import scipy.stats as stats
import statistics
import random
import re
import itertools  
import concurrent.futures
import asyncio



@dataclass
class Matrix:
    """    
    Class to create a matrix from a csv file and export it to a graphml file    
    
    Attributes:
        list_nodes (list): list of nodes
        pos_code_nodes (dict): dictionary with the position of each node
        G (nx.DiGraph): graph
        ady_list (list): list of adyacents nodes
        G_small (nx.DiGraph): small graph
    """

    list_nodes: list
    pos_code_nodes: dict
    G = nx.DiGraph()
    ady_list : list 
    G_small = nx.DiGraph()
    
    def __post_init__(self):        
        #self.list_nodes = pickle.load(open('./dataset/nodes.pkl', 'rb'))
        #self.pos_code_nodes = pickle.load(open('./dataset/pos_code_nodes.pkl', 'rb'))
        #self.ady_list = [[] for _ in range(len(self.list_nodes))]
        pass

    def insert_nodes(self):
        self.G.add_nodes_from(self.list_nodes)
    
    def insert_weighted_edges(self):
        for i in range(len(self.ady_list)):
            for j in range(len(self.ady_list[i])):
                y = self.ady_list[i][j][0]
                weight = self.ady_list[i][j][1]
                self.G.add_edge(self.list_nodes[y], self.list_nodes[i], weight=weight)
    
    def export_graph_to_csv(self, path = "./dataset/graph_19k_3.5m.csv"):
        '''
        Export graph to csv file for Cosmos viewer
        '''
        nx.write_edgelist(self.G, path, delimiter=",", data=['weight']) # type: ignore

    def export_graph_to_csv_size(self, size):

        G = nx.DiGraph()
        G.add_nodes_from(self.list_nodes[:size])

        df = pd.read_csv('./dataset/matrix/0 file.csv')
        

        ady_list = [[] for _ in range(size)]

        for i in range(len(self.list_nodes[:size])):
            row = df.loc[i].to_list()
            for j in range(size):
                if row[j] != 0:
                    ady_list[i].append((j, row[j]))

        for i in range(len(ady_list)):
            for j in range(len(ady_list[i])):
                y = ady_list[i][j][0]
                weight = ady_list[i][j][1]
                G.add_edge(self.list_nodes[y], self.list_nodes[i], weight=weight)

        self.G_small = G

        nx.readwrite.edgelist.write_edgelist(G, f"0set_size{str(size)}.csv", delimiter=",", data=['weight'])        
        nx.gexf.write_gexf(G, f"0set_size{str(size)}.gexf")
        
    def export_graph_to_graphml(self, path = "./dataset/graph_19k_3.5m.gml"):        
        nx.graphml.write_graphml(self.G, path)

    def export_graph_to_adjlist(self, path = "./dataset/graph_19k_3.5m.adyl"):        
        nx.adjlist.write_adjlist(self.G, path)

    def load_ady_matrix(self, count = 0):

        for x in range(count+1):
            print(f'Loading: {x} file.csv')
            df = pd.read_csv(f'./dataset/matrix/{x} file.csv')            
            for i in range(df.shape[0]):
                row = df.loc[i].to_list()
                for j in range(len(row)):
                    if row[j] != 0:
                        self.ady_list[i + x * 642].append((j, row[j]))
            print(f'Finished: {x} file.csv')


        with open(f'./dataset/adym_{count}.pkl', 'wb') as f:
            pickle.dump(self.ady_list, f)
   
    def save_graph_obj(self, path = './dataset/graph-1.x.pkl'):
        
        with open(path, 'wb') as f:
            pickle.dump(self.G, f)

    def load_matrix_obj(self, path = './dataset/graph_19k_3.5m.pkl'):
        self.G = pickle.load(open(path, 'rb'))
    
    def read_adym(self, path = './dataset/adym_30.pkl'):
        self.ady_list = pickle.load(open(path, 'rb'))

    def save_attributed_graph(self, path = './dataset/outputs/attributed_graph.csv'):
        '''
        Save the graph with the attributes of each node in csv file.
        '''

        with open(path, 'w') as f:
            for node in self.G.nodes:
                
                f.write(f'{node},{str(nx.degree(self.G, node))},{str(self.G.in_degree[node])},{str(self.G.out_degree[node])},') # type: ignore
                f.write(str(self.G.degree(node, weight='weight')) + ',') # type: ignore
                f.write(str(self.G.in_degree(node, weight='weigth')) + ',')
                f.write(str(self.G.out_degree(node, weight='weigth')) + ',')
                f.write(str(self.G.nodes[node]['eigenvector_centrality']) + ',')
                f.write(str(self.G.nodes[node]['eigenvector_centrality_weighted']) + ',')
                f.write(str(self.G.nodes[node]['pagerank']) + ',')
                f.write(str(self.G.nodes[node]['degree_centrality']) + ',')
                f.write(str(self.G.nodes[node]['core_number']) + '\n')

    def export_infomap_iterations(self, folder_version = 'NetsType_1.0', init = 0, end = 10):
        '''
        Export the infomap iterations to pickle file.

        Parameters
        ----------
        path : str, optional (default='./dataset/outputs/infomap_iterations.pkl')
            The path where the pickle file will be saved.
        init : int, optional (default=0)
            The initial iteration range.
        end : int, optional (default=10)
            The final iteration range.

        '''

        from cdlib import algorithms

        all_iterations = []
        count = 0
        for j in range(1, 12):

            self.G = pickle.load(open('dataset/' + folder_version + '/network'+ str(j) + '/network'+ str(j) + '.pkl', 'rb'))
            
            for _ in range(init, end):
                result = algorithms.infomap(self.G, flags='--seed ' + str(random.randint(0,1000)))
                result = result.communities
                communities = [list(x) for x in result]

                if count == 0:
                    self.export_Simple(folderpath = folder_version, filepath= '/network'+ str(j) + '_Infomap' + '.txt', result= communities)
                    count += 1

                if len(communities) > 1:
                    all_iterations.append(communities)

            count = 0

            pickle.dump(all_iterations, open('output/' + folder_version + '/network'+ str(j) + '_Infomap' + '.pkl', 'wb'))

            all_iterations = []

           

    # ALGORITMOS DE COMUNIDADES

    def lovain_concurrent(self, weight = 'weight', resolution = 1, threshold = 1e-07, seed = None , n = 10):

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
        seed : list (lenght=n), random_state, or None (default)
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

        if seed:
            with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
                communities = pool.starmap(nx_comm.louvain_communities, [(self.G, weight, resolution, threshold, seed[i]) for i in range(n)])
        else:            
            with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
                communities = pool.starmap(nx_comm.louvain_communities, [(self.G, weight, resolution, threshold) for _ in range(n)])

        return communities

    def lpa_wrapper(self, G, weight = 'weight', seed = 1):

        import networkx.algorithms.community as nx_comm
        return list(nx_comm.asyn_lpa_communities(G, weight, seed)) # type: ignore
        

    def asyn_lpa_concurrent(self, weight = 'weight', seed = None , n = 10):
        
        '''  
        This functiosn is for execute asyn_lpa algorithm in parallel.

        Parameters
        ----------
        G : Graph

        weight : string
            The edge attribute representing the weight of an edge.
            If None, each edge is assumed to have weight one. In this
            algorithm, the weight of an edge is used in determining the
            frequency with which a label appears among the neighbors of a
            node: a higher weight means the label appears more often.

        seed : list(integer) with length = n , random_state, or None (default = 1)
            Indicator of random number generation state.
            See :ref:`Randomness<randomness>`.
        
        n :int, optional (default=10)
            Number of times to execute the algorithm.

        Returns
        -------
        communities : iterable
            Iterable of communities given as sets of nodes.

        Notes
        -----
        Edge weight attributes must be numerical.

        '''

        if seed:
            with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
                communities = pool.starmap(self.lpa_wrapper, [(self.G, weight, seed[i]) for i in range(n)])
        else:
            with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
                communities = pool.starmap(self.lpa_wrapper, [(self.G, weight) for _ in range(n)])
        

        return [com for com in communities]

    def greedy_modularity_concurrent(self, weight=None, resolution=1, cutoff=1, best_n=None, n = 10):
        
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
            communities = pool.starmap(nx_comm.greedy_modularity_communities, [(self.G, weight ,resolution, cutoff,best_n) for _ in range(n)])
        return communities
    
    def infomap_concurrent(self, seed = None, n = 10):

        number_of_pool = n if n < multiprocessing.cpu_count() else multiprocessing.cpu_count()
        '''
        This functiosn is for execute infomap algorithm in parallel.
        Parameters
        ----------

            G (networkx.Graph): Graph to be clustered.

            n (int, optional): Number of times to execute the algorithm. Defaults to 10.

            seed : list(integer) with length = n , random_state, or None (default = 1)
            Indicator of random number generation state.
            See :ref:`Randomness<randomness>`.

        Returns:
        ----
            list (cdlib.classes.node_clustering.NodeClustering): List of communities.
                
        NodeClustering type Properties:
            communities: List of communities
            graph: A networkx/igraph object
            method_name: Community discovery algorithm name
            method_parameters: Configuration for the community discovery algorithm used
            overlap: Boolean, whether the partition is overlapping or not
        '''
        from cdlib import algorithms # type: ignore Only necesary with cdlib environment #TODO

        if seed:
            with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
                communities = pool.starmap(algorithms.infomap, [(self.G, '--seed ' + f'{seed[i]}') for i in range(n)])
        else:
            with multiprocessing.Pool(number_of_pool) as pool:
                communities = pool.map(algorithms.infomap, [self.G for _ in range(n)])

        return communities
    # Tools for Algorithms

    def communities_length(self, communities):

        '''
        This function is for calculate the length of each partition in the community.
        '''

        a = [len(com) for com in communities]
        a.sort(reverse=True)

        return a
    
    def save_communities(self, communities, algorithm : str, params: list = [], seed = None):

        if algorithm == 'louvain':
            
            if params == []:
                params = [1, 1e-07, None]

            for i in range(len(communities)):
                params_name = ''.join(['_' + str(param) for param in params])
                params_name += '_seed_' + str(seed[i]) if seed else ''

                with open('./dataset/outputs/' + algorithm + '/' + algorithm + params_name + '_iter_' + str(i) , 'wb+') as f:
                    pickle.dump(communities[i], f)

        elif algorithm == 'greedy':

            if params == []:
                params = [1, 1, 0]

            for i in range(len(communities)):
                params_name = ''.join(['_' + str(param) for param in params])

                with open('./dataset/outputs/' + algorithm + '/' + algorithm + params_name + '_iter_' + str(i) , 'wb+') as f:
                    pickle.dump(communities[i], f)

        elif algorithm == 'lpa':

            for i in range(len(communities)):

                params_name = ''
                params_name += '_seed_' + str(seed[i]) if seed else ''

                with open('./dataset/outputs/' + algorithm + '/' + algorithm + params_name + '_iter_' + str(i) , 'wb+') as f:
                    pickle.dump(communities[i], f)

        elif algorithm == 'infomap':

            for i in range(len(communities)):

                params_name = ''
                params_name += '_seed_' + str(seed[i]) if seed else ''

                with open('./dataset/outputs/' + algorithm + '/' + algorithm + params_name + '_iter_' + str(i) , 'wb+') as f:
                    pickle.dump(communities[i], f)


    def load_communities(self, algorithm : str, resolution = 1, threshold = 1e-07 , seed = 1, iter = 0) -> list:

        if algorithm == 'louvain':
            path = './dataset/outputs/' + algorithm + '/' + algorithm + '_'+ str(resolution) + '_' + str(threshold) + '_seed_' + str(seed) + '_iter_' + str(iter)
        
            with open(path, 'rb') as f:
                return pickle.load(f)
        
        return []

    def load_all_communities(self, algorithm : str, infomap_flags = False) -> list: 
        '''
        This function is for load all the communities generated by one algorithm.

        Parameters
        ----------
        algorithm : str
            The name of the algorithm to load the communities.
        infomap_flags : bool, optional
            If the algorithm is infomap, this parameter is for load the communities and returned the communities in list of list format. The default is False.
        
        Returns
        -------
        list : list
            A list of communities present in the folder of the algorithm.
        '''

        if algorithm == 'louvain':
            
            paths = os.listdir('./output/' + algorithm + '/')
        
            communities = []
            for path in paths:
                with open('./output/' + algorithm + '/' + path, 'rb') as f:
                    community = pickle.load(f)
                    sorted_community = sorted(community, key = lambda x: len(x), reverse = True)
                    communities.append(sorted_community)

            return communities
        if algorithm == 'lpa':

            paths = os.listdir('./output/' + algorithm + '/')
        
            communities = []
            
            for path in paths:
                with open('./output/' + algorithm + '/' + path, 'rb') as f:
                    communities.append(pickle.load(f))
            
            return communities
        
        if algorithm == 'greedy':

            paths = os.listdir('./output/' + algorithm + '/')

            communities = []

            for path in paths:
                with open('./output/' + algorithm + '/' + path, 'rb') as f:
                    communities.append(pickle.load(f))
        
            
            return communities

        if algorithm == 'infomap':

            paths = os.listdir('./output/' + algorithm + '/')

            communities = []

            for path in paths:
                with open('./output/' + algorithm + '/' + path, 'rb') as f:
                    communities.append(pickle.load(f))
        
            if infomap_flags:
                communities = [community['communities'] for community in communities]

                iteration = []
                for community in communities:
                    iteration.append([set(com) for com in community])

                return iteration
            
            return communities
        
        return []    

    def load_all_algorithm_communities(self, algorithms: list) -> list:

        '''
        This function is for load all the communities generated by a list of algorithms.

        Parameters
        ----------
        algorithms : list
            A list of algorithms's name to load the communities.
        
        Returns
        -------
        list : list
            A list of communities present in the folder of the algorithms.
        '''

        iterations = []
        allowed_algorithms = ['louvain', 'lpa', 'greedy', 'infomap']

        for algorithm in algorithms:
            
            if algorithm in allowed_algorithms:
                if algorithm != 'infomap':
                    iterations.extend(self.load_all_communities(algorithm=algorithm))
                else:
                    iterations.extend(self.load_all_communities(algorithm=algorithm, infomap_flags=True))

            

        return iterations
    
    def compare_communities_limits(self, communities, limit: list = [1,1]):

        '''
        This function is for compare the communities generated by an algorithm with a limit.

        Parameters
        ----------
        communities : list
            A list of communities generated by one algorithm.
        limit : list (default = [1,1])
            A list of two elements. The first element is the cuantity of the biggest communities (by length) must be compared.
            The second element is the cuantity of the smallest communities (by length) must be compared.
        
        Returns
        -------
        list : list with two elements. 
            A list of elemments .
        '''

        return None
    
    def create_sub_graph(self, community):

        '''
        This function is for create a subgraph from one community.

        Parameters
        ----------
        community : list
            A list of nodes that form a community.
        
        Returns
        -------
        Graph : NetworkX DiGraph
            A subgraph of the original graph.
        '''

        return self.G.subgraph(community)
    
    def add_degree_property(self, G):

        '''
        This function is for add a property to the nodes of the graph. The property is the degree of the node.

        Parameters
        ----------
        G : NetworkX DiGraph
            A graph.
        
        Returns
        -------
        Graph : NetworkX DiGraph
            A graph with the property degree.
        '''

        for node in G.nodes():
            G.nodes[node]['degree'] = G.degree(node)

        return G

    def nodes_in_communities(self, communities: list):

        '''
        This function is for return a dict with the name of the nodes is teh key and the values
         is the communities to they belong.

        Parameters
        ----------
        communities : list
            A list of communities generated by one algorithm.
        
        Returns
        -------
        result : dict
            A dict of nodes in their communities.
        '''

        nodes = {}
        for i in range(len(communities)):
            for node in communities[i]:
                nodes[node] = i

        return nodes

    def save_dict_to_csv(self, dict, name):

        '''
        This function is for save a dict to a csv file.

        Parameters
        ----------
        dict : dict
            A dict to save.
        name : str
            The name of the csv file.
        '''

        with open('./dataset/outputs/graph/' + name + '.csv', 'w') as f:
            for key, value in dict.items():
                f.write(str(key) + ',' + str(value) + '\n')

    def add_property(self, measure: list):

        '''
        This function is for add a property to the nodes of the graph.

        Parameters
        ----------
        property : str
            The name of the property to add.
        '''

        if 'eigenvector_centrality' in measure:

            data = nx.eigenvector_centrality(self.G)

            data_weighted = nx.eigenvector_centrality(self.G, weight='weight') 

            for node in self.G.nodes():
                self.G.nodes[node]['eigenvector_centrality'] = data[node] # type: ignore
                self.G.nodes[node]['eigenvector_centrality' + '_weighted'] = data_weighted[node] # type: ignore
            
            print('eigenvector_centrality added')
        
        if 'pagerank' in measure:

            data = nx.pagerank(self.G)           

            for node in self.G.nodes():
                self.G.nodes[node]['pagerank'] = data[node]

            print('pagerank added')
        
        if 'degree_centrality' in measure:

            data = nx.degree_centrality(self.G)           

            for node in self.G.nodes():
                self.G.nodes[node]['degree_centrality'] = data[node]

            print('degree_centrality added')
        
        if 'core_number' in measure:

            data = nx.core_number(self.G)            # type: ignore

            for node in self.G.nodes():
                self.G.nodes[node]['core_number'] = data[node] # type: ignore
            
            print('core_number added')
        
    def participation_coefficient(self, communities: list):

        '''
        This function is for calculate the participation coefficient of a community.

        Parameters
        ----------        
        communities : list
            A list of comunities.
        
        Returns
        -------
        result : dict
            The participation coefficient of every node.
        '''

        data_nodes = {}
        
        start_time = datetime.datetime.now()

        count = 0
        
        for node in self.G.nodes():

            k_i = self.G.degree(nbunch=node) # type: ignore

            suma = 0
            
            neighbors = list(nx.neighbors(self.G, node)) # type: ignore

            for community in communities:

                subgraph = nx.subgraph(self.G, community) # type: ignore
                                
                k_i_s = 0

                for n in neighbors:
                    
                    count += 1
                    if count == 10000000:
                        print('10M operations calculated in participation coefficient')
                        print('Time for 10M operations: ' + str(datetime.datetime.now() - start_time))
                        count = 0
                        start_time = datetime.datetime.now()

                    if n in subgraph.nodes():
                        k_i_s += 1
                        continue

                suma += (k_i_s / k_i) ** 2 # type: ignore

            data_nodes[node] = 1 - suma

        return data_nodes

    def insert_measure_dict(self, measure: str, dict: dict):

        '''
        This function is for insert a measure to the nodes of the graph.

        Parameters
        ----------
        measure : str
            The name of the measure.
        dict : dict
            A dict with the values of the measure.
        '''

        for node in self.G.nodes():
            self.G.nodes[node][measure] = dict[node]
    

    def RoughClustering(self, communities: list):

        '''
        This function is for calculate the rough clustering of a community list.

        Parameters
        ----------        
        communities : list
            A list of list of comunities.
        
        Returns
        -------
        result : dict
            The rough clustering of every node.
        '''
        

         # Create a list of nodes from the graph
        nodes = list(self.G.nodes())
        
        match_array = np.zeros((len(nodes), len(nodes)), dtype=int)

        node_hash = {}

        for i in range(len(nodes)):
            node_hash[nodes[i]] = i
                
        Gb0 = nx.Graph()

        Gb0.add_nodes_from(nodes)
              
        print('Updating ocurances in match array')
        start_time = datetime.datetime.now()
        for community in communities:
            self.update_match_array(community=community, match_array=match_array, hash=node_hash)
            #self.update_match_array_faster(communities=community, match_array=match_array, hash=node_hash)        
        
        end_time = datetime.datetime.now()
        print(f'Ocurances updated in match array in: {end_time - start_time}')

        b0 = len(communities) * 0.75

        edges_list = []

        print('Calculating edges')        

        # Find the indices where match_array is greater than or equal to b0
        i, j = np.where(match_array >= b0)

        # Use those indices to create the edges_list
        edges_list = [(nodes[i[k]], nodes[j[k]]) for k in range(len(i))]

        print('Edges calculated')
        
        
        Gb0.add_edges_from(edges_list)


        #nx.draw(Gb0, with_labels=True, font_weight='bold')
        #plt.show()

        print('Calculating strong components')
        strongcomponents = nx.connected_components(Gb0)

        seeds = []

        for component in strongcomponents:
            seeds.append(nx.function.subgraph(self.G, component))
            
        print('Strong components calculated. n = '+ str(len(seeds)))

        seeds.sort(key=len, reverse=True)

        print('Calculating k')
        k = self.calculate_k(communities)

        print('Real k: ' + str(k + 1))


        # list of set of nodes that represents the  coverage of the graph . The first set is the inferior coverage and the second set is the superior coverage.
        coverage_inferior = [set() for _ in range(k + 1)]
        coverage_superior = [set() for _ in range(k + 1)]


        for i in range(k + 1):
            coverage_inferior[i] = coverage_inferior[i].union(set(seeds[i].nodes())) # type: ignore
            coverage_superior[i] = coverage_superior[i].union(set(seeds[i].nodes())) # type: ignore

        similarity_values = []

        print('Calculating similarity')
        start_time = datetime.datetime.now()

        for j in range( k + 1, len(seeds)):


            grj1 = seeds[j] # Subgraph of the seed that will be calculated its similarity with all k seeds

            #similarity_values = [self.similarity_between_subgraphs(grj1, coverage_inferior[i], match_array, node_hash) for i in range(k + 1)]

            similarity_values = [self.similarity_between_subgraphs_faster(grj1, coverage_inferior[i], match_array, node_hash) for i in range(k + 1)]
                        
            edges_values = [self.edges_between_subgraphs(grj1, seeds[i]) for i in range(k + 1)]
            
            #total_nodes = sum([len(seed.nodes()) for seed in seeds]) + len(grj1.nodes())
            #total_edges = (total_nodes * (total_nodes-1)) / 2
            total_edges = max(edges_values)

            if total_edges == 0:
                edges_values_normalized = [0 for _ in range(k + 1)]
            else:
                edges_values_normalized = [x / total_edges for x in edges_values]

            max_similarity =  max(similarity_values) 
            
            
            
            
            if max_similarity == 0:
                similarity_values_normalized = [0 for _ in range(k + 1)]
            else:
                similarity_values_normalized = [svalues / max_similarity for svalues in similarity_values]

            

            similarity_values_normalized = [(similarity_values_normalized[i] + edges_values_normalized[i])/2 for i in range(len(similarity_values_normalized))]
            
            max_similarity_index = similarity_values_normalized.index(max(similarity_values_normalized)) 

            #print('Max similarity: ' + str(max_similarity))
            #print(similarity_values_normalized)
            #print(grj1.nodes())
            T = []

            gamma = 0.8

            for i in range(len(similarity_values_normalized)):
                if  max_similarity_index != i and similarity_values_normalized[i] >= gamma:
                    T.append(i)
            
            if len(T) > 0:
                for element in T:
                    coverage_superior[element] = coverage_superior[element].union(set(grj1.nodes()))
                    coverage_superior[max_similarity_index] = coverage_superior[max_similarity_index].union(set(grj1.nodes()))
            else:
                coverage_superior[max_similarity_index] = coverage_superior[max_similarity_index].union(set(grj1.nodes()))
                coverage_inferior[max_similarity_index] = coverage_inferior[max_similarity_index].union(set(grj1.nodes()))
        end_time = datetime.datetime.now()
        print(f'Similarity calculated in: {end_time - start_time}')        
        return (coverage_inferior, coverage_superior)

    
    def update_dict_count(self, dict_node: dict ,community: list):
        

        '''
        This function is for update the dict of the count of the communities.

        Parameters
        ----------
        community : list
            A list of nodes of a community.
        '''
        
        dict_nodes = dict_node.copy()

        # For each community
        for com in community:
            # For each node in the community
            com = list(com)
            com.sort()
            for i in range(len(com)):
                nodei = com[i]
                # For each node in the community except the current node
                for j in range(i+1, len(com)):
                    nodej = com[j]
                    # If the current combination of nodes is not in the dict
                    if (nodei, nodej) not in dict_nodes:
                        # Add it with a count of 1
                        dict_nodes[(nodei,nodej)] = 1
                    else:
                        # If it is in the dict, add 1 to the count
                        dict_nodes[(nodei,nodej)] = dict_nodes[(nodei,nodej)] + 1

        return dict_nodes             

    def update_match_array(self, match_array: np.ndarray, community: list, hash: dict):
        

        '''
        This function is for update the match array of the count of the communities.

        Parameters
        ----------
        match_array : np.ndarray
            A numpy array with the count of the ocurrances of the nodes in the communities.
        community : list
            A list of nodes of a community.
        '''
        
        # For each community
        for com in community:
            # For each node in the community
            com = list(com)
            com.sort()
            for i in range(len(com)):
                nodei = com[i]
                # For each node in the community except the current node
                for j in range(i+1, len(com)):
                    nodej = com[j]
                    # Increment the count of the current combination of nodes
                    nodei_hash = hash[nodei]
                    nodej_hash = hash[nodej]
                    match_array[nodei_hash, nodej_hash] += 1

    def update_match_array_faster(self, match_array: np.ndarray, communities: list, hash: dict):
    

        '''
        This function is for update the match array of the count of the communities.

        Parameters
        ----------
        match_array : np.ndarray
            A numpy array with the count of the ocurrances of the nodes in the communities.
        community : list
            A list of nodes of a community.
        '''
     

        for community in communities:
            for combination in itertools.combinations(community, 2):
                nodei_hash = hash[combination[0]]
                nodej_hash = hash[combination[1]]
                match_array[nodei_hash, nodej_hash] += 1


    def calculate_k(self, communities: list, statistic = True, ) -> int:

        '''
        This function is for calculate the k of a community.

        Parameters
        ----------        
        communities : list
            A list of comunities.
        statistic : bool
            If True, the function will calculate the k of all communities with the statistic method.
        
        Returns
        -------
        result : int
            The k of the community.
        '''

        k = 0
        data = []
        data_array = np.array([])

        if statistic:
            for community in communities:
                for i in range(len(community)):
                    partition = community[i]
                    vale = np.full(len(partition), i)
                    data_array = np.concatenate((data_array, vale), axis=None)
                    #data.append(len(partition))

            mean = statistics.mean(data_array)
            
            stddev = statistics.stdev(data_array)    

            print('mean: ' + str(mean))
            print('stddev: ' + str(stddev))
            

                        
            k = int(mean + (2 * stddev))
            return k
        else:
            list_data = []
            for community in communities:
                for i in range(len(community)):
                    partition = community[i]                    
                    list_data.append(len(partition))

            mean = statistics.mean(list_data)
            
            stddev = statistics.stdev(list_data)    

            print('mean: ' + str(mean))
            print('stddev: ' + str(stddev))
            

            count = 0
            for x in list_data:
                if x in range(int(mean + (-2 * stddev)) , int(mean + (2 * stddev))):
                    count += 1
            
            return int(count / len(communities))
    
    def merge_communities_dict(self, general_dict: dict, dict_new: dict):
    
        '''
        Merges two dictionaries into a single dictionary.

        Args:
        -----
            general_dict (dict): General dictionary where the values of the new dictionary will be added.
            dict_new (dict): New dictionary with values that will be added to the general dictionary.
        '''

        # For each key in the new dictionary:
        for key in dict_new:
            # If the key is not in the general dictionary
            if key not in general_dict:
                # Add it to the general dictionary
                general_dict[key] = dict_new[key]
            # If the key is in the general dictionary
            else:
                # Add the value from the new dictionary to the value in the general dictionary
                general_dict[key] = general_dict[key] + dict_new[key]

    def important_nodes(self):

        '''
        This function is for calculate the important nodes of the graph.

        Returns
        -------
        result : dict
            A dict with the important nodes.
        '''

       

        dict_degree = dict(self.G.degree()) # type: ignore

        dict_weighted = dict(self.G.degree(weight='weight')) # type: ignore

        list_degree = list(dict_degree.items())

        list_weighted = list(dict_weighted.items())

        list_degree.sort(key=lambda x: x[1], reverse=True)

        list_weighted.sort(key=lambda x: x[1], reverse=True)


        list_degree = list_degree[:7200]

        list_weighted = list_weighted[:7200]

        list_degree = [x[0] for x in list_degree]

        list_weighted = [x[0] for x in list_weighted]

        intersection = set(list_degree).intersection(set(list_weighted))

        return intersection

    def similarity_between_subgraphs(self, subgraph1: nx.Graph, subgraph2: set, match_array: np.ndarray, hash: dict) -> float:

        '''
        This function is for calculate the similarity between two subgraphs.
        
        Parameters
        ----------
        subgraph1 : nx.Graph
            A subgraph.
        subgraph2 : set
            A set of nodes with represent an inferior coverage.
        match_array : np.ndarray
            A numpy array with the count of the ocurrances of the nodes in the communities.
        hash : dict
            A dict with the hash of the nodes correspondaing to the index of the match array.
        
        Returns
        -------

        result : float
            The similarity between the two subgraphs between 0 and 1.
        '''
        result = 0

        
        
        for node1 in subgraph1.nodes():
            node_hash1 = hash[node1]
            for node2 in subgraph2:
                node_hash2 = hash[node2]
                result += match_array[node_hash1, node_hash2]
                
                

        
        return result / (len(subgraph1.nodes()) * len(subgraph2))
    
    def similarity_between_subgraphs_faster(self, subgraph1: nx.Graph, subgraph2: set, match_array: np.ndarray, hash: dict) -> float:
            '''
            This function is for calculate the similarity between two subgraphs.
            
            Parameters
            ----------
            subgraph1 : nx.Graph
                A subgraph.
            subgraph2 : set
                A set of nodes with represent an inferior coverage.
            match_array : np.ndarray
                A numpy array with the count of the ocurrances of the nodes in the communities.
            hash : dict
                A dict with the hash of the nodes correspondaing to the index of the match array.
            
            Returns
            -------
            result : float
                The similarity between the two subgraphs between 0 and 1.
            '''
            nodes1 = np.array([hash[node1] for node1 in subgraph1.nodes()])
            nodes2 = np.array([hash[node2] for node2 in subgraph2])

            result = np.sum(match_array[nodes1[:, None], nodes2])
            
            return result / (len(subgraph1.nodes()) * len(subgraph2))
    

    def edges_between_subgraphs(self, subgraph1: nx.Graph, subgraph2: nx.Graph) -> int:
            
            '''
            This function is for calculate the number of edges between two subgraphs.
            
            Parameters
            ----------
            subgraph1 : nx.Graph
                A subgraph.
            subgraph2 : nx.Graph
                A subgraph.
            
            Returns
            -------
    
            result : int
                The number of edges between the two subgraphs.
            '''
            
    
            total = len(subgraph1.edges()) + len(subgraph2.edges())


            subgraph3 = self.G.subgraph(list(subgraph1.nodes()) + list(subgraph2.nodes())) # type: ignore


            return len(subgraph3.edges()) - total
    
    def export_RC(self, folderpath, filepath, rc):

        '''
        This function is for export the RC file.

        Parameters
        ----------
        path : str
            The path where the file will be saved.
        rc : list(set)
            A list of sets with the communities. Include solapated nodes.
        '''
        if not os.path.exists(f'output/{folderpath}'):
            os.mkdir(f'output/{folderpath}')
        
        with open(f'output/{folderpath}{filepath}', 'w') as file:
            for i in range(len(rc[0])):
                temp = set()
                temp = temp.union(rc[0][i])
                temp = temp.union(rc[1][i])
                line = ''
                temp = list(temp)
                temp.sort()
                for node in temp:
                    line += str(node) + ' '                   
                file.write(line + '\n')

    def export_Simple(self, folderpath, filepath, result):

        '''
        This function is for export the RC file.

        Parameters
        ----------
        path : str
            The path where the file will be saved.
        rc : list(set)
            A list of sets with the communities. Include solapated nodes.
        '''
        if not os.path.exists('output/' + folderpath):
            os.mkdir('output/' + folderpath)
        with open('output/' + folderpath + filepath, 'w') as file:
            for ci in result:
                ci.sort()
                line = ''
                for node in ci:
                    line += str(node) + ' '                   
                file.write(line + '\n')

    def evaluate_nodesmatch_with_gt(self, nodesmatch: list, gt: list) -> list:
        
        response = []
        total_matches = 0
        with open('LFRBenchamark/community'+ str(j) + '_GT.dat', 'r') as f:
            lines = f.readlines()
            for comm in value[0]:
                maximum = 0
                for line in lines:
                    gt_comm = line.split(' ')
                    int_gt_comm = [int(x) for x in gt_comm]
                    if len(set(comm).intersection(set(int_gt_comm))) > maximum:
                        maximum = len(set(comm).intersection(set(int_gt_comm)))
                total_matches += maximum
        
        response.append(total_matches/1100)
        total_matches = 0

        return response

        
    # Begining of Horacio's Region

    def edgeWithinComm(self, vi, si, w, directed):
        sugSi = nx.subgraph(self.G, si) # type: ignore
        if not directed:
            ki = sugSi.degree(vi, weight = w)
        else:
            ki = sugSi.out_degree(vi, weight = w)
        return ki

    def degMeanByCommunity(self, si, w, directed):
        sum = 0
        for vi in si:
            edgeWith = self.edgeWithinComm(vi, si, w, directed)
            sum += edgeWith
        return sum/len(si)

    def standDesvCommunity(self, si, w, directed):
        dMeanCi = self.degMeanByCommunity(si, w, directed)
        sum = 0 
        for vi in si:
            ki = self.edgeWithinComm(vi, si, w, directed)
            sum += math.pow((ki - dMeanCi), 2)
        desv = math.sqrt(sum/len(si))
        return desv

    def withinCommunityDegree(self, w, commList, algName, directed):
        zi = dict()
        run_community_mean = dict()
        run_community_vi = dict()
        run_community_desv = dict()

        # iterate through all runs
        for run_i in range(0, len(commList)):
            clustering_i = commList[run_i]
            # iterate through all commnunities
            for c_j in range(0, len(clustering_i)):
                sumMean = 0
                key_RunComm = (run_i, c_j)
                for vi in clustering_i[c_j]:
                    key_RunCommVi = (run_i, c_j, vi)
                    if key_RunCommVi not in run_community_vi:
                        ki = self.edgeWithinComm(vi, clustering_i[c_j], w, directed)
                        run_community_vi[key_RunCommVi] = ki
                        sumMean += ki
                if key_RunComm not in run_community_mean:
                    ciMean = sumMean/len(clustering_i[c_j])
                    run_community_mean[key_RunComm] = ciMean
        # ki, mean(si) -> DONE
        print('ki, mean(si) -> DONE')

        # iterate through all runs
        for run_i in range(0, len(commList)):
            clustering_i = commList[run_i]
            # iterate through all commnunities
            for c_j in range(0, len(clustering_i)):
                sumDesv = 0
                key_runC = (run_i, c_j)
                meanCj = run_community_mean[key_runC]
                for vi in clustering_i[c_j]:
                    key_RunCVi = (run_i, c_j, vi)
                    kVi = run_community_vi[key_RunCVi]
                    sumDesv += math.pow((kVi - meanCj), 2)
                if key_runC not in run_community_desv:
                    desvCj = math.sqrt(sumDesv/len(clustering_i[c_j]))
                    run_community_desv[key_runC] = desvCj
        print('desviation -> DONE')
        # print(run_community_vi)

        
        # iterate through all vertex
        for vi in self.G.nodes:
            # iterate through all runs
            for run_i in range(0, len(commList)):
                clustering_i = commList[run_i]
                # iterate through all commnunities
                for c_j in range(0, len(clustering_i)):
                    # identify the Ci 
                    if vi in clustering_i[c_j]:
                        desvSi = run_community_desv[(run_i, c_j)]
                        ziValue = 0
                        if desvSi != 0:
                            meanSi = run_community_mean[(run_i, c_j)]
                            ki = run_community_vi[(run_i, c_j, vi)]
                            ziValue = (ki - meanSi)/desvSi
                        if vi not in zi:
                            listTupleValue = []
                            listTupleValue.append(ziValue)
                            zi[vi] = listTupleValue
                        else:
                            zi[vi].append(ziValue)
                        break
        
            # count+=1
            # print('vertex: ', vi, ' Done', ' numbers of items: ', count)
        print('zi -> DONE')

        direct = 'directed_' if directed else 'notDirected_'
        weg = 'weighted' if w == 'weight' else 'notWeighted'
        nameFile = 'within_' + algName + '_' + direct + weg
        # create a binary pickle file 
        f = open('./dataset/outputs/' + nameFile ,"wb")

        # write the python object (dict) to pickle file
        pickle.dump(zi, f)

        # close file
        f.close()

    def calculateWithin(self, algList = ['louvain', 'greedy', 'lpa', 'infomap'], wegList = ['weight', 'none'], directList = [True, False]):
        for alg_i in algList:
            for weg_j in wegList:
                for direct_k in directList:
                    commList =  m.load_all_communities(alg_i)
                    if alg_i == 'infomap':
                        commListModified = []
                        for ci in commList:
                            commListModified.append(ci['communities'])
                        self.withinCommunityDegree(weg_j, commListModified, alg_i, direct_k)
                    else:
                        self.withinCommunityDegree(weg_j, commList, alg_i, direct_k)
                    print('finish Directed: ', direct_k)
                print('finish Weight: ', weg_j)
            print('finish Algorithm: ', alg_i)

    def edgeMetricsExport(self, pathFile):
        with open(pathFile, 'w') as f:
            for u,v in self.G.edges:
                f.write(u + ', ' + v + ', ' + str(self.G.edges[u,v]['weight']) + ', ' + str(self.G.edges[u,v]['edge_betweenes']))
                f.write('\n')
                
        

        
    def update_graph_with_within_measure(self, communities ,algorithms = ['louvain', 'greedy', 'lpa', 'infomap']):

        for algorithm in algorithms:
                

            wdnw = pickle.load(open('dataset/outputs/FlyCircuitResult/within_'+ algorithm + '_directed_notWeighted', 'rb'))
            wdw = pickle.load(open('dataset/outputs/FlyCircuitResult/within_'+ algorithm + '_directed_weighted', 'rb'))
            wunw = pickle.load(open('dataset/outputs/FlyCircuitResult/within_'+ algorithm + '_notDirected_notWeighted', 'rb'))
            wuw = pickle.load(open('dataset/outputs/FlyCircuitResult/within_'+ algorithm + '_notDirected_weighted', 'rb'))

            files = [(wdnw, 'withing_directed_notWeighted'), (wdw, 'withing_directed_Weighted'), (wunw, 'withing_notDirected_notWeighted'), (wuw,  'withing_notDirected_Weighted')]

            for data, name in files:
                for k, v in data.items():
                    
                    for i in range(len(v)):

                        community = communities[i]
                        community_number = 0

                        if algorithm == 'infomap':
                            for j in range(len(community['communities'])):                        
                                if k in community['communities'][j]:
                                    community_number = j
                                    break
                        else:
                            for j in range(len(community)):                        
                                if k in community[j]:
                                    community_number = j
                                    break
                            
                        key = (algorithm, str(i), str(community_number))
                        if key not in self.G.nodes[k]['data'].keys():
                            self.G.nodes[k]['data'][(algorithm, str(i), str(community_number))] = {name: v[i]} # type: ignore
                        else:
                            self.G.nodes[k]['data'][(algorithm, str(i), str(community_number))][name] = v[i]
            

    # End Horacio's Region

    def apply_measures_to_communities_nodes(self, algorithm: str ,communities: list):
        
        

        for i in range(len(communities)):            
                
            start_time = datetime.datetime.now()
            iter_number = i
            iter_communities = communities[i]['communities'] if algorithm == 'infomap' else communities[i]
            data_participation = {}
            community_number = 0

            data_participation = self.participation_coefficient(iter_communities)

            print('Participation coefficient calculated for community ' + str(i) + ' of ' + str(len(communities)) + ' run of ' + algorithm + ' algorithm')
            print('Time for participation coefficient: ' + str(datetime.datetime.now() - start_time))

            for k, v in data_participation.items():

                for j in range(len(iter_communities)):
                    if k in iter_communities[j]:
                        community_number = j
                        break

                if 'data' not in self.G.nodes[k].keys():
                    self.G.nodes[k]['data'] = {(algorithm, str(iter_number), str(community_number)): {'participation_coefficient': v}} # type: ignore
                else:
                    if (algorithm, str(iter_number), str(community_number)) not in self.G.nodes[k]['data'].keys():
                        self.G.nodes[k]['data'][(algorithm, str(iter_number), str(community_number))] = {'participation_coefficient': v} # type: ignore
                    else:
                        self.G.nodes[k]['data'][(algorithm, str(iter_number), str(community_number))]['participation_coefficient'] = v
        
        # Withing Region

        # start_time = datetime.datetime.now()       
        # data_whiting_degree = self.withinCommunityDegree('None', communities)
        # print('Time for whiting degree: ' + str(datetime.datetime.now() - start_time))
        
        # start_time = datetime.datetime.now()
        # data_withing_degree_weighted = self.withinCommunityDegree('weight', communities)
        # print('Time for whiting degree weighted: ' + str(datetime.datetime.now() - start_time))

        # self.save_attributed_graph()
        # print('Graph saved with whiting degree and participation coefficient')
        # print('*******************************************************')
        
        # community_number = 0

        # for k, v in data_whiting_degree.items():

        #     for j in range(len(v)):
        #         for i in range(len(communities[j])):
        #             if k in communities[j][i]:
        #                 community_number = i
        #                 break
                                
        #         if 'data' not in self.G.nodes[k].keys():
        #             self.G.nodes[k]['data'] = {(algorithm, str(j), str(community_number)): {'whiting_degree': v[j]}} # type: ignore        
        #         else:
        #             if (algorithm, str(j), str(community_number)) not in self.G.nodes[k]['data'].keys():
        #                 self.G.nodes[k]['data'][(algorithm, str(j), str(community_number))] = {'whiting_degree': v[j]}
        #             else:
        #                 self.G.nodes[k]['data'][(algorithm, str(j), str(community_number))]['whiting_degree'] = v[j]

        # for k, v in data_withing_degree_weighted.items():

        #     for j in range(len(v)):
        #         for i in range(len(communities[j])):
        #             if k in communities[j][i]:
        #                 community_number = i
        #                 break
        #         if 'data' not in self.G.nodes[k].keys():
        #             self.G.nodes[k]['data'] = {(algorithm, str(j), str(community_number)): {'whiting_degree_weighted': v[j]}} # type: ignore        
        #         else:
        #             if (algorithm, str(j), str(community_number)) not in self.G.nodes[k]['data'].keys():
        #                 self.G.nodes[k]['data'][(algorithm, str(j), str(community_number))] = {'whiting_degree_weighted': v[j]}
        #             else:
        #                 self.G.nodes[k]['data'][(algorithm, str(j), str(community_number))]['whiting_degree_weighted'] = v[j]

        # end Whiting Region
                    
                        

    def apply_measures_to_communities(self, communities: list):                   
        pass
        # conductance_dict = {}

        #     for ci in range(len(iter_communities)):

        #         for cj in range(ci + 1, len(iter_communities)):

        #             community_i = iter_communities[ci]
        #             community_j = iter_communities[cj]

        #             data_conductance = nx.conductance(self.G, community_i, community_j, weight='weight')
                
        #             conductance_dict[(ci, cj)] = data_conductance
        
    def save_attributed_graph_to_csv(self, path: str = 'dataset/outputs/attributed_graph.csv'):

        columns = ['id', 'degree', 'in_degree', 'out_degree', 'weight', 'in_weight', 'out_weight']

        columns.extend(self.G.nodes['104198-F-000000'].keys())

        columns.remove('data')

        keys = self.G.nodes['104198-F-000000']['data'].keys()

        for alg, iter, ci in keys:
            columns.append(alg + '_' + iter + '_' + 'ci')

        for alg, iter, ci in keys:
            columns.append(alg + '_' + iter + '_' + 'ci' + '_participation_coefficient')
            columns.append(alg + '_' + iter + '_' + 'ci' + '_whiting_directed_weighted')
            columns.append(alg + '_' + iter + '_' + 'ci' + '_whiting_directed_notweighted')
            columns.append(alg + '_' + iter + '_' + 'ci' + '_whiting_notdirected_weighted')
            columns.append(alg + '_' + iter + '_' + 'ci' + '_whiting_notdirected_notweighted')
            

        df = pd.DataFrame(columns=columns)

        for node in self.G.nodes():

            data  = {'id': node,  'degree': self.G.degree(node), 'in_degree': self.G.in_degree(node), 'out_degree': self.G.out_degree(node),  # type: ignore
                            'weight': self.G.degree(node, weight='weight'), 'in_weight': self.G.in_degree(node, weight='weight'), 'out_weight': self.G.out_degree(node, weight='weight'),  # type: ignore
                            'eigenvector_centrality': self.G.nodes[node]['eigenvector_centrality'], 'eigenvector_centrality_weighted': self.G.nodes[node]['eigenvector_centrality_weighted'],
                            'pagerank': self.G.nodes[node]['pagerank'], 'degree_centrality': self.G.nodes[node]['degree_centrality'], 'core_number': self.G.nodes[node]['core_number'],
                            'closeness_centrality': self.G.nodes[node]['closeness_centrality'], 'clustering_coefficient': self.G.nodes[node]['clustering_coefficient'], 
                            'vertex_betweenes': self.G.nodes[node]['vertex_betweenes']}
            
            for alg, iter, ci in self.G.nodes[node]['data'].keys():
                data[alg + '_' + iter + '_' + 'ci'] = ci
                data[alg + '_' + iter + '_' + 'ci' + '_participation_coefficient'] = self.G.nodes[node]['data'][(alg, iter, ci)]['participation_coefficient']
                data[alg + '_' + iter + '_' + 'ci' + '_whiting_directed_weighted'] = self.G.nodes[node]['data'][(alg, iter, ci)]['withing_directed_Weighted']
                data[alg + '_' + iter + '_' + 'ci' + '_whiting_directed_notweighted'] = self.G.nodes[node]['data'][(alg, iter, ci)]['withing_directed_notWeighted']
                data[alg + '_' + iter + '_' + 'ci' + '_whiting_notdirected_weighted'] = self.G.nodes[node]['data'][(alg, iter, ci)]['withing_notDirected_Weighted']
                data[alg + '_' + iter + '_' + 'ci' + '_whiting_notdirected_notweighted'] = self.G.nodes[node]['data'][(alg, iter, ci)]['withing_notDirected_notWeighted']

            

            df.loc[len(df)] = data # type: ignore

        df.to_csv(path, index=False)

        

def save_all_communities_tocsv(algorithm: str, communities: list):

    '''
    This function is for save all the communities generated by one algorithm to a csv file.

    Parameters
    ----------
    algorithm : str
        The name of the algorithm that was used to save the communities.
    '''

    df = pd.read_csv('dataset/outputs/all.csv')
    
    df.set_index('id', inplace=True)

    if algorithm != 'infomap':
    

        for i in range(len(communities)):

            for z in range(len(communities[i])):
                
                for _,value in enumerate(communities[i][z]):
                    df.loc[value, i] = z

        for i in range(len(communities)):
            df[i] = df[i].astype('Int64')
    else:

        for i in range(len(communities)):

            for z in range(len(communities[i].communities)):
                for _,value in enumerate(communities[i].communities[z]):
                    df.loc[value, i] = z

        for i in range(len(communities)):
            df[i] = df[i].astype('Int64')

    df.to_csv('dataset/outputs/all_' + algorithm + '.csv')
    
                
def writter(lis, name):

    with open('./dataset/outputs/' + name, 'w') as f:
        for (id, value) in lis:
            f.write(f'{id}, {value}\n')  

def run_and_save_algorithm(m: Matrix, algorithm, params, n, seed = []) :

    if seed == []:
        seed = [x for x in range(n)]

    if algorithm == 'louvain':

        communities = m.lovain_concurrent(seed= seed,  n=10)

        for com in communities:
            print(m.communities_length(com))

        m.save_communities(communities, 'louvain', params=params, seed= seed)
    
    elif algorithm == 'greedy':

        communities = m.greedy_modularity_concurrent(resolution=2 , n=n)

        for com in communities:
            print(m.communities_length(com))

        m.save_communities(communities, 'greedy', params=params )
    
    elif algorithm == 'lpa':

        communities = m.asyn_lpa_concurrent(weight = 'weight', seed = seed , n = n)

        for com in communities:
            print(m.communities_length(com))

        m.save_communities(communities, 'lpa', params=params, seed = seed )

    elif algorithm == 'infomap':

        communities = m.infomap_concurrent(seed=seed, n = n)

        for com in communities:
            print(m.communities_length(com.communities))

        m.save_communities(communities, 'infomap', params=params, seed = seed )
      
def small(communities: list):
    
    com = []
    data = {}

    for comm, path in communities:

        data['communities'] = comm.communities
        data['method_name'] = comm.method_name
        data['method_parameters'] = comm.method_parameters
        data['overlap'] = comm.overlap

        com.append((data, path))
        data = {}

    return com

def apply_central_dispersion_measures(path: str):

    data = pd.read_csv('output/attributed_graph-1.4.1.csv', header=0)

    #data = data[[' weight', ' edge_betweenness']]

    data = data[['degree', 'in_degree', 'out_degree', 'weight', 'in_weight', 'out_weight', 'eigenvector_centrality', 'eigenvector_centrality_weighted', 'pagerank', 'degree_centrality', 'core_number', 'closeness_centrality', 'clustering_coefficient', 'vertex_betweenes']]

    #Calculate summary statistics
    summary = data.describe()
    

    # Add mean, median, mode, variance, standard deviation, skewness, and kurtosis
    summary.loc['mean'] = data.mean()
    summary.loc['median'] = data.median()
    summary.loc['mode'] = data.mode().iloc[0]
    summary.loc['variance'] = data.var()    
    summary.loc['skewness'] = data.skew()
    summary.loc['kurtosis'] = data.kurtosis()

    sw_test_results = pd.DataFrame({'shapiro_w': []})
    for col in data.columns:
        sw, p = stats.shapiro(data[col])
        sw_test_results.loc[col] = [p]

   
    # Add the Shapiro-Wilk test results to the original DataFrame
    summary.loc['shapiro_w'] = sw_test_results['shapiro_w']

    # Export summary statistics to a CSV file
    summary.to_csv(path)

def nmi_overlapping_evaluate(foldername: str) -> None:

    '''
    Evaluate the overlapping detection methods using NMI

    Parameters
    ----------
    foldername: str
        Path to folder containing  GT communities and detected communities

    Returns
    -------
    None
    '''

    from cdlib import evaluation, NodeClustering

    files = os.listdir('dataset/' + foldername)
    files.remove('GT')
    files.remove('README.txt')

    for file in files:
        nodes = []
        match = re.search(r'(network)(\d+)', file)
        number = '0'
        if match:
            number = str(match.group(2))
            with open('dataset/' + foldername + '/GT/community' + number + '_GT.dat', 'r') as f:
                lines = f.readlines()        
                for line in lines:
                    data = line.split(' ')
                    inter_data = [int(x) for x in data]
                    nodes.append(inter_data)
        
        G = pickle.load(open('dataset/' + foldername + '/' + file + '/' + file + '.pkl', 'rb')) 
    

        nodeClustA = NodeClustering(communities=nodes, graph=G, method_name='GT', method_parameters={}, overlap=True)
    
        nodes = []
        outputs = ['', '_Lpa', '_Louvain', '_Greedy', '_Infomap']
        
        with open('output/' + foldername + '/' + foldername + '_result.txt', 'a') as f:
            f.write('network' + number + '\n')

        for output in outputs:
            with open('output/' + foldername + '/network' + number + output + '.txt', 'r') as f:
                lines = f.readlines()        
                for line in lines:
                    line = line.strip('\n').rstrip()
                    data = line.split(' ')
                    inter_data = [int(x) for x in data]
                    nodes.append(inter_data)

            nodeClustB = NodeClustering(communities=nodes, graph=G, method_name=output, method_parameters={}, overlap=True)
            
            match_resoult = evaluation.overlapping_normalized_mutual_information_MGH(nodeClustA, nodeClustB)

            with open('output/' + foldername + '/' + foldername + '_result.txt', 'a') as f:
                f.write(output + ': ' + str(match_resoult.score) + '\n')

            nodes = []
        with open('output/' + foldername + '/' + foldername + '_result.txt', 'a') as f:
            f.write('------------------------\n')

def nmi_overlapping_evaluateTunning(foldername: str) -> None:
    '''
    Evaluate the overlapping detection methods using NMI

    Parameters
    ----------
    foldername: str
        Path to folder containing  GT communities and detected communities

    Returns
    -------
    None
    '''
    from cdlib import evaluation, NodeClustering

    files = os.listdir('dataset/' + foldername)
    files.remove('GT')
    files.remove('README.txt')

    files = sorted(files, key=lambda x: int("".join([i for i in x if i.isdigit()])))

    dictResult = dict()

    for file in files:
        if '.pkl' in file:
            continue
        nodes = []
        match = re.search(r'(network)(\d+)', file)
        number = '0'
        if match:
            number = str(match.group(2))
            with open('dataset/' + foldername + '/GT/community' + number + '_GT.dat', 'r') as f:
                lines = f.readlines()        
                for line in lines:
                    data = line.split(' ')
                    inter_data = [int(x) for x in data]
                    nodes.append(inter_data)
        
        G = pickle.load(open('dataset/' + foldername + '/' + file + '/' + file + '.pkl', 'rb')) 
    
        # GT created
        nodeClustA = NodeClustering(communities=nodes, graph=G, method_name='GT', method_parameters={}, overlap=True)
    
        nodes = []

        with open('output/' + foldername + '/' + foldername + '_result.txt', 'a') as f:
            f.write('network' + number + '\n')
        
        # read files result
        filesResultAlg = os.listdir('output/' + foldername)

        # remove .txt, .pkl
        if os.path.exists('output/' + foldername + '/' + foldername + '_result.txt'):
            filesResultAlg.remove(foldername + '_result.txt')

        # if os.path.exists('output/' + foldername + '/' + foldername + '_result.pkl'):
        #     filesResultAlg.remove(foldername + '_result.pkl')

        # sorted files
        filesResultAlg = sorted(filesResultAlg, key=lambda x: int("".join([i for i in x if i.isdigit()])))
            
        for file_i in filesResultAlg:
            if '.pkl' not in file_i and file == file_i.split('_')[0]:
                with open(f'output/{foldername}/{file_i}', 'r') as f:
                    lines = f.readlines()        
                    for line in lines:
                        line = line.strip('\n').rstrip()
                        data = line.split(' ')
                        inter_data = [int(x) for x in data]
                        nodes.append(inter_data)

                    # community created
                    nodeClustB = NodeClustering(communities=nodes, graph=G, method_name=file_i, method_parameters={}, overlap=True)
                    
                    # evaluate GT vs community
                    match_resoult = evaluation.overlapping_normalized_mutual_information_MGH(nodeClustA, nodeClustB)

                    algName = file_i.split('_')[1].removesuffix('.txt')
                    fileNameMod = algName

                    with open('output/' + foldername + '/' + foldername + '_result.txt', 'a') as f:
                        f.write(fileNameMod + ': ' + str(match_resoult.score) + '\n')
                    
                    
                    if not fileNameMod in dictResult.keys():
                        dictResult[fileNameMod] = {'Algorithms/Parameters': fileNameMod, file: match_resoult.score}
                    else:
                        dictResult[fileNameMod][file] = match_resoult.score
                        
                    nodes = []
        with open('output/' + foldername + '/' + foldername + '_result.txt', 'a') as f:
            f.write('------------------------\n')

        pickle.dump(dictResult, open('output/' + foldername + '/' + foldername + '_result.pkl', 'wb'))
        

def nmi_overlapping_evaluateTunning(foldername: str) -> None:
    '''
    Evaluate the overlapping detection methods using NMI

    Parameters
    ----------
    foldername: str
        Path to folder containing  GT communities and detected communities

    Returns
    -------
    None
    '''
    from cdlib import evaluation, NodeClustering

    files = os.listdir('dataset/' + foldername)
    files.remove('GT')
    files.remove('README.txt')

    files = sorted(files, key=lambda x: int("".join([i for i in x if i.isdigit()])))

    dictResult = dict()

    for file in files:
        if '.pkl' in file:
            continue
        nodes = []
        match = re.search(r'(network)(\d+)', file)
        number = '0'
        if match:
            number = str(match.group(2))
            with open('dataset/' + foldername + '/GT/community' + number + '_GT.dat', 'r') as f:
                lines = f.readlines()        
                for line in lines:
                    data = line.split(' ')
                    inter_data = [int(x) for x in data]
                    nodes.append(inter_data)
        
        G = pickle.load(open('dataset/' + foldername + '/' + file + '/' + file + '.pkl', 'rb')) 
    
        # GT created
        nodeClustA = NodeClustering(communities=nodes, graph=G, method_name='GT', method_parameters={}, overlap=True)
    
        nodes = []

        with open('output/' + foldername + '/' + foldername + '_result.txt', 'a') as f:
            f.write('network' + number + '\n')
        
        # read files result
        filesResultAlg = os.listdir('output/' + foldername)

        # remove .txt, .pkl
        if os.path.exists('output/' + foldername + '/' + foldername + '_result.txt'):
            filesResultAlg.remove(foldername + '_result.txt')

        # if os.path.exists('output/' + foldername + '/' + foldername + '_result.pkl'):
        #     filesResultAlg.remove(foldername + '_result.pkl')

        # sorted files
        filesResultAlg = sorted(filesResultAlg, key=lambda x: int("".join([i for i in x if i.isdigit()])))
            
        for file_i in filesResultAlg:
            if '.pkl' not in file_i and file == file_i.split('_')[0]:
                with open(f'output/{foldername}/{file_i}', 'r') as f:
                    lines = f.readlines()        
                    for line in lines:
                        line = line.strip('\n').rstrip()
                        data = line.split(' ')
                        inter_data = [int(x) for x in data]
                        nodes.append(inter_data)

                    # community created
                    nodeClustB = NodeClustering(communities=nodes, graph=G, method_name=file_i, method_parameters={}, overlap=True)
                    
                    # evaluate GT vs community
                    match_resoult = evaluation.overlapping_normalized_mutual_information_MGH(nodeClustA, nodeClustB)

                    algName = file_i.split('_')[1].removesuffix('.txt')
                    fileNameMod = algName

                    with open('output/' + foldername + '/' + foldername + '_result.txt', 'a') as f:
                        f.write(fileNameMod + ': ' + str(match_resoult.score) + '\n')
                    
                    
                    if not fileNameMod in dictResult.keys():
                        dictResult[fileNameMod] = {'Algorithms/Parameters': fileNameMod, file: match_resoult.score}
                    else:
                        dictResult[fileNameMod][file] = match_resoult.score
                        
                    nodes = []
        with open('output/' + foldername + '/' + foldername + '_result.txt', 'a') as f:
            f.write('------------------------\n')

        pickle.dump(dictResult, open('output/' + foldername + '/' + foldername + '_result.pkl', 'wb'))
        

def runRoughClustering(m : Matrix, folder_version = 'NetsType_1.1'):

    all_iterations = []
    
    if os.path.exists(f'dataset/{folder_version}'):
        list_dir = os.listdir(f'dataset/{folder_version}')

        list_dir.remove('README.txt')
        list_dir.remove('GT')

        for net in list_dir:
        
            m.G = pickle.load(open(f'dataset/{folder_version}/{net}/{net}.pkl', 'rb'))

            n = 0
            top = 10
            
            print(f'async_lpa Algorithm running ' + str(top) + f' times in {net}')
            for _ in range(n, top):
                result = nx.algorithms.community.label_propagation.asyn_lpa_communities(m.G, seed=random.randint(0, 10000))
                communities = [list(x) for x in result]
                if len(communities) > 1:
                    all_iterations.append(communities) # type: ignore
                #print(all_iterations[-1])
            print('async_lpa Algorithm finished')
            
            # Range of Resolution 3.5 - 5.5
            print(f'Greedy Algorithm running {str(1)} times in {net}')
            for _ in range(0, 1):
                result = nx.algorithms.community.greedy_modularity_communities(m.G, resolution= random.uniform(3.5, 5.5))  # type: ignore
                result = [list(x) for x in result] # type: ignore
                
                for _ in range(0, int(top/1.5)):        
                    all_iterations.append(result) 
                #print(all_iterations[-1])

            print(f'Greedy Algorithm finished')

            # Range of Resolution 2 - 3.5
            print(f'Louvain Algorithm running {str(top)}  times in {net}')
            for _ in range(n, top):
                result = nx.algorithms.community.louvain.louvain_communities(m.G, seed=random.randint(0, 10000), resolution= random.uniform(2, 3.5)) # type: ignore
                #print(result)
                all_iterations.append([list(x) for x in result]) # type: ignore
            print('Louvain Algorithm finished')

            print(f'Infomap Algorithm loading {str(top)} times in {net}')
            infomap_results = pickle.load(open(f'output/{folder_version}/{net}_Infomap.pkl', 'rb'))
            all_iterations.extend(infomap_results) # type: ignore
            print('Infomap Algorithm finished')

            value = m.RoughClustering(communities=all_iterations)

            all_iterations = []


            exportpath_RC = f'/{net}_RC.txt'

            m.export_RC(folder_version, exportpath_RC, value)

def runRoughClustering_on_FlyCircuit(m: Matrix ,version_dataset: str,  iterations: list):



    f'''
    Run RoughClustering algorithm on FlyCircuit dataset

    Parameters
    ----------
    m: Matrix
        Matrix object with graph and tools
    iterations: list
        List of iterations of all algorithms used to commutiies detection
    version_dataset: str
        Version of dataset

    Output
    ------
    The result of RoughClustering algorithm is stored in output/FlyCircuit/FlyCircuit_{version_dataset}_RC.txt
            
    '''

    if not os.path.exists('output/FlyCircuit/'):
        os.makedirs('output/FlyCircuit/')

    consensus = m.RoughClustering(communities=iterations)

    m.export_RC('FlyCircuit', f'/FlyCircuit_{version_dataset}_RC.txt', consensus)

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

def runAlgorithmSimple(m: Matrix, folder_version = 'NetsType_1.3'):

    for j in range(1, 12):

        m.G = pickle.load(open('dataset/' + folder_version + '/network'+ str(j) + '/network'+ str(j) + '.pkl', 'rb'))

        n = 0
        top = 1

        exportpath_Simple = folder_version

        for _ in range(n, top):
            result = nx.algorithms.community.label_propagation.asyn_lpa_communities(m.G, seed=random.randint(0, 10000))
            communities = [list(x) for x in result]
            m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Lpa.txt', communities)
        
        for _ in range(n, int(top)):
            result = nx.algorithms.community.greedy_modularity_communities(m.G, resolution= random.uniform(3.5, 5.5))         # type: ignore
            communities = [list(x) for x in result] # type: ignore
            m.export_Simple(exportpath_Simple, '/network'+ str(j) +'_Greedy.txt', communities)
        
        for _ in range(n, top):
            result = nx.algorithms.community.louvain.louvain_communities(m.G, seed=random.randint(0, 10000), resolution= random.uniform(2, 3.5)) # type: ignore
            communities = [list(x) for x in result] # type: ignore
            m.export_Simple(exportpath_Simple, '/network'+ str(j) + '_Louvain.txt', communities)

    print('done')

def plot_degree_distribution(m: Matrix):

    # Calculate the degree distribution
    degree_sequence = sorted([d for _, d in m.G.degree(weight='weight')], reverse=True) # type: ignore
    
    degree_count = {}
    for degree in degree_sequence:
        if degree in degree_count:
            degree_count[degree] += 1            
        else:
            degree_count[degree] = 1
    deg, cnt = zip(*degree_count.items())

    

    # Plot the degree distribution
    with plt.style.context('tableau-colorblind10'):
        #plt.style.use('seaborn-v0_8-dark')
        #plt.bar(deg, cnt, width=0.80)
        plt.hist(degree_sequence, bins=100, alpha=0.5)
        plt.hist([1], bins=1, alpha=0.5, color='red')
        #plt.gca().invert_xaxis()
        plt.title("Weight Distribution")
        plt.ylabel("Count")
        plt.xlabel("Degree")

        data = []

        for d, c in zip(deg, cnt):
            for _ in range(c):
                data.append(d)
                
        # Add mean line
        mean = np.average(data)
    
        plt.axvline(mean, color='r', linestyle='--', label=f'Mean: {mean:.2f}') # type: ignore
        #plt.text(mean, max(cnt), f'Mean: {mean:.2f}', ha='left', va='top', fontsize=12, color='r')
        
        
        total = sum([ d * c for d, c in zip(deg, cnt)])
        portion = total * 0.8
        
        index = 0
        
        count = 0
        increment = 0
        for d, c in zip(deg, cnt):
            increment += (d * c)
            count += c
            if increment >= portion:
                plt.axvline(d, color='g', linestyle='--', label=f'80%: {d}')
                #plt.text(v, max(cnt), f'80%: {v}', ha='left', va='top', fontsize=12, color='g')
                index = deg.index(d)
                break

        
        print(count)
        print(index)
        plt.legend()
        plt.show()

def save_comunities_summary():

    df = pd.read_csv('output/summary_communities.csv', sep=',', header=0)

    greedycolums = [f'greedy_{i}_ci' for i in range(0, 8)]
    louvaincolums = [f'louvain_{i}_ci' for i in range(0, 10)]
    lpacolums = [f'lpa_{i}_ci' for i in range(0, 9)]
    infomapcolums = [f'infomap_{i}_ci' for i in range(0, 10)]

    gp = [f'greedy_{i}_ci_participation_coefficient' for i in range(0, 8)]
    gwdw = [f'greedy_{i}_ci_whiting_directed_weighted' for i in range(0, 8)]
    gwdnw = [f'greedy_{i}_ci_whiting_directed_notweighted' for i in range(0, 8)]
    gndw = [f'greedy_{i}_ci_whiting_notdirected_weighted' for i in range(0, 8)]
    gndnw = [f'greedy_{i}_ci_whiting_notdirected_notweighted' for i in range(0, 8)]
    
    lpap = [f'lpa_{i}_ci_participation_coefficient' for i in range(0, 9)]
    lpadw = [f'lpa_{i}_ci_whiting_directed_weighted' for i in range(0, 9)]
    lpadnw = [f'lpa_{i}_ci_whiting_directed_notweighted' for i in range(0, 9)]
    lpandw = [f'lpa_{i}_ci_whiting_notdirected_weighted' for i in range(0, 9)]
    lpandnw = [f'lpa_{i}_ci_whiting_notdirected_notweighted' for i in range(0, 9)]

    infomap_p = [f'infomap_{i}_ci_participation_coefficient' for i in range(0, 10)]
    infomap_dw = [f'infomap_{i}_ci_whiting_directed_weighted' for i in range(0, 10)]
    infomap_dnw = [f'infomap_{i}_ci_whiting_directed_notweighted' for i in range(0, 10)]
    infomap_ndw = [f'infomap_{i}_ci_whiting_notdirected_weighted' for i in range(0, 10)]
    infomap_ndnw = [f'infomap_{i}_ci_whiting_notdirected_notweighted' for i in range(0, 10)]
    
    lp = [f'louvain_{i}_ci_participation_coefficient' for i in range(0, 10)]
    ldw = [f'louvain_{i}_ci_whiting_directed_weighted' for i in range(0, 10)]
    ldnw = [f'louvain_{i}_ci_whiting_directed_notweighted' for i in range(0, 10)]
    lndw = [f'louvain_{i}_ci_whiting_notdirected_weighted' for i in range(0, 10)]
    lndnw = [f'louvain_{i}_ci_whiting_notdirected_notweighted' for i in range(0, 10)]

    columns_set = [gp, gwdw, gwdnw, gndw, gndnw, lpap, lpadw, lpadnw, lpandw, lpandnw, infomap_p, infomap_dw, infomap_dnw, infomap_ndw, infomap_ndnw, lp, ldw, ldnw, lndw, lndnw]
    columns_set_name = ['gp', 'gwdw', 'gwdnw', 'gndw', 'gndnw', 'lpap', 'lpadw', 'lpadnw', 'lpandw', 'lpandnw', 'infomap_p', 'infomap_dw', 'infomap_dnw', 'infomap_ndw', 'infomap_ndnw', 'lp', 'ldw', 'ldnw', 'lndw', 'lndnw']
    for i in range(len(columns_set)):
        df = df.assign(**{columns_set_name[i]: df[columns_set[i]].mean(axis=1)})

    
    
    result = pd.concat([df.iloc[:, 0:1], df.iloc[:, -20:]], axis=1)

    result.to_csv('output/summary_communities_measures.csv', index=False)

def evaluate_overlaping(net_path : str):

    folder_path_gt = f'dataset/{net_path}/GT/'
    folder_path_alg = f'output/{net_path}/'

    list_files = os.listdir(folder_path_gt)
    list_files.sort(key=lambda f: int(re.sub(r'\D', '', f)))
    data_final = []
    
    

    for file in list_files:
        nodes_gt = []
        nodes_alg = []
        match = re.search(r'(community)(\d+)(_GT.dat)', file)
        number = '0'
        if match:
            number = str(match.group(2))
            with open(f'{folder_path_gt}community{number}_GT.dat', 'r') as f:
                lines = f.readlines()        
                for line in lines:
                    data = line.split(' ')
                    inter_data = [int(x) for x in data]
                    nodes_gt.append(inter_data)

            list_nodes_overlaping_gt = detect_nodes_with_overlapping(nodes_gt)

            

            with open(f'{folder_path_alg}network{number}_RC.txt', 'r') as f:
                lines = f.readlines()        
                for line in lines:
                    data = line.split(' ')
                    data.remove('\n')
                    inter_data = [int(x) for x in data]
                    nodes_alg.append(inter_data)

            list_nodes_overlaping_rc = detect_nodes_with_overlapping(nodes_alg)


            result = compare_overlapings(list_nodes_overlaping_gt, list_nodes_overlaping_rc)
            
            data_final.append(result)
    
    
    return data_final


def analyze_overlaping(net_type : str):

    gt_files = os.listdir(f'dataset/{net_type}/GT/')
    rc_files = os.listdir(f'output/{net_type}/')

    rc_files = [file for file in rc_files if file.endswith('_RC.txt')]

    gt_files.sort(key=lambda f: int(re.sub(r'\D', '', f))) 
    rc_files.sort(key=lambda f: int(re.sub(r'\D', '', f)))

    

    for i in range(1, 12):

        nodes_gt_overlaping = detect_nodes_with_overlapping(read_communities_from_dat(f'dataset/{net_type}/GT/community{i}_GT.dat'))
        nodes_rc_overlaping = detect_nodes_with_overlapping(read_communities_from_dat(f'output/{net_type}/network{i}_RC.txt'))
        
        nodes_gt_overlaping = dict(sorted(nodes_gt_overlaping.items(), key=lambda x: x[0]))

        nodes_rc_overlaping = dict(sorted(nodes_rc_overlaping.items(), key=lambda x: x[0]))

        

        

        match = set(nodes_gt_overlaping.keys()).intersection(set(nodes_rc_overlaping.keys()))    

        
        #whiting_gt = pickle.load(open(f'output/{net_type}/within_GT_network{i}.pkl', 'rb'))
        #whithing_rc = pickle.load(open(f'output/{net_type}/within_RC_network{i}.pkl', 'rb'))

        pc_gt = pickle.load(open(f'dataset/{net_type}/network{i}/network{i}_GT_PC.pkl', 'rb'))
        pc_rc = pickle.load(open(f'output/{net_type}/network{i}_RC_PC.pkl', 'rb'))

        #nodes_whiting_gt = {key: whiting_gt[key] for key in nodes_gt_overlaping.keys()}
        #nodes_whiting_rc = {key: whithing_rc[key] for key in nodes_rc_overlaping.keys()}
        
        nodes_id = list(range(1000))

        #nodes_whiting_gt_inf = {key: whiting_gt[key] for key in nodes_id if key not in nodes_gt_overlaping.keys()}
        #nodes_whiting_rc_inf = {key: whithing_rc[key] for key in nodes_id if key not in nodes_rc_overlaping.keys()}
        
        

        # # Create a line plot of the first dictionary values
        # plt.scatter(nodes_whiting_gt_inf.keys(), nodes_whiting_gt_inf.values(), label='GT', marker='^') # type: ignore

        # # Create a line plot of the second dictionary values
        # plt.scatter(nodes_whiting_rc_inf.keys(), nodes_whiting_rc_inf.values(), label='RC') # type: ignore

        # # Set the plot title and axis labels
        # plt.title('GT vs RC - Whiting in not overlaping nodes')
        # plt.xlabel('Keys')
        # plt.ylabel('Values')

        # # Add a legend to the plot
        # plt.legend()

      

        # gt_mean = np.mean(list(nodes_whiting_gt_inf.values())) # type: ignore
        # rc_mean = np.mean(list(nodes_whiting_rc_inf.values())) # type: ignore
        
        # plt.axhline(y=gt_mean, color='b', linestyle='--')
        # plt.axhline(y=rc_mean, color='orange')

        # # Show the plot
        # plt.show()
        
        

        nodes_pc_gt = {key: pc_gt[key] for key in nodes_gt_overlaping.keys() if key not in match}
        nodes_pc_rc = {key: pc_rc[key] for key in nodes_rc_overlaping.keys() if key not in match}

        nodes_pc_gt_match = {key: pc_gt[key] for key in nodes_gt_overlaping.keys() if key in match}
        nodes_pc_rc_match = {key: pc_rc[key] for key in nodes_rc_overlaping.keys() if key in match}

        #nodes_pc_gt_inf = {key: pc_gt[key] for key in nodes_id if key not in nodes_gt_overlaping.keys()}
        #nodes_pc_rc_inf = {key: pc_rc[key] for key in nodes_id if key not in nodes_rc_overlaping.keys()}
        
               
       
        # Create a dot plot of the GT overlaping nodes but without match with RC values
        plt.scatter(nodes_pc_gt.keys(), nodes_pc_gt.values(), label='GT', marker='^', c= 'orange') # type: ignore

        # Create a dot plot of the RC overlaping nodes but without match with GT values
        plt.scatter(nodes_pc_rc.keys(), nodes_pc_rc.values(), label='RC', marker='s', c='blue') # type: ignore

        # Create a dot plot of the GT overlaping nodes with only match with RC values
        plt.scatter(nodes_pc_gt_match.keys(), nodes_pc_gt_match.values(), label='GT Match', marker='^', c= 'red') # type: ignore

        # Create a dot plot of the RC overlaping nodes with only match with GT values
        plt.scatter(nodes_pc_rc_match.keys(), nodes_pc_rc_match.values(), label='RC Match', marker='s', c='red') # type: ignore

        #for node in match:
        #    plt.axvline(x=node, ymin=nodes_pc_rc[node] , ymax=nodes_pc_gt[node], color='r')

        # Set the plot title and axis labels
        plt.title('GT vs RC - PC in overlaping nodes')
        plt.xlabel('Nodes')
        plt.ylabel('PC Values')

        

        gt_mean_values = list(nodes_pc_gt.values())
        rc_mean_values = list(nodes_pc_rc.values())
        gt_mean_values.extend(list(nodes_pc_gt_match.values()))
        rc_mean_values.extend(list(nodes_pc_rc_match.values()))

        gt_mean = np.mean(gt_mean_values) # type: ignore
        rc_mean = np.mean(rc_mean_values) # type: ignore
        
        plt.axhline(y=gt_mean, color='b', linestyle='--', label='GT Mean') # type: ignore
        plt.axhline(y=rc_mean, color='orange', label='RC Mean') # type: ignore

        # Add a legend to the plot
        plt.legend()

        # Show the plot
        #plt.show()

        if not os.path.exists(f'output/{net_type}/PC_Overlaping_score.csv'):
            with open(f'output/{net_type}/PC_Overlaping_score.csv', 'a+') as f:
                f.write(f'Net,PC_Mean_GT, PC_Mean_GT,T_Possitive,F_Possitive \n')

        with open(f'output/{net_type}/PC_Overlaping_score.csv', 'a+') as f:            
            f.write(f'Network{i},{gt_mean},{rc_mean},{len(match)},{len(nodes_pc_rc.keys())} \n')

        if not os.path.exists(f'output/{net_type}/img'):
            os.makedirs(f'output/{net_type}/img')
        
        #plt.savefig(f'output/{net_type}/img/PC_network{i}.png', dpi=500)
        plt.close()

def detect_nodes_with_overlapping(node_list: list[list[int]]) -> dict[int, int]:
    from collections import Counter
    nodes_with_overlapping = {}
    node_counts = Counter(node for com in node_list for node in com)
    
    for node, count in node_counts.items():
        if count > 1:
            nodes_with_overlapping[node] = count

    return nodes_with_overlapping

def compare_overlapings(gt : dict, algorithm : dict) -> tuple[list[int], list[int]]:
    true_overlaping = []
    false_overlaping = []

    for node in algorithm.keys():
        if node in gt.keys():
            true_overlaping.append(node)
        else:
            false_overlaping.append(node)

    return (true_overlaping, false_overlaping)

def lpa_wrapper(G, seed = 1):

        import networkx.algorithms.community as nx_comm
        return list(nx_comm.asyn_lpa_communities(G, seed=seed)) # type: ignore
    

def stability(sequence : int, num_run : int, net_path : str):

    folder_path_gt = f'dataset/{net_path}'

    files = os.listdir(folder_path_gt)
    files.remove('GT')
    files.remove('README.txt')

    for file in files:
        G = pickle.load(open(f'{folder_path_gt}/{file}/{file}.pkl', 'rb'))        

        for seq in range(sequence):
            

            # print(f'async_lpa Algorithm running ' + str(seq) + f' times in {file}')
            
            # with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            #     communities = pool.starmap(lpa_wrapper, [(G, i) for i in range(num_run)])
            #     communities = [[list(x) for x in com] for com in communities]                
            #     os.makedirs(f'output/stability/{net_path}/{file}/', exist_ok=True)
            #     pickle.dump(communities, open(f'output/stability/{net_path}/{file}/async_lpa_{num_run}_run_{seq}.pkl', 'wb'))
            # print('async_lpa Algorithm finished')


            # print(f'louvain Algorithm running ' + str(seq) + f' times in {file}')

            # with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            #     communities = pool.starmap(nx_comm.louvain_communities, [(G, 'weight', random.uniform(2, 3.5), 1e-07, i) for i in range(num_run)])
            #     communities = [[list(x) for x in com] for com in communities]
            #     pickle.dump(communities, open(f'output/stability/{net_path}/{file}/louvain_{num_run}_run_{seq}.pkl', 'wb'))
            # print('louvain Algorithm finished')


            print(f'greedy Algorithm running ' + str(seq) + f' times in {file}')

            with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
                communities = pool.starmap(nx_comm.greedy_modularity_communities, [(G, None ,random.uniform(3.5, 5.5), 1, None) for _ in range(int(num_run/1.5))])
                communities = [[list(x) for x in com] for com in communities]
                os.makedirs(f'output/stability/{net_path}/{file}/', exist_ok=True)
                pickle.dump(communities, open(f'output/stability/{net_path}/{file}/greedy_{num_run}_run_{seq}.pkl', 'wb'))
                            
            print(f'Greedy Algorithm finished')



def stability_infomap(sequence : int, num_run : int, net_path : str):
    from cdlib import algorithms # type: ignore Only necesary with cdlib environment #TODO

    folder_path_gt = f'dataset/{net_path}'

    files = os.listdir(folder_path_gt)
    files.remove('GT')
    files.remove('README.txt')

    for file in files:
        if file == 'network1' or file == 'network2':
            continue
        G = pickle.load(open(f'{folder_path_gt}/{file}/{file}.pkl', 'rb'))        

        for seq in range(sequence):

            print(f'infomap Algorithm running ' + str(num_run) + f' times in {file}')
            
            communities = []
            
            with concurrent.futures.ProcessPoolExecutor() as executor:
                futures = [executor.submit(algorithms.infomap, G, flags='--seed ' + str(i+1)) for i in range(num_run)]
                for future in concurrent.futures.as_completed(futures):
                    communities.append(future.result().communities)
            
            pickle.dump(communities, open(f'output/stability/{net_path}/{file}/infomap_{num_run}_run_{seq}.pkl', 'wb'))
            
            print('Infomap Algorithm finished')

def run_RC_sequences(sequence : int, folder_version: str, r: int):
    
    m = Matrix([], {},[])

    folder_path = f'output/stability/{folder_version}'
    
    max_iter = 1000 if folder_version == 'NetType_1.6' else 100

    folder_list = os.listdir(folder_path)
    for net in folder_list:
        
        m.G = pickle.load(open(f'dataset/{folder_version}/{net}/{net}.pkl', 'rb'))

        for i in range(sequence):

            all_communities = []
            
            greedy_communities = pickle.load(open(f'{folder_path}/{net}/greedy_{r}_run_{0}.pkl', 'rb'))
            all_communities.extend(greedy_communities)

            louvain_communities = pickle.load(open(f'{folder_path}/{net}/louvain_{r}_run_{i}.pkl', 'rb'))
            louvain_communities = [list(x) for x in louvain_communities] # type: ignore
            all_communities.extend(louvain_communities)

            async_lpa_communities = pickle.load(open(f'{folder_path}/{net}/async_lpa_{r}_run_{i}.pkl', 'rb'))
            async_lpa_communities = [list(x) for x in async_lpa_communities] # type: ignore
            all_communities.extend(async_lpa_communities)

            infomap_communities = pickle.load(open(f'{folder_path}/{net}/infomap_{r}_run_{i}.pkl', 'rb'))
            all_communities.extend(infomap_communities)

                    
            value = m.RoughClustering(communities=all_communities)

            all_communities = []


            exportpath_RC = f'{net}_RC.txt'

            m.export_RC(f'stability/{folder_version}/{net}/', exportpath_RC, value)

def read_communities_from_dat(path : str) -> list[list[int]]:
    
    communities = []
    with open(path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            data = line.split(' ')
            if data[-1] == '\n':
                data.remove('\n')
            inter_data = [int(x) for x in data]
            communities.append(inter_data)
    return communities

def apply_PC_to_GT(net_version : str):

    folder_path = f'dataset/{net_version}'
    
    
    folder_list = os.listdir(folder_path)
    folder_list.remove('GT')
    folder_list.remove('README.txt')
    m = Matrix([], {},[])

    for net in folder_list:
        
        m.G = pickle.load(open(f'dataset/{net_version}/{net}/{net}.pkl', 'rb'))

        communities = []
        match = re.search(r'(network)(\d+)', net)
        number = '0'
        if match:
            number = str(match.group(2))
            with open(f'dataset/{net_version}/GT/community{number}_GT.dat', 'r') as f:
                lines = f.readlines()        
                for line in lines:
                    data = line.split(' ')
                    inter_data = [int(x) for x in data]
                    communities.append(inter_data)

        dict_node = m.participation_coefficient(communities)

        pickle.dump(dict_node, open(f'dataset/{net_version}/{net}/{net}_GT_PC.pkl', 'wb'))

def apply_PC_to_RC(net_version: str):
    
    folder_path = f'output/{net_version}'

    files_list = os.listdir(folder_path)
    m = Matrix([], {},[])

    for file in files_list:
        if file.endswith('_RC.txt'):
            m.G = pickle.load(open(f'dataset/{net_version}/{file[:-7]}/{file[:-7]}.pkl', 'rb'))
            communities = read_communities_from_dat(f'{folder_path}/{file}')
            dict_node = m.participation_coefficient(communities)
            pickle.dump(dict_node, open(f'{folder_path}/{file[:-7]}_RC_PC.pkl', 'wb'))

def increse_greedy_files(net_version: str):

    folder_path = f'output/stability/{net_version}'

    folder_list = os.listdir(folder_path)

    for folder in folder_list:

        for i in range(4):

            communities = pickle.load(open(f'{folder_path}/{folder}/greedy_10_run_{i}.pkl', 'rb'))
            for j in range(i*5, (i+1)*5):
                file_path = f'{folder_path}/{folder}/greedy_100_run_{j}.pkl'
                if not os.path.exists(file_path):
                    pickle.dump(communities, open(file_path, 'wb'))

def evaluate_stability(net_version: str, num_iter: int):

    algorithms_names = ['louvain', 'infomap', 'greedy', 'async_lpa']

    folder_path = f'output/stability/{net_version}'

    folder_list = os.listdir(folder_path)

    iter_list = [10, 100, 1000] if net_version == 'NetsType_1.4' else [10, 50, 100]
    
    m = Matrix([], {},[])
    
    for folder in folder_list:

        
        for i in range(20):

            all_communities = []

            for algorithm in algorithms_names:
                communities = pickle.load(open(f'{folder_path}/{folder}/{algorithm}_{num_iter}_run_{i}.pkl', 'rb'))
                all_communities.extend(communities)

            
            m.G = pickle.load(open(f'dataset/{net_version}/{folder}/{folder}.pkl', 'rb'))


            for iter in iter_list:
                comunities_subset = []
                for j in range(4):
                    index = j * num_iter
                    comunities_subset.extend(all_communities[index:(index + iter)])
                     
                
                value = m.RoughClustering(communities=comunities_subset)
                m.export_RC(f'stability/{net_version}/{folder}/', f'{folder}_RC_{iter}_run_{i}.txt', value)
            
def process_folder(net_version: str, num_iter: int, folder: str, algorithms_names: list[str], folder_path: str, iter_list: list[int]):

    m = Matrix([], {},[])

    for i in range(20):

            all_communities = []

            for algorithm in algorithms_names:
                communities = pickle.load(open(f'{folder_path}/{folder}/{algorithm}_{num_iter}_run_{i}.pkl', 'rb'))
                all_communities.extend(communities)

            
            m.G = pickle.load(open(f'dataset/{net_version}/{folder}/{folder}.pkl', 'rb'))


            for iter in iter_list:
                comunities_subset = []
                for j in range(4):
                    index = j * num_iter
                    comunities_subset.extend(all_communities[index:(index + iter)])
                     
                
                value = m.RoughClustering(communities=comunities_subset)
                m.export_RC(f'stability/{net_version}/{folder}/', f'{folder}_RC_{iter}_run_{i}.txt', value)

def process_folder_parallel(net_version: str, num_iter: int, folder: str, algorithms_names: list[str], folder_path: str, iter_list: list[int]):

    m = Matrix([], {},[])
    
    params = []

    for i in range(20):

            all_communities = []

            for algorithm in algorithms_names:
                communities = pickle.load(open(f'{folder_path}/{folder}/{algorithm}_{num_iter}_run_{i}.pkl', 'rb'))
                all_communities.extend(communities)

            
            m.G = pickle.load(open(f'dataset/{net_version}/{folder}/{folder}.pkl', 'rb'))


            for iter in iter_list:
                comunities_subset = []
                for j in range(4):
                    index = j * num_iter
                    comunities_subset.extend(all_communities[index:(index + iter)])
                     
                params.append((m, comunities_subset, [net_version, folder, iter, i]))

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.starmap(wrapper_rc, params)        


def evaluate_stability_parallel(net_version: str, num_iter: int):

    algorithms_names = ['louvain', 'infomap', 'greedy', 'async_lpa']

    folder_path = f'output/stability/{net_version}'

    folder_list = os.listdir(folder_path)

    iter_list = [10, 100, 1000] if net_version == 'NetsType_1.4' else [10, 50, 100]

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.starmap(process_folder, [(net_version, num_iter, folder, algorithms_names, folder_path, iter_list) for folder in folder_list])



def wrapper_rc(m: Matrix, iterations: list, params: list):
    
    net_version = params[0]
    folder = params[1]
    iter = params[2]
    i = params[3]

    value = m.RoughClustering(communities=iterations)
    m.export_RC(f'stability/{net_version}/{folder}/', f'{folder}_RC_{iter}_run_{i}.txt', value)

def evaluate_modularity(net_version: str):
    
    algorithms_names = ['louvain', 'infomap', 'greedy', 'async_lpa']

    folder_path = f'output/stability/{net_version}'

    folder_list = os.listdir(folder_path)

    iter_list = [10, 100, 1000] if net_version == 'NetsType_1.4' else [10, 50, 100]
    
    m = Matrix([], {},[])

    modularity_dict = {}

    for folder in folder_list:
        
        print(f'Processing {folder}', datetime.datetime.now())

        m.G = pickle.load(open(f'dataset/{net_version}/{folder}/{folder}.pkl', 'rb'))

        for algorithm in algorithms_names:

            modularity_dict[algorithm] = {}


            for i in range(20):

                communities = pickle.load(open(f'{folder_path}/{folder}/{algorithm}_{100}_run_{i}.pkl', 'rb'))
                
                modularity_list = []
                for com in communities:
                    mod = nx.algorithms.community.quality.modularity(m.G, com)
                    modularity_list.append(mod)

                for iter in iter_list:
                    position = modularity_list.index(max(modularity_list[0:iter]))
                    if i not in modularity_dict[algorithm].keys():
                        modularity_dict[algorithm][i] = {iter: (position, modularity_list[position])}
                    else:
                        modularity_dict[algorithm][i][iter] = (position, modularity_list[position])
        print(f'Finished {folder}', datetime.datetime.now())
    pickle.dump(modularity_dict, open(f'output/stability/{net_version}/modularity_dict.pkl', 'wb'))         

def nmi_stability(foldername: str):
    
    algorithms_names = ['louvain', 'infomap', 'greedy', 'async_lpa']   

    from cdlib import evaluation, NodeClustering

    files = os.listdir('dataset/' + foldername)
    files.remove('GT')
    files.remove('README.txt')

    files = sorted(files, key=lambda x: int("".join([i for i in x if i.isdigit()])))

    iter_list = [10, 100, 1000] if foldername == 'NetsType_1.4' else [10, 50, 100]
    iter_lenght = 1000 if foldername == 'NetsType_1.4' else 100

    modularity_dict = pickle.load(open(f'output/stability/{foldername}/modularity_dict.pkl', 'rb'))

    df = pd.DataFrame(columns=['Network','Algorithm', 'Sequence','NMI', 'Index_Higher' , 'Iterations'])

    for file in files:
        print(f'Processing file: {file}', datetime.datetime.now())

        if '.pkl' in file:
            continue
        nodes = []
        match = re.search(r'(network)(\d+)', file)
        number = '0'
        if match:
            number = str(match.group(2))
            with open('dataset/' + foldername + '/GT/community' + number + '_GT.dat', 'r') as f:
                lines = f.readlines()        
                for line in lines:
                    data = line.split(' ')
                    inter_data = [int(x) for x in data]
                    nodes.append(inter_data)
        
        G = pickle.load(open('dataset/' + foldername + '/' + file + '/' + file + '.pkl', 'rb')) 
    
        # GT created
        nodeClustA = NodeClustering(communities=nodes, graph=G, method_name='GT', method_parameters={}, overlap=True)

        for iter in iter_list:
            for i in range(20):
                nodes = []
                with open(f'output/stability/{foldername}/{file}/{file}_RC_{iter}_run_{i}.txt', 'r') as f:
                    lines = f.readlines()        
                    for line in lines:
                        line = line.strip('\n').rstrip()
                        data = line.split(' ')
                        inter_data = [int(x) for x in data]
                        nodes.append(inter_data)

                nodeClustB = NodeClustering(communities=nodes, graph=G, method_name='RC', method_parameters={}, overlap=True)
                nmi = evaluation.overlapping_normalized_mutual_information_MGH(nodeClustA, nodeClustB)
                new_row = pd.DataFrame({'Network': file, 'Algorithm': 'RC', 'Sequence': i, 'NMI': nmi.score, 'Index_Higher': 0, 'Iterations': iter}, index=[0])
                df = pd.concat([df, new_row], ignore_index=True)
                
                for algorithm in algorithms_names:

                    communities = pickle.load(open(f'output/stability/{foldername}/{file}/{algorithm}_{iter_lenght}_run_{i}.pkl', 'rb'))
                    
                    index = modularity_dict[algorithm][i][iter][0]
                    communities = communities[index]
                    
                    nodeClustB = NodeClustering(communities=communities, graph=G, method_name=algorithm, method_parameters={}, overlap=True)
                    nmi = evaluation.overlapping_normalized_mutual_information_MGH(nodeClustA, nodeClustB)
                    
                   
                    new_row = pd.DataFrame({'Network': file, 'Algorithm': algorithm, 'Sequence': i, 'NMI': nmi.score, 'Index_Higher': index, 'Iterations': iter}, index=[0])
                    df = pd.concat([df, new_row], ignore_index=True)
        nodes = []
        print(f'Finished file: {file}', datetime.datetime.now())        
    df.to_csv(f'output/stability/{foldername}/nmi_stability.csv', index=False)

def save_nmi_pkl_to_cvs(path: str):

    pickle_data = pickle.load(open(path, 'rb'))

    df = pd.DataFrame(columns=['Network','Algorithm', 'Sequence','NMI', 'Index_Higher' , 'Iterations'])

    for data in pickle_data:
        new_row = pd.DataFrame({'Network': data[1], 'Algorithm': data[2], 'Sequence': data[3], 'NMI': data[0], 'Index_Higher': data[4], 'Iterations': data[5]}, index=[0])
        df = pd.concat([df, new_row], ignore_index=True)
    
    df.to_csv(f'{path[:-4]}.csv', index=False)

def nmi_stability_parallel(foldername: str):
    
    algorithms_names = ['louvain', 'infomap', 'greedy', 'async_lpa']   

    from cdlib import evaluation, NodeClustering

    files = os.listdir('dataset/' + foldername)
    files.remove('GT')
    files.remove('README.txt')

    files = sorted(files, key=lambda x: int("".join([i for i in x if i.isdigit()])))

    iter_list = [10, 100, 1000] if foldername == 'NetsType_1.4' else [10, 50, 100]
    iter_lenght = 1000 if foldername == 'NetsType_1.4' else 100

    modularity_dict = pickle.load(open(f'output/stability/{foldername}/modularity_dict.pkl', 'rb'))
   

    data_inpust = []

    for file in files:
        print(f'Processing file: {file}', datetime.datetime.now())

        if '.pkl' in file:
            continue
        nodes = []
        match = re.search(r'(network)(\d+)', file)
        number = '0'
        if match:
            number = str(match.group(2))
            with open('dataset/' + foldername + '/GT/community' + number + '_GT.dat', 'r') as f:
                lines = f.readlines()        
                for line in lines:
                    data = line.split(' ')
                    inter_data = [int(x) for x in data]
                    nodes.append(inter_data)
        
        G = pickle.load(open('dataset/' + foldername + '/' + file + '/' + file + '.pkl', 'rb')) 
    
        # GT created
        nodeClustA = NodeClustering(communities=nodes, graph=G, method_name='GT', method_parameters={}, overlap=True)

        for iter in iter_list:
            for i in range(20):
                nodes = []
                with open(f'output/stability/{foldername}/{file}/{file}_RC_{iter}_run_{i}.txt', 'r') as f:
                    lines = f.readlines()        
                    for line in lines:
                        line = line.strip('\n').rstrip()
                        data = line.split(' ')
                        inter_data = [int(x) for x in data]
                        nodes.append(inter_data)

                nodeClustB = NodeClustering(communities=nodes, graph=G, method_name='RC', method_parameters={}, overlap=True)
                
                                
                data_inpust.append((nodeClustA, nodeClustB, file, 'RC', i, 0, iter))
                
                
                for algorithm in algorithms_names:

                    communities = pickle.load(open(f'output/stability/{foldername}/{file}/{algorithm}_{iter_lenght}_run_{i}.pkl', 'rb'))
                    
                    index = modularity_dict[algorithm][i][iter][0]
                    communities = communities[index]
                    
                    nodeClustB = NodeClustering(communities=communities, graph=G, method_name=algorithm, method_parameters={}, overlap=True)
                    
                   
                    data_inpust.append((nodeClustA, nodeClustB, file, algorithm, i, index, iter))
                    
        nodes = []
        print(f'Finished file: {file}', datetime.datetime.now())        

    pickle.dump(data_inpust, open(f'output/stability/{foldername}/nmi_stability_data_inputs.pkl', 'wb'))

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.starmap(wrapper_nmi_stability, data_inpust)
    
    pickle.dump(results, open(f'output/stability/{foldername}/nmi_stability.pkl', 'wb'))


def wrapper_nmi_stability(nodeClustA, nodeClustB, file: str, algorithm: str, i: int, index: int, iter: int):
    
    from cdlib import evaluation

    value = evaluation.overlapping_normalized_mutual_information_MGH(nodeClustA, nodeClustB)

    return (value.score, file, algorithm, i, index, iter)           


if __name__ == '__main__':

    print(datetime.datetime.now())
    start_time = datetime.datetime.now()

    m = Matrix([], {},[])

    #analyze_overlaping('NetsType_1.4')

    #nmi_stability('NetsType_1.6')

    save_nmi_pkl_to_cvs('output/stability/NetsType_1.6/nmi_stability.pkl')

    datas = pickle.load(open('output/stability/NetsType_1.6/nmi_stability.pkl', 'rb'))

    print(datas[0])
    
    # df = pd.read_csv('output/stability/NetsType_1.4/nmi_stability.csv', header=0)
    # df = df[df['Algorithm'] == 'RC']
    # print(df.describe())
    #evaluate_modularity('NetsType_1.6')

    #evaluate_stability_parallel('NetsType_1.6', 100)
  
    #evaluate_stability('NetsType_1.4', 1000)
        
    #stability(4, 10, 'NetsType_1.6')

    #run_RC_sequences(1, 'NetsType_1.4', 100)
    
    #apply_PC_to_GT('NetsType_1.6')

    
    #apply_PC_to_RC('NetsType_1.6')
    #PC_data = pickle.load(open('output/NetsType_1.4/network10_RC_PC.pkl', 'rb')) 
    #print(sorted(PC_data.items(), key=lambda x: x[1], reverse=True)[-10:-1])
    #stability_infomap(20, 100, 'NetsType_1.6')

    # FlyCircuit Region

    #m.load_matrix_obj(path='dataset/attributed_graph-1.4.fly')
    
    #nx.write_gexf(m.G, 'output/flycircuit.gexf')
    
    #plot_degree_distribution(m)

    

    #iterations = m.load_all_algorithm_communities(algorithms=['louvain', 'greedy', 'gnfomap', 'lpa'])

    #runRoughClustering_on_FlyCircuit(m, '1.4', iterations=iterations)

    #end FlyCircuit Region

    # Benchmark Region

    #generate_pkl('NetsType_1.6')

    #m.export_infomap_iterations(folder_version='NetsType_1.6', end=5)

    #runAlgorithmSimple(m, folder_version='NetsType_1.6')
    
    # runRoughClustering('NetsType_1.6')
    # nmi_overlapping_evaluate('NetsType_1.1')
    #nmi_overlapping_evaluateTunning('NetsType_1.6')

    #end Benchmark Region
    
    #print(len(pickle.load(open('output/NetsType_1.1/network2_Infomap.pkl', 'rb'))))

    end_time = datetime.datetime.now()
    real_time = end_time - start_time    
    print(f'Elapsed time: {real_time}')
    
    
    

   