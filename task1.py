import networkx as nx
from utils import read_graph

def sort_node_centrality_list(centrality_list):
    """
    Sorts list of pairs (node, centrality_measure) 
    by the second element in that pair
    """
    return sorted([(node, c) for node, c in centrality_list], key=lambda x: x[1], reverse=True)

def map_centrality_list(centrality_list, labels):
    """
    Map (node_number, centrality_measure) to (node_label, centrality_measure)
    """
    return list(map(lambda elem: (labels[elem[0]], elem[1]), centrality_list))

def find_sn100(centrality_list):
    """
    Find the almighty SN100 dolphin in a given
    (mapped) centrality list
    """
    for i, elem in enumerate(centrality_list):
        dolphin, importance = elem
        if dolphin == "SN100":
            return i + 1
    return None

#############################
##### Degree centrality #####
#############################

def degree_centrality(G):
    deg_centrality = [(node, G.degree[node] / (G.number_of_nodes()-1)) for node in G.nodes()]
    return sort_node_centrality_list(deg_centrality)


#############################################
##### Clustering coefficient centrality #####
#############################################

def get_link_triads(G, i, j):
    """
    Return link triads containing nodes i and j
    In other words return the number of triangles
    involving edge (i, j)
    """
    return len(set(G[i]).intersection(G[j]))

def get_node_triads(G, i):
    """
    Return triads containing node i
    In other words count number of triangles 
    involving node i
    """
    t = 0
    for j in G[i]:
        if G.degree[i] <= G.degree[j]:
            t += get_link_triads(G, i, j) / 2
        else:
            t += get_link_triads(G, j, i) / 2
    return t

def get_node_clustering_coef(G, i):
    """
    Return node clustering coef for node i
    """
    k_i = G.degree[i]
    if k_i <= 1:
        return 0
    return get_node_triads(G, i) * 2 / (k_i*k_i-k_i)

def clustering_centrality(G):
    """
    Calculates clustering coeficient for each node
    """
    cluster_centrality = [(node, get_node_clustering_coef(G, node)) for node in G.nodes()]
    return sort_node_centrality_list(cluster_centrality)

def mu_clustering_centrality(G):
    """
    Calculates mu-clustering coeficient for each node
    
    C_i^{mu} = C_i * (k_i - 1) (additionally, can divide with mu)
    where mu is the maximum number of triangles over any edge
    """
    mu = max([get_link_triads(G, i, j) for i, j in G.edges()])
    mu_centrality = [(node, get_node_clustering_coef(G, node) * (G.degree[node] - 1) / mu) for node in G.nodes()]
    return sort_node_centrality_list(mu_centrality)


################################
##### Closeness centrality #####
################################

def get_distances(graph, start):
    """
    Return a list of distances starting from
    this node - i.e. how far away every
    other node is from the starting node
    """
    distances = {}
    queue = []
    distances[start] = 0
    queue.append(start)
    while queue:
        node = queue.pop(0)
        for neighbor in graph[node]:
            if neighbor not in distances:
                distances[neighbor] = distances[node] + 1
                queue.append(neighbor)
    return [d for d in distances.values() if d > 0]

def closeness_centrality(G):
    """
    Return closeness centrality for each node using formula

    l_i^{-1} = 1/(n-1) sum_{j =/= i} 1/(d_ij)

    where d_ij is the (shortest) distance between nodes i and j
    """
    constant = G.number_of_nodes() - 1
    close_centrality = [(node, sum([1/d for d in get_distances(G, node)]) / constant) for node in G.nodes()]
    return sort_node_centrality_list(close_centrality)


##################################
##### Betweenness centrality #####
##################################

def betweenness_centrality(G):
    """
    Calculates betweenness centrality for each node
    """ 
    between_centrality = nx.betweenness_centrality(G).items()
    return sort_node_centrality_list(between_centrality)


##################################
##### Eigenvector centrality #####
##################################

def eigenvector_centrality(G, eps=1e-6):
    """
    Calculates eigenvector centrality for each node
    """

    # define list of default centrality values for each node (set value = 1)
    # E = [1] * G.number_of_nodes()
    E = {node: 1 for node in G.nodes()}

    # define default difference between old and new values (set value > eps)
    diff = 1

    while diff > eps:
        # U = sum of neighboring centrality values for each node
        # U = [sum([E[neighbor] for neighbor in G[node]]) for node in G.nodes()]
        U = {node: sum([E[neighbor] for neighbor in G[node]]) for node in G.nodes()}

        # calculate normalizing constant u
        # u = sum(U)
        u = sum(U.values())

        # update all nodes' values in U according to the algo
        # U = [U[node] * G.number_of_nodes()/u for node in G.nodes()]
        U = {node: U[node] * G.number_of_nodes()/u for node in G.nodes()}

        # update diff = sum of old value - new value over each node
        diff = sum([abs(E[node]-U[node]) for node in G.nodes()])

        # update E
        E = U

    # eigen_centrality = [(i, E[i]) for i in range(len(E))]
    eigen_centrality = E.items()
    return sort_node_centrality_list(eigen_centrality) 


###############################
##### PageRank centrality #####
###############################

def pagerank_centrality(G, alpha=0.85, eps=1e-6):
    """
    Calculates pagerank centrality for each node
    """

    # define list of default centrality values for each node (set value = 1/|nodes|)
    # P = [1 / G.number_of_nodes()] * G.number_of_nodes()
    P = {node: 1 / G.number_of_nodes() for node in G.nodes()}

    # define default difference between old and new values (set value > eps)
    diff = 1

    while diff > eps:
        # init U according to the formula in the algo
        # U = [sum([P[j] * alpha / G.degree(j) for j in G[i]]) for i in G.nodes()]
        U = {node: sum([P[j] * alpha / G.degree(j) for j in G[node]]) for node in G.nodes()}

        # calculate normalizing constant u
        # u = sum(U)
        u = sum(U.values())

        # update all nodes' values in U according to the algo
        # U = [U[i] + (1-u) / G.number_of_nodes()  for i in G.nodes()]
        U = {node: U[node] + (1-u) / G.number_of_nodes() for node in G.nodes()}

        # update diff = sum of old value - new value over each node
        diff = sum([abs(P[i]-U[i]) for i in G.nodes()])

        # update P
        P = U

    # build sorted list of pairs (node, centrality)
    # node_centrality = [(i, P[i]) for i in range(len(P))]
    pr_centrality = P.items()
    return sort_node_centrality_list(pr_centrality)


graph, labels = read_graph("dolphins", directed=False)
centrality_functions = [
    degree_centrality,
    clustering_centrality,
    mu_clustering_centrality,
    closeness_centrality,
    betweenness_centrality,
    eigenvector_centrality,
    pagerank_centrality
]

for f in centrality_functions:
    result = map_centrality_list(f(graph), labels)
    sn100_rank = find_sn100(result)
    print(f"Function: {f.name}, most important dolphin: {result[0]} \
        \n\tSN100 rank: {sn100_rank}, value: {result[sn100_rank-1][1]}")
    print()