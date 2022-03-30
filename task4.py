from utils import read_graph, degree_distribution, calculate_gamma
import matplotlib.pyplot as plt
from random import randrange
import networkx as nx
from math import log

def random_walk_induced_graph(graph):
    """
    Do a random walk on the graph, until 
    10% of the nodes is sampled

    Return the induced graph
    """
    n = graph.number_of_nodes()
    current = randrange(0, n) # random starting node
    walk = [current]
    sample_size = round(0.1 * n) # need to sample 10% of the nodes
    visited = set()
    visited.add(current)

    while len(visited) < sample_size:
        neighbors = list(graph[current].keys())
        index = randrange(0, len(neighbors))
        node = neighbors[index]
        if node not in visited:
            visited.add(node)
        walk.append(node)
        current = node

    return graph.subgraph(visited)

def average_degree(graph):
    return 2*graph.number_of_edges() / graph.number_of_nodes()

def plot():
    """
    Plot stuff
    """
    cutoff = 3
    cutoff_induced = 4

    # Original graph distribution, gamma
    graph, labels = read_graph("social", directed=False)
    n = graph.number_of_nodes()
    # print(f"Checking if our original graph is small world.")
    # print(f"\tAverage distance: {nx.average_shortest_path_length(graph)}")
    # print(f"\tln(n) / ln(<k>): {log(n) / log(average_degree(graph))}")
    # print(f"\tAverage clustering: {nx.average_clustering(graph)}")
    degrees = [graph.degree[node] for node in graph.nodes() if graph.degree[node] > 0]
    gamma = calculate_gamma(degrees, cutoff)
    # print(f"Gamma for original graph (cutoff {cutoff}): {gamma}\n")
    distribution = degree_distribution(degrees, n)
    max_degree = max(degrees)

    # Random walk induced graph and the distribution
    induced = random_walk_induced_graph(graph)
    n_induced = induced.number_of_nodes()
    # print(f"Checking if our induced graph is small world.")
    # print(f"\tAverage distance: {nx.average_shortest_path_length(induced)}")
    # print(f"\tln(n) / ln(<k>): {log(n_induced) / log(average_degree(induced))}")
    # print(f"\tAverage clustering: {nx.average_clustering(induced)}")
    degrees_induced = [induced.degree[node] for node in induced.nodes() if induced.degree[node] > 0]
    gamma_induced = calculate_gamma(degrees_induced, cutoff_induced)
    print(f"Gamma for induced graph (cutoff {cutoff_induced}): {gamma_induced}\n")
    distribution_induced = degree_distribution(degrees_induced, n_induced)
    
    # Power laws
    power_law_x = range(1, max_degree + 1)
    power_law_original = [elem ** (-gamma) for elem in power_law_x]
    power_law_induced = [elem ** (-gamma_induced) for elem in power_law_x]
    
    # Plot
    plt.style.use("ggplot")
    plt.figure(figsize=(13, 7))
    plt.plot(distribution.keys(), distribution.values(), 'o', color="firebrick", label=f"Original graph deg. dis.")
    plt.plot(distribution_induced.keys(), distribution_induced.values(), 'X', color="dodgerblue", label=f"Induced graph deg. dis.")
    plt.plot(power_law_x, power_law_induced, "--", color="forestgreen", 
        label=f"Power-law for induced graph, gamma {round(gamma_induced, 2)}, cutoff {cutoff_induced}")
    plt.plot(power_law_x, power_law_original, "--", color="goldenrod", 
        label=f"Power-law for original graph, gamma {round(gamma, 2)}, cutoff {cutoff}")
    plt.xlabel("Induced graph degrees")
    plt.ylabel("Distribution p_k / power-law k^(-gamma)")
    plt.xscale("log")
    plt.yscale("log")
    plt.title("Random walk induced graph degree distribution")
    plt.legend(loc=3)
    plt.subplots_adjust(left=0.085, bottom=0.105, right=0.964, top=0.94)
    plt.show()
        
plot()

"""
Checking if our original graph is ultra small world.
        Average distance: 7.4855400514784
        ln(n) / ln(<k>): 6.119185664133248
        Average clustering: 0.2659452243010437
Gamma for original graph (cutoff 3): 2.0605610597029367
Gamma for original graph (cutoff 4): 2.152434367238218
Gamma for original graph (cutoff 5): 2.232880279297983

Checking if our induced graph is ultra small world.
        Average distance: 5.525734584346886
        ln(n) / ln(<k>): 2.969112093861218
        Average clustering: 0.46075452939143896
Gamma for induced graph (cutoff 3): 1.7318751394190044 (varies between 1.7x and 1.8x)
Gamma for induced graph (cutoff 4): 1.8966716742827932 (varies between 1.8x and 1.9x)
Gamma for induced graph (cutoff 5): 2.06510132214632 (varies between 1.9x and 2.0x)
"""