from utils import read_graph, degree_distribution, calculate_gamma
import matplotlib.pyplot as plt
from random import randrange
# import networkx as nx
# from math import log

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

def plot():
    """
    Plot stuff
    """
    cutoff = 3

    # Original graph distribution, gamma
    graph, labels = read_graph("social", directed=False)
    n = graph.number_of_nodes()
    # print(f"Checking if our original graph is ultra small world.")
    # print(f"\tAverage distance: {nx.average_shortest_path_length(graph)}, ln(ln(n)): {log(log(n))}")
    degrees = [graph.degree[node] for node in graph.nodes() if graph.degree[node] > 0]
    gamma = calculate_gamma(degrees, cutoff)
    print(f"Gamma for original graph: {gamma}")
    distribution = degree_distribution(degrees, n)
    max_degree = max(degrees)

    # Random walk induced graph and the distribution
    induced = random_walk_induced_graph(graph)
    n_induced = induced.number_of_nodes()
    # print(f"Checking if our induced graph is ultra small world.")
    # print(f"\tAverage distance: {nx.average_shortest_path_length(induced)}, ln(ln(n)): {log(log(n_induced))}")
    print(f"Nodes in random walk: {n_induced}")
    degrees_induced = [induced.degree[node] for node in induced.nodes() if induced.degree[node] > 0]
    gamma_induced = calculate_gamma(degrees_induced, cutoff)
    print(f"Gamma for induced graph: {gamma_induced}")
    distribution_induced = degree_distribution(degrees_induced, n_induced)
    
    # Power laws
    power_law_x = range(1, max_degree+1)
    power_law_calc = [elem ** (-gamma) for elem in power_law_x]
    power_law_2 = [elem ** (-2) for elem in power_law_x]
    power_law_1_5 = [elem ** (-1.5) for elem in power_law_x] # this one is mostly for reference
    
    plt.style.use("ggplot")
    plt.figure(figsize=(13, 7))
    plt.plot(distribution.keys(), distribution.values(), 'o', color="firebrick", label=f"Original graph deg. dis.")
    plt.plot(distribution_induced.keys(), distribution_induced.values(), 'X', color="dodgerblue", label=f"Induced graph deg. dis.")
    plt.plot(power_law_x, power_law_calc, "--", color="goldenrod", label=f"Power-law w/ calc gamma {round(gamma, 2)}, cutoff {cutoff}")
    plt.plot(power_law_x, power_law_1_5, "--", color="darkorchid", label=f"Power-law w/ gamma 1.5")
    plt.plot(power_law_x, power_law_2, "--", color="forestgreen", label=f"Power-law w/ gamma 2")
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
Looks like the random walk induced graph is not scale free
while the original graph is, at least looking at the
degree distribution

Seems like both are "ultra small world" (<d> ~ lnln(n))
or even small world

Checking if our original graph is ultra small world.
        Average distance: 7.4855400514784, ln(ln(n)): 2.2274442304788287
Gamma for original graph: 2.0605610597029367
Checking if our induced graph is ultra small world.
        Average distance: 5.105590247429701, ln(ln(n)): 1.942123419797389
Nodes in random walk: 1068
Gamma for induced graph: 1.7318751394190044
"""