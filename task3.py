from utils import read_graph
from random import sample
import matplotlib.pyplot as plt
import networkx as nx

def bfs(graph, start, visited):
    """
    Do BFS to find number of nodes in this connected component
    """
    queue = [start]
    visited.add(start)
    num_visited = 0
    
    while queue:
        current = queue.pop()
        num_visited += 1
        for neighbor in graph[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    
    return num_visited

def lcc(graph):
    """
    Return the largest connected component of the graph by size
    """
    visited = set()
    largest = 0
    for node in graph.nodes():
        if node not in visited:
            largest = max(largest, bfs(graph, node, visited))
    return largest

def prune_random(graph, num_to_remove):
    """
    Remove num_to_remove random nodes from the graph
    """
    while num_to_remove > 0:
        node = sample(list(graph.nodes()), 1)[0]
        graph.remove_node(node)
        num_to_remove -= 1

def prune_hubs(graph, num_to_remove, iteration):
    """
    Remove num_to_remove highest degree nodes from the graph
    """
    degrees = sorted(graph.degree, key=lambda elem: elem[1], reverse=True)
    for i in range(num_to_remove):
        current_max, d = degrees[i]
        graph.remove_node(current_max)
        # print(f"\t({iteration}) At {i} out of {num_to_remove}")

def plot_failure(graph):
    """
    Plot stuff
    """
    graph_hubs = graph.copy()
    n = graph.number_of_nodes()
    ten_percent = round(0.1 * n) # amount to remove each time
    fractions_random = [] # ith element = fraction of nodes in LCC after removing (i*10)% random nodes
    fractions_hubs = []

    for i in range(0, 6):
        print(f"Iteration {i}")
        if i > 0:
            # remove 10% of nodes each time
            prune_random(graph, ten_percent)
            prune_hubs(graph_hubs, ten_percent, i)

        lcc_random_number = lcc(graph)
        lcc_hubs_number = lcc(graph_hubs)

        lcc_random_fraction = lcc_random_number / graph.number_of_nodes()
        lcc_hubs_fraction = lcc_hubs_number / graph_hubs.number_of_nodes()

        fractions_random.append(lcc_random_fraction)
        fractions_hubs.append(lcc_hubs_fraction)

    plt.style.use("ggplot")
    plt.figure(figsize=(10, 6))
    x = range(0, 6)
    plt.plot(x, fractions_random, 'x', color="firebrick", label="Random removal")
    plt.plot(x, fractions_hubs, '.', color="forestgreen", label="Hub removal")
    plt.title(f"Failure simulation for {graph.name}")
    plt.xlabel("Percentage removed (divided by 10)")
    plt.ylabel("Fraction of nodes in LCC")
    # write the actual values on the graph, slightly above/below the point
    for i, j in zip(x, fractions_random):
        text = (-7, 7) if i == 0 else (-14, 7)
        plt.annotate(str(round(j, 4)), xy=(i, j), xytext=text, textcoords='offset points')
    for i, j in zip(x, fractions_hubs):
        text = (-7, -14) if i == 0 else (-14, -14)
        plt.annotate(str(round(j, 4)), xy=(i, j), xytext=text, textcoords='offset points')
    plt.legend(loc=3) # lower left
    plt.show()

graph, labels = read_graph("internet", directed=False)
erdos = nx.gnm_random_graph(graph.number_of_nodes(), graph.number_of_edges())
erdos.name = "ER"
plot_failure(graph)
plot_failure(erdos)