import networkx as nx
from collections import defaultdict
from math import log

def read_graph(filename, directed=False):
    if directed:
        G = nx.DiGraph(name=filename)
    else:
        G = nx.Graph(name=filename)
    labels = {}
    with open(f"graphs/{filename}.net", 'r', encoding='utf8') as f:
        f.readline()
        
        for line in f:
            if line.startswith("*"):
                break
            else:
                node_stuff = line.split('"')
                node = int(node_stuff[0]) - 1
                label = node_stuff[1].strip()
                G.add_node(node)
                labels[node] = label
                
        for line in f:
            node1_str, node2_str = line.split()[:2]
            G.add_edge(int(node1_str) - 1, int(node2_str) - 1)
            if not directed:
                G.add_edge(int(node2_str) - 1, int(node1_str) - 1)
                
    return G, labels

def degree_distribution(sequence, n):
    """
    Return the degree distribution for this degree sequence
    """
    distribution = defaultdict(int)
    for degree in sequence:
        distribution[degree] += 1
    # import json
    # print(json.dumps(dict(sorted(distribution.items(), key=lambda elem: elem[0], reverse=True)), indent=4))
    for degree in distribution:
        distribution[degree] /= n
    return distribution

def calculate_gamma(degree_sequence, cutoff):
    # new_n = 0
    # s = 0
    # for degree in degree_sequence:
    #     if degree >= cutoff:
    #         s += log(degree / (cutoff - 0.5))
    #         new_n += 1
    # return 1 + new_n / s
    cutoff_degrees = [degree for degree in degree_sequence if degree >= cutoff]
    return 1 + len(cutoff_degrees) / (sum([log(deg / (cutoff - 0.5)) for deg in cutoff_degrees]))
