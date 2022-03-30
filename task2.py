from utils import read_graph, degree_distribution, calculate_gamma
import matplotlib.pyplot as plt

def plot(graph_label, cutoff, in_cutoff, out_cutoff):
    graph, labels = read_graph(graph_label, directed=True)
    n = graph.number_of_nodes()

    degrees = [graph.degree[node] for node in graph.nodes() if graph.degree[node] > 0]
    # regular_gamma = calculate_gamma(degrees, cutoff)
    
    in_degrees = [graph.in_degree[node] for node in graph.nodes() if graph.in_degree[node] > 0]
    in_gamma = calculate_gamma(in_degrees, in_cutoff)
    
    out_degrees = [graph.out_degree[node] for node in graph.nodes() if graph.out_degree[node] > 0]
    # out_gamma = calculate_gamma(out_degrees, out_cutoff)
    
    # print(f"Regular gamma: {regular_gamma} w/ cutoff {cutoff}")
    print(f"In gamma: {in_gamma} w/ cutoff {in_cutoff}")
    # print(f"Out gamma: {out_gamma} w/ cutoff {out_cutoff}")
    
    distribution = degree_distribution(degrees, n)
    in_distribution = degree_distribution(in_degrees, n)
    out_distribution = degree_distribution(out_degrees, n)

    plt.style.use("ggplot")
    plt.figure(figsize=(13, 7))
    plt.plot(distribution.keys(), distribution.values(), 'o', color="firebrick", label=f"Degree dis.")
    plt.plot(in_distribution.keys(), in_distribution.values(), 'x', color="forestgreen", label=f"In-degree dis.")
    plt.plot(out_distribution.keys(), out_distribution.values(), '*', color="goldenrod", label=f"Out-degree dis.")
    plt.xlabel("Degrees")
    plt.ylabel("Various degree distributions")
    plt.xscale("log")
    plt.yscale("log")
    plt.title(f"Distributions for {graph_label}.net")
    plt.legend()
    plt.subplots_adjust(left=0.07, bottom=0.105, right=0.964, top=0.94)
    plt.show()

plot("java", cutoff=4, in_cutoff=5, out_cutoff=4)
plot("lucene", cutoff=4, in_cutoff=4, out_cutoff=4)