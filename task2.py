from utils import read_graph, degree_distribution, calculate_gamma
import matplotlib.pyplot as plt

def plot(graph_label):
    graph, labels = read_graph(graph_label, directed=True)
    n = graph.number_of_nodes()
    cutoff = 3

    degrees = [graph.degree[node] for node in graph.nodes() if graph.degree[node] > 0]
    regular_gamma = calculate_gamma(degrees, cutoff)
    # regular_est_max = min(degrees) * (n ** (1 / (regular_gamma - 1)))
    # regular_est_max = cutoff * (n ** (1 / (regular_gamma - 1)))
    # regular_real_max = max(degrees)
    # print(f"(Regular) Max degree: {regular_real_max}, estimated max: {regular_est_max}")

    in_degrees = [graph.in_degree[node] for node in graph.nodes() if graph.in_degree[node] > 0]
    in_gamma = calculate_gamma(in_degrees, cutoff)
    # in_est_max = min(in_degrees) * (n ** (1 / (in_gamma - 1)))
    # in_est_max = cutoff * (n ** (1 / (in_gamma - 1)))
    # in_real_max = max(in_degrees)
    # print(f"(In) Max degree: {in_real_max}, estimated max: {in_est_max}")

    out_degrees = [graph.out_degree[node] for node in graph.nodes() if graph.out_degree[node] > 0]
    out_gamma = calculate_gamma(out_degrees, cutoff)
    # out_est_max = min(out_degrees) * (n ** (1 / (out_gamma - 1)))
    # out_est_max = cutoff * (n ** (1 / (out_gamma - 1)))
    # out_real_max = max(out_degrees)
    # print(f"(Out) Max degree: {out_real_max}, estimated max: {out_est_max}")

    print(f"Regular gamma: {regular_gamma}, in gamma: {in_gamma}, out gamma: {out_gamma}")
    max_degree = max(degrees + in_degrees + out_degrees)
    power_law_x = range(1, max_degree + 1)

    distribution = degree_distribution(degrees, n)
    in_distribution = degree_distribution(in_degrees, n)
    out_distribution = degree_distribution(out_degrees, n)
    power_law_distribution = [elem ** (-2) for elem in power_law_x]

    plt.style.use("ggplot")
    plt.figure(figsize=(13, 7))
    plt.plot(distribution.keys(), distribution.values(), 'o', color="firebrick", label=f"Degree dis.")
    plt.plot(in_distribution.keys(), in_distribution.values(), 'x', color="forestgreen", label=f"In-degree dis.")
    plt.plot(out_distribution.keys(), out_distribution.values(), '*', color="goldenrod", label=f"Out-degree dis.")
    plt.plot(power_law_x, power_law_distribution, "--", color="darkorchid", label=f"Power-law w/ gamma -2")
    plt.xlabel("Degrees")
    plt.ylabel("Various degree distributions")
    plt.xscale("log")
    plt.yscale("log")
    plt.title(f"Distributions for {graph_label}.net")
    plt.legend()
    plt.subplots_adjust(left=0.07, bottom=0.105, right=0.964, top=0.94)
    plt.show()

plot("java")
plot("lucene")

"""
java.net
    Regular gamma: 1.96, in gamma: 1.78, out gamma: 2.21 za cutoff 3
    regular degree - 
    in degree - 
    out degree -
    I seriously have no idea

lucene.net
    Regular gamma: 2.12, in gamma: 1.75, out gamma: 2.46 za cutoff 3

    regular degree - seems like power law?
    in degree - doesn't seem like power law? 
    out degree - seems like power law?

"""