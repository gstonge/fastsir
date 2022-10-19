from fastsir import DiscreteSIR as sir
import numpy as np
import networkx as nx

#Network parameters
N = 1000
mean_degree = 10
p = mean_degree/(N-1)
G = nx.erdos_renyi_graph(N,p)
max_degree = np.max([G.degree(n) for n in G])
edgelist = list(G.edges())

#contagion parameters
recovery_probability = 1
infection_func = lambda q,k: 1-(1-q)**k
qlist = np.linspace(0.05,0.2,20)

#contagion parameters
recovery_probability = 1
infection_func = lambda q,k: 1-(1-q)**k
qlist = np.linspace(0.05,0.2,20)

#test
infection_probability = infection_func(qlist[2],np.arange(max_degree+1))
process = sir(edgelist,recovery_probability,infection_probability)
initial_infected_nodes = {np.random.randint(N)}
process.infect_node_set(initial_infected_nodes)

process.evolve(2)
process.get_current_macro_state()
