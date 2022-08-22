from graph_tool.all import *
import numpy as np
import pandas as pd
import matplotlib.cm
import pickle
from collections import defaultdict

class seastar_hsbm():
    '''
    Class for selfunity detection using hSBMs
    '''

    def __init__(self):
        self.g = None # network

        self.samples = [] # list of seastar sample indices
        self.otus = [] # list of taxa labels

        self.state = None # inference state from graphtool
        self.dls = None # list of description lengths throughout optimization

    def make_graph(self, f, index_col="index"):
        '''
        Build the graph from csv file f, setting sample and otu indices as well
        '''
        df = pd.read_csv(f, index_col=index_col)
        otu_all = [c for c in df.columns if "k__" in c]
        otu_filter = [col for col in otu_all if sum(df[col]) > 100]
        df = df[otu_filter]
        samples = np.array(df.index)
        otus = np.array(otu_filter)

        g = Graph(directed=False)
        g.add_vertex(len(samples) + len(otus))

        # add necessary "graph properties"
        g.vp.name = g.new_vertex_property("string", np.append(samples, otus))
        g.vp.part = g.new_vertex_property("int", np.append(np.ones(len(samples)), np.zeros(len(otus))))
        g.ep.counts = g.new_edge_property("int")

        for si in range(len(samples)):
            for ti in range(len(samples), len(samples) + len(otus)):
                c = int(df.loc[g.vp.name[si], g.vp.name[ti]])
                if c > 0:
                    e = g.add_edge(si, ti)
                    g.ep.counts[e] = np.log(c) + 1

        self.g = g
        self.samples = samples
        self.otus = otus

    def fit(self, deg_corr=True, niter=1000):
        '''
        Find clusters in the network based on the maximum likelihood for the hSBM 
        '''
        part = self.g.vp.part
        counts = self.g.ep.counts
        state_args = {"clabel": part, "pclabel": part, "eweight": counts, "deg_corr": deg_corr}
        state_tmp = minimize_nested_blockmodel_dl(self.g, state_args=state_args)
        dl_curr = state_tmp.entropy()
        dls = [dl_curr]

        for _ in range(niter): # this should be sufficiently large
            res = state_tmp.multiflip_mcmc_sweep(beta=np.inf, niter=10)
            dl_curr += res[0]
            dls.append(dl_curr)

        # try to remove some of the empty levels
        L = 0
        for s in state_tmp.levels:
            L += 1
            if s.get_nonempty_B() == 2:
                break

        # save the final fit state
        self.state = state_tmp.copy(bs=state_tmp.get_bs()[:L] + [np.zeros(1)])
        self.dls = dls

    def save_model(self, path):
        '''
        Save the trained model in the specified path as a pickle
        '''
        if '.pickle' not in path:
            path += '.pickle'
        with open(path, 'wb') as f:
            pickle.dump(self, f)

    def save_model_dict(self, path):
        '''
        Save the model in an R-friendly format (a dictionary)
        '''
        if '.pickle' not in path:
            path += '.pickle'

        membership = dict()
        for lvl in range(len(self.state.levels) - 1):
            state_proj = self.state.project_level(lvl)
            mem = defaultdict(lambda: list())
            for v, i in enumerate(state_proj.get_blocks().a):
                mem[str(i)].append(self.g.vp.name[v])
            membership[lvl] = dict(mem)

        saver = {"membership": membership, "dls": self.dls}
        with open(path, 'wb') as f:
            pickle.dump(saver, f)

    def load_model(self, path):
        '''
        Load the trained model from the specified path to the pickle file
        '''
        if '.pickle' not in path:
            path += '.pickle'
        with open(path, 'rb') as f:
            obj = pickle.load(f)
            self.__dict__.update(obj.__dict__)

    def plot(self, path):
        part = self.g.vp.part
        counts = self.g.ep.counts
        self.state.draw(
            layout="bipartite", output=path, subsample_edges=1000, hshortcuts=1, 
            # edge_pen_width=prop_to_size(g.ep.counts, ma=4, power=1.2, log=False), 
            edge_pen_width=1.5,
            vertex_size=1.6,
            vertex_fill_color=part,
            vertex_color=part,
            edge_color=prop_to_size(counts, ma=10, power=1),
            ecmap=(matplotlib.cm.inferno, .4), edge_gradient=[], eorder=counts
        )


