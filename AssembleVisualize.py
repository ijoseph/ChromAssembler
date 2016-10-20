import graphviz
import os
import seaborn as sns
import matplotlib.pyplot as plt


class DeBruijnGraphVisualize:
    @staticmethod
    def to_dot(de_bruijn_graph, weights=False):
        """ Return string with graphviz representation.  If 'weights'
            is true, label edges corresponding to distinct k-1-mers
            with weights, instead of drawing separate edges for
            k-1-mer copies. """
        g = graphviz.Digraph(comment='DeBruijn graph')
        for node in iter(de_bruijn_graph.graph.keys()):
            g.node(node.k_minus_1_mer, node.k_minus_1_mer)
        for src, dsts in iter(de_bruijn_graph.graph.items()):
            for dst in dsts:
                g.edge(src.k_minus_1_mer, dst.k_minus_1_mer)
        return g

    @staticmethod
    def show_graph(de_bruijn_graph, temp_folder = "."):
        """
        Show graphviz plot from graph
        :param de_bruijn_graph:
        :param temp_folder:
        :return:
        """
        dot_string = DeBruijnGraphVisualize.to_dot(de_bruijn_graph)
        dot_string.render(filename=os.path.join(temp_folder, str(de_bruijn_graph) + ".dot"), view=True)


    @staticmethod
    def show_fragment_length_distribution(de_bruijn_graph):
        """
        Show fragment length distribution
        :return:
        """
        lengths = map(len, de_bruijn_graph.fragment_list)
        ax = sns.distplot(lengths)
        ax.set_title(label= "Length distribution of {0} fragments".format(len(lengths)))
        plt.show()


