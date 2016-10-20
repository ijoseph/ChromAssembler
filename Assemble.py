# -*- coding: utf-8 -*-
"""
@author: Isaac C. Joseph
"""
import argparse, sys, multiprocessing, os
import graphviz
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import math
import logging

class Assembler(object):
    """
    Performs file handing and uses an assembler_name to assemble reads
    """

    def __init__(self, fragments_fasta, temp_folder = ".", assembler = "DeBrujinGraph", output = sys.stdout):
        """
        Constructor
        """
        self.assembler_name = assembler
        self.output = output
        self.temp_folder = temp_folder

        self.sequence_list = self.read_fasta(fragments_fasta)

    def assemble(self):
        """ Chooses the approrpiate assembler and tells it to assemble"""
        if self.assembler_name == "DeBrujinGraph":
            assembler= DeBruijnGraph2(sequence_list=self.sequence_list)
        else:
            raise NotImplementedError("Assemblers beyond DeBruijinGraph")

        return assembler.assemble()

    def read_fasta(self, fragments_fasta):
        """
        Parses FASTA text file
        :return: list of strings with each sequence
        """
        all_sequences = []
        this_sequence = []
        for line in fragments_fasta:
            if line.startswith(">"): # Header line
                if len(this_sequence): # If we have been reading sequence, finish doing so
                    all_sequences += ["".join(this_sequence)]
                    this_sequence = []
                continue

            this_sequence += [line.strip()] # Add to current sequence

        return all_sequences + ["".join(this_sequence)]


        # ax = sns.distplot(map(len, self.sequence_list))
        # ax.set_title( label= "{0} is min".format(min(map(len,self.sequence_list))))
        #
        # plt.show()


class DeBruijnGraph:
    """ De Bruijn directed multigraph built from a collection of
        strings. User supplies strings and k-mer length k.  Nodes
        are k-1-mers.  An Edge corresponds to the k-mer that joins
        a left k-1-mer to a right k-1-mer. """

    @staticmethod
    def chop(sequence, k):
        """ Chop sequence into k-mers of given length """
        for i in range(len(sequence) - (k - 1)):
            yield (sequence[i:i + k], sequence[i:i + k - 1], sequence[i + 1:i + k])

    class Node:
        """ Node representing a k-1 mer.  Keep track of # of
            incoming/outgoing edges so it's easy to check for
            balanced, semi-balanced. """

        def __init__(self, k_minus_1_mer):
            self.k_minus_1_mer = k_minus_1_mer
            self.in_degree = 0
            self.out_degree = 0

        def is_semi_balanced(self):
            return abs(self.in_degree - self.out_degree) == 1

        def is_balanced(self):
            return self.in_degree == self.out_degree

        def __hash__(self):
            return hash(self.k_minus_1_mer)

        def __str__(self):
            return self.k_minus_1_mer

    def __init__(self, sequence_list):
        """
        Build de Bruijn multigraph given string iterator and k-mer
            length k
        """

        self.sequence_list = sequence_list

        self.graph = {}  # multimap from nodes to neighbors
        self.nodes = {}  # maps k-1-mers to Node objects
        chosen_k = self.choose_k()
        self.build_graph(k= chosen_k, sequence_list=sequence_list)



    def choose_k(self):
        """
        Choose k for k-mer De Bruijin Graph assembly.
        Chooses k as ⌊1/2 minimum read length⌋
        :return:
        """

        min_read_length = min(map(len, self.sequence_list))

        chosen_k = int(math.floor(.5 * min_read_length))

        logging.info("Chosen k: {0}".format(chosen_k))

        return chosen_k


    def build_graph(self, k, sequence_list):

        seen_k_mers = set()

        for sequence in sequence_list:
            for k_mer, k_minus_1_mer_left, k_minus_1_mer_right in self.chop(sequence, k):
                if k_mer in seen_k_mers:
                    continue
                seen_k_mers.add(k_mer)
                if k_minus_1_mer_left in self.nodes:
                    nodeL = self.nodes[k_minus_1_mer_left]
                else:
                    nodeL = self.nodes[k_minus_1_mer_left] = self.Node(k_minus_1_mer_left)
                if k_minus_1_mer_right in self.nodes:
                    nodeR = self.nodes[k_minus_1_mer_right]
                else:
                    nodeR = self.nodes[k_minus_1_mer_right] = self.Node(k_minus_1_mer_right)
                nodeL.out_degree += 1
                nodeR.in_degree += 1
                self.graph.setdefault(nodeL, []).append(nodeR)

        self.assess_graph()

    def assess_graph(self):
        """
        Keep track of head and tail nodes in the case of a graph with Eularian walk (not cycle); set tail and head
        :return:
        """
        # Iterate over nodes; tally # balanced, semi-balanced, neither
        self.number_semi_balanced_nodes, self.number_balanced_nodes, self.number_unbalanced_nodes = 0, 0, 0

        self.head, self.tail = None, None
        for node in iter(self.nodes.values()):
            if node.is_balanced():
                self.number_balanced_nodes += 1
            elif node.is_semi_balanced():
                if node.in_degree == node.out_degree + 1:
                    self.tail = node
                if node.in_degree == node.out_degree - 1:
                    self.head = node
                self.number_semi_balanced_nodes += 1
            else:
                self.number_unbalanced_nodes += 1

    def number_of_nodes(self):
        """ Return # nodes """
        return len(self.nodes)

    def number_of_edges(self):
        """ Return # edges """
        return len(self.graph)

    def has_eulerian_walk(self):
        """ Return true iff graph has Eulerian walk. """
        return self.number_unbalanced_nodes == 0 and self.number_semi_balanced_nodes == 2

    def has_eulerian_cycle(self):
        """ Return true iff graph has Eulerian cycle. """
        return self.number_unbalanced_nodes == 0 and self.number_semi_balanced_nodes == 0

    def is_eulerian(self):
        """ Return true iff graph has Eulerian walk or cycle """
        # technically, if it has an Eulerian walk
        return self.has_eulerian_walk() or self.has_eulerian_cycle()

    def get_eulerian_walk_or_cycle(self):
        """
        Find and return sequence of nodes (represented by
            their k-1-mer labels) corresponding to Eulerian walk
            or cycle
        """
        assert self.is_eulerian() # graph g has an Eulerian cycle
        graph_copy = self.graph

        if self.has_eulerian_walk():
            graph_copy = graph_copy.copy()
            graph_copy.setdefault(self.tail, []).append(self.head)

        tour = []
        source_node = next(iter(graph_copy.keys()))

        node = source_node
        while len(graph_copy[node]) > 0 :
            node = graph_copy[node].pop() # move to neighboring node and continue
            tour.append(node)


        if self.has_eulerian_walk():
            # Adjust node list so that it starts at head and ends at tail
            sti = tour.index(self.head)
            tour = tour[sti:] + tour[:sti]

        # Return node list as string
        return list(map(str, tour))

    def assemble(self):
        # dot_file = dbg.to_dot()
        # dot_file.render(filename=os.path.join(self.temp_folder,
        #                                       "Figures/{0}_smol_example.dot").format(k), view=True)

        assert self.is_eulerian(), \
            "Cannot assemble these sequences with k = {0}; " \
            "try specifying k to be different than the above?".format(self.choose_k())

        # return "".join(dbg.eulerianWalkOrCycle())
        walk = self.get_eulerian_walk_or_cycle()
        return walk[0] + ''.join(map(lambda x: x[-1], walk[1:]))

class DeBruijnGraph2(DeBruijnGraph):
    def to_dot(self, weights=False):
        """ Return string with graphviz representation.  If 'weights'
            is true, label edges corresponding to distinct k-1-mers
            with weights, instead of drawing separate edges for
            k-1-mer copies. """
        g = graphviz.Digraph(comment='DeBruijn graph')
        for node in iter(self.graph.keys()):
            g.node(node.k_minus_1_mer, node.k_minus_1_mer)
        for src, dsts in iter(self.graph.items()):
            if weights:
                weightmap = {}
                if weights:
                    for dst in dsts:
                        weightmap[dst] = weightmap.get(dst, 0) + 1
                for dst, v in weightmap.iteritems():
                    g.edge(src.k_minus_1_mer, dst.k_minus_1_mer, label=str(v))
            else:
                for dst in dsts:
                    g.edge(src.k_minus_1_mer, dst.k_minus_1_mer)
        return g



def main():
    parser = argparse.ArgumentParser(description="description")
    parser.add_argument("--fragments", type=file, help="Input fragments in FASTA format")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output assembled sequence to file [standard out]")
    parser.add_argument('-p', '--numThreads', type=int, default=multiprocessing.cpu_count() / 2,
                        help="Number of threads to run [(# of CPUs)/2]")
    parser.add_argument('--tempFolder', default="temp/", help="temporary folder for files [./temp]")
    namespace = parser.parse_args()

    # assembler_name = Assembler(fragments_fasta=namespace.fragments, output=namespace.output,
    #                       assembler_name="ShortestCommonSuperstring")

    assembler = Assembler(fragments_fasta=namespace.fragments)
    namespace.output.write(assembler.assemble())
    print len(assembler.assemble())






if __name__ == '__main__':
    main()