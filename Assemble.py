# -*- coding: utf-8 -*-
"""
@author: Isaac C. Joseph ijoseph@berkeley.edu

Implements Assembler, an assembly class runner, and
DeBruijnGraph, an assembler adapted from Dr. Ben Langmead:
http://www.langmead-lab.org/teaching-materials/.
"""
import argparse, sys, math, logging


class Assembler(object):
    """
    Performs file handing and uses an `assembler_name`-type assembler to assemble reads
    """



    def __init__(self, fragments_fasta, assembler_name ="DeBrujinGraph"):
        """
        Constructor: set assembler based on name, read in fragments
        """
        self.assembler_name = assembler_name
        self.fragment_list = self.read_fasta(fragments_fasta)

        if self.assembler_name == "DeBrujinGraph":
            self.assembler = DeBruijnGraph(fragment_list=self.fragment_list)
        else:
            raise NotImplementedError("Assemblers other than DeBruijinGraph not implemented")

    def assemble(self):
        """ Simply tells assembler to do its work."""
        return self.assembler.assemble()

    def read_fasta(self, fragments_fasta):
        """
        Parses FASTA text file
        :return: list of fragment strings
        """
        all_fragments = []
        this_fragment = []
        max_line_len = 0
        for line in fragments_fasta:
            if line.startswith(">"): # Header line
                if len(this_fragment): # If we have been reading fragment, finish doing so
                    all_fragments += ["".join(this_fragment)]
                    this_fragment = []
                continue

            this_fragment += [line.strip()] # Add to current fragment

            if len(line.strip()) > max_line_len:
                self.line_length = len(line.strip()) # record line length for later printing
                max_line_len = len(line.strip())
        return all_fragments + ["".join(this_fragment)]


    def write_fasta(self, assembled_sequence, output):
        """
        Write FASTA format of sequence, using line lengths learned during reading

        :param assembled_sequence: sequence to write as fasta format
        :param output: where to write fasta format; object with .write() method.
        """

        output.write(">assembled_sequence\n")

        seq_len = len(assembled_sequence)
        i = 0
        while i < seq_len:
            output.write(assembled_sequence[i:min(i+self.line_length, seq_len)] + "\n")
            i += self.line_length



class DeBruijnGraph:
    """ De Bruijn directed multigraph built from a collection of
        strings. Nodes are k-1-mers.  An Edge corresponds to the k-mer that joins
        a left k-1-mer to a right k-1-mer. """

    @staticmethod
    def chop(fragment, k):
        """ Chop fragment into k-mers of given length """
        for i in range(len(fragment) - (k - 1)):
            yield (fragment[i:i + k], fragment[i:i + k - 1], fragment[i + 1:i + k])

    class Node:
        """
        Node representing a k-1 mer.  Keep track of # of
        incoming/outgoing edges so it's easy to check for
        balanced, semi-balanced.
        """

        def __init__(self, k_minus_1_mer):
            self.k_minus_1_mer = k_minus_1_mer
            self.in_degree = 0
            self.out_degree = 0

        def is_semi_balanced(self):
            return abs(self.in_degree - self.out_degree) == 1

        def is_balanced(self):
            return self.in_degree == self.out_degree

        def __hash__(self):
            """ Hash used for mapping nodes to neighbors"""
            return hash(self.k_minus_1_mer)

        def __str__(self):
            return self.k_minus_1_mer

    def __init__(self, fragment_list):
        """
        :param fragment_list: list of fragment strings for assembly
        Initialize De Bruijn Graph, choosing k and building graph
        """

        self.fragment_list = fragment_list

        self.graph = {}  # multimap from nodes to neighbors
        self.nodes = {}  # maps k-1-mers to Node objects
        self.chosen_k = self.choose_k()
        self.build_graph(k= self.chosen_k, fragment_list=fragment_list)

    def choose_k(self):
        """
        Choose k for k-mer De Bruijin Graph assembly.
        Chooses k as ⌊1/2 minimum read length⌋
        :return:
        """

        min_read_length = min(map(len, self.fragment_list))

        chosen_k = int(math.floor(.5 * min_read_length))

        logging.info("Chosen k: {0}".format(chosen_k))

        return chosen_k

    def build_graph(self, k, fragment_list):
        """
        Builds graph by chopping each fragment in fragment_list and linking
        left k-1-mer to right k-1-mer
        """

        seen_k_mers = set()

        for fragment in fragment_list:
            for k_mer, k_minus_1_mer_left, k_minus_1_mer_right in self.chop(fragment, k):
                # if we've already seen this k-mer, continue because we don't want multi-edges. @ijoseph modification
                if k_mer in seen_k_mers:
                    continue
                seen_k_mers.add(k_mer)

                # Find left and right k-1-mers if already part of the graph
                if k_minus_1_mer_left in self.nodes:
                    nodeL = self.nodes[k_minus_1_mer_left]
                else:
                    nodeL = self.nodes[k_minus_1_mer_left] = self.Node(k_minus_1_mer_left)
                if k_minus_1_mer_right in self.nodes:
                    nodeR = self.nodes[k_minus_1_mer_right]
                else:
                    nodeR = self.nodes[k_minus_1_mer_right] = self.Node(k_minus_1_mer_right)

                # Record increased degrees
                nodeL.out_degree += 1
                nodeR.in_degree += 1

                # Add this edge from the left and right k-m-mer
                self.graph.setdefault(nodeL, []).append(nodeR)

        self.assess_graph()

    def assess_graph(self):
        """
        Keep track of head and tail nodes in the case of a graph with Eularian walk (not cycle); set tail and head;
        tally # of nodes balanced, semi-balanced, or unbalanced for later Eulerian check
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
        Find and return fragment of nodes (represented by
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
        """ Assemble by finding the Eulerian walk if possible"""
        assert self.is_eulerian(), \
            "Cannot assemble these fragments with k = {0}; " \
            "try specifying k to be different than the above?".format(self.choose_k())

        # return "".join(dbg.eulerianWalkOrCycle())
        walk = self.get_eulerian_walk_or_cycle() # walk is a list of Nodes

        # Each node contributes one character to the output fragment
        return walk[0] + ''.join(map(lambda x: x[-1], walk[1:]))

    def __str__(self):
        """ String representation"""
        output = "De Bruijn Graph; k= {0} with {1} nodes, {2} edges".format(self.chosen_k,
                                                                            self.number_of_nodes(),
                                                                            self.number_of_edges())
        if self.is_eulerian():
            output += "; Eulerian"
        else:
            output += "; NOT Eulerian"

        return output


def main():
    parser = argparse.ArgumentParser(description="Assemble fragments using de Bruijn graph")
    parser.add_argument("--fragments", type=file, help="Input fragments in FASTA format")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output assembled fragment to file [standard out]")
    namespace = parser.parse_args()

    assembler = Assembler(fragments_fasta=namespace.fragments)
    assembled_sequence = assembler.assemble()
    assembler.write_fasta(assembled_sequence, output=namespace.output)

if __name__ == '__main__':
    main()