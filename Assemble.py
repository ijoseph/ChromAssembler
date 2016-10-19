# -*- coding: utf-8 -*-
"""
@author: Isaac C. Joseph
"""
import argparse, sys, multiprocessing
import graphviz
import itertools


class Assembler(object):
    """
    Performs file handing and uses an assembler to assemble reads
    """

    def __init__(self, fragments_fasta, assembler = "DeBrujinGraph", output = sys.stdout):
        """
        Constructor
        """
        self.assembler = assembler
        self.output = output

        self.sequence_list = self.read_fasta(fragments_fasta)

    def assemble(self):
        if self.assembler == "DeBrujinGraph":
            for k in range(8, 9):
                dbg = DeBruijnGraph2(strIter=self.sequence_list, k=k)

                self.output.write("k ={0}: {1}\n".format(k, dbg.isEulerian()))

                dot_file = dbg.to_dot()
                dot_file.render(filename="Figures/{0}_smol_example.dot".format(k), view=True)


                if dbg.isEulerian():
                    # return "".join(dbg.eulerianWalkOrCycle())
                    walk = dbg.eulerianWalkOrCycle()
                    return walk[0] + ''.join(map(lambda x: x[-1], walk[1:]))




        elif self.assembler == "ShortestCommonSuperstring":
            print ShortestCommonSuperstring.scs(self.sequence_list)

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

    def __init__(self, strIter, k):
        """
        Build de Bruijn multigraph given string iterator and k-mer
            length k
        """
        self.G = {}  # multimap from nodes to neighbors
        self.nodes = {}  # maps k-1-mers to Node objects
        seen_k_mers = set()

        for sequence in strIter:
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
                self.G.setdefault(nodeL, []).append(nodeR)
        # Iterate over nodes; tally # balanced, semi-balanced, neither
        self.nsemi, self.nbal, self.nneither = 0, 0, 0
        # Keep track of head and tail nodes in the case of a graph with
        # Eularian walk (not cycle)
        self.head, self.tail = None, None
        for node in iter(self.nodes.values()):
            if node.is_balanced():
                self.nbal += 1
            elif node.is_semi_balanced():
                if node.in_degree == node.out_degree + 1:
                    self.tail = node
                if node.in_degree == node.out_degree - 1:
                    self.head = node
                self.nsemi += 1
            else:
                self.nneither += 1

    def nnodes(self):
        """ Return # nodes """
        return len(self.nodes)

    def nedges(self):
        """ Return # edges """
        return len(self.G)

    def hasEulerianWalk(self):
        """ Return true iff graph has Eulerian walk. """
        return self.nneither == 0 and self.nsemi == 2

    def hasEulerianCycle(self):
        """ Return true iff graph has Eulerian cycle. """
        return self.nneither == 0 and self.nsemi == 0

    def isEulerian(self):
        """ Return true iff graph has Eulerian walk or cycle """
        # technically, if it has an Eulerian walk
        return self.hasEulerianWalk() or self.hasEulerianCycle()

    def eulerianWalkOrCycle(self):
        """ Find and return sequence of nodes (represented by
            their k-1-mer labels) corresponding to Eulerian walk
            or cycle """
        assert self.isEulerian()
        g = self.G
        if self.hasEulerianWalk():
            g = g.copy()
            g.setdefault(self.tail, []).append(self.head)
        # graph g has an Eulerian cycle
        tour = []
        src = next(iter(g.keys()))  # pick arbitrary starting node

        def __visit(n):
            while len(g[n]) > 0:
                dst = g[n].pop()
                __visit(dst)
            tour.append(n)

        __visit(src)
        tour = tour[::-1][:-1]  # reverse and then take all but last node

        if self.hasEulerianWalk():
            # Adjust node list so that it starts at head and ends at tail
            sti = tour.index(self.head)
            tour = tour[sti:] + tour[:sti]

        # Return node list
        return list(map(str, tour))



class DeBruijnGraph2(DeBruijnGraph):
    def to_dot(self, weights=False):
        """ Return string with graphviz representation.  If 'weights'
            is true, label edges corresponding to distinct k-1-mers
            with weights, instead of drawing separate edges for
            k-1-mer copies. """
        g = graphviz.Digraph(comment='DeBruijn graph')
        for node in iter(self.G.keys()):
            g.node(node.k_minus_1_mer, node.k_minus_1_mer)
        for src, dsts in iter(self.G.items()):
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


class ShortestCommonSuperstring:
    @staticmethod
    def overlap(a, b, min_length=3):
        """ Return length of longest suffix of 'a' matching
            a prefix of 'b' that is at least 'min_length'
            characters long.  If no such overlap exists,
            return 0. """
        start = 0  # start all the way at the left
        while True:
            start = a.find(b[:min_length], start)  # look for b's suffx in a
            if start == -1:  # no more occurrences to right
                return 0
            # found occurrence; check for full suffix/prefix match
            if b.startswith(a[start:]):
                return len(a) - start
            start += 1  # move just past previous match


    @staticmethod
    def scs(ss):
        """ Returns shortest common superstring of given
            strings, which must be the same length """
        shortest_sup = None

        import operator
        from collections import Counter
        from math import factorial
        def npermutations(l):
            num = factorial(len(l))
            mults = Counter(l).values()
            den = reduce(operator.mul, (factorial(v) for v in mults), 1)
            return num / den

        total_pairs = npermutations(ss)

        for (i, ssperm )in enumerate(itertools.permutations(ss)):
            print("Pair {0} of {1} ({2:.70f} %)".format(i+1, total_pairs, (float(i+1)/total_pairs) * 100))
            sup = ssperm[0]  # superstring starts as first string
            for i in range(len(ss) - 1):
                # overlap adjacent strings A and B in the permutation
                olen = ShortestCommonSuperstring.overlap(ssperm[i], ssperm[i + 1], min_length=1)
                # add non-overlapping portion of B to superstring
                sup += ssperm[i + 1][olen:]
            if shortest_sup is None or len(sup) < len(shortest_sup):
                shortest_sup = sup  # found shorter superstring
        return shortest_sup  # return shortest



def main():
    parser = argparse.ArgumentParser(description="description")
    parser.add_argument("--fragments", type=file, help="Input fragments in FASTA format")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output assembled sequence to file [standard out]")
    parser.add_argument('-p', '--numThreads', type=int, default=multiprocessing.cpu_count() / 2,
                        help="Number of threads to run [(# of CPUs)/2]")
    parser.add_argument('--tempFolder', default="temp/", help="temporary folder for files [./temp]")
    namespace = parser.parse_args()

    # assembler = Assembler(fragments_fasta=namespace.fragments, output=namespace.output,
    #                       assembler="ShortestCommonSuperstring")

    assembler = Assembler(fragments_fasta=namespace.fragments)






if __name__ == '__main__':
    main()