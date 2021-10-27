#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Charlotte des Mares de Trébons"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Charlotte des Mares de Trébons"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Charlotte des Mares de Trébons"
__email__ = "chadmdt@outlook.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file) as file:
        i = 0
        for line in file:
            i += 1
            if i%2==0 and i%4!=0:
                yield line.rstrip()

def cut_kmer(read, kmer_size):
    for i in range(0,len(read)-kmer_size+1):
        yield read[i:i+kmer_size]

def build_kmer_dict(fastq_file, kmer_size):
    dict = {}
    for line in read_fastq(fastq_file):
        for kmer in cut_kmer(line, kmer_size):
            occ = dict.pop(kmer, None)
            if occ != None:
                dict[kmer] = occ + 1
            else:
                dict[kmer] = 1
    return dict

def build_graph(kmer_dict):
    digraph = nx.DiGraph()
    for kmer in kmer_dict.keys():
        digraph.add_edge(kmer[:-1],kmer[1:],weight=kmer_dict[kmer])
    return digraph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node is True and delete_sink_node is True:
            graph.remove_nodes_from(path)
        elif delete_entry_node is True and delete_sink_node is False:
            graph.remove_nodes_from(path[:-1])
        elif delete_entry_node is False and delete_sink_node is True:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph

def std(data):
    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    if std(weight_avg_list) > 0:
        best_index = weight_avg_list.index(max(weight_avg_list))
        del path_list[best_index]
        return remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
    elif std(path_length) > 0:
        best_index = path_length.index(max(path_length))
        del path_list[best_index]
        return remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
    else:
        best_index = randint(0,len(path_list)-1)
        del path_list[best_index]
        return remove_paths(graph,path_list,delete_entry_node,delete_sink_node)


def path_average_weight(graph, path):
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = nx.all_simple_paths(graphe,ancestor_node,descendant_node)
    path_length = []
    weight_avg_list = []
    for path in path_list:
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight)
    return select_best_path(graph, path_list, weight_avg_list)

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    node_list = []
    for node in graph.nodes():
        if(len([p for p in graph.predecessors(node)]) == 0):
            node_list.append(node)
    return node_list

def get_sink_nodes(graph):
    node_list = []
    for node in graph.nodes():
        if(len([p for p in graph.successors(node)]) == 0):
            node_list.append(node)
    return node_list

def get_contigs(graph, starting_nodes, ending_nodes):
    result = []
    for st_node in starting_nodes:
        for en_node in ending_nodes:
            if(nx.has_path(graph,st_node,en_node) is False):
                break
            for n in nx.all_simple_paths(graph,st_node,en_node):
                contig = n[0]
                for i in range(1,len(n)):
                    kmer = n[i]
                    contig += kmer[-1]
                result.append([contig,len(contig)])
    return result

def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as file:
        for i,element in enumerate(contigs_list):
            file.write('>contig_'+str(i)+' len='+str(element[1])+'\n')
            file.write(fill(element[0]))
            file.write('\n')


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
