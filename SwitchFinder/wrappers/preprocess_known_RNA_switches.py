import argparse
import IO as IO
import numpy as np

def handler(raw_args = None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="the fasta file with the target sequences", type=str)
    parser.add_argument("-o", help="output filename", type=str)
    parser.add_argument("--length", help="fragment length", type=int)


    parser.set_defaults(
                        f = "",
                        o = "",
                        length = 186
                        )
    args = parser.parse_args(raw_args)
    return args


# altschulEriksonDinuclShuffle.py
# P. Clote, Oct 2003
# NOTE: One cannot use function "count(s,word)" to count the number
# of occurrences of dinucleotide word in string s, since the built-in
# function counts only nonoverlapping words, presumably in a left to
# right fashion.
# code from here: https://github.com/wassermanlab/BiasAway/blob/master/altschulEriksonDinuclShuffle.py

import sys, string, random


def computeCountAndLists(s):
    # WARNING: Use of function count(s,'UU') returns 1 on word UUU
    # since it apparently counts only nonoverlapping words UU
    # For this reason, we work with the indices.

    # Initialize lists and mono- and dinucleotide dictionaries
    List = {}  # List is a dictionary of lists
    List['A'] = [];
    List['C'] = [];
    List['G'] = [];
    List['T'] = [];
    nuclList = ["A", "C", "G", "T"]
    s = s.upper()
    s = s.replace("T", "T")
    nuclCnt = {}  # empty dictionary
    dinuclCnt = {}  # empty dictionary
    for x in nuclList:
        nuclCnt[x] = 0
        dinuclCnt[x] = {}
        for y in nuclList:
            dinuclCnt[x][y] = 0

    # Compute count and lists
    nuclCnt[s[0]] = 1
    nuclTotal = 1
    dinuclTotal = 0
    for i in range(len(s) - 1):
        x = s[i];
        y = s[i + 1]
        List[x].append(y)
        nuclCnt[y] += 1;
        nuclTotal += 1
        dinuclCnt[x][y] += 1;
        dinuclTotal += 1
    assert (nuclTotal == len(s))
    assert (dinuclTotal == len(s) - 1)
    return nuclCnt, dinuclCnt, List


def chooseEdge(x, dinuclCnt):
    numInList = 0
    for y in ['A', 'C', 'G', 'T']:
        numInList += dinuclCnt[x][y]
    z = random.random()
    denom = dinuclCnt[x]['A'] + dinuclCnt[x]['C'] + dinuclCnt[x]['G'] + dinuclCnt[x]['T']
    numerator = dinuclCnt[x]['A']
    if z < float(numerator) / float(denom):
        dinuclCnt[x]['A'] -= 1
        return 'A'
    numerator += dinuclCnt[x]['C']
    if z < float(numerator) / float(denom):
        dinuclCnt[x]['C'] -= 1
        return 'C'
    numerator += dinuclCnt[x]['G']
    if z < float(numerator) / float(denom):
        dinuclCnt[x]['G'] -= 1
        return 'G'
    dinuclCnt[x]['T'] -= 1
    return 'T'


def connectedToLast(edgeList, nuclList, lastCh):
    D = {}
    for x in nuclList: D[x] = 0
    for edge in edgeList:
        a = edge[0];
        b = edge[1]
        if b == lastCh: D[a] = 1
    for i in range(2):
        for edge in edgeList:
            a = edge[0];
            b = edge[1]
            if D[b] == 1: D[a] = 1
    ok = 0
    for x in nuclList:
        if x != lastCh and D[x] == 0: return 0
    return 1


def eulerian(s):
    nuclCnt, dinuclCnt, List = computeCountAndLists(s)
    # compute nucleotides appearing in s
    nuclList = []
    for x in ["A", "C", "G", "T"]:
        if x in s: nuclList.append(x)
    # compute numInList[x] = number of dinucleotides beginning with x
    numInList = {}
    for x in nuclList:
        numInList[x] = 0
        for y in nuclList:
            numInList[x] += dinuclCnt[x][y]
    # create dinucleotide shuffle L
    firstCh = s[0]  # start with first letter of s
    lastCh = s[-1]
    edgeList = []
    for x in nuclList:
        if x != lastCh: edgeList.append([x, chooseEdge(x, dinuclCnt)])
    ok = connectedToLast(edgeList, nuclList, lastCh)
    return ok, edgeList, nuclList, lastCh


def shuffleEdgeList(L):
    n = len(L);
    barrier = n
    for i in range(n - 1):
        z = int(random.random() * barrier)
        tmp = L[z]
        L[z] = L[barrier - 1]
        L[barrier - 1] = tmp
        barrier -= 1
    return L


def dinuclShuffle(s):
    ok = 0
    while not ok:
        ok, edgeList, nuclList, lastCh = eulerian(s)
    nuclCnt, dinuclCnt, List = computeCountAndLists(s)

    # remove last edges from each vertex list, shuffle, then add back
    # the removed edges at end of vertex lists.
    for [x, y] in edgeList: List[x].remove(y)
    for x in nuclList: shuffleEdgeList(List[x])
    for [x, y] in edgeList: List[x].append(y)

    # construct the eulerian path
    L = [s[0]];
    prevCh = s[0]
    for i in range(len(s) - 2):
        ch = List[prevCh][0]
        L.append(ch)
        del List[prevCh][0]
        prevCh = ch
    L.append(s[-1])
    t = "".join(L)
    return t


def shuffle_by_rolling_window(fasta_filename_loc, output_fasta_filename_loc, shuffling_number):
    with open(fasta_filename_loc, 'r') as rf:
        with open(output_fasta_filename_loc, 'w') as wf:
            bigline = rf.read()
            bigarray = bigline.split('\n>')
            bigarray = [x for x in bigarray if x != '']
            for entry in bigarray:
                index = entry.find('\n')
                annotation = entry[0:index]
                annotation = annotation.replace('>', '')
                sequence = entry[index:].replace('\n', '')

                info_string = '>%s_original\n%s\n' % (annotation, sequence)
                wf.write(info_string)

                for i in range(shuffling_number):
                    new_sequence = dinuclShuffle(sequence)
                    info_string = '>%s_shuffling_%d\n%s\n' % (annotation, i, new_sequence)
                    wf.write(info_string)

def shuffle_and_write_sequences(fasta_dict, output_file):
    with open(output_file, 'w') as wf:
        for name in fasta_dict:
            info_string = '>%s_original\n%s\n' % (name, fasta_dict[name])
            wf.write(info_string)
            new_sequence = dinuclShuffle(fasta_dict[name])
            info_string = '>%s_shuffling_1\n%s\n' % (name, new_sequence)
            wf.write(info_string)


def main(raw_args = None):
    args = handler(raw_args)
    fasta_dict = IO.read_fasta(args.f)
    # filter by length
    print(len(fasta_dict))
    fasta_dict = {x : fasta_dict[x] for x in fasta_dict if len(fasta_dict[x]) <= args.length }
    print(len(fasta_dict))
    # shuffle
    shuffle_and_write_sequences(fasta_dict, args.o)


if __name__ == "__main__":
    main()