#!/usr/bin/env python3
import sys
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def gfa_to_circular_fasta(gfa_path, fasta_output_prefix):
    segments = {}
    with open(gfa_path) as f:
        for line in f:
            if line.startswith('S\t'):
                _, seg_id, seq, *_ = line.strip().split('\t')
                segments[seg_id] = seq
    if not segments:
        raise RuntimeError("No segments (S lines) found in GFA.")

    UG = nx.MultiGraph()
    UG.add_nodes_from(segments)
    links = []
    with open(gfa_path) as f:
        for idx, line in enumerate(f):
            if line.startswith('L\t'):
                _, a, ao, b, bo, cigar, *_ = line.strip().split('\t')
                UG.add_edge(a, b, key=idx, a_ori=ao, b_ori=bo, cigar=cigar)
                links.append((a, ao, b, bo, cigar))

    comp_id = 0
    for comp in nx.connected_components(UG):
        comp_id += 1
        sub_nodes = list(comp)
        subG = UG.subgraph(sub_nodes).copy()
        if subG.number_of_edges() == 0:
            print(f"[INFO] Component {comp_id} has no edges, skipping.")
            continue

        longest = max(sub_nodes, key=lambda x: len(segments[x]))
        if not nx.is_eulerian(subG):
            print(f"[WARN] Component {comp_id} is not Eulerian, skipping.")
            continue

        try:
            circuit = list(nx.eulerian_circuit(subG, source=longest, keys=True))
        except TypeError:
            circuit = list(nx.eulerian_circuit(subG, source=longest))

        path = []
        first = circuit[0]
        if len(first) == 3:
            u0, v0, k0 = first
        else:
            u0, v0 = first
            k0 = list(subG[u0][v0].keys())[0]
        data0 = subG[u0][v0][k0]
        ori0 = data0['a_ori'] if u0 == data0.get('a') else data0['b_ori']
        path.append(f"{u0}{ori0}")

        for edge in circuit:
            if len(edge) == 3:
                u, v, key = edge
            else:
                u, v = edge
                key = list(subG[u][v].keys())[0]
            data = subG[u][v][key]
            ori = data['b_ori'] if v == data.get('b') else data['a_ori']
            if v == longest and edge == circuit[-1]:
                break
            path.append(f"{v}{ori}")

        seq_pieces = []
        for entry in path:
            seg_id, ori = entry[:-1], entry[-1]
            seq = Seq(segments[seg_id])
            if ori == '-':
                seq = seq.reverse_complement()
            seq_pieces.append(str(seq))
        full_seq = ''.join(seq_pieces)

        header = f"comp{comp_id} {';'.join(path)}"
        rec = SeqRecord(Seq(full_seq), id=header, description="")
        out_path = f"{fasta_output_prefix}_comp{comp_id}.fa"
        SeqIO.write(rec, out_path, "fasta")
        print(f"Component {comp_id} written to {out_path}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python3 a.py assembly_graph.gfa output_prefix\n")
        sys.exit(1)
    gfa_to_circular_fasta(sys.argv[1], sys.argv[2])

