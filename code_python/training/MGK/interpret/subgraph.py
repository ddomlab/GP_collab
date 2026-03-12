#!/usr/bin/env python
# -*- coding: utf-8 -*-
import itertools
import rdkit.Chem.AllChem as Chem
import networkx as nx


def inter(a, b):
    return list(set(a) & set(b))


def get_node_name(idx, g, clusters, non_cluster_atoms):
    if idx in non_cluster_atoms:
        return idx
    else:
        for i, cluster in enumerate(clusters):
            if idx in cluster:
                return 'cluster_%i' % i


def cluster_fusion(idx_clusters):
    for i in range(len(idx_clusters)):
        for j in range(i, len(idx_clusters)):
            if idx_clusters[i] & idx_clusters[j]:
                idx_clusters[i] = idx_clusters[i] | idx_clusters[j]
                del idx_clusters[j]
                return cluster_fusion(idx_clusters)
    return idx_clusters


def get_strict_atom_smarts(atom):
    ''' Copy from Rdchiral.
    For an RDkit atom object, generate a SMARTS pattern that
    matches the atom as strictly as possible
    '''
    ChiralType = Chem.rdchem.ChiralType

    symbol = atom.GetSmarts()
    if atom.GetSymbol() == 'H':
        symbol = '[#1]'

    if '[' not in symbol:
        symbol = '[' + symbol + ']'

    # Explicit stereochemistry - *before* H
    if atom.GetChiralTag() != ChiralType.CHI_UNSPECIFIED:
        if '@' not in symbol:
            # Be explicit when there is a tetrahedral chiral tag
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                tag = '@'
            elif atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                tag = '@@'
            if ':' in symbol:
                symbol = symbol.replace(':', ';{}:'.format(tag))
            else:
                symbol = symbol.replace(']', ';{}]'.format(tag))

    if 'H' not in symbol:
        H_symbol = 'H{}'.format(atom.GetTotalNumHs())
        # Explicit number of hydrogens: include "H0" when no hydrogens present
        if ':' in symbol:  # stick H0 before label
            symbol = symbol.replace(':', ';{}:'.format(H_symbol))
        else:
            symbol = symbol.replace(']', ';{}]'.format(H_symbol))

    # Explicit degree
    if ':' in symbol:
        symbol = symbol.replace(':', ';D{}:'.format(atom.GetDegree()))
    else:
        symbol = symbol.replace(']', ';D{}]'.format(atom.GetDegree()))

    # Explicit formal charge
    if '+' not in symbol and '-' not in symbol:
        charge = atom.GetFormalCharge()
        charge_symbol = '+' if (charge >= 0) else '-'
        charge_symbol += '{}'.format(abs(charge))
        if ':' in symbol:
            symbol = symbol.replace(':', ';{}:'.format(charge_symbol))
        else:
            symbol = symbol.replace(']', ';{}]'.format(charge_symbol))

    # Explicit ring infomation
    if 'R' not in symbol:
        ring_symbol = 'R' if atom.IsInRing() else '!R'
        if ':' in symbol:
            symbol = symbol.replace(':', ';{}:'.format(ring_symbol))
        else:
            symbol = symbol.replace(']', ';{}]'.format(ring_symbol))

    return symbol


def get_fragments_smarts(mol: Chem.Mol,
                         nb_nodes_min: int = 1,
                         nb_nodes_max: int = 5,
                         StrictSmarts: bool = True):
    # get clusters information, here we group each ring as a cluster.

    idx_clusters = [set(idx_cluster) for idx_cluster in mol.GetRingInfo().AtomRings()]
    idx_clusters = cluster_fusion(idx_clusters)

    idx_non_cluster_atoms = list(range(mol.GetNumAtoms()))
    for cluster in idx_clusters:
        for idx in cluster:
            idx_non_cluster_atoms.remove(idx)

    # Create networkx graph, the clusters will be coarse grained into nodes.
    g = nx.Graph()
    # create nodes
    for i, cluster in enumerate(idx_clusters):
        g.add_node('cluster_%i' % i)
        g.nodes['cluster_%i' % i]['atomContribution'] = 0.
        g.nodes['cluster_%i' % i]['atoms_idx'] = list(cluster)
        for idx in cluster:
            atom = mol.GetAtomWithIdx(idx)
            g.nodes['cluster_%i' % i]['atomContribution'] += float(atom.GetProp('atomNote'))
    for idx in idx_non_cluster_atoms:
        atom = mol.GetAtomWithIdx(idx)
        g.add_node(idx)
        g.nodes[idx]['atomContribution'] = float(atom.GetProp('atomNote'))
        g.nodes[idx]['atoms_idx'] = [idx]
    # create edges

    for bond in mol.GetBonds():
        bidx = bond.GetBeginAtomIdx()
        bnode = get_node_name(bidx, g, idx_clusters, idx_non_cluster_atoms)
        eidx = bond.GetEndAtomIdx()
        enode = get_node_name(eidx, g, idx_clusters, idx_non_cluster_atoms)
        if bnode != enode:
            g.add_edge(bnode, enode)
    # find all possible subgraph
    idx_subgraph = []
    for nb_nodes in range(nb_nodes_min, nb_nodes_max):
        for selected_nodes in itertools.combinations(g, nb_nodes):
            subgraph = g.subgraph(selected_nodes)
            if nx.is_connected(subgraph):
                idx = []
                for value in list(subgraph.nodes.values()):
                    idx.extend(value['atoms_idx'])
                idx_subgraph.append(idx)
    smarts_fragments = []
    for idx in idx_subgraph:
        if StrictSmarts:
            symbols = [atom.GetSmarts() for atom in mol.GetAtoms()]
            for i in idx:
                symbols[i] = get_strict_atom_smarts(mol.GetAtomWithIdx(i))
            smarts = Chem.MolFragmentToSmiles(mol,
                                              atomsToUse=idx,
                                              atomSymbols=symbols,
                                              allHsExplicit=True,
                                              isomericSmiles=True,
                                              allBondsExplicit=True)
        else:
            smarts = Chem.MolFragmentToSmiles(mol,
                                              atomsToUse=idx,
                                              allHsExplicit=True,
                                              isomericSmiles=True,
                                              allBondsExplicit=True)
        if smarts not in smarts_fragments:
            smarts_fragments.append(smarts)
    return smarts_fragments
