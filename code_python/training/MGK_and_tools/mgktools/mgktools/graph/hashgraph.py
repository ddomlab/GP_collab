#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import networkx as nx
from graphdot import Graph
from graphdot.graph._from_networkx import _from_networkx
from rxntools.reaction import *
from mgktools.graph.from_rdkit import _from_rdkit, default_config, AtomBondFeaturesConfig


CWD = os.path.dirname(os.path.abspath(__file__))


class HashGraph(Graph):
    """
    Class HashGraph is hashable Graph, where the Class Graph in GraphDot is unhashable.

    For each HashGraph, a string is needed that defined by the user as its hash value. Usually I use SMILES string.

    HashGraph is sortable, which is important if you want to calculate the graph kernel matrix only once at the
    beginning.
    """

    def __init__(self, hash=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hash = hash

    def __eq__(self, other):
        if self.hash == other.hash:
            return True
        else:
            return False

    def __lt__(self, other):
        if self.hash < other.hash:
            return True
        else:
            return False

    def __gt__(self, other):
        if self.hash > other.hash:
            return True
        else:
            return False

    def __hash__(self):
        return hash(self.hash)

    @classmethod
    def from_smiles(self, smiles, hash=None, atom_bond_features_config=default_config):
        mol = Chem.MolFromSmiles(smiles)
        g = self.from_rdkit(mol, hash or smiles, atom_bond_features_config)
        return g

    @classmethod
    def from_rdkit(cls, mol, hash=None, atom_bond_features_config=default_config):
        atom_bond_features_config.preprocess(mol)
        g = _from_rdkit(cls, mol, atom_bond_features_config)
        g.normalize_concentration()
        g.hash = hash
        # g = g.permute(rcm(g))
        return g

    @classmethod
    def agent_from_reaction_smarts(cls, reaction_smarts, hash,
                                   atom_bond_features_config=default_config):
        cr = ChemicalReaction(reaction_smarts)
        return cls.agent_from_cr(cr, hash, atom_bond_features_config=atom_bond_features_config)

    @classmethod
    def agent_from_cr(cls, cr, hash=None, atom_bond_features_config=default_config):
        if len(cr.agents) == 0:
            return HashGraph.from_smiles(smiles='[He]',
                                         hash=hash,
                                         atom_bond_features_config=atom_bond_features_config)

        agents = HashGraph.from_rdkit(mol=cr.agents[0],
                                      hash=hash,
                                      atom_bond_features_config=atom_bond_features_config).to_networkx()
        for mol in cr.agents[1:]:
            g = HashGraph.from_rdkit(mol=mol,
                                     hash=hash,
                                     atom_bond_features_config=atom_bond_features_config).to_networkx()
            agents = nx.disjoint_union(agents, g)
        g = _from_networkx(cls, agents)
        g.hash = hash
        return g

    @classmethod
    def from_reaction_smarts(cls, reaction_smarts, hash):
        cr = ChemicalReaction(reaction_smarts)
        return cls.from_cr(cr, hash)

    @classmethod
    def from_cr(cls, cr, hash):
        atom_bond_features_config = AtomBondFeaturesConfig(reaction_center=cr.ReactingAtomsMN,
                                                           reactant_or_product='reactant')
        reaction = HashGraph.from_rdkit(mol=cr.reactants[0],
                                        atom_bond_features_config=atom_bond_features_config).to_networkx()
        for reactant in cr.reactants[1:]:
            g = HashGraph.from_rdkit(mol=reactant,
                                     atom_bond_features_config=atom_bond_features_config).to_networkx()
            reaction = nx.disjoint_union(reaction, g)

        atom_bond_features_config = AtomBondFeaturesConfig(reaction_center=cr.ReactingAtomsMN,
                                                           reactant_or_product='product')
        for product in cr.products:
            g = HashGraph.from_rdkit(mol=product,
                                     atom_bond_features_config=atom_bond_features_config).to_networkx()
            reaction = nx.disjoint_union(reaction, g)

        g = _from_networkx(cls, reaction)
        g.hash = hash
        if g.nodes.to_pandas()['ReactingCenter'].max() <= 0:
            raise RuntimeError(f'No reacting atoms are found in reactants:　'
                               f'{cr.reaction_smarts}')
        if g.nodes.to_pandas()['ReactingCenter'].min() >= 0:
            raise RuntimeError(f'No reacting atoms are found in products:　'
                               f'{cr.reaction_smarts}')
        return g

    @classmethod
    def from_reaction_template(cls, template_smarts):
        template = ReactionTemplate(template_smarts)
        _rdkit_config = AtomBondFeaturesConfig(reaction_center=template.ReactingAtomsMN,
                                               reactant_or_product='reactant',
                                               IsSanitized=False,
                                               set_morgan_identifier=False)
        reaction = Graph.from_rdkit(
            template.reactants[0], _rdkit_config).to_networkx()
        for reactant in template.reactants[1:]:
            g = Graph.from_rdkit(reactant, _rdkit_config).to_networkx()
            reaction = nx.disjoint_union(reaction, g)

        _rdkit_config = AtomBondFeaturesConfig(reaction_center=template.ReactingAtomsMN,
                                               reactant_or_product='product',
                                               IsSanitized=False,
                                               set_morgan_identifier=False)
        for product in template.products:
            g = Graph.from_rdkit(product, _rdkit_config).to_networkx()
            reaction = nx.disjoint_union(reaction, g)

        g = _from_networkx(cls, reaction)
        if g.nodes.to_pandas()['ReactingCenter'].max() <= 0:
            raise RuntimeError(f'No reacting atoms are found in reactants:　'
                               f'{template_smarts}')
        if g.nodes.to_pandas()['ReactingCenter'].min() >= 0:
            raise RuntimeError(f'No reacting atoms are found in products:　'
                               f'{template_smarts}')
        return g

    @classmethod
    def reactant_from_reaction_smarts(cls, reaction_smarts, hash):
        cr = ChemicalReaction(reaction_smarts)

        _rdkit_config = AtomBondFeaturesConfig(reaction_center=cr.ReactingAtomsMN,
                                               reactant_or_product='reactant')
        reaction = HashGraph.from_rdkit(
            cr.reactants[0], '1', _rdkit_config).to_networkx()
        for reactant in cr.reactants[1:]:
            g = HashGraph.from_rdkit(reactant, '1', _rdkit_config). \
                to_networkx()
            reaction = nx.disjoint_union(reaction, g)

        g = _from_networkx(cls, reaction)
        g.hash = hash
        return g

    @classmethod
    def product_from_reaction_smarts(cls, reaction_smarts, hash):
        cr = ChemicalReaction(reaction_smarts)

        _rdkit_config = AtomBondFeaturesConfig(reaction_center=cr.ReactingAtomsMN,
                                               reactant_or_product='reactant')
        reaction = HashGraph.from_rdkit(
            cr.products[0], '1', _rdkit_config).to_networkx()
        for product in cr.products[1:]:
            g = HashGraph.from_rdkit(product, '1', _rdkit_config). \
                to_networkx()
            reaction = nx.disjoint_union(reaction, g)

        g = _from_networkx(cls, reaction)
        g.hash = hash
        return g

    def update_concentration(self, concentration: float):
        self.nodes['Concentration'] *= concentration

    def normalize_concentration(self):
        sum_c = self.nodes['Concentration'].sum()
        self.nodes['Concentration_norm'] = self.nodes['Concentration'] / sum_c
