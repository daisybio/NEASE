"""Residue contacts between the protein chains of a structure (PDB entries and AlphaFold models), with gemmi."""
from collections import defaultdict

import gemmi

PROTEIN = (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)


def protein_subchains(st):
    return {s for e in st.entities if e.entity_type == gemmi.EntityType.Polymer and e.polymer_type in PROTEIN
            for s in e.subchains}


def build_assembly(st, max_chains):
    '''
    Protein chains of the first biological assembly (the asymmetric unit if there is none or it has more than
    max_chains chains). Chains are named '<author chain>~<copy>' so that copies can be traced back.
    Returns (model, 'assembly' or 'asu').
    '''
    proteins = protein_subchains(st)
    operations = []
    if st.assemblies:
        for generator in st.assemblies[0].generators:
            operations += [(set(generator.subchains) & proteins, op.transform) for op in generator.operators]
    n_chains = sum(len({ch.name for ch in st[0] for r in ch if r.subchain in subs}) for subs, _ in operations)
    kind = 'assembly'
    if not operations or n_chains > max_chains:
        operations, kind = [(proteins, None)], 'asu'
    model = gemmi.Model(1)
    for copy, (subchains, transform) in enumerate(operations, 1):
        part = gemmi.Model(1)
        for chain in st[0]:
            new = gemmi.Chain(f'{chain.name}~{copy}')
            for residue in chain:
                if residue.subchain in subchains:
                    new.add_residue(residue)
            if len(new):
                part.add_chain(new)
        if transform is not None:
            part.transform_pos_and_adp(transform)
        for chain in part:
            model.add_chain(chain)
    return model, kind


def drop_low_confidence(model, min_plddt):
    '''Remove residues whose C-alpha B-factor (pLDDT in AlphaFold models) is below min_plddt.'''
    for chain in model:
        for i in reversed(range(len(chain))):
            ca = chain[i].find_atom('CA', '*')
            if ca is None or ca.b_iso < min_plddt:
                del chain[i]


def residue_contacts(model, distance):
    '''
    Residue pairs of different chains with heavy atoms within distance.
    Returns {(chain_1, chain_2): {(label_seq_1, label_seq_2), ...}} with chain_1 < chain_2.
    '''
    search = gemmi.NeighborSearch(model, gemmi.UnitCell(), max(distance, 5)).populate(include_h=False)
    contact_search = gemmi.ContactSearch(distance)
    contact_search.ignore = gemmi.ContactSearch.Ignore.SameChain
    pairs = defaultdict(set)
    for contact in contact_search.find_contacts(search):
        a, b = contact.partner1, contact.partner2
        sa, sb = a.residue.label_seq, b.residue.label_seq
        if sa is None or sb is None:
            continue
        if a.chain.name < b.chain.name:
            pairs[(a.chain.name, b.chain.name)].add((sa, sb))
        else:
            pairs[(b.chain.name, a.chain.name)].add((sb, sa))
    return pairs
