#!/usr/bin/env python
'''
Takes a phylogeny and a taxonomy and labels nodes representing higher-order 
taxa in the phylogeny.

Both the phylogeny and the taxonomy should be trees in a format supported by
BioPython.

Usage:

    python phylolabel.py tree_file taxonomy_file 
                         [tree_format] [taxonomy_format] [output_format]

Output will be printed to stdout.
'''
import Bio.Phylo as bp
from time import time
import sys

tree_file, tax_file = sys.argv[1:3]
try: tree_format = sys.argv[3]
except: tree_format = 'newick'
try: tax_format = sys.argv[4]
except: tax_format = 'newick'
try: output_format = sys.argv[5]
except: output_format = 'newick'

def convert_labels(t):
    '''Replace underscores with spaces in taxonomy labels for consistency.'''
    for x in t.find_elements():
        if x.name:
            x.name = x.name.replace('_', ' ')

# read in the tree and taxonomy
tree = bp.read(tree_file, tree_format)
taxonomy = bp.read(tax_file, tax_format)

# standardize labels
convert_labels(tree)
convert_labels(taxonomy)

# cache labels
tree.cache_labels()
taxonomy.cache_labels()

# get all named terminal nodes from phylogeny
tree_species = [sp for sp in tree.get_terminals() if sp.name]

# create two dictionaries, mapping the nodes from the phylogeny to the taxonomy
# and vice versa
tax_to_tree = {}
tree_to_tax = {}
for sp in tree_species:
    x = taxonomy.find_any(sp.name)
    if x:
        tree_to_tax[sp] = x
        tax_to_tree[x] = sp

# walk through species of phylogeny, marking the common ancestors of each
# taxonomic grouping they're a member of
done = set()
for sp in tree_species:
    if not sp in tree_to_tax: continue
    if sp.name in done: continue
    
    tax_sp = tree_to_tax[sp]
    tax_parents = tax_sp.get_parents(False)
    for parent in tax_parents:
        if parent.name in done: continue
        
        fellows = parent.find_elements()
        tree_fellows = (tax_to_tree[x] for x in fellows if x in tax_to_tree)
        group_root = tree.common_ancestor(tree_fellows)
        if not group_root.name:
            group_root.name = parent.name
        
        done.add(parent.name)
        
    done.add(sp.name)

# write output to stdout
print tree.format(output_format)
