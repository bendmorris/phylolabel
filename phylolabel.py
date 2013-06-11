#!/usr/bin/env python
'''
Takes a phylogeny and a taxonomy and labels nodes representing higher-order 
taxa in the phylogeny.

Both the phylogeny and the taxonomy should be trees in a format supported by
BioPython.

For help on usage, run `python phylolabel.py -h`
'''
import Bio.Phylo as bp


def convert_labels(t):
    '''Replace underscores with spaces in taxonomy labels for consistency.'''
    for x in t.find_elements():
        if x.name:
            x.name = x.name.replace('_', ' ')

def label_tree(phylogeny, taxonomy, tax_root=None):
    '''Add taxonomic labels to phylogeny. Operates on the phylogeny in-place.
    Returns a set of labels common to both the tree and taxonomy.
    
    tax_root, if provided, should be the name of a node in the taxonomy; this 
    will be used to take a subset of the taxonomy, avoiding taxonomic homonym 
    issues.
    '''

    # standardize labels
    convert_labels(phylogeny)
    convert_labels(taxonomy)
    
    # subset the taxnoomy if necessary
    if tax_root:
        top_node = taxonomy.find_any(tax_root)
        if top_node:
            if hasattr(top_node, '_parent'):
                top_node._parent.clades.remove(top_node)
            taxonomy = bp.BaseTree.Tree(root=top_node)

    # index labels
    phylogeny.index_labels()
    taxonomy.index_labels()

    # get all named terminal nodes from phylogeny
    phylogeny_species = [sp for sp in phylogeny.get_terminals() if sp.name]

    # create two dictionaries, mapping the nodes from the phylogeny to the taxonomy
    # and vice versa
    taxonomy_to_phylogeny = {}
    phylogeny_to_taxonomy = {}
    for sp in phylogeny_species:
        x = taxonomy.find_any(sp.name)
        if x:
            phylogeny_to_taxonomy[sp] = x
            taxonomy_to_phylogeny[x] = sp
    
    # walk through species of phylogeny, marking the common ancestors of each
    # taxonomic grouping they're a member of
    done = set()
    for sp in phylogeny_species:
        if not sp in phylogeny_to_taxonomy: continue
        if sp.name in done: continue
        
        tax_sp = phylogeny_to_taxonomy[sp]
        tax_parents = tax_sp.get_parents(False)
        for parent in tax_parents:
            if parent.name in done: continue
            
            fellows = parent.find_elements()
            fellows = (taxonomy_to_phylogeny[x] for x in fellows if x in taxonomy_to_phylogeny)
            group_root = phylogeny.common_ancestor(fellows)
            
            if not group_root.name:
                # the node is currently unlabeled, so label it
                group_root.name = parent.name

            else:
                # the node was already labeled, so split it into two nodes
                new_clade = bp.BaseTree.Clade(name=parent.name, branch_length=0)
                
                # how many nodes are there in the tree already separated by
                # branches of length 0?
                old_otus = [taxonomy.find_any(group_root.name)]
                for x in group_root.get_parents(False):
                    if x.branch_length == 0 and x.name:
                        old_otus.append(taxonomy.find_any(x.name))
                    else:
                        break
                new_otu = parent
                
                placed = False
                for (n, old_otu) in enumerate(old_otus):
                    if old_otu in new_otu.get_parents(False):
                        # existing node should be parent of new node
                        for _ in range(n):
                            group_root = group_root._parent
                        
                        for child in [x for x in group_root.clades]:
                            group_root.clades.remove(child)
                            new_clade.clades.append(child)
                            
                        group_root.clades.append(new_clade)
                        placed = True
                        break
                
                if not placed:
                    # new node should be parent of existing nodes
                    for x in group_root.get_parents(False):
                        if x.branch_length == 0 and x.name:
                            group_root = x
                        else: break

                    if group_root is phylogeny.root:
                        phylogeny.root = new_clade
                    else:
                        old_parent = group_root._parent
                        old_parent.clades.remove(group_root)
                        old_parent.clades.append(new_clade)

                    new_clade.clades.append(group_root)
            
            done.add(parent.name)
            
        done.add(sp.name)

    return done
    
    
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('phylogeny_file', help='path to the phylogeny')
    parser.add_argument('taxonomy_file',  help='path to the taxonomy')
    parser.add_argument('-p', '--phylogeny_format', nargs='?', default='newick',
                        help='phylogeny format (%s)' % 
                        (','.join(bp._io.supported_formats.keys())))
    parser.add_argument('-t', '--taxonomy_format', nargs='?', default='newick', help='taxonomy format')
    parser.add_argument('-o', '--output_format', nargs='?', default='newick', help='output format')
    parser.add_argument('-r', '--root', nargs='?', default=None, help='name of OTU to use as root of taxonomy')

    args = parser.parse_args()
    
    # read in the tree and taxonomy
    phylogeny = bp.read(args.phylogeny_file, args.phylogeny_format)
    taxonomy = bp.read(args.taxonomy_file, args.taxonomy_format)
    
    label_tree(phylogeny, taxonomy, tax_root=args.root)
    
    # write output to stdout
    print phylogeny.format(args.output_format)
