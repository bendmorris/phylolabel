Takes a phylogeny and a taxonomy and labels nodes representing higher-order 
taxa in the phylogeny.

Both the phylogeny and the taxonomy should be trees in a format supported by
BioPython.

Usage:

    python phylolabel.py tree_file taxonomy_file 
                         [tree_format] [taxonomy_format] [output_format]

Output will be printed to stdout.

Note: this currently relies on functionality in my development fork of BioPython 
which hasn't been pulled into the main branch yet. You can get my fork here:

    https://github.com/bendmorris/biopython