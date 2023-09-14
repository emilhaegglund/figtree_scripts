# FigTree-scripts
Collection of scripts for preparing tree-files for visualisation in FigTree.


## Getting started
Clone this repo
```
git clone https://github.com/emilhaegglund/figtree_scripts.git
cd figtree_scripts
```
Create new environment and install required packages using conda.
```
conda create env -f environment.yaml
```

## Usage
Using the collapse_phylogeny.py script it is possible to collapse monophyletic nodes and when the output file is open in FigTree collapsed nodes will be displayed with triangles. Also, leaf and triangles will be annotated based on the given annotation file.
```
python collapse_phylogeny.py \
    --tree examples/collapse_phylogeny.nwk \
    --annotation examples/collapse_annotation.tsv \
    --settings examples/figtree_settings.txt \
    --column  class \
    --output collapased_phylogeny.nxs
```

This is an example of an annotation file with default column names. The first column is required to contain the accessions of the taxa in
the phylogeny.

| taxa         | domain   | phylum          | class          | order            | family            | genus    | species              |
|--------------|----------|-----------------|----------------|------------------|-------------------|----------|----------------------|
| HGE70390.1   | Bacteria | Poribacteria    | DTPJ01         | DTPJ01           | DTPJ01            | DTPJ01   | DTPJ01 sp011334435   |
| MBS4020802.1 | Bacteria | Firmicutes_D    | Dethiobacteria | Dethiobacterales | Dethiobacteraceae | JAGXRN01 | JAGXRN01 sp018335875 |
| GAB61033.1   | Bacteria | Planctomycetota | Brocadiae      | Brocadiales      | Brocadiaceae      | Jettenia | Jettenia caeni       |



If there are additional columns in the annotation, they will be omitted.

A custom set of columns can be defined using the `--header` argument with a comma-separated list of column values to use. The annotation line will be in the same order as the columns are given to the header option.
```
--header genus,species,strain
```

### Coming
Scripts to color labels of a phylogeny
