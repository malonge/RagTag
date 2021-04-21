![logo_text](logo/logo_banner.svg)

## RagTag 

[![DOI](https://zenodo.org/badge/242898323.svg)](https://zenodo.org/badge/latestdoi/242898323) [![RELEASE](https://img.shields.io/github/v/release/malonge/RagTag?color=EE7733)](https://github.com/malonge/RagTag/releases/tag/v2.0.0) [![CONDA](https://img.shields.io/conda/dn/bioconda/ragtag?color=009988&label=conda)](https://anaconda.org/bioconda/ragtag) [![GitHub](https://img.shields.io/github/license/malonge/RagTag?color=CC3311)](https://github.com/malonge/RagTag/blob/master/LICENSE)

RagTag is a collection of command-line utilities for improving modern genome assemblies. Tasks include:

- Homology-based sequence [correction](https://github.com/malonge/RagTag/wiki/correct)
- Homology-based sequence [scaffolding](https://github.com/malonge/RagTag/wiki/scaffold)
- Homology-based continuous scaffolding and gap-filling ([patching](https://github.com/malonge/RagTag/wiki/patch))
- Scaffold [merging](https://github.com/malonge/RagTag/wiki/merge)
  
Ragtag also provides a [collection of command line utilities](https://github.com/malonge/RagTag/wiki/Usage) for working with common genome assembly file formats.

## Getting Started

```bash
# install with conda
conda install -c bioconda ragtag

# correct contigs
ragtag.py correct ref.fasta query.fasta

# scaffold contigs
ragtag.py scaffold ref.fa ragtag_output/query.corrected.fasta

# scaffold with multiple references
ragtag.py scaffold -o out_1 ref1.fasta query.fasta
ragtag.py scaffold -o out_2 ref2.fasta query.fasta
ragtag.py merge query.fasta out_*/*.agp

# use Hi-C to resolve conflicts
ragtag.py merge -b hic.bam query.fasta out_*/*.agp

# make joins and fill gaps in target.fa using sequences from query.fa
ragtag.py patch target.fa query.fa
```

## Docs
Please see the [Wiki](https://github.com/malonge/RagTag/wiki) for detailed documentation.

## Dependencies
- [Minimap2](https://github.com/lh3/minimap2), [Unimap](https://github.com/lh3/unimap), or [Nucmer](http://mummer.sourceforge.net/)
- Python 3 (with the following auto-installed packages)
    - numpy
    - intervaltree
    - pysam
    - networkx
    
## Citation

Alonge, Michael, et al. ["RaGOO: fast and accurate reference-guided scaffolding of draft genomes."](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1829-6) Genome biology 20.1 (2019): 1-17.

## Acknowledgments

Many of the major algorithmic improvements relative to RaGOO's first release were provided by Aleksey Zimin, lead developer of the [MaSuRCA assembler](https://github.com/alekseyzimin/masurca). [Luca Venturini](https://github.com/lucventurini) suggested and initially implemented many feature enhancments, such as pysam integration. RagTag "merge" was inspired by [CAMSA](https://doi.org/10.1186/s12859-017-1919-y). The developer of CAMSA, [Sergey Aganezov](https://github.com/aganezov), helped review relevant RagTag code. RagTag "patch" was inspired by [Grafter](https://github.com/mkirsche/Grafter), a scaffolding tool written by [Melanie Kirsche](https://github.com/mkirsche). Melanie provided guidance for the RagTag implementation. [Michael Schatz](http://schatz-lab.org/) has provided guidance for the whole project.   
