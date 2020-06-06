# RagTag

#### Reference-guided genome assembly correction and scaffolding.

## Getting Started

```bash
git clone https://github.com/malonge/RagTag
cd RagTag
python3 setup.py install

# correct contigs
ragtag.py correct ref.fasta query.fasta

# scaffold contigs
ragtag.py scaffold ref.fa ragtag_output/query.corrected.fasta
```

## Dependencies
- [Minimap2](https://github.com/lh3/minimap2) or [Nucmer](http://mummer.sourceforge.net/)
- Python 3 (with the following auto-installed packages)
    - numpy
    - intervaltree
    - pysam
    
## Citation

Alonge, Michael, et al. ["RaGOO: fast and accurate reference-guided scaffolding of draft genomes."](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1829-6) Genome biology 20.1 (2019): 1-17.


## \*More Docs Soon\*