from ragoo_utilities.Aligner import Minimap2Aligner, NucmerAligner

f = "/Users/malonge/Projects/RaGOO2/ref.fa"
q = "/Users/malonge/Projects/RaGOO2/query.fa"

#x = Minimap2Aligner(f, q, "/usr/local/bin/minimap2", "-x asm5", "/Users/malonge/Projects/RaGOO2/alns", in_overwrite=False)
x = NucmerAligner(f, q, "nucmer", "-l 500 -c 1000", "/Users/malonge/Projects/RaGOO2/alns", in_overwrite=False)
x.run_aligner()