from ragoo_utilities.Aligner import Minimap2Aligner

f = "/Users/malonge/Projects/RaGOO2/ref.fa"
q = "/Users/malonge/Projects/RaGOO2/query.fa"

x = Minimap2Aligner(f, q, "minimap2", "-cx asm5", "/Users/malonge/Projects/RaGOO2/alns", in_overwrite=True)
x.run_aligner()