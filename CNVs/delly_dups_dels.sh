# Use for loops to call duplications and deletions for every sample, and create bcf files (input files for delly)

for f in *.mkdup.bam; do delly call -t DUP -g Anopheles_gambiae.AgamP4.dna.toplevel.fa "$f" -o "${f%.*}.dup.bcf"; done
 
for f in *.mkdup.bam; do delly call -t DEL -g Anopheles_gambiae.AgamP4.dna.toplevel.fa "$f" -o "${f%.*}.del.bcf"; done

# and insertions?

# use Jody's script to merge and create a multi-sample vcf of deletions and duplications
# https://jodyphelan.gitbook.io/tutorials/ngs/fst-with-delly

# filtering
