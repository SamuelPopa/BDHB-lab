#Lab Activity 1: Working with DNA Sequences

#Importing Biopython Modules
#Open your Python environment and import the necessary modules:
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#Creating a DNA Sequence Object
#Let's create a DNA sequence:

dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
print("DNA Sequence:", dna_seq)

#Transcription
#Transcribe the DNA to RNA:
rna_seq = dna_seq.transcribe()
print("RNA Sequence:", rna_seq)

#Translation
#Translate the RNA to a protein sequence:
protein_seq = rna_seq.translate()
print("Protein Sequence:", protein_seq)

#Reverse Complement
#Get the reverse complement of the DNA sequence:
rev_comp_seq = dna_seq.reverse_complement()
print("Reverse Complement:", rev_comp_seq)

#Analyzing the Sequence
#Calculate the GC content:
from Bio.SeqUtils import GC
gc_content = GC(dna_seq)
print("GC Content:", gc_content)


#Lab Activity 2: Retrieving and Analyzing Real Genomic Data
#Accessing NCBI Databases

We'll retrieve the mRNA sequence of the BRCA1 gene.
#Setting Up Entrez
from Bio import Entrez
Entrez.email = "your.email@example.com"  # Replace with your email

#Fetching the Sequence
handle = Entrez.efetch(db="nucleotide", id="NM_007294.3", rettype="fasta", retmode="text")
record = handle.read()
print(record)

#Parsing the Sequence
from Bio import SeqIO
from io import StringIO
seq_record = SeqIO.read(StringIO(record), "fasta")
print("Sequence ID:", seq_record.id)
print("Sequence Length:", len(seq_record.seq))

#Analyzing the Sequence
#Identify Open Reading Frames (ORFs)
#An ORF is a sequence of DNA that could potentially encode a protein.
#In Python, you can write a function to find start and stop codons.
#Calculate Nucleotide Frequencies
from collections import Counter
nucleotide_counts = Counter(seq_record.seq)
print("Nucleotide Counts:", nucleotide_counts)
#Use matplotlib in Python or ggplot2 in R to create a bar chart of nucleotide frequencies.

#Lab Activity 3: Exploring Genetic Variation
#Let's investigate known mutations in the BRCA1 gene.
#Accessing SNP Data
#We'll use dbSNP to retrieve information about single nucleotide polymorphisms.
    #Fetching SNP Data
handle = Entrez.esearch(db="snp", term="BRCA1[gene] AND Homo sapiens[organism]")
record = Entrez.read(handle)
snp_ids = record["IdList"]
print(f"Found {len(snp_ids)} SNPs associated with BRCA1.")
    #Fetching SNP Details
snp_record <- entrez_fetch(db="snp", id=snp_ids[1], rettype="docset", retmode="text")
cat(snp_record)


