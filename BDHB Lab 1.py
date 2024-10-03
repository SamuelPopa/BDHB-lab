#BDHB Lab 1
#Objectives:

#Access and retrieve data from online biological databases.
#Parse and manipulate sequence data.
#Perform basic sequence analyses.
#Lab Activity 1: Accessing Biological Databases
#Instructions:

#"First, we'll learn how to access data from NCBI databases programmatically.

#Set Up Entrez:

from Bio import Entrez
Entrez.email = "your.email@example.com"  # Replace with your email
#Note: NCBI requires an email address for usage tracking.

#Search the Database:
#Let's search for the BRCA1 gene in Homo sapiens.

handle = Entrez.esearch(db="nucleotide", term="BRCA1[Gene] AND Homo sapiens[Organism]")
record = Entrez.read(handle)
id_list = record["IdList"]
print(f"Found {len(id_list)} records.")

#Fetch a Record:
handle = Entrez.efetch(db="nucleotide", id=id_list[0], rettype="gb", retmode="text")
from Bio import SeqIO
record = SeqIO.read(handle, "genbank")
print(f"Sequence ID: {record.id}")
print(f"Description: {record.description}")

#Lab Activity 2: Parsing and Analyzing Sequence Data
#Instructions:

#Next, we'll parse the sequence data and perform basic analyses.
#Exploring the SeqRecord Object:
print(f"Sequence Length: {len(record.seq)}")
print(f"Features: {len(record.features)}")

#Calculating GC Content:
from Bio.SeqUtils import GC
gc_content = GC(record.seq)
print(f"GC Content: {gc_content:.2f}%")

#Transcription and Translation:
rna_seq = record.seq.transcribe()
protein_seq = record.seq.translate(to_stop=True)
print(f"Protein Sequence: {protein_seq}")

#Finding Motifs:
motif = "ATG"  # Start codon
positions = [pos for pos in range(len(record.seq)) if record.seq.startswith(motif, pos)]
print(f"Motif '{motif}' found at positions: {positions}")

#Lab Activity 3: Exenting analysis to multiple sequences
#Instructions:

#Download a multi-FASTA file containing mitochondrial sequences from different species.
#Parse the file and calculate GC content for each sequence.
#Compare the GC content across species and discuss any patterns observed.
