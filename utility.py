import re
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
import requests
import numpy as np


# Basic Functions
def parse_fasta(file):
    """
    Parse a FASTA file and return a dictionary of sequences.
    Handles both text and binary file inputs.
    """
    sequences = {}
    description = None
    sequence_parts = []

    # Decode the file content if it's in bytes
    for line in file:
        if isinstance(line, bytes):  # Handle binary mode
            line = line.decode("utf-8").strip()
        else:  # Handle text mode
            line = line.strip()

        if line.startswith(">"):
            if description and sequence_parts:
                sequences[description] = "".join(sequence_parts)
                sequence_parts = []
            description = line[1:]  # Exclude ">" in description
        else:
            sequence_parts.append(line)

    if description and sequence_parts:
        sequences[description] = "".join(sequence_parts)

    return sequences





def clean_sequence(sequence):
    """
    Removes unwanted characters (like newline and carriage return) from the DNA sequence.
    """
    return sequence.replace("\n", "").replace("\r", "")



def transcribe_dna_to_rna(dna_sequence):
    return dna_sequence.upper().replace("T", "U")


def translate_orf_with_ambiguity(sequence):
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }
    protein = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        if len(codon) != 3:
            break
        if "N" in codon:
            protein.append("X")
        else:
            protein.append(codon_table.get(codon, "X"))
    return "".join(protein)


def find_orfs(dna_sequence):
    dna_sequence = clean_sequence(dna_sequence)
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    for frame in range(3):
        i = frame
        while i < len(dna_sequence) - 2:
            codon = dna_sequence[i:i + 3]
            if codon == start_codon:
                for j in range(i + 3, len(dna_sequence) - 2, 3):
                    stop_codon = dna_sequence[j:j + 3]
                    if stop_codon in stop_codons:
                        orfs.append(dna_sequence[i:j + 3])
                        i = j
                        break
            i += 3
    return orfs



def calculate_gc_content(sequence):
    gc_count = sum(1 for base in sequence if base in "GC")
    return (gc_count / len(sequence)) * 100 if sequence else 0


def find_cpg_islands(sequence, window_size=200, gc_threshold=50, cpg_ratio_threshold=0.6):
    cpg_islands = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        gc_content = calculate_gc_content(window)
        cpg_count = window.count("CG")
        expected_cpg = (window.count("C") * window.count("G")) / len(window)
        cpg_ratio = cpg_count / expected_cpg if expected_cpg > 0 else 0
        if gc_content > gc_threshold and cpg_ratio > cpg_ratio_threshold:
            cpg_islands.append((i, i + window_size))
    return cpg_islands



def codon_usage(sequence):
    codons = [sequence[i:i + 3] for i in range(0, len(sequence) - 2, 3) if len(sequence[i:i + 3]) == 3]
    return Counter(codons)



def find_restriction_sites(sequence):
    """
    Find restriction enzyme recognition sites in a DNA sequence.

    Args:
        sequence (str): The DNA sequence.

    Returns:
        dict: A dictionary of enzymes and their cut positions.
    """
    # Convert the DNA sequence to a Seq object
    seq_obj = Seq(sequence)
    
    # Define the enzymes to search for
    enzymes = RestrictionBatch(["EcoRI", "BamHI", "HindIII", "NotI"])
    
    # Perform the restriction analysis
    analysis = Analysis(enzymes, seq_obj)
    
    # Return the full analysis results
    return analysis.full()




# JASPAR Database API Functions
def fetch_matrices_for_species(taxonomy_id):
    url = f"https://jaspar.elixir.no/api/v1/matrix/?species={taxonomy_id}&format=json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return [matrix["matrix_id"] for matrix in data.get("results", [])]
    else:
        raise Exception(f"Failed to fetch matrices for species ID {taxonomy_id}: {response.status_code}")




def fetch_pfm_from_jaspar(matrix_id):
    base_url = "https://jaspar.elixir.no/api/v1/matrix"
    url = f"{base_url}/{matrix_id}/?format=json"
    response = requests.get(url)

    # Log the API response
    print(f"Fetching PFM for Matrix ID: {matrix_id}")
    print(f"Response Code: {response.status_code}")
    print(f"Response Content: {response.text}")

    if response.status_code == 200:
        try:
            pfm_data = response.json().get("pfm", None) 
        except ValueError:
            raise ValueError(f"Invalid JSON response for matrix ID {matrix_id}")
        
        if not pfm_data:
            raise ValueError(f"PFM data missing in the response for matrix ID {matrix_id}")
        
        print("Successfully fetched PFM inside function")
        return {
            "A": pfm_data["A"],
            "C": pfm_data["C"],
            "G": pfm_data["G"],
            "T": pfm_data["T"],
        }
    else:
        raise ConnectionError(f"Failed to fetch PFM for {matrix_id}: HTTP {response.status_code}")




def convert_pfm_to_pwm(pfm, pseudocount=0.1):
    pwm = {}
    total_counts = {}
    for pos in range(len(pfm["A"])):
        total_counts[pos] = sum([pfm[base][pos] for base in "ACGT"])
        if total_counts[pos] == 0:
            raise ValueError(f"Zero total count at position {pos} in PFM.")

    for base in "ACGT":
        pwm[base] = []
        for pos in range(len(pfm[base])):
            total = total_counts[pos] + 4 * pseudocount
            pwm_score = np.log2((pfm[base][pos] + pseudocount) / total / 0.25)
            pwm[base].append(pwm_score)

    print(f"PWM for PFM: {pwm}") 
    return pwm





def scan_sequence_with_pwm(sequence, pwm):
    sequence = sequence.upper()
    motif_length = len(pwm["A"])
    scores = []

    if len(sequence) < motif_length:
        raise ValueError(f"DNA sequence is too short for the PWM motif length ({motif_length}).")

    for i in range(len(sequence) - motif_length + 1):
        subseq = sequence[i : i + motif_length]
        if any(base not in "ACGT" for base in subseq):
            print(f"Skipping invalid subsequence at position {i}: {subseq}")
            continue
        score = sum(pwm[base][pos] for pos, base in enumerate(subseq))
        scores.append((i, score))
        print(f"Position: {i}, Subsequence: {subseq}, Score: {score}")  # Debug log

    print(f"Binding sites: {scores}")  
    return scores





def apply_jaspar_to_sequence(sequence, matrix_id):
    pfm = fetch_pfm_from_jaspar(matrix_id)
    pwm = convert_pfm_to_pwm(pfm)
    binding_sites = scan_sequence_with_pwm(sequence, pwm)
    return binding_sites

# Example usage
if __name__ == "__main__":
    user_sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    jaspar_id = "MA0001.1" 

    try:
        sites = apply_jaspar_to_sequence(user_sequence, jaspar_id)
        print("Binding sites found:")
        for position, score in sites:
            print(f"Position: {position}, Score: {score:.2f}")
    except Exception as e:
        print(f"Error: {e}")
