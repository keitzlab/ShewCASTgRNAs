import streamlit as st
from Bio.KEGG import REST
import re

# Function to load genes from the 'shewgenesparse.txt' file
def load_gene_list():
    genes = []
    file_path = 'shewgenesparse.txt'  # Assuming the file is in the same directory as the script
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split(" ", 1)
                if len(parts) == 2:
                    genes.append({"id": parts[0], "name": parts[1]})
    except Exception as e:
        st.error(f"Error loading gene list: {e}")
    
    return genes

# Function to fetch the first 500 bp of a gene's DNA sequence from KEGG
def get_gene_dna_sequence_500bp(organism, gene_name):
    try:
        search_result = REST.kegg_find("genes", gene_name).read()
        gene_id = None
        for line in search_result.strip().split("\n"):
            if line.startswith(f"{organism}:"):
                gene_id = line.split("\t")[0]
                break
        if not gene_id:
            return None
        gene_data = REST.kegg_get(gene_id).read()
        ntseq_match = re.search(r"NTSEQ\s+\d+\n(.+?)(?:\n//|AASEQ)", gene_data, re.DOTALL)
        if ntseq_match:
            dna_sequence = ntseq_match.group(1).replace(" ", "").replace("\n", "")
            return dna_sequence[:500]
    except Exception as e:
        st.error(f"An error occurred: {e}")
    return None

# Function to find NCC motifs with context
def find_ncc_motifs_with_context(sequence, prefix, suffix):
    sequence = sequence.upper()
    contexts = []
    pattern = re.compile(r"(.CC)")
    i = 0
    while i < len(sequence) - 2 and len(contexts) < 3:
        match = pattern.search(sequence, i)
        if match:
            start = match.start()
            context_start = start + len(match.group(1))
            context_end = context_start + 32
            if context_end <= len(sequence):
                context = sequence[context_start:context_end]
                formatted_context = f"{prefix}{context}{suffix}"
                contexts.append(formatted_context)
            else:
                break
            i = context_end
        else:
            break
    return contexts

# Golden Gate Assembly Instructions
def golden_gate_assembly():
    st.markdown("""
    ## Golden Gate Assembly Instructions:
    To perform Golden Gate Assembly, use the following protocol for each of the generated gRNA gBlocks:
    
    | Reagent                          | Volume         |
    |----------------------------------|----------------|
    | pCAST GFP Dropout Backbone (100 ng/µL) | 3.41 µL        |
    | gRNA gBlock (10 ng/µL)           | 3.36 µL        |
    | Water                            | 3.36 µL        |
    | T4 DNA Ligase Buffer            | 1 µL           |
    | T4 DNA Ligase                   | 0.5 µL         |
    | BsaI                             | 0.5 µL         |
    
    ### Thermal Cycler Program:
    1. **37°C** for 1.5 minutes
    2. **16°C** for 3 minutes
    Repeat steps 1 and 2 **25 times**.
    
    After 25 cycles:
    1. **37°C** for 5 minutes
    2. **80°C** for 10 minutes
    3. Hold at **4°C** indefinitely.
    """, unsafe_allow_html=True)

# Streamlit app UI
def app():
    # Header with UT Austin colors
    st.markdown("""
    <style>
    .title {
        color: #bf5700;
        font-size: 2.5em;
        font-weight: bold;
    }
    </style>
    """, unsafe_allow_html=True)

    # Title and instructions
    st.markdown("<h1 class='title'>Shewanella CAST gRNA Generator</h1>", unsafe_allow_html=True)

    # Gene dropdown with search functionality
    genes = load_gene_list()
    gene_names = [f"{gene['id']}: {gene['name']}" for gene in genes]
    
    selected_gene = st.selectbox("Select Shewanella Gene:", gene_names)

    # Custom DNA input
    custom_dna = st.text_area("Or Paste DNA Sequence (up to 500 bp):", height=200)

    # Prefix and suffix for the gRNA
    prefix = "tttttaaatacgggctgacaatcagtgatcgttacgatcatcgatcgatcgtacgatcgcggctcacgatcatcttactgtgtcagtggtctcacatagtgaactgccgagtaggtagctgataac"
    suffix = "gtgaactgccgagtaggtagctgataaccgttagagacctagctagctagcatcgatcgggatcgaatcgactagtatatcgatcgatcgatcggcatcataactgactgtgtcacgtagcatcgtaggtagcgtagctacgatcattttggcggcgaaatgtcgacgttagtactgcagtactgcagtac"

    # Button to generate gRNA gBlocks
    if st.button("Generate gRNA gBlocks", key="generate"):
        # Handle custom DNA or selected gene
        dna_sequence = None
        if custom_dna:
            dna_sequence = custom_dna[:500]  # Take only the first 500 bases
        elif selected_gene:
            gene_name = selected_gene.split(":")[1]
            dna_sequence = get_gene_dna_sequence_500bp("son", gene_name)

        # Find NCC motifs if DNA sequence is provided
        if dna_sequence:
            contexts = find_ncc_motifs_with_context(dna_sequence, prefix, suffix)

            st.markdown("### gRNA gBlocks Generated:")
            for i, context in enumerate(contexts, 1):
                # Display entire gRNA sequence in a text area
                st.text_area(f"gRNA {i} Sequence", value=context, height=100, key=f"gRNA_{i}", disabled=True)

        else:
            st.error("No DNA sequence found. Please select a gene or paste a DNA sequence.")
    
    # Golden Gate Assembly instructions
    golden_gate_assembly()

if __name__ == "__main__":
    app()
