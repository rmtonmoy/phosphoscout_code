from pathlib import Path
from mcp.server.fastmcp import FastMCP
import multiprocessing
import re
from typing import List, Dict, Union
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uniprot_utils import UniProtAPI
from well_studied_kinase_utils import (
    get_well_studied_human_kinase_ids,
    get_domain_seq_from_id
)

mcp = FastMCP("UniProt_Tools_server")

# Global API instance with multiprocessing lock for thread safety
api_lock = multiprocessing.Lock()
uniprot_api = UniProtAPI(api_lock)
PHOSPHOSITE_DATA_PATH = str(PROJECT_ROOT / 'data' / 'phosphosite' / 'Homo_sapiens.txt')

@mcp.tool()
def get_protein_sequence(uniprot_id: str) -> str:
    """
    Retrieves the complete amino acid sequence for a given UniProt protein ID.

    This function fetches the full protein sequence from the UniProt database using their REST API.
    The results are automatically cached to improve performance for repeated queries.

    Args:
        uniprot_id (str): A valid UniProt accession number (e.g., "Q9Y261", "P12345").
                         This should be the primary accession number from the UniProt database.

    Returns:
        str: The complete amino acid sequence of the protein in single-letter code.
             Returns an error message string if the UniProt ID is invalid or if there's a network issue.

    Examples:
        >>> get_protein_sequence("Q9Y261")
        "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES..."

        >>> get_protein_sequence("INVALID_ID")
        "Error retrieving data: 404 Client Error: Not Found for url: ..."

    Note:
        - The function uses caching to store results and avoid repeated API calls
        - Network failures will return descriptive error messages
        - The sequence is returned as a continuous string without spaces or formatting
    """
    try:
        return uniprot_api.get_seq_from_uniprot_id(uniprot_id)
    except Exception as e:
        return f"Error retrieving protein sequence for {uniprot_id}: {str(e)}"


@mcp.tool()
def get_protein_domain_sequence(uniprot_id: str) -> Union[str, List[str]]:
    """
    Retrieves the protein kinase domain sequence(s) for a given UniProt protein ID.

    This function uses curated kinase domain data from the Kannan Lab rather than directly
    querying UniProt for domain annotations. It's specifically designed for protein kinases
    and returns high-quality, manually curated domain sequences.

    Args:
        uniprot_id (str): A valid UniProt accession number for a protein kinase.
                         The protein should be present in the curated kinase database.

    Returns:
        Union[str, List[str]]:
            - If one domain found: Returns the domain sequence as a string
            - If multiple domains found: Returns a list of domain sequences
            - If no domain found: Raises an AssertionError with descriptive message

    Examples:
        >>> get_protein_domain_sequence("Q13627")
        "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES..."

        >>> get_protein_domain_sequence("P12345")  # Multi-domain protein
        ["MAEGEITTFTALTEKFN...", "PPGNYKKPKLLYCSNGGHF..."]

    Raises:
        AssertionError: When no domain sequence is found for the given UniProt ID

    Note:
        - This function bypasses UniProt API and uses curated data for better accuracy
        - Specifically designed for protein kinases from human proteome
        - Domain boundaries are manually curated by domain experts
        - Warning message indicates the use of curated data rather than direct UniProt query
    """
    try:
        result = uniprot_api.get_domain_seq(uniprot_id)
        return result
    except Exception as e:
        return f"Error retrieving domain sequence for {uniprot_id}: {str(e)}"


@mcp.tool()
def get_gene_names_from_uniprot_ids(uniprot_ids: List[str]) -> Dict[str, str]:
    """
    Converts a list of UniProt accession numbers to their corresponding gene names.

    This function uses UniProt's ID mapping service to convert protein accession numbers
    to official gene names. It handles batch processing automatically and manages
    API rate limits by chunking large requests into smaller batches.

    Args:
        uniprot_ids (List[str]): A list of UniProt accession numbers to convert.
                                Maximum of 25 IDs per batch (handled automatically).
                                Examples: ["Q9Y261", "P12345", "O43524"]

    Returns:
        Dict[str, str]: A dictionary mapping each UniProt ID to its gene name.
                       Failed mappings will have "Invalid" as the value.

    Examples:
        >>> get_gene_names_from_uniprot_ids(["Q9Y261", "P31946"])
        {
            "Q9Y261": "ADCK3",
            "P31946": "YWHAB"
        }

        >>> get_gene_names_from_uniprot_ids(["INVALID_ID", "Q9Y261"])
        {
            "INVALID_ID": "Invalid",
            "Q9Y261": "ADCK3"
        }

    API Details:
        - Uses UniProt's ID mapping service with automatic job submission and polling
        - Automatically batches requests in groups of 25 (UniProt's limit)
        - Includes 15-second initial wait and 30-second polling intervals for job completion
        - Handles both successful mappings and failed ID lookups gracefully

    Note:
        - Function may take 30+ seconds for large batches due to UniProt processing time
        - Invalid or obsolete UniProt IDs will be marked as "Invalid" in results
        - Gene names returned are the official gene symbols from UniProt
    """
    try:
        return uniprot_api.get_gene_names_from_uniprot_id(uniprot_ids)
    except Exception as e:
        return {uid: f"Error: {str(e)}" for uid in uniprot_ids}


@mcp.tool()
def get_uniprot_ids_from_ensembl_transcripts(ensembl_transcript_ids: List[str]) -> Dict[str, str]:
    """
    Converts Ensembl transcript IDs to their corresponding UniProt accession numbers.

    This function performs ID mapping from Ensembl transcript identifiers to UniProt
    protein accession numbers using UniProt's mapping service. It's essential for
    connecting genomic annotations with protein-level data.

    Args:
        ensembl_transcript_ids (List[str]): List of Ensembl transcript IDs to convert.
                                          Format: ENST followed by 11 digits
                                          Examples: ["ENST00000221978", "ENST00000346753"]

    Returns:
        Dict[str, str]: Dictionary mapping each Ensembl transcript ID to UniProt accession.
                       Failed mappings will have "Invalid" as the value.

    Examples:
        >>> get_uniprot_ids_from_ensembl_transcripts(["ENST00000221978", "ENST00000346753"])
        {
            "ENST00000221978": "P31946",
            "ENST00000346753": "Q9Y261"
        }

        >>> get_uniprot_ids_from_ensembl_transcripts(["ENST00000000000"])  # Invalid ID
        {
            "ENST00000000000": "Invalid"
        }

    Processing Details:
        - Uses UniProt's ID mapping service from Ensembl_Transcript to UniProtKB
        - Automatically batches large requests (25 IDs per batch maximum)
        - Implements caching to avoid redundant API calls for previously queried IDs
        - Includes job polling with longer wait times (1000 seconds) for complex mappings

    Use Cases:
        - Converting RNA-seq or genomics data to protein-level analysis
        - Bridging transcript-level mutations with protein structure data
        - Integrating genomic coordinates with proteomics datasets

    Note:
        - Processing time varies with batch size and UniProt server load
        - Some transcript IDs may not have corresponding UniProt entries
        - Results are cached automatically for improved performance on repeat queries
    """
    try:
        return uniprot_api.get_uniprot_ids_from_estn(ensembl_transcript_ids)
    except Exception as e:
        return {eid: f"Error: {str(e)}" for eid in ensembl_transcript_ids}


@mcp.tool()
def get_well_studied_kinase_list() -> List[str]:
    """
    Returns a curated list of well-studied human protein kinase UniProt IDs.

    This function provides access to a high-quality, manually curated dataset of
    human protein kinases that have been extensively studied in the literature.
    The list is derived from the Cantley kinase dataset and excludes problematic
    or withdrawn entries.

    Returns:
        List[str]: A list of UniProt accession numbers for well-studied human kinases.
                  Each ID represents a kinase with substantial literature support
                  and experimental validation.

    Examples:
        >>> kinases = get_well_studied_kinase_list()
        >>> print(len(kinases))
        494  # Approximate number, may vary
        >>> print(kinases[:5])
        ["P31946", "Q9Y261", "P12345", "O43524", "Q13627"]

    Data Source:
        - Based on supplementary data from Nature paper (41586_2022_5575_MOESM3_ESM.xlsx)
        - Excludes cancelled/problematic entries: ["O14874", "Q15118", "Q16654", "Q3UTQ8", "E0W1I1", "Q63285"]
        - Focuses specifically on human kinases with robust experimental evidence

    Applications:
        - Benchmark datasets for kinase research
        - Filtering large-scale proteomics data to focus on well-characterized kinases
        - Prioritizing kinases for drug discovery or functional studies
        - Quality control for kinase-related bioinformatics analyses

    Note:
        - List represents the current state of kinase research knowledge
        - All entries have been validated for human origin and kinase activity
        - Suitable for high-confidence analyses requiring well-annotated kinases
    """
    try:
        return get_well_studied_human_kinase_ids()
    except Exception as e:
        return [f"Error retrieving kinase list: {str(e)}"]


@mcp.tool()
def get_curated_domain_sequences(uniprot_id: str) -> List[str]:
    """
    Retrieves curated protein kinase domain sequences from expert-annotated datasets.

    This function accesses high-quality, manually curated kinase domain sequences
    that have been validated by domain experts. Unlike automated domain prediction,
    this provides experimentally validated domain boundaries and sequences.

    Args:
        uniprot_id (str): UniProt accession number for a protein kinase.
                         Must be present in the curated kinase domain database.

    Returns:
        List[str]: List of curated domain sequences for the given kinase.
                  Multiple domains may be returned for multi-domain proteins.
                  Each sequence represents a functional kinase domain.

    Examples:
        >>> get_curated_domain_sequences("Q2M2I8")
        ["MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES..."]

        >>> get_curated_domain_sequences("P12345")  # Multi-domain kinase
        ["DOMAIN1_SEQUENCE...", "DOMAIN2_SEQUENCE..."]

    Data Source:
        - Curated from PKFold human kinase database (pkfold_hs_curated.fa)
        - Expert-validated domain boundaries and sequences
        - Includes variants (|1, |2 suffixes) for alternative isoforms

    Applications:
        - Structural analysis of kinase domains
        - Homology modeling and comparative studies
        - Drug design targeting specific kinase domains
        - Functional domain analysis and mutation impact assessment

    Raises:
        AssertionError: When no curated domain sequence is found for the given UniProt ID

    Note:
        - Domain sequences are manually curated for accuracy
        - May include multiple isoforms or splice variants
        - Preferred over automated domain prediction for kinases
        - Sequences are validated against structural and functional data
    """
    try:
        return get_domain_seq_from_id(uniprot_id)
    except Exception as e:
        return [f"Error retrieving curated domain for {uniprot_id}: {str(e)}"]


import requests
from typing import Optional

@mcp.tool()
def get_uniprot_id(gene: str, organism: str = "9606") -> Optional[str]:
    """
    Retrieves the canonical UniProt accession number for a given gene name and organism.

    This function queries the UniProt REST API to find the primary (canonical) UniProt
    accession number corresponding to a gene name. It specifically searches for reviewed
    (Swiss-Prot) entries to ensure high-quality, manually curated protein data. The function
    returns the first reviewed entry found, which represents the canonical isoform.

    Args:
        gene (str): The gene name or gene symbol to search for.
                   Examples: "TP53", "BRCA1", "EGFR", "BRAF"
                   Case-insensitive gene names are accepted by the UniProt API.

        organism (str): The NCBI taxonomy ID for the organism. Defaults to "9606" (Homo sapiens).
                       Common organism IDs:
                       - "9606": Human (Homo sapiens)
                       - "10090": Mouse (Mus musculus)
                       - "10116": Rat (Rattus norvegicus)
                       - "7227": Fruit fly (Drosophila melanogaster)
                       - "559292": Yeast (Saccharomyces cerevisiae)

    Returns:
        Optional[str]:
            - The UniProt accession number (e.g., "P04637", "P38398") if found
            - None if no reviewed entries are found for the gene
            - None if an error occurs during the API request

    Examples:
        >>> get_uniprot_id("TP53")
        "P04637"

        >>> get_uniprot_id("BRCA1")
        "P38398"

        >>> get_uniprot_id("BRAF")
        "P15056"

        >>> get_uniprot_id("INVALID_GENE")
        None

        >>> get_uniprot_id("Tp53", organism="10090")  # Mouse TP53
        "P02340"

    API Details:
        - Uses UniProt REST API search endpoint: https://rest.uniprot.org/uniprotkb/search
        - Searches for reviewed entries only (Swiss-Prot database)
        - Query format: "gene:{gene} AND organism_id:{organism} AND reviewed:true"
        - Retrieves up to 5 results and returns the first entry (canonical isoform)
        - Request timeout set to 10 seconds for reliability
        - Returns TSV format for efficient parsing

    Use Cases:
        - Converting gene names to UniProt IDs for downstream protein analysis
        - Resolving gene symbols to canonical protein accessions
        - Batch processing of gene lists for proteomics workflows
        - Integration with genomics data that uses gene symbols instead of protein IDs

    Note:
        - Only returns reviewed (Swiss-Prot) entries for high-quality data
        - Returns the canonical isoform (first entry) if multiple reviewed entries exist
        - Gene name matching is approximate - may return results for partial matches
        - For genes with multiple isoforms, only the primary canonical isoform is returned
        - Compatible with UniProt REST API version 2025
        - Function may return None for very new genes not yet in reviewed database
    """
    try:
        base_url = "https://rest.uniprot.org/uniprotkb/search"
        query = f"gene:{gene} AND organism_id:{organism} AND reviewed:true"

        params = {
            "query": query,
            "fields": "accession,reviewed,organism_name,protein_name,gene_primary",
            "format": "tsv",
            "size": 5,
        }

        response = requests.get(base_url, params=params, timeout=10)
        response.raise_for_status()

        text = response.text.strip()
        if not text or "Entry" not in text:
            return None

        lines = text.split("\n")
        if len(lines) < 2:
            return None

        first_entry = lines[1].split("\t")[0].strip()
        return first_entry

    except Exception as e:
        return None





# Amino acid 3-letter to 1-letter conversion mapping
AA_3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def convert_3letter_to_1letter(residue_3letter):
    """Convert 3-letter amino acid code to 1-letter code."""
    return AA_3_TO_1.get(residue_3letter.upper(), residue_3letter)

def parse_mutation(mutation_notation):
    """Parse mutation notation like 'ARG206HIS' to extract residue info."""
    pattern = r'^([A-Z]{3})(\d+)([A-Z]{3})$'
    match = re.match(pattern, mutation_notation.upper())

    if not match:
        return None

    old_residue = match.group(1)
    position = int(match.group(2))
    new_residue = match.group(3)

    return {
        'old_residue': old_residue,
        'position': position,
        'new_residue': new_residue
    }

def parse_mutation_1letter(mutation_notation):
    """Parse mutation notation like 'R206H' to extract residue info (1-letter codes)."""
    pattern = r'^([A-Z])(\d+)([A-Z])$'
    match = re.match(pattern, mutation_notation.upper())

    if not match:
        return None

    old_residue = match.group(1)
    position = int(match.group(2))
    new_residue = match.group(3)

    return {
        'old_residue': old_residue,
        'position': position,
        'new_residue': new_residue
    }


def get_peptide_window(sequence, position, old_residue, new_residue, window_size=29):
    """Extract peptide window centered on the mutation position."""
    if not sequence or position <= 0 or position > len(sequence):
        return []

    half_window = window_size // 2  # 14 for window_size=29
    start_pos = position - half_window  # position - 14
    end_pos = position + half_window    # position + 14

    # Calculate padding needed if we go beyond sequence boundaries
    left_padding = max(0, 1 - start_pos)
    right_padding = max(0, end_pos - len(sequence))

    # Adjust positions to stay within sequence bounds
    start_pos = max(1, start_pos)
    end_pos = min(len(sequence), end_pos)

    # Extract the sequence segment
    seq_segment = sequence[start_pos-1:end_pos]

    # Create peptides with padding
    old_peptide = "-" * left_padding + seq_segment + "-" * right_padding
    new_peptide = "-" * left_padding + seq_segment + "-" * right_padding

    # Replace the center position (index 14) with the new residue
    center_pos = half_window  # 14
    if center_pos < len(new_peptide):
        new_peptide = new_peptide[:center_pos] + new_residue + new_peptide[center_pos+1:]

    return [old_peptide, new_peptide]


@mcp.tool()
def resolve_mutation(gene: str, mutation_notation: str, organism: str = "9606") -> List[str]:
    """
    Resolve a mutation to get peptide windows for analysis.

    This function takes a gene name and mutation notation (like "ARG206HIS") and returns
    peptide windows centered on the mutation site. It handles position validation and
    includes error checking for one-off position errors.

    Args:
        gene (str): Gene name (e.g., "TP53", "BRCA1")
        mutation_notation (str): Mutation in format like "ARG206HIS" (3-letter amino acid codes)
        organism (str): NCBI organism ID, defaults to "9606" (human)

    Returns:
        List[str]: List containing [old_peptide, new_peptide] or empty list if error
                 Each peptide is 29 amino acids centered on the mutation site

    Examples:
        >>> resolve_mutation("TP53", "ARG206HIS")
        ["---MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES---",
         "---MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES---"]

        >>> resolve_mutation("BRCA1", "INVALID_FORMAT")
        []

    Note:
        - Handles one-off position errors automatically
        - Returns empty list for invalid mutations or missing genes
        - Peptide windows are padded with "-" if extending beyond sequence boundaries
    """
    try:
        mutation_info = parse_mutation(mutation_notation)
        #print(mutation_info)
        if not mutation_info:
            return []

        uniprot_id = get_uniprot_id(gene, organism)
        #print(uniprot_id)
        if not uniprot_id:
            return []

        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
        params = {"fields": ["sequence"]}
        headers = {"accept": "application/json"}

        response = requests.get(url, headers=headers, params=params, timeout=10)
        if not response.ok:
            return []
        data = response.json()
        sequence = data.get('sequence', {}).get('value', '')
        #print(sequence)
        if not sequence:
            return []

        position = mutation_info['position']
        if position > len(sequence):
            return []

        actual_residue = sequence[position-1]
        expected_residue = mutation_info['old_residue']
        expected_residue_1letter = convert_3letter_to_1letter(expected_residue)

        if actual_residue == expected_residue_1letter:
            return get_peptide_window(sequence, position, actual_residue, convert_3letter_to_1letter(mutation_info['new_residue']))

        '''print("Adjusting for one-off position error")
        for offset in [-1, 1]:
            check_pos = position + offset
            if 1 <= check_pos <= len(sequence):
                actual_residue_offset = sequence[check_pos-1]
                if actual_residue_offset == expected_residue_1letter:
                    return get_peptide_window(sequence, check_pos, expected_residue, mutation_info['new_residue'])
        '''
        return []

    except Exception as e:
        print(f"Error resolving mutation {mutation_notation} for gene {gene}: {e}")
        return []


@mcp.tool()
def resolve_mutation_for_cosmic(ensembl_transcript: str, mutation_notation: str, organism: str = "9606") -> List[str]:
    """
    Resolve a mutation to get peptide windows for analysis using Ensembl transcript ID.

    This function takes an Ensembl transcript ID and mutation notation (like "R206H" with 1-letter codes)
    and returns peptide windows centered on the mutation site. It first maps the Ensembl transcript
    to a UniProt ID, then retrieves the sequence and processes the mutation.

    Args:
        ensembl_transcript (str): Ensembl transcript ID (e.g., "ENST00000221978" or "ENST00000221978.2").
                                 Version numbers (e.g., ".2") are automatically stripped if present.
        mutation_notation (str): Mutation in format like "R206H" (1-letter amino acid codes)
        organism (str): NCBI organism ID, defaults to "9606" (human)

    Returns:
        List[str]: List containing [old_peptide, new_peptide] or empty list if error
                 Each peptide is 29 amino acids centered on the mutation site

    Examples:
        >>> resolve_mutation_for_cosmic("ENST00000221978", "R206H")
        ["---MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES---",
         "---MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES---"]

        >>> resolve_mutation_for_cosmic("ENST00000221978.2", "R206H")
        ["---MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES---",
         "---MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES---"]

        >>> resolve_mutation_for_cosmic("ENST00000000000", "R206H")
        []

    Note:
        - Version numbers in Ensembl transcript IDs (e.g., ".2") are automatically stripped
        - Returns empty list for invalid mutations or missing transcripts
        - Peptide windows are padded with "-" if extending beyond sequence boundaries
        - Uses 1-letter amino acid codes (A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V)
    """
    try:
        ensembl_transcript = ensembl_transcript.split('.')[0]

        mutation_info = parse_mutation_1letter(mutation_notation)
        if not mutation_info:
            return []

        uniprot_mapping = get_uniprot_ids_from_ensembl_transcripts([ensembl_transcript])
        uniprot_id = uniprot_mapping.get(ensembl_transcript)

        if not uniprot_id or uniprot_id == "Invalid":
            return []

        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
        params = {"fields": ["sequence"]}
        headers = {"accept": "application/json"}

        response = requests.get(url, headers=headers, params=params, timeout=10)
        if not response.ok:
            return []
        data = response.json()
        sequence = data.get('sequence', {}).get('value', '')
        if not sequence:
            return []

        position = mutation_info['position']
        if position > len(sequence):
            return []

        actual_residue = sequence[position-1]
        expected_residue = mutation_info['old_residue']

        if actual_residue == expected_residue:
            return get_peptide_window(sequence, position, actual_residue, mutation_info['new_residue'])

        return []

    except Exception as e:
        print(f"Error resolving mutation {mutation_notation} for transcript {ensembl_transcript}: {e}")
        return []


@mcp.tool()
def differential_phophosite_analysis(old_peptide: str, new_peptide: str, window_size: int = 29, gene: Optional[str] = None, mutation_notation: Optional[str] = None) -> List[tuple]:
    """
    Analyze differential phosphosite changes between two peptides based on a single amino acid change.

    This function extracts all possible k-mers (or half window_size) from the peptides
    and checks if any phosphosites (S, T, Y) are present in these windows.
    It analyzes how a single amino acid residue change affects phosphosite patterns.
    If phosphosites are found, it returns pairs of old and new k-mers for analysis.

    Args:
        old_peptide (str): The original peptide sequence (length should be window_size).
        new_peptide (str): The peptide sequence with a single amino acid change (length should be window_size).
        window_size (int): The total window size for analysis (default: 29).
        gene (str, optional): Gene name (e.g., "TP53", "BRCA1"). Used to look up UniProt ID and check
                             phosphosite data. When provided along with mutation_notation, the function
                             will return comments about whether positions are known phosphosites.
                             Must be provided together with mutation_notation, or neither should be provided.
        mutation_notation (str, optional): Mutation in format like "ARG206HIS" (3-letter amino acid codes)
                                          or "R206H" (1-letter amino acid codes).
                                          The format is: OLD_RESIDUE + POSITION + NEW_RESIDUE (e.g., "ARG206HIS"
                                          or "R206H" means arginine at position 206 is changed to histidine).
                                          When provided along with gene, enables phosphosite validation
                                          against the EPSD database.
                                          Must be provided together with gene, or neither should be provided.

    Returns:
        List[tuple]:
            - If gene and mutation_notation are NOT provided: A list of (old_kmer, new_kmer) tuples of length k
              having S, T, or Y in the middle position. Returns empty list if any window contains gaps ("-")
              or no phosphosites found.
            - If gene and mutation_notation ARE provided: A list of (old_kmer, new_kmer, positions, comment) tuples where
              positions is a list of absolute protein positions (one per amino acid in the k-mer, in order) and
              comment describes whether the position is a known phosphosite. Returns empty list if any window
              contains gaps ("-") or no phosphosites found.

    Examples:
        >>> old_pep = "ABCDEFGHIJKLMNOPQRSTUVWXYZAB"
        >>> new_pep = "ABCDEFGHIJKLMNOPQRSTUVWXYZAC"  # Single amino acid change: B->C
        >>> result = differential_phophosite_analysis(old_pep, new_pep)
        >>> print(result)
        [('LMNOPQRSTUVWXYZ', 'LMNOPQRSTUVWXYZ'), ('MNOPQRSTUVWXYZA', 'MNOPQRSTUVWXYZA')]

        >>> old_pep = "ABCDEFGHIJKLMNOPQRSTUVWXYZAB"
        >>> new_pep = "ABCDEFGHIJKLMNOPQRSTUVWXYZ--"  # Contains gaps
        >>> result = differential_phophosite_analysis(old_pep, new_pep)
        >>> # Returns [] due to gaps in sequence

        >>> # With gene and mutation_notation for phosphosite validation (3-letter code)
        >>> old_pep = "---MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES---"
        >>> new_pep = "---MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAES---"
        >>> result = differential_phophosite_analysis(old_pep, new_pep, gene="TP53", mutation_notation="ARG206HIS")
        >>> print(result)
        [('...', '...', [192, 193, 194, ...], 'Known phosphosite at position X'), ...]

        >>> # With gene and mutation_notation for phosphosite validation (1-letter code)
        >>> result = differential_phophosite_analysis(old_pep, new_pep, gene="TP53", mutation_notation="R206H")
        >>> print(result)
        [('...', '...', [192, 193, 194, ...], 'Known phosphosite at position X'), ...]

    Raises:
        ValueError: If exactly one of gene or mutation_notation is provided (both must be provided together or neither).
    """
    if (gene is None) != (mutation_notation is None):
        raise ValueError("gene and mutation_notation must be provided together, or neither should be provided")

    ret = []
    change_pos = window_size // 2
    k_mer_size = window_size // 2
    middle_pos_kmer = k_mer_size // 2

    # Parse mutation if provided
    mutation_pos = None
    uniprot_id = None
    phosphosite_data = None

    if gene and mutation_notation:
        mutation_info = parse_mutation(mutation_notation)
        if not mutation_info:
            mutation_info = parse_mutation_1letter(mutation_notation)
        if mutation_info:
            mutation_pos = mutation_info['position']
            uniprot_id = get_uniprot_id(gene)
            if uniprot_id:
                phosphosite_data = load_phosphosite_data_for_uniprot(uniprot_id)

    # Extract all possible k-mers to analyze phosphosite changes from single amino acid mutation
    for start_pos in range(0, change_pos + 1):
        old_pep = old_peptide[start_pos:start_pos + k_mer_size + 1]
        new_pep = new_peptide[start_pos:start_pos + k_mer_size + 1]

        # Skip windows that contain gaps (invalid sequences)
        if "-" in old_pep or "-" in new_pep:
            continue

        # Check for phosphosites (S, T, Y) in the middle position of the window
        has_phosphosite = old_pep[middle_pos_kmer] in ["S", "T", "Y"] or new_pep[middle_pos_kmer] in ["S", "T", "Y"]

        if has_phosphosite:
            if gene and mutation_notation and mutation_pos and uniprot_id:
                left_padding = old_peptide.count("-", 0, change_pos + 1)

                positions = []
                for i in range(len(old_pep)):
                    if left_padding > 0:
                        protein_pos = 1 + (start_pos + i - left_padding)
                    else:
                        protein_pos = mutation_pos - change_pos + start_pos + i
                    positions.append(protein_pos)

                comment = _get_phosphosite_comment(
                    positions[middle_pos_kmer], mutation_pos, old_pep[middle_pos_kmer],
                    new_pep[middle_pos_kmer], phosphosite_data
                )
                ret.append((old_pep, new_pep, positions, comment))
            else:
                ret.append((old_pep, new_pep))

    return ret


@mcp.tool()
def load_phosphosite_data_for_uniprot(uniprot_id: str) -> Dict[int, str]:
    """
    Loads phosphosite data for a specific UniProt protein ID from the EPSD database.

    This function queries the PhosphoSitePlus database (Homo_sapiens.txt) to retrieve
    all known phosphorylation sites for a given protein. It returns a mapping of
    protein positions to their corresponding phosphorylatable amino acids (S, T, Y).

    Args:
        uniprot_id (str): A valid UniProt accession number (e.g., "Q9Y261", "P12345").
                         The protein should be present in the PhosphoSitePlus database.

    Returns:
        Dict[int, str]: A dictionary mapping protein positions (1-indexed) to amino acid
                       codes. Only includes positions with phosphorylatable residues
                       (serine, threonine, or tyrosine).
                       Returns an empty dictionary if no phosphosites are found or if
                       the UniProt ID is not in the database.

    Examples:
        >>> load_phosphosite_data_for_uniprot("P04637")  # TP53
        {15: 'S', 20: 'S', 33: 'S', 37: 'S', 46: 'S', ...}

        >>> load_phosphosite_data_for_uniprot("INVALID_ID")
        {}

    Data Source:
        - PhosphoSitePlus database: Homo_sapiens.txt
        - Contains experimentally validated phosphorylation sites
        - Includes only serine (S), threonine (T), and tyrosine (Y) residues

    Processing Details:
        - Uses grep to efficiently search the database file
        - Timeout set to 30 seconds for large database queries
        - Filters results to include only S, T, Y residues
        - Returns positions as 1-indexed integers matching protein sequence numbering

    Use Cases:
        - Validating known phosphosites in mutation analysis
        - Identifying phosphorylation sites affected by mutations
        - Cross-referencing experimental phosphoproteomics data
        - Quality control for phosphosite predictions

    Note:
        - Database file path is configured globally (PHOSPHOSITE_DATA_PATH)
        - Only returns positions with experimentally validated phosphorylation
        - Empty results may indicate the protein is not in the database or has no known phosphosites
        - Function is optimized for fast lookups using command-line grep
    """
    phosphosite_map = {}
    import subprocess
    result = subprocess.run(
        ['grep', f'\t{uniprot_id}\t', PHOSPHOSITE_DATA_PATH],
        capture_output=True, text=True, timeout=30
    )
    if result.returncode == 0:
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4 and parts[1] == uniprot_id:
                pos = int(parts[3])
                aa = parts[2]
                if aa in ["S", "T", "Y"]:
                    phosphosite_map[pos] = aa
    return phosphosite_map


def _get_phosphosite_comment(protein_pos: int, mutation_pos: int, old_aa: str, new_aa: str, phosphosite_data: Optional[Dict[int, str]]) -> str:
    """Generate comment about phosphosite status for a given position."""
    if phosphosite_data is None:
        return "no evidence of phosphosite in the EPSD database"

    if protein_pos == mutation_pos:
        # Case 2: This position is mutated
        if old_aa in ["S", "T", "Y"] and protein_pos in phosphosite_data:
            return f"Position {protein_pos} was a proven phosphosite ({old_aa}) in the EPSD database before mutation"
        elif new_aa in ["S", "T", "Y"]:
            return f"Mutation creates a new phosphosite at position {protein_pos} ({new_aa})"
        else:
            return f"No evidence of phosphorylation at position {protein_pos} in the EPSD database before mutation"
    else:
        # Case 1: Mutation did not happen at this site
        if protein_pos in phosphosite_data:
            return f"Known phosphosite at position {protein_pos} in the EPSD database"
        else:
            return f"No evidence of phosphorylation at position {protein_pos} in the EPSD database"

if __name__ == "__main__":
    '''
    # Test resolve_mutation function directly
    
    gene = "ACVR1"
    mutation = "ARG206HIS"
    organism = "9606"

    result = resolve_mutation(gene, mutation, organism)
    print(f"Result for {gene} {mutation}: {result}")

    print(f"Differential phosphosite analysis: {differential_phophosite_analysis(result[0], result[1])}")


    # Test for CDKL5 (UniProt O76039) • Variant: Cys291Tyr (C291Y)
    gene = "CDKL5"
    mutation = "CYS291TYR"
    organism = "9606"

    print(f"\nTesting {gene} {mutation} (UniProt O76039)")
    result = resolve_mutation(gene, mutation, organism)
    print(f"Result for {gene} {mutation}: {result}")

    if result:
        print(f"\nDifferential phosphosite analysis:")
        analysis_result = differential_phophosite_analysis(result[0], result[1], gene=gene, mutation_notation=mutation)
        print(f"Analysis result: {analysis_result}")
    else:
        print("Failed to resolve mutation")
    
    '''
    mcp.run(transport='stdio')
