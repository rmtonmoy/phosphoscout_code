from pathlib import Path
from mcp.server.fastmcp import FastMCP
import requests
from pydantic import BaseModel, field_validator
from typing import List, Dict
import pandas as pd
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[2]

mcp = FastMCP("Phosphorylation_Tools")

import os
import sys
sys.path.append(str(PROJECT_ROOT / 'tools' / 'kinase-library' / 'src'))
import kinase_library as kl
import requests
sys.path.insert(0, str(Path(__file__).resolve().parent))
from kinase_proximity import check_proximity

def get_gene_name(uniprot_id):
    """
    Retrieve the primary gene name for a given UniProt accession ID.
    Works with the 2025 UniProt REST API.
    Example: get_gene_name("Q96GX5") -> "MASTL"
    """
    try:
        base_url = f"https://rest.uniprot.org/uniprotkb/search"
        query = f"accession:{uniprot_id}"

        params = {
            "query": query,
            "fields": "accession,gene_primary,organism_name,reviewed",
            "format": "tsv",
            "size": 1,
        }

        response = requests.get(base_url, params=params, timeout=10)
        response.raise_for_status()

        text = response.text.strip()
        if not text or "Entry" not in text:
            return None

        lines = text.split("\n")
        if len(lines) < 2:
            return None

        # Split the second line (first data entry)
        fields = lines[1].split("\t")
        gene_name = fields[1].strip() if len(fields) > 1 else None

        return gene_name or None

    except Exception as e:
        return None
def convert_uniprot_ids_to_gene_names(uniprot_change_pairs: list[(str, float)]) -> list[(str, float)]:
    """
    Convert a list of (UniProt ID, change) pairs to (gene name, change) pairs.

    Args:
        uniprot_change_pairs: List of tuples containing (UniProt ID, change value)

    Returns:
        List of tuples containing (gene name, change value). If gene name lookup fails,
        falls back to using the UniProt ID as the identifier.

    Example:
        convert_uniprot_ids_to_gene_names([("P31749", 0.15), ("Q96GX5", -0.23)])
        -> [("AKT1", 0.15), ("MASTL", -0.23)]
    """
    result_with_gene_names = []
    for uniprot_id, change in uniprot_change_pairs:
        gene_name = get_gene_name(uniprot_id)
        if gene_name:
            result_with_gene_names.append((gene_name, change))
        else:
            # Fallback to UniProt ID if gene name lookup fails
            result_with_gene_names.append((uniprot_id, change))

    return result_with_gene_names
def _get_domain_seq_from_id(uniprot_id):
    curated_dom_seq_path = str(PROJECT_ROOT / 'data' / 'kinase_data' / 'pkfold_hs_curated.fa')
    domain_seqs = []

    for record in SeqIO.parse(curated_dom_seq_path, "fasta"):
        if record.id == uniprot_id or record.id == f"{uniprot_id}|1" or record.id == f"{uniprot_id}|2":
            domain_seqs.append(str(record.seq))

    if(len(domain_seqs) == 0):
        raise AssertionError("Domain sequence not found")

    return domain_seqs
def get_Y_kinase_list() -> List[str]:
    """Get list of Y-kinase names.

    Returns:
        List[str]: List of Y-kinase names.
    """
    return kl.get_kinase_list('tyrosine')
def get_ST_kinase_list() -> List[str]:
    """Get list of serine/threonine kinase UniProt IDs.

    Returns:
        List[str]: List of serine/threonine kinase UniProt IDs.
    """
    cantley_kinase_path = str(PROJECT_ROOT / 'data' / 'kinase_data' / '41586_2022_5575_MOESM3_ESM.xlsx')
    data = pd.read_excel(cantley_kinase_path, sheet_name = 1, engine='openpyxl')
    ids = data['Uniprot id'].to_list()
    cancelled = ["O14874", "Q15118", "Q16654", "Q3UTQ8", "E0W1I1", "Q63285"]
    return [id for id in ids if id not in cancelled]
class PhosformerRequest(BaseModel):
    kinase: str
    protein_segment: str

    @field_validator('kinase')
    @classmethod
    def validate_kinase(cls, v: str) -> str:
        # Check if it's a valid amino acid sequence
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        if not all(char.upper() in valid_aa for char in v):
            raise ValueError(
                'Kinase must be a valid amino acid sequence containing only standard amino acid letters. '
                'Do not use UniProt IDs (e.g., P12345) or protein names. '
                'Provide the actual amino acid sequence of the kinase domain.'
            )

        # Check minimum length for kinase domain
        if len(v) < 20:
            raise ValueError(
                f'Kinase sequence is too short ({len(v)} amino acids). '
                'Kinase domains are typically 200-400 amino acids long. '
                'You may have provided a UniProt ID or protein name instead of the actual domain sequence. '
                'Please provide the full amino acid sequence of the kinase domain.'
            )

        # Check maximum reasonable length
        if len(v) > 900:
            raise ValueError(
                f'Kinase sequence is unusually long ({len(v)} amino acids). '
                'Kinase domains are typically 200-900 amino acids long. '
                'Please provide only the kinase domain sequence, not the full protein sequence.'
            )

        return v.upper()

    @field_validator('protein_segment')
    @classmethod
    def validate_protein_segment(cls, v: str) -> str:
        if len(v) != 15:
            raise ValueError('Protein segment must be exactly 15 characters long')
        middle_char = v[7]  # 0-indexed, so 7 is the middle position
        if middle_char not in ['S', 'T']:
            raise ValueError('Middle character of protein segment must be S or T')
        return v
class KlRequest(BaseModel):
    kinase: str
    protein_segment: str

    @field_validator('kinase')
    @classmethod
    def validate_kinase(cls, v: str) -> str:
        if not v or v not in kl.get_kinase_list('tyrosine'):
            raise ValueError('Kinase name must be in the list of tyrosine kinases')
        return v.strip()

    @field_validator('protein_segment')
    @classmethod
    def validate_protein_segment(cls, v: str) -> str:
        if len(v) != 15:
            raise ValueError('Protein segment must be exactly 15 characters long')

        middle_char = v[7]
        if middle_char.upper() != 'Y':
            raise ValueError('Middle character of protein segment must be Y (case insensitive)')

        # Convert Y to lowercase for processing
        protein_segment_list = list(v)
        if protein_segment_list[7].upper() == 'Y':
            protein_segment_list[7] = 'y'

        return ''.join(protein_segment_list)

@mcp.tool()
def get_phosphorylation_confidence_for_ST_kinase(kinase: str, protein_segments: List[str]) -> List[float]:
    """Get phosphorylation probabilities for a kinase against one or more protein segments.

    IMPORTANT: This function requires actual amino acid sequences, NOT UniProt IDs or protein names.

    Args:
        kinase: Single kinase amino acid sequence (200-500 amino acids).
                Do NOT provide UniProt IDs (e.g., P12345), gene names, or protein names.
                Example: "MKTSLSQGQKPKGTKGTGKDDDKPTFYISTKGTDPTKPKLTSADLKKIGGDLL..."
        protein_segments: List of protein segment amino acid sequences (must be exactly 15 chars with S or T in middle position).
                   Can contain a single protein segment or multiple protein segments.
                   Example: ["LARKRRNSRDGDPLP"] for single protein segment
                   Example: ["LARKRRNSRDGDPLP", "MKKKRRRSTAAAAAA"] for multiple protein segments

    Returns:
        List[float]: List of phosphorylation probability scores between 0 and 1,
                    in the same order as the input protein segments. For single protein segment, returns a list with one element.

    Raises:
        ValueError: If kinase is not a valid amino acid sequence or wrong length,
                   or if any protein segment format is incorrect
    """

    if not protein_segments:
        raise ValueError('At least one protein segment must be provided')

    # Validate the kinase once
    try:
        validated_kinase = PhosformerRequest(kinase=kinase, protein_segment=protein_segments[0]).kinase
    except ValueError as e:
        raise ValueError(f'Kinase validation error: {str(e)}')

    # Validate each protein segment and create kinase-protein_segment pairs
    validated_kinases = []
    validated_substrates = []

    for i, protein_segment in enumerate(protein_segments):
        try:
            request = PhosformerRequest(kinase=kinase, protein_segment=protein_segment)
            validated_kinases.append(validated_kinase)
            validated_substrates.append(request.protein_segment)
        except ValueError as e:
            raise ValueError(f'Validation error for protein segment {i+1}: {str(e)}')

    server_url = "http://localhost:5000"
    payload = {
        "kinases": validated_kinases,
        "substrates": validated_substrates
    }
    response = requests.post(
        f"{server_url}/predict",
        json=payload,
        headers={'Content-Type': 'application/json'}
    )

    if response.status_code == 200:
        return response.json()['probabilities']
    else:
        return f"Error: {response.json()}"

@mcp.tool()
def get_phosphorylation_confidence_for_Y_kinase(kinase: str, protein_segments: List[str]) -> List[float]:
    """Get kinase scores and percentiles for one or more protein segments with Y in the middle position.

    Args:
        kinase: Kinase name (e.g., 'AKT1', 'CDK2')
        protein_segments: List of 15-character protein segment sequences with Y in the middle position (7th character).
                   Can contain a single protein segment or multiple protein segments.
                   Y can be uppercase or lowercase in input, will be converted to lowercase for processing

    Returns:
        List[float]: List of kinase scores for each protein segment. For single protein segment, returns a list with one element.
    """
    if not protein_segments:
        raise ValueError('At least one protein segment must be provided')

    # Validate the kinase once using the first protein segment
    try:
        validated_kinase = KlRequest(kinase=kinase, protein_segment=protein_segments[0]).kinase
    except ValueError as e:
        raise ValueError(f'Kinase validation error: {str(e)}')

    # Validate each protein segment and create kinase-protein_segment pairs
    validated_kinases = []
    validated_substrates = []

    for i, protein_segment in enumerate(protein_segments):
        try:
            request = KlRequest(kinase=kinase, protein_segment=protein_segment)
            validated_kinases.append(validated_kinase)
            validated_substrates.append(request.protein_segment)
        except ValueError as e:
            raise ValueError(f'Validation error for protein segment {i+1}: {str(e)}')

    # Get kinase and score protein segments
    kin = kl.get_kinase(validated_kinase)
    kin.score(validated_substrates)
    # Convert to list to avoid pandas FutureWarning about Series.__getitem__
    return kin.percentile(validated_substrates).tolist()

@mcp.tool()
def get_predicted_change_in_substrate_specifity_for_Y_kinase(wildtype_protein_segment: str, mutant_protein_segment: str) -> List[dict]:
    """Get predicted change in protein segment specificity for all Y-kinases between wildtype and mutant protein segments.

    This function compares phosphorylation confidence scores for all tyrosine kinases
    between a wildtype protein segment and a mutant protein segment, returning the kinases
    with the largest absolute changes in specificity.

    Args:
        wildtype_protein_segment: Wildtype protein segment sequence (15 characters with Y in middle position).
                      Example: "RRRLYSLRASTSKSN"
        mutant_protein_segment: Mutant protein segment sequence (15 characters with Y in middle position).
                       Example: "RRRLYSLRASTSKSN" (same as wildtype) or "RRRLYSLRASTSKSN" (different)

    Returns:
        List[dict]: List of dictionaries containing kinase_name and change_value for the
                   top 10 kinases with largest absolute changes in protein segment specificity.
                   Each dictionary has keys "kinase", "change", "wildtype_value", and "mutant_value".
                   The change value represents: mutant_confidence - wildtype_confidence
                   - Positive values indicate increased phosphorylation likelihood
                   - Negative values indicate decreased phosphorylation likelihood
                   - Results are sorted by absolute change magnitude (descending)

        Example return:
        [
            {"kinase": "SRC", "change": 0.234, "wildtype_value": 0.456, "mutant_value": 0.690},
            {"kinase": "ABL1", "change": -0.189, "wildtype_value": 0.723, "mutant_value": 0.534},
            {"kinase": "EGFR", "change": 0.156, "wildtype_value": 0.234, "mutant_value": 0.390},
            ...
        ]

    Raises:
        ValueError: If protein segments are not exactly 15 characters or don't have Y in middle position
    """
    if len(wildtype_protein_segment) != 15:
        raise ValueError('Wildtype protein segment must be exactly 15 characters long')
    if len(mutant_protein_segment) != 15:
        raise ValueError('Mutant protein segment must be exactly 15 characters long')

    valid_wildtype = wildtype_protein_segment[7].upper() == 'Y'
    valid_mutant = mutant_protein_segment[7].upper() == 'Y'

    data: Dict[str, Dict[str, float]] = {}
    for kin in get_Y_kinase_list():
        wildtype_confidence = get_phosphorylation_confidence_for_Y_kinase(kin, [wildtype_protein_segment])[0] if valid_wildtype else 0.0
        mutant_confidence = get_phosphorylation_confidence_for_Y_kinase(kin, [mutant_protein_segment])[0] if valid_mutant else 0.0
        change = mutant_confidence - wildtype_confidence
        data[kin] = {"change": change, "wildtype_value": wildtype_confidence, "mutant_value": mutant_confidence}

    return [{"kinase": kinase, "change": info["change"], "wildtype_value": info["wildtype_value"], "mutant_value": info["mutant_value"]}
            for kinase, info in data.items()]

@mcp.tool()
def get_predicted_change_in_substrate_specifity_for_ST_kinase(wildtype_protein_segment: str, mutant_protein_segment: str) -> List[dict]:
    """Get predicted change in protein segment specificity for all ST-kinases between wildtype and mutant protein segments.

    This function compares phosphorylation confidence scores for all serine/threonine kinases
    between a wildtype protein segment and a mutant protein segment, returning the kinases
    with the largest absolute changes in specificity. For kinases with multiple domains,
    the domain showing the maximum absolute change is selected.

    Args:
        wildtype_protein_segment: Wildtype protein segment sequence (15 characters; middle residue unrestricted).
                      Example: "RRRLSSLRASTSKSN"
        mutant_protein_segment: Mutant protein segment sequence (15 characters; middle residue unrestricted).
                       Example: "RRRLSSLRASTSKSN" (same as wildtype) or "RRRLSSLRASTSKSN" (different)

    Returns:
        List[dict]: List of dictionaries containing gene_name and change_value for the
                   top 10 kinases with largest absolute changes in protein segment specificity.
                   Each dictionary has keys "kinase", "change", "wildtype_value", and "mutant_value".
                   Gene names are retrieved from UniProt IDs; if lookup fails,
                   UniProt ID is used as fallback.
                   The change value represents: mutant_confidence - wildtype_confidence
                   - Positive values indicate increased phosphorylation likelihood
                   - Negative values indicate decreased phosphorylation likelihood
                   - Results are sorted by absolute change magnitude (descending)

        Example return:
        [
            {"kinase": "AKT1", "change": 0.234, "wildtype_value": 0.456, "mutant_value": 0.690},
            {"kinase": "CDK2", "change": -0.189, "wildtype_value": 0.723, "mutant_value": 0.534},
            {"kinase": "MAPK1", "change": 0.156, "wildtype_value": 0.234, "mutant_value": 0.390},
            ...
        ]

    Raises:
        ValueError: If protein segments are not exactly 15 characters
        AssertionError: If domain sequence cannot be found for a kinase
    """
    if len(wildtype_protein_segment) != 15:
        raise ValueError('Wildtype protein segment must be exactly 15 characters long')
    if len(mutant_protein_segment) != 15:
        raise ValueError('Mutant protein segment must be exactly 15 characters long')

    valid_wildtype = wildtype_protein_segment[7].upper() in ['S', 'T']
    valid_mutant = mutant_protein_segment[7].upper() in ['S', 'T']

    data: Dict[str, Dict[str, float]] = {}
    for kin in get_ST_kinase_list():
        for dom_seq in _get_domain_seq_from_id(kin):
            wildtype_confidence = get_phosphorylation_confidence_for_ST_kinase(dom_seq, [wildtype_protein_segment])[0] if valid_wildtype else 0.0
            mutant_confidence = get_phosphorylation_confidence_for_ST_kinase(dom_seq, [mutant_protein_segment])[0] if valid_mutant else 0.0
            change = mutant_confidence - wildtype_confidence

            if kin not in data or abs(change) > abs(data[kin]["change"]):
                data[kin] = {"change": change, "wildtype_value": wildtype_confidence, "mutant_value": mutant_confidence}

    # Convert UniProt IDs to gene names; fallback to UniProt ID on failure
    result: List[Dict[str, float]] = []
    for uniprot_id, info in data.items():
        gene_name = get_gene_name(uniprot_id) or uniprot_id
        result.append({"kinase": gene_name, "change": info["change"], "wildtype_value": info["wildtype_value"], "mutant_value": info["mutant_value"]})
    return result

@mcp.tool()
def get_substantial_change_in_specificity_with_proximity_check(gene_name: str, wildtype_protein_segment: str, mutant_protein_segment: str) -> List[dict]:
    """Get substantial specificity changes for proximal kinases between wildtype and mutant protein segments.

    This function first identifies kinases that are genomically proximal to a target
    gene using a proximity check. It then computes phosphorylation confidence for those
    proximal kinases on a wildtype and a mutant protein segment, returning only the cases
    with substantial changes.

    The central residue of the protein segments determines the kinase class and scoring:
    - If either protein segment has Y in the middle (position 8), tyrosine kinases are used.
      Scores are percentiles, and substantial change is defined as crossing 50 ↔ 85.
    - Otherwise, serine/threonine kinases are used. Scores are probabilities in [0, 1],
      and substantial change is defined as crossing 0.40 ↔ 0.60.

    Args:
        gene_name: Target gene symbol used for proximity checks.
        wildtype_protein_segment: Wildtype 15-character protein segment (middle position index 7).
                       Must contain Y for Y-kinase mode; otherwise ST-kinase mode is used.
        mutant_protein_segment: Mutant 15-character protein segment (middle position index 7).
                       Must contain Y for Y-kinase mode; otherwise ST-kinase mode is used.

    Returns:
        List[dict]: For each proximal kinase with substantial specificity change, a dict with:
            - "kinase": Kinase identifier (gene name for ST mode; kinase name for Y mode).
            - "change": mutant_confidence - wildtype_confidence.
            - "wildtype_value": Confidence on the wildtype protein segment.
            - "mutant_value": Confidence on the mutant protein segment.
        For Y protein segments, confidences are percentiles (0–100). For ST protein segments, confidences
        are probabilities (0.0–1.0). Results are filtered by the thresholds above and are
        not sorted.

    Raises:
        AssertionError: If a kinase domain sequence cannot be found for ST kinases.
    """
    if wildtype_protein_segment[7].upper() == 'Y' or mutant_protein_segment[7].upper() == 'Y':
        proximal_kinase: List[str] = []
        for kin in get_Y_kinase_list():
            proximity, _ = check_proximity(kin, gene_name)
            if proximity:
                proximal_kinase.append(kin)

        valid_wildtype = wildtype_protein_segment[7].upper() == 'Y'
        valid_mutant = mutant_protein_segment[7].upper() == 'Y'
        data: Dict[str, Dict[str, float]] = {}
        for kin in proximal_kinase:
            wildtype_conf = get_phosphorylation_confidence_for_Y_kinase(kin, [wildtype_protein_segment])[0] if valid_wildtype else 0.0
            mutant_conf = get_phosphorylation_confidence_for_Y_kinase(kin, [mutant_protein_segment])[0] if valid_mutant else 0.0
            change = mutant_conf - wildtype_conf
            data[kin] = {"change": change, "wildtype_value": wildtype_conf, "mutant_value": mutant_conf}
        # Apply specificity-change filter for Y protein segments (percentile scale)
        return [
            {"kinase": kinase, "change": info["change"], "wildtype_value": info["wildtype_value"], "mutant_value": info["mutant_value"]}
            for kinase, info in data.items()
            if (
                (info["wildtype_value"] < 50 and info["mutant_value"] > 85) or
                (info["wildtype_value"] > 85 and info["mutant_value"] < 50)
            )
        ]
    else:
        proximal_kinase: List[str] = []
        for kin in get_ST_kinase_list():
            kin_gene = get_gene_name(kin)
            if not kin_gene:
                continue
            proximity, _ = check_proximity(kin_gene, gene_name)
            if proximity:
                proximal_kinase.append(kin)

        valid_wildtype = wildtype_protein_segment[7].upper() in ['S', 'T']
        valid_mutant = mutant_protein_segment[7].upper() in ['S', 'T']
        data: Dict[str, Dict[str, float]] = {}
        for kin in proximal_kinase:
            for dom_seq in _get_domain_seq_from_id(kin):
                wildtype_conf = get_phosphorylation_confidence_for_ST_kinase(dom_seq, [wildtype_protein_segment])[0] if valid_wildtype else 0.0
                mutant_conf = get_phosphorylation_confidence_for_ST_kinase(dom_seq, [mutant_protein_segment])[0] if valid_mutant else 0.0
                change = mutant_conf - wildtype_conf
                if kin not in data or abs(change) > abs(data[kin]["change"]):
                    data[kin] = {"change": change, "wildtype_value": wildtype_conf, "mutant_value": mutant_conf}

        # Convert UniProt IDs to gene names and apply specificity-change filter for ST protein segments (probability scale)
        result: List[Dict[str, float]] = []
        for uniprot_id, info in data.items():
            if (
                (info["wildtype_value"] < 0.40 and info["mutant_value"] > 0.60) or
                (info["wildtype_value"] > 0.60 and info["mutant_value"] < 0.40)
            ):
                gene = get_gene_name(uniprot_id) or uniprot_id
                result.append({"kinase": gene, "change": info["change"], "wildtype_value": info["wildtype_value"], "mutant_value": info["mutant_value"]})
        return result

def _test_get_phosphorylation_confidence_for_ST_kinase() -> None:
    kinase_id = "P31749"
    kinase_sequence = _get_domain_seq_from_id(kinase_id)[0]
    protein_segments = ["RRRLSSLSASTSKSN"]
    try:
        probabilities = get_phosphorylation_confidence_for_ST_kinase(kinase_sequence, protein_segments)
        print(f"Probabilities for {kinase_id}: {probabilities}")
    except ValueError as exc:
        print(f"Validation failed for docstring example: {exc}")
    except Exception as exc:
        print(f"Request failed: {exc}")

def _test_get_substantial_change_in_specificity_with_proximity_check() -> None:
    payload = {
        "gene_name": "CDKL5",
        "wildtype_protein_segment": "DRYLTEQCLNHPTFQ",
        "mutant_protein_segment": "DRYLTEQYLNHPTFQ",
    }
    try:
        result = get_substantial_change_in_specificity_with_proximity_check(**payload)
        print(f"Specificity changes for {payload['gene_name']}: {result}")
    except ValueError as exc:
        print(f"Validation failed: {exc}")
    except Exception as exc:
        print(f"Request failed: {exc}")

def _test_get_substantial_change_in_specificity_with_proximity_check_for_tyrosine_kinases() -> None:
    payload = {
        "gene_name": "EGFR",
        "wildtype_protein_segment": "AQQCSGRYRGKSPSD",
        "mutant_protein_segment": "AQQCSGRCRGKSPSD",
    }
    try:
        result = get_substantial_change_in_specificity_with_proximity_check(**payload)
        print(f"Specificity changes for {payload['gene_name']}: {result}")
    except ValueError as exc:
        print(f"Validation failed: {exc}")
    except Exception as exc:
        print(f"Request failed: {exc}")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test-st":
        _test_get_phosphorylation_confidence_for_ST_kinase()
    elif len(sys.argv) > 1 and sys.argv[1] == "test-gene-protein-segment":
        _test_get_substantial_change_in_specificity_with_proximity_check()
    elif len(sys.argv) > 1 and sys.argv[1] == "test-gene-protein-segment-for-tyrosine-kinases":
        _test_get_substantial_change_in_specificity_with_proximity_check_for_tyrosine_kinases()
    else:
        mcp.run(transport='stdio')
