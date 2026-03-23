import requests
import os
import json
import re
from pathlib import Path
from typing import Set

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# Point to the cache directory within the project
CACHE_DIR = str(PROJECT_ROOT / 'data' / 'cache' / 'kinase_proximity')
os.makedirs(CACHE_DIR, exist_ok=True)

UNIPROT_CACHE_FILE = os.path.join(CACHE_DIR, "uniprot_id_cache.json")
SUBCELL_CACHE_FILE = os.path.join(CACHE_DIR, "subcellular_locations_cache.json")
PATHWAYS_CACHE_FILE = os.path.join(CACHE_DIR, "pathways_cache.json")
PROXIMITY_CACHE_FILE = os.path.join(CACHE_DIR, "proximity_cache.json")
HPA_CACHE_FILE = os.path.join(CACHE_DIR, "hpa_subcellular_locations_cache.json")


def _load_cache(file_path):
    try:
        if os.path.exists(file_path):
            with open(file_path, "r") as f:
                return json.load(f)
    except Exception:
        return {}
    return {}


def _save_cache(file_path, cache_obj):
    try:
        with open(file_path, "w") as f:
            json.dump(cache_obj, f)
    except Exception:
        # Cache failures should not kill the program.
        pass


_uniprot_cache = _load_cache(UNIPROT_CACHE_FILE)
_subcell_cache = _load_cache(SUBCELL_CACHE_FILE)
_pathways_cache = _load_cache(PATHWAYS_CACHE_FILE)
_proximity_cache = _load_cache(PROXIMITY_CACHE_FILE)
_hpa_cache = _load_cache(HPA_CACHE_FILE)


def get_uniprot_id(gene_name, organism="9606"):
    """
    Find the UniProt accession for a given gene name (symbol) and organism.
    Defaults to human (NCBI taxonomy ID 9606). Returns the primary accession
    or None if not found.

    Network / HTTP / JSON errors are allowed to raise
    (e.g. if IP is blocked or UniProt is unreachable).
    """
    # Check cache first
    cache_key = f"{gene_name}|{organism}"
    if cache_key in _uniprot_cache:
        return _uniprot_cache[cache_key]

    # Construct the query URL for UniProt search API (reviewed Swiss-Prot entries only for accuracy)
    query = f"gene:{gene_name}+AND+organism_id:{organism}+AND+reviewed:true"
    url = (
        "https://rest.uniprot.org/uniprotkb/search"
        f"?query={query}&fields=accession&format=json"
    )

    # Let RequestException / HTTPError / JSONDecodeError bubble up
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    data = response.json()

    # Parse the JSON to get the first accession
    results = data.get("results")
    if results:
        acc = results[0].get("primaryAccession")
        _uniprot_cache[cache_key] = acc
        _save_cache(UNIPROT_CACHE_FILE, _uniprot_cache)
        return acc
    return None


def get_subcellular_locations(uniprot_id):
    """
    Retrieve a set of subcellular compartment names for the protein
    with given UniProt ID. Uses UniProt API to get the
    "Subcellular location" annotations.

    Network / HTTP / JSON errors are allowed to raise.
    """
    # Check cache first
    if uniprot_id in _subcell_cache:
        return set(_subcell_cache[uniprot_id])

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    locations = set()

    # Let RequestException / HTTPError / JSONDecodeError bubble up
    r = requests.get(url, timeout=10)
    r.raise_for_status()
    entry = r.json()

    # UniProt JSON 'comments' may contain subcellular location info
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            # Each such comment can have one or multiple subcellular locations listed
            if "subcellularLocations" in comment:
                # Newer UniProt JSON structure might use 'subcellularLocations'
                for loc in comment["subcellularLocations"]:
                    # Each loc might have a "location" dict with a 'value'
                    if "location" in loc and "value" in loc["location"]:
                        locations.add(loc["location"]["value"])
            elif "locations" in comment:
                # Some JSON structures use 'locations'
                for loc in comment["locations"]:
                    if "location" in loc and "value" in loc["location"]:
                        locations.add(loc["location"]["value"])
            # Also check text fields (for older format or additional notes)
            if "text" in comment:
                text = comment["text"].get("value", "")
                # Sometimes location info might be embedded in text (e.g. isoform notes)
                # Simple approach: split by ';' and look for possible location-like tokens.
                for part in text.split(';'):
                    part = part.strip()
                    if part:
                        words = part.split()
                        if words:
                            locations.add(words[0])

    # Persist to cache as a list
    _subcell_cache[uniprot_id] = sorted(list(locations))
    _save_cache(SUBCELL_CACHE_FILE, _subcell_cache)
    return locations


def get_pathways(uniprot_id, species="9606"):
    """
    Retrieve a set of pathway names that the protein (UniProt ID) is involved in.
    Uses Reactome Content Service API to map UniProt ID to pathway list
    for the given species.

    All Reactome request / HTTP / JSON errors are ignored and treated
    as empty results, since upstream connectivity can be unstable.
    """
    # Check cache first
    cache_key = f"{uniprot_id}|{species}"
    if cache_key in _pathways_cache:
        return set(_pathways_cache[cache_key])

    pathways = set()
    url = (
        "https://reactome.org/ContentService/data/mapping/UniProt/"
        f"{uniprot_id}/pathways?species={species}"
    )
    headers = {"Accept": "application/json"}

    try:
        res = requests.get(url, headers=headers, timeout=10)

        # If Reactome has no pathways for this UniProt (404), treat as empty
        # and do not cache so that future requests can re-check.
        if res.status_code == 404:
            return pathways

        res.raise_for_status()
        data = res.json()
    except Exception:
        return pathways

    # The returned data is likely a list of pathway objects with names
    for pathway in data:
        name = pathway.get("displayName") or pathway.get("name") or pathway.get("stId")
        if name:
            pathways.add(name)

    # Persist to cache as a list
    _pathways_cache[cache_key] = sorted(list(pathways))
    _save_cache(PATHWAYS_CACHE_FILE, _pathways_cache)
    return pathways


HPA_SEARCH_URL = "https://www.proteinatlas.org/api/search_download.php"


def get_hpa_subcellular_locations(gene_name: str, species: str = "9606") -> Set[str]:
    """
    Fetch subcellular location annotations for a gene from the Human Protein Atlas.

    This uses the search_download.php API with:
        search=<gene_name>, format=json, columns=g,scl, compress=no

    Because `search` is a free-text query, the response may contain many genes.
    We therefore:
      - Filter rows where row["Gene"] == gene_name (exact match).
      - Parse the "Subcellular location" field, which may be:
          * a list of strings, or
          * a single string with ';' or ',' separators.

    Returns:
        A set of subcellular location strings (may be empty if no matching row is found).

    Network / HTTP / JSON errors are allowed to raise
    (e.g. if HPA is unreachable or your IP is blocked).
    """
    # HPA subcellular is human-centric; for non-human species, just skip.
    if str(species) not in ("9606", "Homo sapiens", "human"):
        return set()

    # Check cache first
    cache_key = f"{gene_name}|{species}"
    if cache_key in _hpa_cache:
        return set(_hpa_cache[cache_key])

    params = {
        "search": gene_name,
        "format": "json",
        "columns": "g,scl",
        "compress": "no",
    }

    # Let RequestException / HTTPError / JSONDecodeError bubble up
    resp = requests.get(HPA_SEARCH_URL, params=params, timeout=10)
    resp.raise_for_status()
    data = resp.json()

    locs: Set[str] = set()
    for row in data:
        if row.get("Gene") != gene_name:
            continue

        scl = row.get("Subcellular location")
        if not scl:
            continue

        if isinstance(scl, list):
            for entry in scl:
                if entry:
                    locs.add(entry.strip())
        elif isinstance(scl, str):
            # Handle formats like "Nucleus; Cytosol" or "Nucleus, Cytosol"
            for entry in re.split(r"[;,]", scl):
                entry = entry.strip()
                if entry:
                    locs.add(entry)

    # Persist to cache as a list
    _hpa_cache[cache_key] = sorted(list(locs))
    _save_cache(HPA_CACHE_FILE, _hpa_cache)
    return locs


def check_proximity(kinase_name, substrate_name, species="9606"):
    """
    Determine if kinase and substrate are in the same pathway or subcellular compartment.
    Returns (proximity_boolean, evidence_string).

    - kinase_name: gene name of the kinase (e.g. 'AKT1')
    - substrate_name: gene name of the substrate (e.g. 'GSK3B')
    - species: taxonomy ID or name (default 9606 for human)

    Logic:
      1) Fetch subcellular locations from HPA using gene symbols (human only).
      2) Map gene symbols to UniProt accessions (if found).
      3) Fetch subcellular locations from UniProt (if IDs exist).
      4) Union HPA + UniProt locations.
      5) Fetch pathway memberships (Reactome) when UniProt IDs exist.
      6) Check for overlaps in compartments/pathways.
      7) Cache and return proximity + textual evidence.

    Any network / HTTP / JSON error in the upstream calls is allowed to raise
    so that repeated IP blocking / request failures are visible.
    """
    # Check cache first (positions are not used, so exclude them from the key)
    cache_key = f"{kinase_name}|{substrate_name}|{species}"
    if cache_key in _proximity_cache:
        cached = _proximity_cache[cache_key]
        return bool(cached.get("proximity", False)), cached.get("evidence", "")

    # Step 1: Get subcellular locations (HPA; using gene symbols, then filter by gene)
    kinase_locs_hpa = get_hpa_subcellular_locations(kinase_name, species=species)
    substrate_locs_hpa = get_hpa_subcellular_locations(substrate_name, species=species)

    # Step 2: Map gene names to UniProt accessions
    kinase_id = get_uniprot_id(kinase_name, organism=species)
    substrate_id = get_uniprot_id(substrate_name, organism=species)

    # Step 3: Get subcellular locations (UniProt), treating missing IDs as empty sets
    if kinase_id:
        kinase_locs_uniprot = get_subcellular_locations(kinase_id)
    else:
        kinase_locs_uniprot = set()

    if substrate_id:
        substrate_locs_uniprot = get_subcellular_locations(substrate_id)
    else:
        substrate_locs_uniprot = set()

    # Step 4: Combine UniProt + HPA locations
    kinase_locs = kinase_locs_uniprot.union(kinase_locs_hpa)
    substrate_locs = substrate_locs_uniprot.union(substrate_locs_hpa)

    # Step 5: Get pathway memberships for each protein (only if UniProt IDs exist)
    # If Reactome has no pathways or returns a not-found response, keep going.
    if kinase_id:
        try:
            kinase_paths = get_pathways(kinase_id, species=species)
        except Exception:
            kinase_paths = set()
    else:
        kinase_paths = set()

    if substrate_id:
        try:
            substrate_paths = get_pathways(substrate_id, species=species)
        except Exception:
            substrate_paths = set()
    else:
        substrate_paths = set()

    # Step 6: Check overlaps
    shared_compartments = kinase_locs.intersection(substrate_locs)
    shared_pathways = kinase_paths.intersection(substrate_paths)

    # Step 7: Determine proximity and compile evidence
    if shared_compartments or shared_pathways:
        proximity = True
        evidence_parts = []

        if shared_compartments:
            comp_list = ", ".join(sorted(shared_compartments))
            evidence_parts.append(
                "common compartment(s): "
                f"{comp_list} (UniProt + HPA subcellular)"
            )

        if shared_pathways:
            path_list = ", ".join(sorted(shared_pathways))
            evidence_parts.append(
                f"shared pathway(s): {path_list} (Reactome)"
            )

        evidence = "Kinase and substrate co-localize: " + "; ".join(evidence_parts)
    else:
        proximity = False
        evidence = (
            "No shared pathways or compartments found for the given "
            "kinase-substrate pair."
        )

    # Persist computed result
    _proximity_cache[cache_key] = {"proximity": proximity, "evidence": evidence}
    _save_cache(PROXIMITY_CACHE_FILE, _proximity_cache)
    return proximity, evidence


def _test_check_proximity() -> None:
    """
    Quick manual test: run check_proximity with known kinase-substrate pairs.

    Usage:
        python kinase_proximity.py test-proximity

    This test checks proximity for several kinase-substrate pairs:
      - AKT1 and GSK3B (known interaction)
      - MAPK1 and MAPK3 (related kinases)
      - A pair that likely doesn't share pathways/compartments
    """
    import sys

    test_cases = [
        {
            "name": "AKT1 and GSK3B (known interaction)",
            "kinase": "AKT1",
            "substrate": "GSK3B",
            "species": "9606",
        },
        {
            "name": "MAPK1 and MAPK3 (related kinases)",
            "kinase": "MAPK1",
            "substrate": "MAPK3",
            "species": "9606",
        },
        {
            "name": "CDK2 and TP53 (cell cycle regulation)",
            "kinase": "CDK2",
            "substrate": "TP53",
            "species": "9606",
        },
    ]

    print("\n=== Testing check_proximity ===\n")
    sys.stdout.flush()

    for i, test_case in enumerate(test_cases, 1):
        print(f"\n--- Test {i}: {test_case['name']} ---")
        print(f"Kinase: {test_case['kinase']}")
        print(f"Substrate: {test_case['substrate']}")
        print(f"Species: {test_case['species']}")
        sys.stdout.flush()

        try:
            proximity, evidence = check_proximity(
                kinase_name=test_case["kinase"],
                substrate_name=test_case["substrate"],
                species=test_case["species"],
            )

            print(f"\nProximity: {proximity}")
            print(f"Evidence: {evidence}")
            sys.stdout.flush()

        except Exception as e:
            print(f"ERROR: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.stdout.flush()

    print("\n=== All tests completed ===\n")
    sys.stdout.flush()


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "test-proximity":
        _test_check_proximity()
    else:
        print("Usage: python kinase_proximity.py test-proximity")
        sys.exit(1)
