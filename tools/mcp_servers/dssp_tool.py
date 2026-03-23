import os
import re
import requests
import warnings
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from mcp.server.fastmcp import FastMCP

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# Path to the mkdssp binary. Override via MKDSSP_PATH env var if the conda-installed
# mkdssp has a Boost version mismatch (common when libboost is upgraded separately).
# See README.md "Install DSSP" section for details.
MKDSSP_PATH = os.environ.get("MKDSSP_PATH", "mkdssp")

mcp = FastMCP("DSSP_Tool")

# Suppress libcifpp / DSSP warnings
warnings.filterwarnings("ignore", message=".*compound.*")
warnings.filterwarnings("ignore", message=".*components.cif.*")
warnings.filterwarnings("ignore", category=UserWarning)

def _get_dssp_for_uniprot(uniprot_id: str, save_path: str = str(PROJECT_ROOT / "structures")) -> DSSP:
    """
    Helper function to download AlphaFold structure, parse it, and run DSSP.
    Returns the DSSP object for further processing.
    """
    os.makedirs(save_path, exist_ok=True)

    pdb_url = None

    # Try UniProt ID directly first (some entries might work)
    meta_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    r = requests.get(meta_url, timeout=30)

    if r.status_code == 200:
        info = r.json()
        if info and len(info) > 0:
            pdb_url = info[0].get("pdbUrl")

    # If direct UniProt ID doesn't work, try AlphaFold model ID format: AF-{uniprot_id}-F1
    if not pdb_url:
        af_model_id = f"AF-{uniprot_id}-F1"
        meta_url = f"https://alphafold.ebi.ac.uk/api/prediction/{af_model_id}"
        r = requests.get(meta_url, timeout=30)
        if r.status_code == 200:
            info = r.json()
            if info and len(info) > 0:
                pdb_url = info[0].get("pdbUrl")

    # If API methods fail, try fallback URLs
    if not pdb_url:
        pdb_url = _try_fallback_urls(uniprot_id)

    if not pdb_url:
        raise ValueError(f"Could not find AlphaFold entry for UniProt ID: {uniprot_id}")

    pdb_file = os.path.join(save_path, os.path.basename(pdb_url))
    if not os.path.exists(pdb_file):
        resp = requests.get(pdb_url)
        if resp.status_code != 200:
            raise ValueError(f"Failed to download PDB file from {pdb_url}")
        with open(pdb_file, "wb") as f:
            f.write(resp.content)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(uniprot_id, pdb_file)
    model = list(structure)[0]

    dssp = DSSP(model, pdb_file, dssp=MKDSSP_PATH)
    if len(dssp) == 0:
        raise RuntimeError(f"DSSP returned no residues for {uniprot_id}")

    return dssp

def _try_fallback_urls(uniprot_id: str) -> str:
    """
    Fallback function to try alternative methods for finding AlphaFold structure URLs.
    Only called if the primary API endpoint fails.
    """
    af_model_id = f"AF-{uniprot_id}-F1"

    # Try direct file URLs for different model versions (v6, v5, v4, v3, v2)
    model_versions = ["v6", "v5", "v4", "v3", "v2"]
    for version in model_versions:
        direct_url = f"https://alphafold.ebi.ac.uk/files/{af_model_id}-model_{version}.pdb"
        try:
            test_resp = requests.head(direct_url, timeout=10, allow_redirects=True)
            if test_resp.status_code == 200:
                return direct_url
        except Exception:
            continue

    # Try parsing the AlphaFold entry page HTML
    entry_page_url = f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}"
    try:
        resp = requests.get(entry_page_url, timeout=30)
        if resp.status_code == 200:
            content = resp.text
            # Look for PDB URLs in various formats
            pdb_patterns = [
                rf'https://alphafold\.ebi\.ac\.uk/files/{af_model_id}[^"\'<>\s]*\.pdb',
                rf'https://alphafold\.ebi\.ac\.uk/files/AF-{uniprot_id}-F1[^"\'<>\s]*\.pdb',
                rf'/files/{af_model_id}[^"\'<>\s]*\.pdb',
            ]
            for pattern in pdb_patterns:
                matches = re.findall(pattern, content)
                if matches:
                    url = matches[0]
                    if url.startswith('/'):
                        url = 'https://alphafold.ebi.ac.uk' + url
                    return url
    except Exception:
        pass

    return None

@mcp.tool()
def get_solvent_accessible_residues_from_uniprot(
    uniprot_id: str,
    save_path: str = str(PROJECT_ROOT / "structures"),
    rsa_threshold: float = 0.25
) -> dict:
    """
    Given a UniProt ID, download its AlphaFold structure, run DSSP,
    and return residues whose normalized solvent accessibility >= rsa_threshold.

    In your environment DSSP output is a flat tuple:
        (index, aa, ss, ASA, phi, psi, ...)
    where ASA is already normalized to [0-1].
    """
    dssp = _get_dssp_for_uniprot(uniprot_id, save_path)

    exposed_positions = []
    for record in dssp:
        pos = record[0]
        asa = record[3]

        if asa >= rsa_threshold:
            exposed_positions.append(pos)

    return {
        "uniprot_id": uniprot_id,
        "rsa_threshold": rsa_threshold,
        "total_residues": len(exposed_positions),
        "solvent_accessible_residues": exposed_positions,
    }

@mcp.tool()
def get_solvent_accesibility(
    uniprot_id: str,
    central_residue: int,
    save_path: str = str(PROJECT_ROOT / "structures"),
    rsa_threshold: float = 0.25
) -> list:
    """
    Given a UniProt ID and a central residue position, download its AlphaFold structure,
    run DSSP, and return a list of tuples for nearby residues (+/-10 from central_residue).

    Each tuple contains: (position, rsa_value, flag)
    where flag is True if rsa_value >= rsa_threshold, False otherwise.

    In your environment DSSP output is a flat tuple:
        (index, aa, ss, ASA, phi, psi, ...)
    where ASA is already normalized to [0-1].
    """
    dssp = _get_dssp_for_uniprot(uniprot_id, save_path)

    position_to_rsa = {}
    for record in dssp:
        pos = record[0]
        asa = record[3]
        position_to_rsa[pos] = asa

    result = []
    start_pos = central_residue - 10
    end_pos = central_residue + 10

    for pos in range(start_pos, end_pos + 1):
        if pos in position_to_rsa:
            rsa_value = position_to_rsa[pos]
            flag = rsa_value >= rsa_threshold
            result.append((pos, rsa_value, flag))

    return result

def _test_dssp_tool() -> None:
    """
    Test function for dssp_tool functionality.

    Tests:
    1. Primary API endpoint with a UniProt ID that should work
    2. Fallback mechanism with P00533 (EGFR) that was failing
    3. get_solvent_accesibility tool function
    4. get_solvent_accessible_residues_from_uniprot tool function

    Usage:
        python dssp_tool.py test
    """
    import sys

    print("\n=== Testing DSSP Tool ===\n")
    sys.stdout.flush()

    test_cases = [
        {
            "name": "Test fallback with P00533 (EGFR) - the failing case",
            "uniprot_id": "P00533",
            "test_type": "fallback",
        },
        {
            "name": "Test primary API with P04637 (TP53)",
            "uniprot_id": "P04637",
            "test_type": "primary",
        },
    ]

    print("--- Test 1: Testing _get_dssp_for_uniprot function ---")
    for i, test_case in enumerate(test_cases, 1):
        print(f"\n  Test {i}.{len(test_cases)}: {test_case['name']}")
        print(f"  UniProt ID: {test_case['uniprot_id']}")
        print(f"  Expected method: {test_case['test_type']}")
        sys.stdout.flush()

        try:
            save_path = "/tmp/test_structures"
            os.makedirs(save_path, exist_ok=True)

            dssp = _get_dssp_for_uniprot(test_case['uniprot_id'], save_path)

            residue_count = len(dssp) if dssp else 0
            print(f"  Success! DSSP returned {residue_count} residues")

            if residue_count > 0:
                first_residue = list(dssp)[0]
                print(f"  First residue: position={first_residue[0]}, AA={first_residue[1]}, ASA={first_residue[3]:.3f}")

            sys.stdout.flush()

        except Exception as e:
            print(f"  ERROR: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.stdout.flush()

    print("\n--- Test 2: Testing fallback function _try_fallback_urls directly ---")
    sys.stdout.flush()

    test_fallback_cases = [
        {"uniprot_id": "P00533", "description": "EGFR (was failing)"},
        {"uniprot_id": "P84243", "description": "H3F3A (from logs)"},
    ]

    for i, test_case in enumerate(test_fallback_cases, 1):
        print(f"\n  Fallback test {i}.{len(test_fallback_cases)}: {test_case['uniprot_id']} ({test_case['description']})")
        sys.stdout.flush()

        try:
            pdb_url = _try_fallback_urls(test_case['uniprot_id'])
            if pdb_url:
                print(f"  Found URL: {pdb_url}")
            else:
                print(f"  No URL found via fallback")
            sys.stdout.flush()
        except Exception as e:
            print(f"  ERROR: {str(e)}")
            sys.stdout.flush()

    print("\n--- Test 3: Testing get_solvent_accesibility tool ---")
    sys.stdout.flush()

    tool_test_cases = [
        {
            "uniprot_id": "P00533",
            "central_residue": 719,
            "description": "EGFR at residue 719 (G719C mutation site)",
        },
        {
            "uniprot_id": "Q6IQ55",
            "central_residue": 1182,
            "description": "Q6IQ55 at residue 1182",
        },
    ]

    for i, test_case in enumerate(tool_test_cases, 1):
        print(f"\n  Tool test {i}.{len(tool_test_cases)}: {test_case['description']}")
        print(f"  UniProt ID: {test_case['uniprot_id']}, Central residue: {test_case['central_residue']}")
        sys.stdout.flush()

        try:
            save_path = "/tmp/test_structures"
            result = get_solvent_accesibility(
                uniprot_id=test_case['uniprot_id'],
                central_residue=test_case['central_residue'],
                save_path=save_path,
                rsa_threshold=test_case.get("rsa_threshold", 0.1)
            )

            print(f"  Success! Returned {len(result)} residues")
            if result:
                print(f"  Sample results (first 5):")
                for pos, rsa, flag in result[:5]:
                    print(f"    Position {pos}: RSA={rsa:.3f}, Accessible={flag}")

            sys.stdout.flush()

        except Exception as e:
            print(f"  ERROR: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.stdout.flush()

    print("\n--- Test 4: Testing get_solvent_accessible_residues_from_uniprot tool ---")
    sys.stdout.flush()

    accessible_test_cases = [
        {
            "uniprot_id": "P00533",
            "description": "EGFR",
        },
    ]

    for i, test_case in enumerate(accessible_test_cases, 1):
        print(f"\n  Accessible residues test {i}.{len(accessible_test_cases)}: {test_case['uniprot_id']} ({test_case['description']})")
        sys.stdout.flush()

        try:
            save_path = "/tmp/test_structures"
            result = get_solvent_accessible_residues_from_uniprot(
                uniprot_id=test_case['uniprot_id'],
                save_path=save_path,
                rsa_threshold=0.25
            )

            print(f"  Success!")
            print(f"  UniProt ID: {result['uniprot_id']}")
            print(f"  RSA threshold: {result['rsa_threshold']}")
            print(f"  Total accessible residues: {result['total_residues']}")
            if result['solvent_accessible_residues']:
                print(f"  Sample accessible positions: {result['solvent_accessible_residues'][:10]}")

            sys.stdout.flush()

        except Exception as e:
            print(f"  ERROR: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.stdout.flush()

    print("\n=== All tests completed ===\n")
    sys.stdout.flush()

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "test":
        _test_dssp_tool()
    else:
        mcp.run(transport="stdio")
