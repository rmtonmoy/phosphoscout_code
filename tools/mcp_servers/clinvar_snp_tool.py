from mcp.server.fastmcp import FastMCP
import dotenv
import json
import subprocess
from typing import Union
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

mcp = FastMCP("ClinVar_SNP_Tool")
dotenv.load_dotenv(str(PROJECT_ROOT / '.env'))
VARIANT_SUMMARY_PATH = str(PROJECT_ROOT / 'data' / 'clinvar' / 'variant_summary.txt')


@mcp.tool()
def extract_variant_record(rsid : Union[str, int]):
    """
    Retrieve structured ClinVar variant information for a given rsID
    from the variant_summary.txt dataset using a high-speed grep search.

    Parameters
    ----------
    rsid : str or int
        The numeric portion of the dbSNP identifier to search for.

    Returns
    -------
    str
        A JSON-formatted string representing a list of variant records.
    """
    rsid = str(rsid)
    grep_output = subprocess.run(
        ["grep", rsid, VARIANT_SUMMARY_PATH],
        capture_output=True,
        text=True,
        check=False
    ).stdout.strip()

    if not grep_output:
        return json.dumps([])

    matched_lines = grep_output.splitlines()

    with open(VARIANT_SUMMARY_PATH, "r", encoding="utf-8", errors="ignore") as f:
        header = f.readline().strip().split("\t")

    exclude_cols = {
        "#AlleleID", "GeneID", "RS# (dbSNP)", "nsv/esv (dbVar)", "RCVaccession",
        "PhenotypeIDS", "OtherIDs", "VariationID", "PositionVCF",
        "SCVsForAggregateGermlineClassification",
        "SCVsForAggregateSomaticClinicalImpact",
        "SCVsForAggregateOncogenicityClassification"
    }

    results = []
    for line in matched_lines:
        fields = line.split("\t")
        row_dict = dict(zip(header, fields))
        filtered = {k: v for k, v in row_dict.items() if k not in exclude_cols and v}
        results.append(filtered)

    return json.dumps(results, indent=4)


if __name__ == "__main__":
    mcp.run()
