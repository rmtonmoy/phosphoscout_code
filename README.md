# PhosphoScout

An agentic system for analyzing protein-coding mutations and their potential impact on phosphorylation post-translational modification (PTM) patterns. The system orchestrates multiple AI agents to collect mutation data, generate mechanistic hypotheses, and validate them against gene expression data.

## Architecture

The system uses a manager-worker agent pattern built on LangGraph:

1. **agent_manager** — Orchestrates the entire pipeline:
   - Calls `mutation_data_collector` to gather variant context, phosphosite windows, and kinase specificity changes
   - Calls `mutation_impact_expert` to generate novel, testable mechanistic hypotheses
   - Calls `gene_expression_expert` to validate hypotheses against single-cell gene expression data (CELLxGENE Census)
   - Performs iterative refinement (up to 3 iterations) if gene expression analysis does not strongly support hypotheses
   - Writes structured JSON reports and HTML reports
   - Optionally sends email notifications with hypothesis summaries

2. **mutation_data_collector** — Collects mutation data using tools for:
   - Protein segment extraction and phosphosite window identification
   - Solvent accessibility analysis (DSSP on AlphaFold structures)
   - Kinase-substrate specificity prediction (Phosformer + Kinase Library)
   - Literature search (PubMed, Google Scholar, arXiv, web)
   - ClinVar variant lookup

3. **mutation_impact_expert** — Generates mechanistic hypotheses connecting mutations to disease phenotypes through phosphorylation rewiring, grounded in literature evidence

4. **gene_expression_expert** — Evaluates hypotheses against differential gene expression data from CELLxGENE Census (single-cell RNA-seq)

## Setup

### 1. Create Conda Environment

```bash
conda create -n phosphoscout python=3.11 -y
conda activate phosphoscout
```

### 2. Install Dependencies

```bash
pip install -r requirements.txt
```

**Note:** The `kinase_library` package is included locally under `tools/kinase-library/`. If it requires separate installation:
```bash
cd tools/kinase-library
pip install -e .
cd ../..
```

### 3. Install DSSP

DSSP (`mkdssp`) is required for solvent accessibility calculations on AlphaFold structures. There are several ways to install it:

**Option A — conda (recommended):**
```bash
conda install -c salilab dssp
```

**Option B — Build from source:**
```bash
git clone https://github.com/PDB-REDO/dssp.git
cd dssp
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
cmake --build build
cmake --install build
```
This requires `libcifpp`, `libmcfp`, and Boost (thread, filesystem, program_options, iostreams) as build dependencies. See https://github.com/PDB-REDO/dssp for full instructions.

**Option C — System package manager (Ubuntu/Debian):**
```bash
sudo apt-get install dssp
```

#### Troubleshooting: Boost version mismatch

A common issue is that `mkdssp` was compiled against a specific Boost version (e.g., `libboost_thread.so.1.73.0`) but a different version is installed in your conda environment. This manifests as:
```
mkdssp: error while loading shared libraries: libboost_thread.so.1.73.0: cannot open shared object file
```

To fix this, set the `MKDSSP_PATH` environment variable to point to a working `mkdssp` binary. For example, if you have a working installation in another conda environment:
```bash
export MKDSSP_PATH=/path/to/working/conda/envs/myenv/bin/mkdssp
```

Or add it to your `.env` file:
```
MKDSSP_PATH=/path/to/working/conda/envs/myenv/bin/mkdssp
```

You can verify your `mkdssp` works by running:
```bash
mkdssp --version
```
If it prints a version without errors, it is correctly linked.

### 4. Install uv (MCP Server Runner)

The MCP servers are launched using `uv`. Install it:
```bash
pip install uv
```

### 5. Configure Environment Variables

```bash
cp .env.example .env
```

Edit `.env` and fill in:
- `OPENAI_API_KEY` — Required. Used by all agents for LLM calls.
- `TAVILY_API_KEY` — Required. Used by the web_search tool.
- `GMAIL_TOKEN` — Optional. JSON token for Gmail email notifications.
- `GMAIL_SENDER_EMAIL` — Optional. Sender email address for notifications.
- `MKDSSP_PATH` — Optional. Path to a working `mkdssp` binary. Defaults to `mkdssp` (found on PATH). Set this if you encounter Boost library version mismatches (see [Install DSSP](#3-install-dssp)).

### 6. Phosformer Server

The phosphorylation tools for S/T kinases require a running Phosformer prediction server at `http://10.2.4.15:5000`. Ensure the server is accessible from your environment, or update the server URL in `tools/mcp_servers/phosphorylation_tools.py`.

## Usage

### Sequential Processing (Default)

Process all mutations in `data/mutations_to_run.txt` one at a time:

```bash
conda activate phosphoscout
python invoke_agent_manager.py
```

### Parallel Processing

Process mutations in parallel using multiprocessing:

```bash
python invoke_agent_manager.py --parallel
```

Optionally specify the number of worker processes:

```bash
python invoke_agent_manager.py --parallel 4
```

### Input Format

The mutations input file (`data/mutations_to_run.txt`) is a JSON array of mutation objects. Each object can contain:

**Type A — Full COSMIC-style input:**
```json
{
  "gene_name": "EGFR",
  "accession_number": "ENST00000275493.2",
  "mutation_cds": "c.2573T>G",
  "mutation_aa": "p.L858R",
  "mutation_description_aa": "Substitution - Missense",
  "aa_mut_start": 858,
  "aa_mut_stop": 858,
  "mutation_description_cds": "Substitution"
}
```

**Type B — Gene + mutation only:**
```json
{
  "gene_name": "ACVR1",
  "mutation_aa": "p.R206H"
}
```

### Output

Reports are generated in:
- `generated_artifacts/reports/raw_reports/` — Structured JSON reports
- `generated_artifacts/reports/html_reports/` — Formatted HTML reports

## Data Files (Not Included in Repository)

The `data/` directory contains several large data files that are **not tracked in git**. Below is a description of each and how to obtain them.

### `data/clinvar/variant_summary.txt` (~3.2 GB) — NCBI ClinVar
The full ClinVar variant summary table (tab-separated). Used by the `clinvar_snp_tool.py` to look up variant records by rsID.

**Download from:** https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
**Steps:**
```bash
mkdir -p data/clinvar
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz -O data/clinvar/variant_summary.txt.gz
gunzip data/clinvar/variant_summary.txt.gz
```
**Location:** `data/clinvar/variant_summary.txt`

### `data/phosphosite/Homo_sapiens.txt` (~221 MB) — PhosphoSitePlus (EPSD)
Experimentally validated phosphorylation sites for human proteins. Tab-separated file with columns: EPSD ID, UniProt ID, AA, Position, Source, Reference. Used by `mcp_server_for_cms_and_uniprot.py` to retrieve known phosphosites for a given UniProt ID.

**Download from:** https://epsd.biocuckoo.cn/ (EPSD database — Eukaryotic Phosphorylation Site Database)
**Location:** `data/phosphosite/Homo_sapiens.txt`

### `data/kinase_data/41586_2022_5575_MOESM3_ESM.xlsx` (~129 KB) — Cantley Lab Kinase Atlas
Supplementary Table 3 from: Johnson, J. L., et al. (2023). "An atlas of substrate specificities for the human serine/threonine kinome." *Nature*, 613(7945), 759-766. Contains curated serine/threonine kinase data used for building the well-studied kinase list.

**Download from:** Supplementary materials of https://doi.org/10.1038/s41586-022-05575-3
**Location:** `data/kinase_data/41586_2022_5575_MOESM3_ESM.xlsx`

### `data/kinase_data/pkfold_hs_curated.fa` (~320 KB) — Curated Kinase Domain Sequences (INCLUDED IN REPO)
A FASTA file of curated human kinase catalytic domain sequences, indexed by UniProt ID. Used by `phosphorylation_tools.py` and `well_studied_kinase_utils.py` for kinase domain sequence lookups. This file is **included in the repository** as it is a curated dataset from the Kannan Lab and is not available for download elsewhere.

**Location:** `data/kinase_data/pkfold_hs_curated.fa`

### `data/cache/` (~14 MB) — Auto-Generated Runtime Caches
These JSON cache files are **automatically generated at runtime** to avoid repeated API calls. They do not need to be provided manually — they will be created on first run.

- `cache/uniprot/estn_to_uniprot.json` — Ensembl transcript ID to UniProt accession mapping
- `cache/uniprot/uniprot_to_seq.json` — UniProt ID to protein sequence mapping
- `cache/uniprot/uniprot_to_domain_seq.json` — UniProt ID to kinase domain sequence mapping
- `cache/kinase_proximity/uniprot_id_cache.json` — Gene name to UniProt accession mapping
- `cache/kinase_proximity/subcellular_locations_cache.json` — UniProt subcellular location data
- `cache/kinase_proximity/hpa_subcellular_locations_cache.json` — Human Protein Atlas subcellular locations
- `cache/kinase_proximity/pathways_cache.json` — Reactome pathway data
- `cache/kinase_proximity/proximity_cache.json` — Kinase-substrate proximity scores

**Location:** `data/cache/`

### `data/mutations_to_run.txt` — User Input File
JSON array of mutation objects to analyze. This is your input file — create it with the mutations you want to process. See the [Input Format](#input-format) section above for the expected schema. A sample file with common cancer mutations is not included in the repository.

**Location:** `data/mutations_to_run.txt`

## Project Structure

```
phosphoscout/
├── invoke_agent_manager.py      # Main entry point
├── constants.py                 # MCP server configuration
├── requirements.txt             # Python dependencies
├── .env                         # API keys (create from .env.example)
├── mini_graphs/
│   └── agent_manager.py         # LangGraph agent graph
├── implementations/
│   ├── helpers.py               # Agent factory and tool loading
│   ├── agent_blueprint.py       # Agent config model
│   ├── generate_docstring.py    # Docstring generator for toolified agents
│   └── agent_toolifying_prompt.txt
├── configs/
│   ├── agents/                  # Agent YAML configurations
│   │   ├── agent_manager.yaml
│   │   ├── mutation_data_collector.yaml
│   │   ├── mutation_impact_expert.yaml
│   │   └── gene_expression_expert.yaml
│   └── docstrings/              # Cached generated docstrings
├── tools/
│   ├── mcp_servers/             # MCP tool servers
│   │   ├── generic_tool_server.py        # File ops, code execution, API query
│   │   ├── literature_tools.py           # Web search (Tavily) + PubMed/Scholar/arXiv/URL/PDF
│   │   ├── phosphorylation_tools.py      # Kinase specificity prediction
│   │   ├── mcp_server_for_cms_and_uniprot.py  # UniProt tools
│   │   ├── dssp_tool.py                  # Solvent accessibility (DSSP)
│   │   ├── clinvar_snp_tool.py           # ClinVar variant lookup
│   │   ├── gene_expression_tools.py      # CELLxGENE Census tools
│   │   ├── report_generation_tools.py    # HTML report generation
│   │   ├── kinase_proximity.py           # Kinase-substrate proximity
│   │   ├── uniprot_utils.py              # UniProt API utilities
│   │   ├── well_studied_kinase_utils.py  # Curated kinase data
│   │   └── cache_manager.py              # JSON cache management
│   └── kinase-library/          # Kinase Library package
├── data/
│   ├── mutations_to_run.txt     # Input mutation list
│   ├── clinvar/
│   │   └── variant_summary.txt  # ClinVar variant database
│   ├── phosphosite/
│   │   └── Homo_sapiens.txt     # PhosphoSitePlus data
│   ├── kinase_data/
│   │   ├── 41586_2022_5575_MOESM3_ESM.xlsx  # Cantley kinase data
│   │   └── pkfold_hs_curated.fa              # Curated kinase domain sequences
│   └── cache/                   # Runtime caches
│       ├── kinase_proximity/
│       └── uniprot/
├── structures/                  # Downloaded PDB structures (auto-populated)
├── generated_artifacts/         # Output reports
│   └── reports/
│       ├── raw_reports/
│       └── html_reports/
└── logs/
```