from langchain_mcp_adapters.client import MultiServerMCPClient
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent

LOG_PATH = str(PROJECT_ROOT / "logs")
CONFIGS_PATH = str(PROJECT_ROOT / "configs")

# MCP Server paths (relative to project root)
GENERIC_TOOL_SERVER_PATH = str(PROJECT_ROOT / "tools" / "mcp_servers" / "generic_tool_server.py")
PHOSPHORYLATION_TOOLS_SERVER_PATH = str(PROJECT_ROOT / "tools" / "mcp_servers" / "phosphorylation_tools.py")
UNIPROT_TOOLS_SERVER_PATH = str(PROJECT_ROOT / "tools" / "mcp_servers" / "mcp_server_for_cms_and_uniprot.py")
BIOMNI_TOOLS_SERVER_PATH = str(PROJECT_ROOT / "tools" / "biomni_tools" / "run_mcp_server.py")
CLINVAR_SNP_TOOL_PATH = str(PROJECT_ROOT / "tools" / "mcp_servers" / "clinvar_snp_tool.py")
DSSP_TOOL_PATH = str(PROJECT_ROOT / "tools" / "mcp_servers" / "dssp_tool.py")
GENE_EXPRESSION_TOOLS_PATH = str(PROJECT_ROOT / "tools" / "mcp_servers" / "gene_expression_tools.py")
REPORT_GENERATION_TOOLS_PATH = str(PROJECT_ROOT / "tools" / "mcp_servers" / "report_generation_tools.py")

# To make biomni tools work
import builtins
import sys
_real_print = builtins.print
def _print_to_stderr(*args, **kwargs):
    kwargs.setdefault("file", sys.stderr)
    return _real_print(*args, **kwargs)
builtins.print = _print_to_stderr

MCP_CLIENT = MultiServerMCPClient(
        {
            "DSSP_Tool": {
                "command": "uv",
                "args": ["run", DSSP_TOOL_PATH],
                "transport": "stdio",
            },
            "BiomniTools": {
                "command": "uv",
                "args": ["run", BIOMNI_TOOLS_SERVER_PATH],
                "transport": "stdio",
            },
            "generic_tool_server": {
                "command": "uv",
                "args": ["run", GENERIC_TOOL_SERVER_PATH],
                "transport": "stdio",
            },
            "Phosphorylation_Tools": {
                "command": "uv",
                "args": ["run", PHOSPHORYLATION_TOOLS_SERVER_PATH],
                "transport": "stdio",
            },
            "UniProt_Tools_server": {
                "command": "uv",
                "args": ["run", UNIPROT_TOOLS_SERVER_PATH],
                "transport": "stdio",
            },
            "ClinVar_SNP_Tool": {
                "command": "uv",
                "args": ["run", CLINVAR_SNP_TOOL_PATH],
                "transport": "stdio",
            },
            "Gene_Expression_Tools": {
                "command": "uv",
                "args": ["run", GENE_EXPRESSION_TOOLS_PATH],
                "transport": "stdio",
            },
            "Report_Generation_Tools": {
                "command": "uv",
                "args": ["run", REPORT_GENERATION_TOOLS_PATH],
                "transport": "stdio",
            }
        }
    )
