# run_mcp_server.py
import sys
import builtins
import logging
import dotenv
import os
from pathlib import Path
from contextlib import redirect_stdout
from io import StringIO

# Protect stdout for MCP protocol - temporarily suppress all stdout output during imports
_original_stdout = sys.stdout
sys.stdout = StringIO()

ENV_PATH = str(Path(__file__).resolve().parents[2] / '.env')
dotenv.load_dotenv(ENV_PATH)

from biomni.agent.a1 import A1

# Restore stdout after imports
sys.stdout = _original_stdout

# Route all prints to stderr, but leave sys.stdout alone for MCP frames
_real_print = builtins.print
def _print_to_stderr(*args, **kwargs):
    kwargs.setdefault("file", sys.stderr)
    return _real_print(*args, **kwargs)
builtins.print = _print_to_stderr

# Prefer proper logging as well
logging.basicConfig(stream=sys.stderr, level=logging.INFO)


fields = [
        "literature",
        "biochemistry",
        "bioengineering",
        "biophysics",
        "glycoengineering",
        "cancer_biology",
        "cell_biology",
        "molecular_biology",
        "genetics",
        "genomics",
        "immunology",
        "microbiology",
        "pathology",
        "pharmacology",
        "physiology",
        "synthetic_biology",
        "systems_biology",
        "support_tools"
    ]

related_fields = [
    "literature"
]

# Protect stdout during initialization to prevent MCP protocol corruption
sys.stdout = StringIO()
agent = A1()
mcp = agent.create_mcp_server(tool_modules=["biomni.tool." + field for field in related_fields])
sys.stdout = _original_stdout

if __name__ == "__main__":
    # MCP reads stdin, writes stdout — untouched
    mcp.run(transport="stdio")
