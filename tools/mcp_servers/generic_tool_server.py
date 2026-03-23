"""
MCP server with web search, file operations, and code execution tools.

Available Tools:
- web_search: Search the public web for up-to-date information and sources
- code_executor: Execute Python code with timeout and security restrictions
- file_writer: Write content to a file with optional directory creation
- ls: List files and directories in a specified path
- file_reader: Read content from a file with configurable encoding
- get_all_agent_data: Get all file contents from a directory concatenated into a single string
- api_query: Query an API endpoint and return JSON response
"""

import os
from pathlib import Path
from mcp.server.fastmcp import FastMCP
from langchain_tavily import TavilySearch
import dotenv
from typing import Union, Optional
import io
import contextlib
import time
import traceback
import asyncio
import builtins
import requests
import json

PROJECT_ROOT = Path(__file__).resolve().parents[2]

dotenv.load_dotenv(str(PROJECT_ROOT / '.env'))
mcp = FastMCP("GenericToolServer")

@mcp.tool()
def web_search(query: Union[str, dict]) -> str:
    """Search the public web for up-to-date information and sources

    Args:
        query: A concise query string or a dictionary with a key-value pair where the key is the query

    Returns:
        Summarized results with citations
    """
    # Extract query string from input
    if isinstance(query, dict):
        if not query:
            raise ValueError("Dictionary input cannot be empty")
        # Use the first key as the query
        search_query = list(query.keys())[0]
    else:
        search_query = query

    tavily_search = TavilySearch(
        max_results=5,
        topic="general",
    )

    # Get the search results
    results = tavily_search.invoke(search_query)

    # Handle different response formats
    if isinstance(results, str):
        return results
    elif isinstance(results, dict):
        # Format the results as a readable string
        formatted_results = f"Search results for: {search_query}\n\n"

        if "results" in results:
            for i, result in enumerate(results["results"], 1):
                formatted_results += f"{i}. "
                if "title" in result:
                    formatted_results += f"{result['title']}\n"
                if "content" in result:
                    formatted_results += f"   {result['content']}\n"
                if "url" in result:
                    formatted_results += f"   Source: {result['url']}\n"
                formatted_results += "\n"
        else:
            # If the structure is different, convert the dict to string
            formatted_results += str(results)

        return formatted_results
    else:
        # Fallback for any other type
        return str(results)




@mcp.tool()
def file_writer(path: str, content: str, create_dirs: bool = True, append: bool = False) -> dict:
    """Write content to a file

    Args:
        path: File path to write to
        content: Content to write to the file
        create_dirs: Create parent directories if they don't exist (default: True)
        append: Append content to existing file instead of overwriting (default: False)

    Returns:
        Dictionary with status and file path
    """
    file_path = Path(path)

    if create_dirs:
        file_path.parent.mkdir(parents=True, exist_ok=True)

    if append:
        with file_path.open('a') as f:
            f.write(content)
    else:
        file_path.write_text(content)

    return {
        "status": "success",
        "path": str(file_path.absolute()),
        "bytes_written": len(content)
    }


@mcp.tool()
def ls(path: str, show_hidden: bool = False) -> dict:
    """List files and directories in a directory

    Args:
        path: Directory path to list
        show_hidden: Include hidden files/directories (default: False)

    Returns:
        Dictionary with status, path, and list of files/directories
    """
    dir_path = Path(path)

    if not dir_path.exists():
        return {
            "status": "error",
            "message": f"Path does not exist: {path}"
        }

    if not dir_path.is_dir():
        return {
            "status": "error",
            "message": f"Path is not a directory: {path}"
        }

    items = []
    for item in dir_path.iterdir():
        if not show_hidden and item.name.startswith('.'):
            continue
        items.append({
            "name": item.name,
            "type": "directory" if item.is_dir() else "file",
            "path": str(item.absolute())
        })

    items.sort(key=lambda x: (x["type"] != "directory", x["name"]))

    return {
        "status": "success",
        "path": str(dir_path.absolute()),
        "count": len(items),
        "items": items
    }


@mcp.tool()
def file_reader(path: str, encoding: str = "utf-8") -> dict:
    """Read content from a file

    Args:
        path: File path to read from
        encoding: File encoding (default: utf-8)

    Returns:
        Dictionary with status, file path, and content
    """
    file_path = Path(path)

    if not file_path.exists():
        return {
            "status": "error",
            "message": f"File does not exist: {path}"
        }

    if not file_path.is_file():
        return {
            "status": "error",
            "message": f"Path is not a file: {path}"
        }

    try:
        content = file_path.read_text(encoding=encoding)
        return {
            "status": "success",
            "path": str(file_path.absolute()),
            "content": content,
            "bytes_read": len(content)
        }
    except UnicodeDecodeError:
        return {
            "status": "error",
            "message": f"Failed to decode file with encoding '{encoding}'"
        }

@mcp.tool()
def get_all_agent_data(directory: str) -> dict:
    """Get all file contents from a directory concatenated into a single string

    Args:
        directory: Directory path to read all files from

    Returns:
        Dictionary with status, directory path, and concatenated content
    """
    # First list all files in the directory
    ls_result = ls(directory)

    if ls_result["status"] != "success":
        return {
            "status": "error",
            "message": f"Failed to list directory: {ls_result.get('message', 'Unknown error')}"
        }

    # Filter only files (not directories)
    files = [item for item in ls_result["items"] if item["type"] == "file"]

    if not files:
        return {
            "status": "success",
            "directory": str(Path(directory).absolute()),
            "files_processed": 0,
            "content": ""
        }

    # Read and concatenate all files
    all_content = []
    files_processed = 0

    for file_item in files:
        file_path = file_item["path"]
        read_result = file_reader(file_path)

        if read_result["status"] == "success":
            all_content.append(f"=== {file_item['name']} ===\n")
            all_content.append(read_result["content"])
            all_content.append("\n\n")
            files_processed += 1
        else:
            # Continue with other files even if one fails
            all_content.append(f"=== {file_item['name']} ===\n")
            all_content.append(f"Error reading file: {read_result.get('message', 'Unknown error')}\n\n")

    return {
        "status": "success",
        "directory": str(Path(directory).absolute()),
        "files_processed": files_processed,
        "total_files": len(files),
        "content": "".join(all_content)
    }


@mcp.tool()
def api_query(url: str, query_params: Union[dict, str, None] = None, method: str = "GET", headers: Union[dict, None] = None, timeout: int = 30) -> dict:
    """Query an API endpoint and return JSON response

    Args:
        url: The API endpoint URL to query
        query_params: Query parameters as a dictionary (e.g., {"key": "value"}) or as a query string (e.g., "key=value&key2=value2"). If None, no query parameters are added.
        method: HTTP method to use (default: "GET"). Supported: GET, POST, PUT, DELETE, PATCH
        headers: Optional dictionary of HTTP headers to include in the request
        timeout: Request timeout in seconds (default: 30)

    Returns:
        Dictionary with status, response data (parsed JSON), status code, and response headers
    """
    if not url or not isinstance(url, str):
        return {
            "status": "error",
            "message": "URL must be a non-empty string",
            "response_data": None
        }

    try:
        # Parse query parameters
        params = None
        if query_params is not None:
            if isinstance(query_params, str):
                # Parse query string manually or let requests handle it
                params = {}
                for pair in query_params.split('&'):
                    if '=' in pair:
                        key, value = pair.split('=', 1)
                        params[key] = value
            elif isinstance(query_params, dict):
                params = query_params

        # Prepare request arguments
        request_kwargs = {
            "url": url,
            "params": params,
            "headers": headers,
            "timeout": timeout
        }

        # Make the request based on method
        method = method.upper()
        if method == "GET":
            response = requests.get(**request_kwargs)
        elif method == "POST":
            response = requests.post(**request_kwargs)
        elif method == "PUT":
            response = requests.put(**request_kwargs)
        elif method == "DELETE":
            response = requests.delete(**request_kwargs)
        elif method == "PATCH":
            response = requests.patch(**request_kwargs)
        else:
            return {
                "status": "error",
                "message": f"Unsupported HTTP method: {method}. Supported methods: GET, POST, PUT, DELETE, PATCH",
                "response_data": None
            }

        # Try to parse JSON response
        try:
            response_data = response.json()
        except json.JSONDecodeError:
            # If response is not JSON, return as text
            response_data = response.text

        return {
            "status": "success",
            "response_data": response_data,
            "status_code": response.status_code,
            "headers": dict(response.headers)
        }

    except requests.exceptions.Timeout:
        return {
            "status": "error",
            "message": f"Request timed out after {timeout} seconds",
            "response_data": None
        }
    except requests.exceptions.ConnectionError:
        return {
            "status": "error",
            "message": "Failed to connect to the API endpoint",
            "response_data": None
        }
    except requests.exceptions.RequestException as e:
        return {
            "status": "error",
            "message": f"Request failed: {str(e)}",
            "response_data": None
        }
    except Exception as e:
        return {
            "status": "error",
            "message": f"Unexpected error: {str(e)}",
            "response_data": None
        }


@mcp.tool()
async def code_executor(code: str, timeout: int = 8) -> dict:
    """Execute Python code with timeout and security restrictions

    Args:
        code: Python code to execute
        timeout: Max execution time in seconds (default 8)

    Returns:
        Dictionary with execution results and metadata
    """
    if not isinstance(code, str) or not code.strip():
        return {
            "status": "error",
            "message": "Code must be a non-empty string",
            "execution_result": None
        }

    try:
        execution_result = await _execute_python_code(code, timeout)
        return {
            "status": "success",
            "execution_result": execution_result
        }
    except Exception as e:
        return {
            "status": "error",
            "message": f"Execution failed: {str(e)}",
            "execution_result": None
        }


async def _execute_python_code(code: str, timeout: int = 8) -> dict:
    """Helper function to execute Python code with minimal security restrictions"""
    stdout = io.StringIO()
    stderr = io.StringIO()

    # Block only the most dangerous modules that could compromise the system
    # Most other modules are safe to use in this controlled environment
    unsafe_modules = {
        # Direct system execution and file system manipulation
        "os", "subprocess", "shutil",
        # Code execution risks (but allow normal imports)
        "exec", "eval",
        # Unpickling risks
        "pickle",
        # Direct socket access (though requests/http libraries are fine)
        "socket",
        # Low-level system access
        "ctypes", "multiprocessing"
    }

    # Only restrict the most dangerous imports, allow everything else
    real_import = builtins.__import__
    def safe_import(name, globals=None, locals=None, fromlist=(), level=0):
        base_name = name.split(".")[0]
        if base_name in unsafe_modules:
            raise ImportError(f"Import of '{name}' is blocked for security reasons")
        return real_import(name, globals, locals, fromlist, level)

    # Use standard builtins with just the import restriction
    safe_globals = {"__builtins__": {**builtins.__dict__, "__import__": safe_import}}
    safe_locals = {}

    # Async timeout wrapper
    async def _run():
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            try:
                exec(code, safe_globals, safe_locals)
                return True
            except Exception:
                traceback.print_exc(file=stderr)
                return False

    start = time.time()
    try:
        success = await asyncio.wait_for(_run(), timeout=timeout)
    except asyncio.TimeoutError:
        stderr.write(f"Execution timed out after {timeout} seconds\n")
        success = False
    except Exception:
        traceback.print_exc(file=stderr)
        success = False

    exec_time = round(time.time() - start, 4)

    return {
        "stdout": stdout.getvalue(),
        "stderr": stderr.getvalue(),
        "exec_time": exec_time,
        "success": success,
    }


if __name__ == "__main__":
    mcp.run(transport='stdio')
