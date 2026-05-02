#!/usr/bin/env python3
"""
Script to invoke the agent_manager agent with mutation data.
Supports multiprocessing for parallel mutation processing.
"""

import asyncio
import json
import os
import shutil
import subprocess
from pathlib import Path
from dotenv import load_dotenv
from langchain_core.messages import HumanMessage
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

PROJECT_ROOT = Path(__file__).resolve().parent
ENV_PATH = str(PROJECT_ROOT / '.env')


def process_mutation_worker(args):
    """
    Worker function for multiprocessing. Must be at module level for pickling.
    Recreates the graph in each worker process.
    """
    mutation_data, mutation_index, total_mutations = args

    # Load environment variables in each worker process
    load_dotenv(ENV_PATH)

    # Run the async mutation processing
    return asyncio.run(process_single_mutation_async(
        mutation_data, mutation_index, total_mutations
    ))


async def process_single_mutation_async(input_data, mutation_index, total_mutations):
    """Async function to process a single mutation. Used by worker processes."""

    # Import and recreate graph in this process
    from mini_graphs.agent_manager import build_agent_manager_graph

    print("\n" + "=" * 80)
    print(f"[Process {os.getpid()}] Processing mutation {mutation_index + 1}/{total_mutations}")
    print("=" * 80)
    print(f"Gene: {input_data['gene_name']}, Mutation: {input_data['mutation_aa']}")
    print(f"Full Input Data: {input_data}")
    print("-" * 80)

    # Recreate graph in this process
    graph = await build_agent_manager_graph()

    # Create a message with all required inputs for agent_manager
    fields = [
        ('gene_name', 'gene_name'),
        ('accession_number', 'accession_number'),
        ('mutation_cds', 'mutation_cds'),
        ('mutation_aa', 'mutation_aa'),
        ('mutation_description_aa', 'mutation_description_aa'),
        ('aa_mut_start', 'aa_mut_start'),
        ('aa_mut_stop', 'aa_mut_stop'),
        ('mutation_description_cds', 'mutation_description_cds')
    ]

    content_lines = ["Analyze mutation with the following inputs:"]
    for field_name, key in fields:
        if key in input_data:
            content_lines.append(f"- {field_name}: {input_data[key]}")

    human_message = HumanMessage(content="\n".join(content_lines))

    # Create the initial state with the message
    initial_state = {"messages": [human_message]}

    # Configure recursion limit to prevent GraphRecursionError
    config = {"recursion_limit": 100}

    try:
        print(f"[Process {os.getpid()}] Invoking agent_manager graph...")
        # Invoke the graph
        result = await graph.ainvoke(initial_state, config=config)

        print(f"\n[Process {os.getpid()}] Mutation {mutation_index + 1} completed successfully!")

        # Print all messages in the conversation
        for j, message in enumerate(result["messages"]):
            print(f"\n[Process {os.getpid()}] Message {j+1}:")
            print(f"Type: {type(message).__name__}")
            print(f"Content: {message.content}")
            if hasattr(message, 'tool_calls') and message.tool_calls:
                print(f"Tool calls: {message.tool_calls}")
            print("-" * 80)

        return True

    except Exception as e:
        print(f"\n[Process {os.getpid()}] Error processing mutation {mutation_index + 1}: {e}")
        import traceback
        traceback.print_exc()
        return False


async def invoke_mutation(input_data, mutation_index, total_mutations):
    """Invoke the agent_manager graph with a single mutation data entry (sequential mode)."""

    # Load environment variables from .env file
    load_dotenv(ENV_PATH)
    from mini_graphs.agent_manager import graph, build_agent_manager_graph

    # Build graph if it wasn't created at module level
    active_graph = graph
    if active_graph is None:
        active_graph = await build_agent_manager_graph()

    print("\n" + "=" * 80)
    print(f"Processing mutation {mutation_index + 1}/{total_mutations}")
    print("=" * 80)
    print(f"Gene: {input_data['gene_name']}, Mutation: {input_data['mutation_aa']}")
    print(f"Full Input Data: {input_data}")
    print("-" * 80)

    # Create a message with all required inputs for agent_manager
    fields = [
        ('gene_name', 'gene_name'),
        ('accession_number', 'accession_number'),
        ('mutation_cds', 'mutation_cds'),
        ('mutation_aa', 'mutation_aa'),
        ('mutation_description_aa', 'mutation_description_aa'),
        ('aa_mut_start', 'aa_mut_start'),
        ('aa_mut_stop', 'aa_mut_stop'),
        ('mutation_description_cds', 'mutation_description_cds')
    ]

    content_lines = ["Analyze mutation with the following inputs:"]
    for field_name, key in fields:
        if key in input_data:
            content_lines.append(f"- {field_name}: {input_data[key]}")

    human_message = HumanMessage(content="\n".join(content_lines))

    # Create the initial state with the message
    initial_state = {"messages": [human_message]}

    # Configure recursion limit to prevent GraphRecursionError
    config = {"recursion_limit": 100}

    try:
        print("Invoking agent_manager graph...")
        result = await active_graph.ainvoke(initial_state, config=config)

        print(f"\nMutation {mutation_index + 1} completed successfully!")

        for j, message in enumerate(result["messages"]):
            print(f"\nMessage {j+1}:")
            print(f"Type: {type(message).__name__}")
            print(f"Content: {message.content}")
            if hasattr(message, 'tool_calls') and message.tool_calls:
                print(f"Tool calls: {message.tool_calls}")
            print("-" * 80)

        return True

    except Exception as e:
        print(f"\nError processing mutation {mutation_index + 1}: {e}")
        import traceback
        traceback.print_exc()
        return False


WEBSITE_DIR = Path("/home/phosphoscout/website")
WEBSITE_REPORTS_DIR = WEBSITE_DIR / "docs" / "reports"
PHOSPHOSCOUT_PYTHON = "/home/phosphoscout/miniconda3/envs/phosphoscout/bin/python"


def publish_reports_to_website(project_root: Path) -> None:
    """Copy generated HTML reports to the website, rebuild, commit, and push."""
    html_reports_dir = project_root / "generated_artifacts" / "reports" / "html_reports"
    if not html_reports_dir.is_dir():
        print("No html_reports directory found — skipping website publish.")
        return

    WEBSITE_REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    copied = []
    for html_file in html_reports_dir.glob("*.html"):
        dest = WEBSITE_REPORTS_DIR / html_file.name
        shutil.copy2(html_file, dest)
        copied.append(html_file.name)

    if not copied:
        print("No HTML reports to publish.")
        return

    print(f"Publishing {len(copied)} report(s) to website: {copied}")

    subprocess.run(
        [PHOSPHOSCOUT_PYTHON, str(WEBSITE_DIR / "build.py")],
        cwd=str(WEBSITE_DIR),
        check=True,
        timeout=120,
    )

    subprocess.run(["git", "add", "-A"], cwd=str(WEBSITE_DIR), check=True, timeout=30)
    result = subprocess.run(
        ["git", "commit", "-m", f"Add report(s): {', '.join(copied)}"],
        cwd=str(WEBSITE_DIR),
        timeout=30,
    )
    if result.returncode == 0:
        subprocess.run(["git", "push"], cwd=str(WEBSITE_DIR), check=True, timeout=60)
        print(f"Website updated and pushed.")
    else:
        print("Nothing new to commit to website (already published).")


def invoke_all_mutations_parallel(max_workers=None):
    """Load mutations from file and process them in parallel using multiprocessing."""

    mutations_file = PROJECT_ROOT / "data" / "mutations_to_run.txt"

    print(f"Loading mutations from: {mutations_file}")
    with open(mutations_file, 'r') as f:
        mutations = json.load(f)

    total_mutations = len(mutations)
    print(f"Found {total_mutations} mutations to process")

    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    print(f"Using {max_workers} worker processes for parallel processing")

    mutation_args = [
        (mutation_data, i, total_mutations)
        for i, mutation_data in enumerate(mutations)
    ]

    successful = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_index = {
            executor.submit(process_mutation_worker, args): i
            for i, args in enumerate(mutation_args)
        }

        for future in as_completed(future_to_index):
            mutation_index = future_to_index[future]
            try:
                success = future.result()
                if success:
                    successful += 1
                else:
                    failed += 1
            except Exception as e:
                print(f"\nUnexpected error processing mutation {mutation_index + 1}: {e}")
                import traceback
                traceback.print_exc()
                failed += 1

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total mutations processed: {total_mutations}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print("=" * 80)

    if successful > 0:
        publish_reports_to_website(PROJECT_ROOT)


async def invoke_all_mutations():
    """Load mutations from file and process each one sequentially."""

    mutations_file = PROJECT_ROOT / "data" / "mutations_to_run.txt"

    print(f"Loading mutations from: {mutations_file}")
    with open(mutations_file, 'r') as f:
        mutations = json.load(f)

    total_mutations = len(mutations)
    print(f"Found {total_mutations} mutations to process")

    successful = 0
    failed = 0

    for i, mutation_data in enumerate(mutations):
        success = await invoke_mutation(mutation_data, i, total_mutations)
        if success:
            successful += 1
        else:
            failed += 1

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total mutations processed: {total_mutations}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print("=" * 80)

    if successful > 0:
        publish_reports_to_website(PROJECT_ROOT)


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--parallel":
        max_workers = None
        if len(sys.argv) > 2:
            try:
                max_workers = int(sys.argv[2])
            except ValueError:
                print(f"Warning: Invalid max_workers value '{sys.argv[2]}', using default")
        invoke_all_mutations_parallel(max_workers=max_workers)
    else:
        asyncio.run(invoke_all_mutations())
