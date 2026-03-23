import os
from pathlib import Path
import dotenv

PROJECT_ROOT = Path(__file__).resolve().parents[1]
dotenv.load_dotenv(str(PROJECT_ROOT / '.env'))

import json
import yaml
import hashlib
import asyncio
from openai import OpenAI

def read_prompt_from_file(filepath: str) -> str:
    """Read the AGENT_TOOLIFYING_PROMPT variable from the text file."""
    with open(filepath, 'r') as f:
        content = f.read()

    start_marker = 'AGENT_TOOLIFYING_PROMPT = """'
    if start_marker in content:
        start_idx = content.index(start_marker) + len(start_marker)
        end_idx = content.rindex('"""')
        prompt = content[start_idx:end_idx].strip()
    else:
        raise ValueError(f"Could not find AGENT_TOOLIFYING_PROMPT in {filepath}")

    return prompt

def read_yaml_file(filepath: str) -> str:
    """Read YAML file and return as string."""
    with open(filepath, 'r') as f:
        return f.read()

async def generate_docstring(blueprint_text: str, agent_name: str) -> str:
    """Generate docstring using OpenAI chat completion with caching."""
    blueprint_hash = hashlib.sha256(blueprint_text.encode('utf-8')).hexdigest()

    docstrings_dir = Path("configs/docstrings")
    docstrings_dir.mkdir(exist_ok=True)

    cache_file = docstrings_dir / f"{agent_name}_{blueprint_hash}.txt"

    if cache_file.exists():
        with open(cache_file, 'r', encoding='utf-8') as f:
            return f.read()

    prompt_file = os.path.join(os.path.dirname(__file__), "agent_toolifying_prompt.txt")
    prompt = read_prompt_from_file(prompt_file)

    client = OpenAI()

    response = client.chat.completions.create(
        model="gpt-5",
        messages=[
            {
                "role": "system",
                "content": prompt,
            },
            {
                "role": "user",
                "content": blueprint_text,
            },
        ],
        response_format={"type": "json_object"},
    )

    if not response.choices or not response.choices[0].message.content:
        raise RuntimeError("Model returned empty content for docstring generation.")

    raw_json = response.choices[0].message.content
    parsed = json.loads(raw_json)

    docstring = parsed.get("docstring", "")
    if not docstring:
        raise ValueError("Response JSON did not contain 'docstring' field")

    with open(cache_file, 'w', encoding='utf-8') as f:
        f.write(docstring)

    return docstring

async def main():
    yaml_file = "../configs/agents/mutation_data_collector.yaml"
    agent_name = Path(yaml_file).stem

    print("Reading YAML file...")
    blueprint_text = read_yaml_file(yaml_file)

    print("Generating docstring...")
    docstring = await generate_docstring(blueprint_text, agent_name)

    print("\n" + "="*80)
    print("GENERATED DOCSTRING:")
    print("="*80 + "\n")
    print(docstring)
    print("\n" + "="*80)

if __name__ == "__main__":
    asyncio.run(main())
