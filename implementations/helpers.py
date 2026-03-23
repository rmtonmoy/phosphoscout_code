from implementations.agent_blueprint import AgentBlueprint
import yaml
import aiofiles
from langgraph.graph import MessagesState
from langchain_openai import ChatOpenAI
from langchain.agents import create_agent
from constants import MCP_CLIENT, CONFIGS_PATH
from implementations.generate_docstring import generate_docstring
import asyncio
import os
_tools_cache = None

from typing import Callable, Awaitable, Any, Dict
from langgraph.graph import StateGraph, START, END
from langchain_core.messages import HumanMessage
from langchain.agents.middleware import TodoListMiddleware
from langchain.chat_models import init_chat_model



def make_agent_runner(
    func_name: str,
    docstring: str,
    agent_node_name: str,
    agent_from_name: Callable[[str], Awaitable[Any]],
):
    """
    Build an async function with a dynamic name + docstring that:
      - builds a StateGraph
      - inserts the agent node
      - runs graph.ainvoke(query)
    """
    async def _runner(query: str) -> Dict[str, Any]:
        builder = StateGraph(MessagesState)
        agent_instance = await agent_from_name(agent_node_name)

        builder.add_node(agent_node_name, agent_instance)

        builder.add_edge(START, agent_node_name)
        builder.add_edge(agent_node_name, END)

        graph = builder.compile()

        human_message = HumanMessage(content=query)
        initial_state = {"messages": [human_message]}
        config = {"recursion_limit": 100}

        result = await graph.ainvoke(initial_state, config=config)
        return result

    _runner.__name__ = func_name
    _runner.__doc__ = docstring
    return _runner


async def toolify_an_agent(agent_name: str):
    assert agent_name != "agent_creator", "Critical Requirement: agent_creator cannot be toolified"
    bp = await get_agent_bp(agent_name)
    yaml_path = f"{CONFIGS_PATH}/agents/{agent_name}.yaml"
    with open(yaml_path, "r") as f:
        blueprint_text = f.read()
    docstring = await generate_docstring(blueprint_text, agent_name)
    
    runner = make_agent_runner(
        func_name=f"{agent_name}_as_tool",
        docstring=docstring,
        agent_node_name=agent_name,
        agent_from_name=agent_from_name,
    )
    return runner

def get_all_agents():
    all_agents = []
    for file in os.listdir(CONFIGS_PATH + "/agents"):
        if file.endswith(".yaml"):
            all_agents.append(file.replace(".yaml", ""))
    return all_agents

async def get_all_tools():
    global _tools_cache
    if _tools_cache is None:
        all_tools = []
        for server_name in MCP_CLIENT.connections:
            try:
                server_tools = await MCP_CLIENT.get_tools(server_name=server_name)
                all_tools.extend(server_tools)
                print(f"Loaded {len(server_tools)} tools from {server_name}")
            except Exception as e:
                print(f"Warning: Failed to load tools from {server_name}: {e}")
        _tools_cache = all_tools
        if not _tools_cache:
            print("WARNING: No MCP tools were loaded from any server!")
    return _tools_cache

async def get_tools_from_names(tools_names: list[str]):    
    all_tools = await get_all_tools()
    if tools_names == ["all_tools"]:
        tools_copy = list(all_tools)
        all_agents = get_all_agents()
        for agent in all_agents:
            if agent != "agent_creator":
                tool = await toolify_an_agent(agent)
                tools_copy.append(tool)
        return tools_copy
    
    allowed_tools = [tool for tool in all_tools if tool.name in tools_names and not ("_as_tool" in tool.name)]
    for tool_name in tools_names:
        if ("_as_tool" in tool_name):
            allowed_tools.append(await toolify_an_agent(tool_name.replace("_as_tool", "")))

    # Warn about missing tools
    found_names = {t.name if hasattr(t, 'name') else t.__name__ for t in allowed_tools}
    missing = [n for n in tools_names if n not in found_names and f"{n}_as_tool" not in found_names]
    if missing:
        print(f"WARNING: Expected tools not found: {missing}")

    return allowed_tools

async def get_agent_bp(agent_name: str)->AgentBlueprint:
    path = f"{CONFIGS_PATH}/agents/{agent_name}.yaml"
    async with aiofiles.open(path, "r") as f:
        content = await f.read()
        config = yaml.safe_load(content)
    return AgentBlueprint(**config)

async def agent_from_name(agent_name: str):    
    bp = await get_agent_bp(agent_name)
    if bp.name != "agent_creator":
        assert "all_tools" not in bp.tools, "Critical Requirement: Only agent_creator can use all_tools"
    if not bp.agents:
        _model = init_chat_model(
            bp.model,
            temperature=1,
            parallel_tool_calls=False
        )
        agent_kwargs = {
            "model": _model,
            "tools": await get_tools_from_names(bp.tools),
            "system_prompt": bp.prompt,
            "middleware": [TodoListMiddleware()],
        }
        if bp.response_format:
            agent_kwargs["response_format"] = bp.response_format
        ret = create_agent(**agent_kwargs)
        return ret
    else:
        assert False, "Agents with agents are not supported yet"
        pass

async def main():
    print(await MCP_CLIENT.get_tools())

if __name__ == "__main__":
    asyncio.run(main())