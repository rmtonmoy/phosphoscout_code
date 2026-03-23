"""
Individual graph for agent_manager agent.
Simple pattern: START -> agent_manager -> END
"""

from langgraph.graph import StateGraph, START, END, MessagesState
from implementations.helpers import agent_from_name
import asyncio

async def build_agent_manager_graph():
    """Build the graph for agent_manager."""
    builder = StateGraph(MessagesState)

    # Get the agent instance
    agent_instance = await agent_from_name("agent_manager")

    # Add the agent as a node
    builder.add_node("agent_manager", agent_instance)

    # Add simple edges: START -> agent -> END
    builder.add_edge(START, "agent_manager")
    builder.add_edge("agent_manager", END)

    return builder.compile()

# Create the graph lazily - only if there's no event loop running
try:
    loop = asyncio.get_running_loop()
    graph = None
except RuntimeError:
    graph = asyncio.run(build_agent_manager_graph())
