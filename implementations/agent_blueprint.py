from pydantic import BaseModel
from langchain_openai import ChatOpenAI
from langgraph.graph import MessagesState
from langchain_core.messages import SystemMessage, BaseMessage
from typing import Any, Union, Sequence, Callable, Optional
from langchain_core.tools import BaseTool

class AgentBlueprint(BaseModel):
    """Configuration blueprint for creating a LangGraph agent.
    
    This class defines the essential components needed to construct an agent,
    including its language model, tools, prompt, and response format.
    
    Attributes:
        name: The identifier for the agent. Used for tracking and logging purposes.
        model: The language model to use for the agent. Can be a model identifier
            string (e.g., "gpt-4", "claude-3-sonnet") or a model instance identifier.
        tools: Optional list of tool names that the agent can use.
            If provided, the agent will have access to these tools for function calling.
            Defaults to None for agents without tool usage.
        agents: Optional list of agent names that this agent can delegate to.
            Enables multi-agent collaboration where this agent can invoke other agents.
            Defaults to None for agents that don't delegate to other agents.
        prompt: The system prompt or instruction for the agent. Defines the agent's
            behavior, role, and task-specific guidance. This will be used to configure
            the agent's base instructions.
        response_format: Optional dictionary schema for structuring the agent's JSON output.
            The model must support `.with_structured_output` to use response_format.
        unit_test: whether to run unit tests for the agent
    """    
    name: str
    model: str 
    tools: Optional[list[str]] = None
    agents: Optional[list[str]] = None
    prompt: str 
    response_format: Optional[dict] = None
    unit_test: bool = False

