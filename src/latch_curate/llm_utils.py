import inspect
from typing import Callable
import json

import openai
from openai.types.chat import ChatCompletionMessageToolCall
import tiktoken

from latch_curate.config import user_config

tiktoken_enc = tiktoken.get_encoding("cl100k_base")

def get_token_count(text: str) -> int:
    return len(tiktoken_enc.encode(text))

def function_to_schema(func_tool: Callable) -> dict:
    assert func_tool.__name__
    assert func_tool.__doc__

    python_to_openai_type_map = {
        str: "string",
        int: "integer",
        float: "number",
        bool: "boolean",
        list: "array",
        dict: "object",
        type(None): "null",
    }

    signature = inspect.signature(func_tool)

    parameters = {}
    for param in signature.parameters.values():
        if param.name == "self":
            continue
        param_type = python_to_openai_type_map[param.annotation]
        parameters[param.name] = {"type": param_type}

    return {
        "type": "function",
        "function": {
            "name": func_tool.__name__,
            "strict": True,
            "description": func_tool.__doc__.strip(),
            "parameters": {
                "properties": parameters,
                "required": list(parameters.keys()),
                "type": "object",
                "additionalProperties": False
            },
        },
    }

def execute_tool_calls(func_tools: list[Callable], tool_calls:
                       ChatCompletionMessageToolCall) -> (any, str, list[any]):
    print(f">>> parsing tool call response: {tool_calls}" )

    discovered_tool = False
    # todo(kenny): handle multiple
    tool_call = tool_calls[0]
    for tool in func_tools:
        if tool.__name__ == tool_call.function.name:
            discovered_tool = True
        fn_name = tool.__name__
        args = json.loads(tool_call.function.arguments)
        params = inspect.signature(tool).parameters
        for param in params.values():
            if param.name == "self":
                continue
            assert param.name in args, f"args: {args}; param_name: {param.name}"
        if not discovered_tool:
            raise ValueError(f"No valid tools available for {tool_calls}")
        try:
            print(f">>> executing tool call: {fn_name}")
            print(f">>> with args: {args}")
            return tool(**args), fn_name, args
        except Exception as e:
            print(f"unable to call tool {tool}")
            print(f">> args {args}")
            raise e.with_traceback()

    if not discovered_tool:
        raise ValueError("Unable to discover func tool in tool call result {tool_call}")

def prompt_model(messages: list[dict[str, any]], model: str = "o4-mini", tools: list[Callable] = []) -> [str, list[ChatCompletionMessageToolCall] | None]:
    openai.api_key = user_config.openai_api_key

    response = openai.chat.completions.create(
        model=model,
        messages=messages,
        tools=[function_to_schema(tool) for tool in tools]
    )

    message = response.choices[0].message
    message_content = message.content.strip() if message.content else None
    tool_calls = message.tool_calls

    return message_content, tool_calls
