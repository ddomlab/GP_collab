"""Source package for the project.

This package contains Python 3.14 compatibility patches for dependencies.
"""

import typing

import torch_geometric.inspector

# Patch torch_geometric.inspector.type_repr for Python 3.14 compatibility
_original_type_repr = torch_geometric.inspector.type_repr


def _patched_type_repr(obj: typing.Any, _globals: dict[str, typing.Any]) -> str:
    if hasattr(obj, "__module__") and obj.__module__ == "typing":
        name = getattr(obj, "_name", None)
        if name is None:
            origin = getattr(obj, "__origin__", None)
            if origin is not None:
                origin_str = str(origin)
                if "Union" in origin_str:
                    name = "Union"
                else:
                    name = origin_str.split(".")[-1]
            else:
                name = str(obj).split(".")[-1]

        args = getattr(obj, "__args__", None)
        if args is None or len(args) == 0:
            return name if name in _globals else f"typing.{name}"

        if (
            name == "Union"
            and len(args) == 2
            and any(arg is type(None) for arg in args)
        ):
            name = "Optional"

        if name == "Optional":
            args = [arg for arg in args if arg is not type(None)]

        args_repr = ", ".join([_patched_type_repr(arg, _globals) for arg in args])
        fullname = name if name in _globals else f"typing.{name}"
        return f"{fullname}[{args_repr}]"

    return _original_type_repr(obj, _globals)


torch_geometric.inspector.type_repr = _patched_type_repr
