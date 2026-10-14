"""KROME-X: the code generator and its configuration.

    from kromex import Preprocessor

Used from the repository, without installation (docs/manual.md). The
modules of this folder import each other by their plain names (codegen,
knet, constants, ...), as when run as scripts, so the folder itself is put
on the module search path."""

import os as _os
import sys as _sys

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

from preprocessor import Preprocessor, ConfigError, GenerationError  # noqa: E402,F401

__all__ = ["Preprocessor", "ConfigError", "GenerationError"]
