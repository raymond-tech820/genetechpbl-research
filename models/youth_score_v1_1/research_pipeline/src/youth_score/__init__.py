"""Donor-aware Youth Score modeling package.

PyTorch must be loaded before pandas/pyarrow on this Windows CUDA environment.
Both stacks ship native runtime DLLs, and the reverse order can make c10.dll fail
to initialize. Importing it here makes every package entry point deterministic.
"""

import os as _os

_os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import torch as _torch  # noqa: F401

__version__ = "0.1.0"
