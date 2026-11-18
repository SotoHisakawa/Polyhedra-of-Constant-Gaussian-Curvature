#!/usr/bin/env python3
"""Compatibility entry point for the paper's corrected nine-vertex torus.

The original script is preserved in archive/legacy_2026-11-18.zip.
Use generate_models.py for the full collection and CLI options.
"""
from generate_models import main


if __name__ == "__main__":
    main(["--model", "embedded_orientable_g1_v9_flat_torus"])
