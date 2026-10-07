#!/usr/bin/env python3
"""Launcher: run this file to start CoSegVCF without installing the package.

    python run_cosegvcf.py            # opens the graphical interface
    python run_cosegvcf.py --help     # command-line options

Keep it in the repository root, next to the `cosegvcf` folder.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cosegvcf.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
