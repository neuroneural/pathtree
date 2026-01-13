"""
Pytest configuration for pathtree tests.
"""

import sys
from pathlib import Path

# Add pathtree package to path for imports
project_root = Path(__file__).parent.parent
pathtree_dir = project_root / "pathtree"
tools_dir = project_root / "tools"

sys.path.insert(0, str(pathtree_dir))
sys.path.insert(0, str(tools_dir))
sys.path.insert(0, str(project_root))
