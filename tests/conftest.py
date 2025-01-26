import sys
from pathlib import Path

# Get the absolute path of the 'src' directory
src_path = Path(__file__).resolve().parent.parent / "src"

print(src_path)

# Add the 'src' directory to the Python path
sys.path.insert(0, str(src_path))
