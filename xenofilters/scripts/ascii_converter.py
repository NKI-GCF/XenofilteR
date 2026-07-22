import os
import sys
from pathlib import Path

# Mapping of non-ASCII characters to their ASCII counterparts
CHAR_MAP = {
    '§': 'S.',
    '±': '+/-',
    '²': '^2',
    '·': '*',
    '×': '*',
    'ä': 'a',
    'ç': 'c',
    'é': 'e',
    'î': 'i',
    'ñ': 'n',
    'ø': 'o',
    'ü': 'u',
    'Σ': 'Sigma',
    '–': '-',   # en-dash
    '—': '--',  # em-dash
    '’': "'",   # right single quote
    '…': '...',
    '′': "'",   # prime
    '←': '<-',
    '→': '->',
    '↔': '<->',
    '∈': ' in ',
    '−': '-',   # minus sign
    '∧': ' and ', 
    '≈': '~=',
    '≠': '!=',
    '≡': '==',
    '≤': '<=',
    '≥': '>=',
    '⊘': '/',
    '│': '|',
    '┐': '+',
    '└': '+',
    '┘': '+',
    '├': '+',
    '┤': '+',
    '┴': '+',
    '═': '=',
    '█': '#',
    '░': '.',
    '▶': '>',
    '►': '>',
    '▼': 'v',
    '◄': '<',
    '○': 'o',
    '◑': 'o',
    '✓': 'ok',  # or 'v' if you prefer a visual tick
    '✗': 'x'
}

# Create a highly efficient translation table
TRANSLATION_TABLE = str.maketrans(CHAR_MAP)

def process_file(filepath: Path) -> None:
    """Reads a file, replaces non-ASCII characters, and writes it back."""
    try:
        # Read the file content
        content = filepath.read_text(encoding='utf-8')
    except UnicodeDecodeError:
        # Skip binary files or files with unknown encodings
        return
    
    # Translate the content using our mapping
    new_content = content.translate(TRANSLATION_TABLE)
    
    # Only write back to disk if changes were actually made
    if content != new_content:
        filepath.write_text(new_content, encoding='utf-8')
        print(f"Updated: {filepath}")

def main(directory: str = "."):
    """Walks through the directory and processes files."""
    print(f"Scanning directory: {os.path.abspath(directory)}...")
    
    target_dir = Path(directory)
    
    # Feel free to add or remove directories you want to ignore (e.g., .git, venv)
    ignore_dirs = {'.git', '__pycache__', 'venv', '.env', 'node_modules'}
    
    for root, dirs, files in os.walk(target_dir):
        # Modify dirs in-place to skip ignored directories
        dirs[:] = [d for d in dirs if d not in ignore_dirs]
        
        for file in files:
            filepath = Path(root) / file
            
            # Optional: skip certain file extensions if needed
            # if filepath.suffix in {'.jpg', '.png', '.pyc', '.pdf'}:
            #     continue
                
            process_file(filepath)
            
    print("Replacement complete.")

if __name__ == "__main__":
    # You can pass a directory path as a command-line argument, defaults to current dir
    target_path = sys.argv[1] if len(sys.argv) > 1 else "."
    main(target_path)
