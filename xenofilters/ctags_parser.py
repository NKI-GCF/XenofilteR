#!/usr/bin/env python3
"""
CTags Parser - Converts ctags output to a token-optimized hierarchical tree.

Usage:
    python3 ctags_parser.py <ctags_file> [options]
    python3 ctags_parser.py tags.txt
    python3 ctags_parser.py tags.txt --skip-tests
    python3 ctags_parser.py tags.txt --skip-test-functions --skip-test-modules
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple


class CTagsParser:
    """Parses ctags output into a hierarchical file and symbol tree."""

    def __init__(self, skip_test_functions: bool = False, skip_test_modules: bool = False):
        """
        Initialize the parser with configuration options.
        
        Args:
            skip_test_functions: If True, exclude symbols starting with 'test' 
                                 or containing test-related patterns in their name.
            skip_test_modules: If True, exclude files/directories related to testing
                              (e.g., 'tests.rs', 'mod tests').
        """
        self.skip_test_functions = skip_test_functions
        self.skip_test_modules = skip_test_modules
        
        # Compile regex patterns for test detection
        self.test_function_pattern = re.compile(r'^test', re.IGNORECASE)
        self.test_module_pattern = re.compile(r'(^|[/\\])tests?\.\w+$|mod\s+tests', re.IGNORECASE)
        
        # Data structure: nested dicts representing directory -> file -> symbols
        self.tree: Dict[str, any] = defaultdict(lambda: defaultdict(list))

    def is_test_function(self, tag_name: str, kind: str) -> bool:
        """
        Determine if a symbol is a test function based on naming conventions.
        
        Args:
            tag_name: The name of the symbol.
            kind: The kind of symbol (e.g., 'f' for function, 'm' for method).
            
        Returns:
            True if the symbol should be considered a test function.
        """
        return bool(self.test_function_pattern.match(tag_name))

    def is_test_module(self, file_path: str) -> bool:
        """
        Determine if a file path represents a test module.
        
        Args:
            file_path: The relative path to the file.
            
        Returns:
            True if the file path indicates a test module.
        """
        return bool(self.test_module_pattern.search(file_path))

    def parse_ctags_line(self, line: str) -> Optional[Tuple[str, str, int, str]]:
        """
        Parse a single ctags data line into its components.
        
        Args:
            line: A raw line from the ctags output file.
            
        Returns:
            A tuple of (tag_name, file_path, line_number, kind) or None if invalid.
        """
        # Remove trailing newline
        line = line.strip()
        
        # Skip metadata lines (those starting with '!_')
        if line.startswith('!') or not line:
            return None
        
        # Split by tabs - ctags format is tab-separated
        parts = line.split('\t')
        if len(parts) < 3:
            return None
        
        tag_name = parts[0]
        file_path = parts[1]
        
        # Extract line number from the third field (format: "LINE_NUMBER;")
        line_info = parts[2]
        line_match = re.match(r'(\d+);"', line_info)
        if not line_match:
            return None
        line_number = int(line_match.group(1))
        
        # Extract kind if available (fourth field)
        kind = parts[3] if len(parts) >= 4 else 'unknown'
        kind = kind.strip()
        
        return tag_name, file_path, line_number, kind

    def add_symbol(self, file_path: str, line_number: int, kind: str, tag_name: str):
        """
        Add a symbol to the tree structure.
        
        Args:
            file_path: Relative path to the file containing the symbol.
            line_number: Line number where the symbol is defined.
            kind: Type of symbol (e.g., 'f', 'c', 's', 'e').
            tag_name: Name of the symbol.
        """
        # Normalize path separators
        file_path = file_path.replace('\\', '/')
        
        # Split path into components
        path_parts = file_path.split('/')
        
        # Navigate/create directory structure
        current = self.tree
        for i, part in enumerate(path_parts):
            if i == len(path_parts) - 1:
                # This is the file
                if part not in current:
                    current[part] = []
                current[part].append((line_number, kind, tag_name))
            else:
                # This is a directory
                if part not in current:
                    current[part] = defaultdict(list)
                current = current[part]

    def parse_file(self, ctags_file_path: str):
        """
        Parse the entire ctags output file and build the tree.
        
        Args:
            ctags_file_path: Path to the ctags output file.
        """
        try:
            with open(ctags_file_path, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    parsed = self.parse_ctags_line(line)
                    if parsed is None:
                        continue
                    
                    tag_name, file_path, line_number, kind = parsed
                    
                    # Apply filters
                    if self.skip_test_modules and self.is_test_module(file_path):
                        continue
                    
                    if self.skip_test_functions and self.is_test_function(tag_name, kind):
                        continue
                    
                    self.add_symbol(file_path, line_number, kind, tag_name)
                    
        except FileNotFoundError:
            print(f"Error: File '{ctags_file_path}' not found.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Error parsing file: {e}", file=sys.stderr)
            sys.exit(1)

    def render_tree(self, node=None, indent_level: int = 0, is_file: bool = False) -> str:
        """
        Render the tree structure to the token-optimized format.
        
        Args:
            node: Current node in the tree (dict for directories, list for files).
            indent_level: Current indentation level.
            is_file: Whether the current node represents a file.
            
        Returns:
            String representation of the tree.
        """
        if node is None:
            node = self.tree
        
        output = []
        indent = "  " * indent_level
        
        if isinstance(node, dict):
            # Sort directories and files alphabetically
            items = sorted(node.items())
            
            for name, contents in items:
                if isinstance(contents, dict):
                    # Directory
                    output.append(f"{indent}{name}/")
                    output.append(self.render_tree(contents, indent_level + 1))
                else:
                    # File
                    output.append(f"{indent}{name}")
                    self._render_symbols(contents, indent_level + 1, output)
        
        return "\n".join(output)

    def _render_symbols(self, symbols: List[Tuple[int, str, str]], indent_level: int, output: List[str]):
        """
        Render symbols for a file, sorted by line number.
        
        Args:
            symbols: List of (line_number, kind, tag_name) tuples.
            indent_level: Current indentation level.
            output: List to append rendered lines to.
        """
        # Sort by line number, then by kind, then by name
        sorted_symbols = sorted(symbols, key=lambda x: (x[0], x[1], x[2]))
        
        indent = "  " * indent_level
        for line_num, kind, tag_name in sorted_symbols:
            output.append(f"{indent}{line_num}: {kind} {tag_name}")

    def print_tree(self):
        """Print the rendered tree to stdout."""
        print(self.render_tree())


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Parse ctags output into a token-optimized hierarchical tree.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s tags.txt
  %(prog)s tags.txt --skip-test-functions
  %(prog)s tags.txt --skip-test-modules
  %(prog)s tags.txt --skip-test-functions --skip-test-modules
        """
    )
    
    parser.add_argument(
        'ctags_file',
        help="Path to the ctags output file (generated with --excmd=nr -R src)"
    )
    
    parser.add_argument(
        '--skip-test-functions',
        action='store_true',
        default=False,
        help="Exclude symbols where the tag name starts with 'test' (case-insensitive)"
    )
    
    parser.add_argument(
        '--skip-test-modules',
        action='store_true',
        default=False,
        help="Exclude file paths containing 'tests.rs' or 'mod tests'"
    )
    
    parser.add_argument(
        '--skip-tests',
        action='store_true',
        default=False,
        help="Shortcut for enabling both --skip-test-functions and --skip-test-modules"
    )
    
    args = parser.parse_args()
    
    # Handle the --skip-tests shortcut
    skip_functions = args.skip_test_functions or args.skip_tests
    skip_modules = args.skip_test_modules or args.skip_tests
    
    # Create parser with specified configuration
    ctags_parser = CTagsParser(
        skip_test_functions=skip_functions,
        skip_test_modules=skip_modules
    )
    
    # Parse the ctags file
    ctags_parser.parse_file(args.ctags_file)
    
    # Print the rendered tree
    ctags_parser.print_tree()


if __name__ == "__main__":
    main()
