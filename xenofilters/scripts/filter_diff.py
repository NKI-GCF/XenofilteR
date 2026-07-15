#!/usr/bin/env python3
import sys
import re
import argparse
from typing import List

# Matches standard ANSI escape sequences outputted by git diff --color
ANSI_ESCAPE = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')

def strip_ansi(text: str) -> str:
    """Removes ANSI color codes from a string for safe parsing."""
    return ANSI_ESCAPE.sub('', text)

def is_rust_file(diff_git_line: str) -> bool:
    """Checks if the diff target is a Rust file (ends in .rs)."""
    clean_line = strip_ansi(diff_git_line).strip()
    # git diff format is: diff --git a/foo.rs b/foo.rs
    # Paths with spaces are quoted: diff --git "a/f.rs" "b/f.rs"
    return clean_line.endswith('.rs') or clean_line.endswith('.rs"')

def check_hunk_for_matches(hunk_lines: List[str], names: List[str], include_comments: bool) -> bool:
    """Determines if the current hunk contains any of the target names."""
    clean_lines = []

    for line in hunk_lines:
        clean_line = strip_ansi(line)
        if clean_line.startswith('@@ '):
            continue

        # Strip diff structural prefixes (+, -, space) to reconstruct clean code
        if len(clean_line) > 0 and clean_line[0] in ('+', '-', ' '):
            clean_lines.append(clean_line[1:])
        else:
            clean_lines.append(clean_line)

    full_text = "\n".join(clean_lines)

    if not include_comments:
        # Regex to match String literals OR Comments.
        # We match strings so we don't accidentally rip out `//` that sits safely inside a string.
        pattern = re.compile(
            r'("(?:\\.|[^"\\])*")'  # Group 1: String literals
            r'|'
            r'(//.*|/\*.*?\*/)',    # Group 2: Line (//) and Block (/* */) comments
            re.DOTALL
        )

        def replacer(match: re.Match) -> str:
            if match.group(1):
                return match.group(1) # Keep the string intact
            return ''                 # Remove the comment

        full_text = pattern.sub(replacer, full_text)

    for name in names:
        if name in full_text:
            return True

    return False

def main():
    parser = argparse.ArgumentParser(description="Filter git diff output by names in Rust files.")
    parser.add_argument('--matches', required=True, help="File containing a list of names to match, one per line.")
    parser.add_argument('--in-comments', action='store_true', help="If provided, also matches names found inside comments.")
    args = parser.parse_args()

    # Load names, ignoring empty lines
    with open(args.matches, 'r', encoding='utf-8') as f:
        names = [line.strip() for line in f if line.strip()]

    file_header: List[str] = []
    hunk_lines: List[str] = []
    is_rs_file = False
    file_printed = False

    def flush_hunk():
        """Evaluates and conditionally prints the active hunk."""
        nonlocal file_printed
        if hunk_lines and is_rs_file:
            if check_hunk_for_matches(hunk_lines, names, args.in_comments):
                if not file_printed:
                    for hl in file_header:
                        sys.stdout.write(hl)
                    file_printed = True
                for hl in hunk_lines:
                    sys.stdout.write(hl)

    # Process stdin dynamically as a stream
    for line in sys.stdin:
        clean_line = strip_ansi(line)

        # 1. Start of a new file definition
        if clean_line.startswith('diff --git '):
            flush_hunk()
            hunk_lines = []
            file_header = [line]
            is_rs_file = is_rust_file(clean_line)
            file_printed = False
            continue

        # 2. Skip Checksum index
        if clean_line.startswith('index '):
            continue

        # 3. Handle structural file header bounds
        if clean_line.startswith('--- ') or clean_line.startswith('+++ '):
            if not hunk_lines:
                file_header.append(line)
            continue

        # 4. Start of a new Hunk
        if clean_line.startswith('@@ '):
            flush_hunk()
            hunk_lines = [line]
            continue

        # 5. Accumulate payload
        if not hunk_lines:
            file_header.append(line)  # Extended file header info (e.g. file modes)
        else:
            hunk_lines.append(line)

    # Process the final hunk trailing at EOF
    flush_hunk()

if __name__ == '__main__':
    # Force stdin/stdout to handle utf-8 safely, mitigating standard Pipe errors
    sys.stdin.reconfigure(encoding='utf-8')
    sys.stdout.reconfigure(encoding='utf-8')
    main()
