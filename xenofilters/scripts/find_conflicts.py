#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

# ANSI color codes for better readability in the terminal
CYAN = '\033[96m'
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RESET = '\033[0m'

def highlight_line(line: str) -> str:
    """Applies color to git conflict markers."""
    if line.startswith('<<<<<<<'):
        return f"{RED}{line}{RESET}"
    elif line.startswith('|||||||'):
        return f"{YELLOW}{line}{RESET}"
    elif line.startswith('======='):
        return f"{CYAN}{line}{RESET}"
    elif line.startswith('>>>>>>>'):
        return f"{GREEN}{line}{RESET}"
    return line

def process_file(filepath: Path, context_lines: int):
    """Reads a file and prints its conflict blocks with optional context."""
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return

    conflict_blocks = []
    in_conflict = False
    start_idx = -1

    # 1. Identify all conflict blocks
    for i, line in enumerate(lines):
        if line.startswith('<<<<<<< '):
            in_conflict = True
            start_idx = i
        elif line.startswith('>>>>>>> '):
            if in_conflict:
                conflict_blocks.append((start_idx, i))
                in_conflict = False

    if not conflict_blocks:
        return

    # 2. Calculate and merge overlapping context windows
    merged_windows = []
    for start, end in conflict_blocks:
        c_start = max(0, start - context_lines)
        c_end = min(len(lines) - 1, end + context_lines)

        if not merged_windows:
            merged_windows.append([c_start, c_end])
        else:
            last_window = merged_windows[-1]
            # If the new window overlaps or touches the previous one, merge them
            if c_start <= last_window[1] + 1:
                last_window[1] = max(last_window[1], c_end)
            else:
                merged_windows.append([c_start, c_end])

    # 3. Print the results
    print(f"\n{CYAN}{'='*60}{RESET}")
    print(f"{CYAN}File: {filepath}{RESET}")
    print(f"{CYAN}{'='*60}{RESET}")

    for idx, (w_start, w_end) in enumerate(merged_windows):
        if idx > 0:
            print(f"{YELLOW}--- (snip) ---{RESET}")
            
        for i in range(w_start, w_end + 1):
            line_content = highlight_line(lines[i].rstrip('\n'))
            print(f"{i + 1:4d} | {line_content}")

def main():
    parser = argparse.ArgumentParser(
        description="Find .rs.orig files and display their merge conflicts."
    )
    # Using -U to mimic git's unified context flag
    parser.add_argument(
        '-U', 
        type=int, 
        default=0, 
        metavar='<nr>',
        help="Number of lines of context to show around conflicts (e.g. -U3)."
    )
    parser.add_argument(
        'path', 
        nargs='?', 
        default='.', 
        type=Path,
        help="Directory to search (default: current directory)."
    )

    args = parser.parse_args()

    # Find all .rs.orig files recursively
    orig_files = list(args.path.rglob('*.rs.orig'))
    
    if not orig_files:
        print(f"No .rs.orig files found in {args.path.resolve()}")
        return

    for filepath in orig_files:
        process_file(filepath, args.U)

if __name__ == '__main__':
    main()
