#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

# ANSI color codes for better readability in the terminal
# ANSI color codes, disabled when piped to a file or another command
USE_COLORS = sys.stdout.isatty()
CYAN = '\033[96m' if USE_COLORS else ''
RED = '\033[91m' if USE_COLORS else ''
GREEN = '\033[92m' if USE_COLORS else ''
YELLOW = '\033[93m' if USE_COLORS else ''
RESET = '\033[0m' if USE_COLORS else ''

def strip_rust_comments(content: str) -> str:
    """Strips // and /* */ comments (even nested) while preserving newlines."""
    out = []
    i, n = 0, len(content)
    in_string, string_escape = False, False
    comment_depth, in_line_comment = 0, False

    while i < n:
        c = content[i]

        if in_line_comment:
            if c == '\n':
                in_line_comment = False
                out.append(c)
            i += 1
            continue

        if comment_depth > 0:
            if c == '/' and i + 1 < n and content[i+1] == '*':
                comment_depth += 1
                i += 2
            elif c == '*' and i + 1 < n and content[i+1] == '/':
                comment_depth -= 1
                i += 2
            else:
                if c == '\n':
                    out.append(c)
                i += 1
            continue

        if in_string:
            out.append(c)
            if string_escape:
                string_escape = False
            elif c == '\\':
                string_escape = True
            elif c == '"':
                in_string = False
            i += 1
            continue

        if c == '"':
            in_string = True
            out.append(c)
            i += 1
        elif c == '/' and i + 1 < n and content[i+1] == '/':
            in_line_comment = True
            i += 2
        elif c == '/' and i + 1 < n and content[i+1] == '*':
            comment_depth += 1
            i += 2
        else:
            out.append(c)
            i += 1

    return "".join(out)

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

def process_file(filepath: Path, context_lines: int, keep_comments: bool, keep_space: bool):
    """Reads a file and prints its conflict blocks in an AI-optimized format."""
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return

    orig_lines = content.split('\n')
    if orig_lines and orig_lines[-1] == '':
        orig_lines.pop()

    if not keep_comments:
        content = strip_rust_comments(content)

    processed_lines = content.split('\n')

    conflict_blocks = []
    in_conflict = False
    start_idx = -1

    # 1. Identify all conflict blocks using original lines
    for i, line in enumerate(orig_lines):
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
        c_end = min(len(orig_lines) - 1, end + context_lines)

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
    out_filepath = str(filepath).replace('.orig', '')
    print(f"{CYAN}# {out_filepath}{RESET}")

    for idx, (w_start, w_end) in enumerate(merged_windows):
        if idx > 0:
            print(f"{YELLOW}--{RESET}")

        for i in range(w_start, w_end + 1):
            orig_line = orig_lines[i]
            proc_line = processed_lines[i] if i < len(processed_lines) else ""

            if orig_line.startswith('<<<<<<<'):
                branch = orig_line[7:].strip()
                print(f"{RED}< {branch} +{i + 1}{RESET}")
            elif orig_line.startswith('|||||||'):
                print(f"{YELLOW}|{RESET}")
            elif orig_line.startswith('======='):
                print(f"{CYAN}={RESET}")
            elif orig_line.startswith('>>>>>>>'):
                print(f"{GREEN}>{RESET}")
            else:
                if not keep_space and not proc_line.strip():
                    continue
                print(proc_line)

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
        '--comments',
        action='store_true',
        help="Do not strip comments (stripped by default)"
    )
    parser.add_argument(
        '--space',
        action='store_true',
        help="Do not remove empty lines (removed by default)"
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
        process_file(filepath, args.U, args.comments, args.space)

if __name__ == '__main__':
    main()
