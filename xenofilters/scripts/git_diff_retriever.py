#!/usr/bin/env python3
import sys
import os
import subprocess
import argparse
import re

def get_git_root():
    try:
        res = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True,
            text=True,
            check=True
        )
        return res.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error: Not in a git repository or git command failed: {e.stderr.strip()}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: git command not found. Please ensure Git is installed and in your PATH.", file=sys.stderr)
        sys.exit(1)

def get_current_branch():
    res = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
        check=True
    )
    return res.stdout.strip()

def get_git_diff(main_branch, feature_branch, extra_args=None):
    if extra_args is None:
        extra_args = []
    try:
        cmd = ["git", "diff", "-p"] + extra_args + [f"{main_branch}...{feature_branch}"]
        res = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        return res.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running git diff: {e.stderr.strip()}", file=sys.stderr)
        sys.exit(1)

def parse_diff(diff_text):
    files = {}
    current_file = None
    file_header_lines = []
    hunk_re = re.compile(r'^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@(.*)')
    
    lines = diff_text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('diff --git '):
            file_header_lines = [line]
            current_file = None
            i += 1
            while i < len(lines) and not lines[i].startswith('@@ ') and not lines[i].startswith('diff --git '):
                file_header_lines.append(lines[i])
                i += 1
            
            for fhl in file_header_lines:
                path = None
                if fhl.startswith('+++ b/'):
                    path = fhl[6:]
                elif fhl.startswith('--- a/'):
                    path = fhl[6:]
                elif fhl.startswith('+++ '):
                    path = fhl[4:]
                    if path.startswith('b/'): path = path[2:]
                elif fhl.startswith('--- '):
                    path = fhl[4:]
                    if path.startswith('a/'): path = path[2:]
                
                if path:
                    if path.startswith('"') and path.endswith('"'):
                        path = path[1:-1]
                    if path != '/dev/null':
                        current_file = path
            continue
        
        m = hunk_re.match(line)
        if m and current_file:
            old_start = int(m.group(1))
            old_len = int(m.group(2)) if m.group(2) is not None else 1
            new_start = int(m.group(3))
            new_len = int(m.group(4)) if m.group(4) is not None else 1
            
            hunk_lines = []
            i += 1
            while i < len(lines) and not lines[i].startswith('@@ ') and not lines[i].startswith('diff --git '):
                hunk_lines.append(lines[i])
                i += 1
            
            hunk = {
                'old_start': old_start,
                'old_len': old_len,
                'new_start': new_start,
                'new_len': new_len,
                'header': line,
                'lines': hunk_lines,
                'file_header': list(file_header_lines)
            }
            if current_file not in files:
                files[current_file] = []
            files[current_file].append(hunk)
        else:
            i += 1
    return files

def main():
    parser = argparse.ArgumentParser(description="Extract git diff hunks covering specific file:line numbers.")
    parser.add_argument("-b", "--feature-branch", help="Feature branch to compare (defaults to current branch)")
    parser.add_argument("-m", "--main-branch", default="master", help="Main branch to compare against (defaults to master)")
    args, unknown_args = parser.parse_known_args()

    git_root = get_git_root()

    feature_branch = args.feature_branch
    if not feature_branch:
        try:
            feature_branch = get_current_branch()
        except subprocess.CalledProcessError as e:
            print(f"Error getting current branch: {e.stderr.strip()}", file=sys.stderr)
            sys.exit(1)

    diff_text = get_git_diff(args.main_branch, feature_branch, unknown_args)
    parsed_diff = parse_diff(diff_text)

    targets = []
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        if ':' in line:
            parts = line.rsplit(':', 1)
            if len(parts) == 2:
                filepath, line_str = parts
                try:
                    line_num = int(line_str)
                    abs_path = os.path.abspath(filepath)
                    rel_path = os.path.relpath(abs_path, git_root)
                    git_style_path = rel_path.replace(os.sep, '/')
                    targets.append((git_style_path, line_num))
                except ValueError:
                    continue

    printed_headers = set()
    printed_hunks = set()

    for target_file, target_line in targets:
        if target_file in parsed_diff:
            for hunk in parsed_diff[target_file]:
                start = hunk['new_start']
                length = hunk['new_len']
                if length > 0 and start <= target_line < start + length:
                    hunk_key = (target_file, hunk['header'])
                    if hunk_key not in printed_hunks:
                        header_key = tuple(hunk['file_header'])
                        if header_key not in printed_headers:
                            print('\n'.join(hunk['file_header']))
                            printed_headers.add(header_key)
                        print(hunk['header'])
                        if hunk['lines']:
                            print('\n'.join(hunk['lines']))
                        printed_hunks.add(hunk_key)

if __name__ == '__main__':
    main()
