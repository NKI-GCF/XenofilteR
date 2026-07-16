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
                if lines[i].startswith('+++ b/'):
                    current_file = lines[i][6:]
                    if current_file.startswith('"') and current_file.endswith('"'):
                        current_file = current_file[1:-1]
                i += 1
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
    parser.add_argument("--original", action="store_true", help="Keep original diff format")
    parser.add_argument("--comments", action="store_true", help="Keep comment lines")
    parser.add_argument("--space", action="store_true", help="Keep empty lines")
    args, unknown_args = parser.parse_known_args()

    git_root = get_git_root()

    feature_branch = args.feature_branch
    if not feature_branch:
        try:
            feature_branch = get_current_branch()
        except subprocess.CalledProcessError as e:
            print(f"Error getting current branch: {e.stderr.strip()}", file=sys.stderr)
            sys.exit(1)

    diff_text = get_git_diff(args.main_branch, feature_branch)
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
                            if not args.original:
                                a_file, b_file = None, None
                                for fhl in hunk['file_header']:
                                    if fhl.startswith('--- '):
                                        path = fhl[4:].strip()
                                        if path.startswith('a/'):
                                            path = path[2:]
                                        if path.startswith('"') and path.endswith('"'):
                                            path = path[1:-1]
                                        a_file = path
                                    elif fhl.startswith('+++ '):
                                        path = fhl[4:].strip()
                                        if path.startswith('b/'):
                                            path = path[2:]
                                        if path.startswith('"') and path.endswith('"'):
                                            path = path[1:-1]
                                        b_file = path
                                if a_file and b_file and a_file == b_file and a_file != '/dev/null':
                                    print(f"-/+ {a_file}")
                                else:
                                    for fhl in hunk['file_header']:
                                        if fhl.startswith('--- ') or fhl.startswith('+++ '):
                                            print(fhl)
                            else:
                                print('\n'.join(hunk['file_header']))
                            printed_headers.add(header_key)

                        lines_to_print = hunk['lines']
                        if not args.original:
                            lines_to_print = list(lines_to_print)
                            while lines_to_print and lines_to_print[0].startswith(' '):
                                lines_to_print.pop(0)

                            processed = []
                            in_tracing = False
                            for line in lines_to_print:
                                if not line: continue
                                prefix = line[0]
                                text = line[1:]
                                text_stripped = text.strip()

                                if not args.space and not text_stripped:
                                    continue
                                if not args.comments and (text_stripped.startswith('//') or text_stripped.startswith('#')):
                                    continue

                                if in_tracing:
                                    if ');' in text_stripped:
                                        in_tracing = False
                                    continue
                                elif 'tracing::' in text_stripped:
                                    if ');' not in text_stripped:
                                        in_tracing = True
                                    continue

                                processed.append(line)
                            lines_to_print = processed

                        if not args.original:
                            print('@')
                        else:
                            print(hunk['header'])

                        if lines_to_print:
                            print('\n'.join(lines_to_print))
                        printed_hunks.add(hunk_key)

if __name__ == '__main__':
    main()
