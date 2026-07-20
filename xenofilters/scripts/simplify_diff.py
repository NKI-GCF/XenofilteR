#!/usr/bin/env python3
import sys
import argparse
import re

def parse_diff(diff_text):
    files_list = []
    files_dict = {}
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
                if lines[i].startswith('--- '):
                    path = lines[i][4:].strip()
                    if path.startswith('a/'):
                        path = path[2:]
                    if path != '/dev/null':
                        current_file = path
                elif lines[i].startswith('+++ '):
                    path = lines[i][4:].strip()
                    if path.startswith('b/'):
                        path = path[2:]
                    if path != '/dev/null':
                        current_file = path
                
                if current_file and current_file.startswith('"') and current_file.endswith('"'):
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
            if current_file not in files_dict:
                files_dict[current_file] = []
                files_list.append(current_file)
            files_dict[current_file].append(hunk)
        else:
            i += 1
    return files_list, files_dict

def main():
    parser = argparse.ArgumentParser(description="Convert git diff output into simplified format.")
    parser.add_argument("--original", action="store_true", help="Keep original diff format")
    parser.add_argument("--comments", action="store_true", help="Keep comment lines")
    parser.add_argument("--space", action="store_true", help="Keep empty lines")
    args = parser.parse_args()

    diff_text = sys.stdin.read()
    files_list, parsed_diff = parse_diff(diff_text)

    for target_file in files_list:
        printed_headers = False
        for hunk in parsed_diff[target_file]:
            if not printed_headers:
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
                                prefix = fhl[:4]
                                path = fhl[4:]
                                if path.startswith('a/') or path.startswith('b/'):
                                    path = path[2:]
                                print(f"{prefix}{path}")
                else:
                    print('\n'.join(hunk['file_header']))
                printed_headers = True

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

if __name__ == '__main__':
    main()
