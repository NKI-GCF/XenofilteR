#!/usr/bin/env python3
import sys
import re
from collections import defaultdict

def coalesce_lines(line_nums):
    """Combines a list of line numbers into unique, sorted ranges (e.g., '144-145, 166')."""
    if not line_nums:
        return ""
    lines = sorted(set(int(l) for l in line_nums))
    ranges = []
    start = end = lines[0]
    
    for n in lines[1:]:
        if n == end + 1:
            end = n
        else:
            ranges.append(f"{start}-{end}" if start != end else str(start))
            start = end = n
    ranges.append(f"{start}-{end}" if start != end else str(start))
    
    return ",".join(ranges)

def summarize_cargo_output(input_text):
    # Map: file_path -> (error_code, message) -> list of line numbers
    summary = defaultdict(lambda: defaultdict(list))
    
    # Matches standard cargo check formats like: 
    # src/aln_stream.rs:169:33: error[E0308]: mismatched types: expected type...
    pattern = re.compile(r'^(?P<file>[^:]+):(?P<line>\d+):(?P<col>\d+):\s+(?P<lvl>\w+)(?:\[(?P<code>[^\]]+)\])?:\s+(?P<msg>.+)$')
    
    for line_raw in input_text.strip().splitlines():
        match = pattern.match(line_raw)
        if not match:
            continue
        
        file_path, line_num, _col, lvl, code, msg = match.groups()
        error_id = code if code else lvl
        
        # Clean up excessive rustc suggestion artifacts in messages if needed
        msg = msg.strip()
        
        summary[file_path][(error_id, msg)].append(line_num)
        
    output = []
    for file_path, errors in summary.items():
        output.append(f"### `{file_path}`")
        for (error_id, msg), lines in errors.items():
            line_str = coalesce_lines(lines)
            output.append(f"- `[{error_id}]` L{line_str}: {msg}")
            
    return "\n".join(output)

if __name__ == "__main__":
    # Reads from standard input (e.g., `cargo check 2>&1 | python3 summarize.py`)
    raw_input = sys.stdin.read()
    if raw_input.strip():
        print(summarize_cargo_output(raw_input))
