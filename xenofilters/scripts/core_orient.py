#!/usr/bin/env python3
import argparse
import sys
import subprocess
from pathlib import Path

def parse_tags_file(tags_path="./tags"):
    """Parses a standard ctags file into a lookup dictionary of {tag_name: [tag_entries]}."""
    tags = {}
    if not Path(tags_path).exists():
        return tags
    
    try:
        with open(tags_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("!_TAG_"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    tag_name = parts[0]
                    file_path = parts[1]
                    ex_cmd = parts[2] # Typically a line number or a search pattern
                    
                    # Extract line number if explicitly provided in ctags fields (usually "line:N")
                    line_num = None
                    for part in parts[3:]:
                        if part.startswith("line:"):
                            line_num = part.split(":")[1]
                            break
                    
                    # Fallback to parsing ex_cmd if no explicit field exists
                    if not line_num:
                        if ex_cmd.isdigit():
                            line_num = ex_cmd
                        else:
                            # Clean up search pattern, e.g., /^fn my_func/;"
                            clean_pattern = ex_cmd.lstrip("/^").rstrip("/;\"")
                            line_num = f"pattern:{clean_pattern}"
                            
                    tags.setdefault(tag_name, []).append({
                        "file": file_path,
                        "line": line_num,
                        "raw_cmd": ex_cmd
                    })
    except Exception as e:
        sys.stderr.write(f"Warning: Failed to parse tags file: {e}\n")
    return tags

def find_line_by_pattern(file_path, pattern_str):
    """Fallback helper to locate a regex/string pattern inside a file and return line number."""
    if pattern_str.startswith("pattern:"):
        pattern = pattern_str[len("pattern:"):]
    else:
        pattern = pattern_str

    pattern = pattern.replace("$", "").replace("^", "")
    pattern_norm = " ".join(pattern.split())
    
    try:
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            for idx, line in enumerate(f, 1):
                line_norm = " ".join(line.split())
                if pattern_norm in line_norm:
                    return str(idx)
    except Exception:
        pass
    return "1"

def get_ripgrep_matches(identifier):
    """Runs ripgrep to find occurrences of the identifier using strict word boundaries."""
    try:
        # Runs: rg -n -H --no-heading -w <identifier>
        result = subprocess.run(
            ["rg", "-n", "-H", "--no-heading", "-w", identifier],
            capture_output=True, text=True, check=False
        )
        matches = []
        if result.returncode == 0 and result.stdout:
            for line in result.stdout.strip().split("\n"):
                parts = line.split(":", 2)
                if len(parts) == 3:
                    matches.append({
                        "file": parts[0],
                        "line": parts[1],
                        "content": parts[2].strip()
                    })
        return matches
    except FileNotFoundError:
        sys.stderr.write("Error: 'rg' (ripgrep) is not installed or not in system PATH.\n")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Orient an AI agent by compiling call sites and declarations for codebase identifiers."
    )
    parser.add_argument(
        "-l", "--file-line",
        action="store_true",
        help="Only output unique file:line pairs."
    )
    parser.add_argument(
        "-t", "--tags-file",
        default="tags",
        help="Path to the ctags file (default: ./tags)."
    )
    args = parser.parse_args()

    # Process names from stdin
    identifiers = [line.strip() for line in sys.stdin if line.strip()]
    if not identifiers:
        sys.stderr.write("Error: No identifiers provided on stdin.\n")
        sys.exit(1)

    tags_db = parse_tags_file(args.tags_file)
    file_line_set = set()

    for ident in identifiers:
        # 1. Resolve Declarations
        decls = []
        tag_entries = tags_db.get(ident, [])
        for entry in tag_entries:
            file_path = entry["file"]
            line_val = entry["line"]
            if line_val and line_val.startswith("pattern:"):
                line_num = find_line_by_pattern(file_path, line_val)
            else:
                line_num = line_val or "1"
            
            snippet = ""
            try:
                with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                    lines = f.readlines()
                    idx = int(line_num) - 1
                    if 0 <= idx < len(lines):
                        snippet = lines[idx].strip()
            except Exception:
                pass

            decls.append({
                "file": file_path,
                "line": line_num,
                "content": snippet
            })
            file_line_set.add(f"{file_path}:{line_num}")

        # 2. Resolve Usages
        usages = get_ripgrep_matches(ident)
        
        # Filter call sites to exclude declaration lines so we don't display duplicates
        decl_keys = {f"{d['file']}:{d['line']}" for d in decls}
        unique_usages = []
        for u in usages:
            key = f"{u['file']}:{u['line']}"
            file_line_set.add(key)
            if key not in decl_keys:
                unique_usages.append(u)

        # Output formatting depending on flags
        if not args.file_line:
            print(f"# IDENTIFIER: {ident}")
            
            if decls:
                print("## DECLARATIONS:")
                for d in decls:
                    print(f"  {d['file']}:{d['line']}  ->  {d['content']}")
            
            if unique_usages:
                print("\n## CALL SITES / USAGES (from ripgrep):")
                for u in unique_usages:
                    print(f"  {u['file']}:{u['line']}  ->  {u['content']}")
            print()

    if args.file_line:
        # Output strictly file:line (sorted for deterministic AI consumption)
        for item in sorted(file_line_set):
            print(item)

if __name__ == "__main__":
    main()
