#!/usr/bin/env python3
import argparse
import subprocess
import json
import re
import sys
import os

KEYWORDS = {"mut", "dyn", "impl", "for", "where", "fn", "struct", "enum", "trait", "const", "static", "ref", "in"}
RUST_IDENT_RE = re.compile(r'\b[a-zA-Z_][a-zA-Z0-9_:]*\b')

def parse_args():
    parser = argparse.ArgumentParser(description="Extract unique identifiers or file:line from Cargo diagnostics.")
    parser.add_argument("--warnings", action="store_true", help="Include warnings alongside errors")
    parser.add_argument("--clippy", action="store_true", help="Run cargo clippy instead of cargo check")
    parser.add_argument("--test", action="store_true", help="Include --tests flag")
    parser.add_argument("--bench", action="store_true", help="Include --benches flag")
    parser.add_argument("--file-line", action="store_true", help="Output ONLY file:line for the errors/warnings")
    parser.add_argument("--tags", type=str, default="tags", help="Path to ctags file (default: 'tags')")
    return parser.parse_args()

def load_tags(filepath):
    """Loads tag names from a standard ctags file into a set."""
    tags = set()
    if os.path.isfile(filepath):
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('!'): 
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 4:
                    # parts[3] is the tag kind. Skip 'c' (implementations of external types)
                    if parts[3] != 'c':
                        tags.add(parts[0])
                elif parts and parts[0]:
                    tags.add(parts[0])
    return tags

def extract_identifiers(text, valid_tags):
    """Extracts valid Rust paths from text, stripping generics, refs, and filtering by tags."""
    found = set()
    # RUST_IDENT_RE natively drops `<>`, `&`, `[]`, splitting generics into separate words
    for match in RUST_IDENT_RE.finditer(text):
        ident = match.group(0)
        
        if ident in KEYWORDS:
            continue
            
        # If it's a path like a::b::C, check if the base name 'C' is in our tags file
        base_name = ident.split("::")[-1]
        
        if valid_tags:
            if base_name in valid_tags or ident in valid_tags:
                found.add(ident)
        else:
            # If no tags file was found/provided, default to accepting it
            found.add(ident)
            
    return found

def main():
    args = parse_args()
    valid_tags = load_tags(args.tags)
    
    if not valid_tags and not args.file_line:
        print(f"Warning: Tags file '{args.tags}' not found. Cannot filter external crate identifiers.", file=sys.stderr)

    cmd = ["cargo", "clippy" if args.clippy else "check", "--message-format=json"]
    if args.test:
        cmd.append("--tests")
    if args.bench:
        cmd.append("--benches")

    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    except FileNotFoundError:
        print("Error: 'cargo' is not in PATH.", file=sys.stderr)
        sys.exit(1)

    unique_findings = set()

    for line in proc.stdout:
        line = line.strip()
        if not line or not line.startswith("{"):
            continue
            
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            continue

        if msg.get("reason") != "compiler-message":
            continue

        diag = msg.get("message", {})
        level = diag.get("level")

        is_error = level == "error"
        is_warning = level == "warning"
        
        if is_error or (args.warnings and is_warning):
            if args.file_line:
                # ONLY output file:line
                for span in diag.get("spans", []):
                    if span.get("is_primary"):
                        file_name = span.get("file_name", "unknown")
                        line_start = span.get("line_start", "?")
                        unique_findings.add(f"{file_name}:{line_start}")
                        break
            else:
                # Extract and filter identifiers
                text = diag.get("message", "")
                backticked_strings = re.findall(r'`([^`]+)`', text)
                
                for bt_str in backticked_strings:
                    unique_findings.update(extract_identifiers(bt_str, valid_tags))

    proc.wait()

    for finding in sorted(unique_findings):
        print(finding)

if __name__ == "__main__":
    main()
