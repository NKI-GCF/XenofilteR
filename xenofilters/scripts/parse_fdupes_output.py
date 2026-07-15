import sys
import argparse
import re
from collections import defaultdict
from dataclasses import dataclass

@dataclass
class Member:
    desc: str
    filepath: str
    start_line: int
    end_line: int
    code: str = ""

def read_code(filepath: str, start: int, end: int) -> str:
    """Reads the specified lines from a file, returning the stripped code."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            # 1-based indexing for start and end
            code_lines = lines[max(0, start - 1):end]
            return "".join(code_lines).strip()
    except Exception as e:
        return f"// Error reading {filepath}: {e}"

def strip_rust_comments(code: str) -> str:
    """
    Strip Rust single-line (`// ...`) and multi-line (`/* ... */`) comments.
    This is a simple heuristic implementation; it does not handle all edge cases
    (like comments inside raw strings or macros), but works for typical source.
    """
    import re

    # Remove multi-line comments /* ... */
    # Use a pattern that avoids nested comments (Rust doesn't allow them anyway)
    code = re.sub(r'/\*.*?\*/', '', code, flags=re.DOTALL)

    # Remove single-line comments // ... (but not inside strings)
    # This is a simplified approach: remove // to end of line unless inside quotes.
    # For most practical cases, this is sufficient.
    def strip_line_comment(line: str) -> str:
        in_string = False
        string_char = None
        i = 0
        while i < len(line):
            c = line[i]
            if not in_string:
                if c == '"' or c == "'":
                    # Check for raw string r#"..."# etc.
                    if c == '"' and i + 1 < len(line) and line[i + 1] == '#':
                        # raw string start, skip to matching #"..."#
                        i += 2
                        while i < len(line):
                            if line[i] == '"' and i + 1 < len(line) and line[i + 1] == '#':
                                i += 2
                                break
                            i += 1
                        in_string = False
                        continue
                    in_string = True
                    string_char = c
                elif c == '/' and i + 1 < len(line) and line[i + 1] == '/':
                    # comment starts here
                    return line[:i].rstrip()
            else:
                if c == string_char and (i == 0 or line[i - 1] != '\\'):
                    in_string = False
                    string_char = None
            i += 1
        return line

    lines = code.splitlines()
    stripped_lines = [strip_line_comment(line) for line in lines]
    return "\n".join(stripped_lines)

def parse_and_process():
    parser = argparse.ArgumentParser(description="Parse cargo dupes output and format for AI processing.")
    parser.add_argument("-i", "--input", type=argparse.FileType('r'), default=sys.stdin,
                        help="Input file containing cargo dupes output (defaults to stdin).")
    parser.add_argument("--with-fingerprint", action="store_true",
                        help="Include the fingerprint in the output.")
    parser.add_argument("--code", action="store_true",
                        help="Fetch and include the relevant code blocks from the file system.")
    parser.add_argument("--compact", action="store_true",
                        help="When using --code, group identical code blocks across groups and print them once.")
    parser.add_argument("--no-comments", action="store_true",
                        help="When using --code, strip Rust comments from the fetched code blocks.")
    parser.add_argument("--prompt", type=str, nargs="?", const="default",
                        help="Prepend a prompt instructing AI to remove duplication. Provide an optional custom prefix.")

    args = parser.parse_args()

    if args.prompt is not None:
        custom_prefix = f"{args.prompt}\n" if args.prompt != "default" else ""
        print(f"{custom_prefix}Please analyze the following code duplications and refactor to remove them if possible. "
              "Use table driven tests if possible. use functions, generics or traits, and only if no other solution is possible: macros."
              "Only apply changes if they do not introduce any performance penalty and successfully decrease overall code size or complexity.\n")

    sections = []
    current_section = None
    current_group = None

    group_header_re = re.compile(r"^Group (\d+) \(fingerprint: ([a-f0-9]+)(?:, similarity: (\d+)%)?, (\d+) members\):")

    # Read and Parse Data
    for line in args.input:
        line_stripped = line.strip()
        if not line_stripped:
            continue

        # Detect new sections by looking for "Duplicates" keywords
        if line_stripped.endswith("Duplicates") and not line_stripped.startswith("Duplication"):
            current_section = {"name": line_stripped, "groups": []}
            sections.append(current_section)
            continue

        if not current_section:
            continue

        # Parse Group Header
        group_match = group_header_re.match(line_stripped)
        if group_match:
            current_group = {
                "id": group_match.group(1),
                "fingerprint": group_match.group(2),
                "similarity": group_match.group(3),
                "members": []
            }
            current_section["groups"].append(current_group)
            continue

        # Parse Members
        if line_stripped.startswith("- ") and current_group is not None:
            content = line_stripped[2:]
            if " at " in content:
                desc_part, loc_part = content.rsplit(" at ", 1)
                if ":" in loc_part:
                    file_part, lines_part = loc_part.rsplit(":", 1)
                    try:
                        if "-" in lines_part:
                            start_l, end_l = map(int, lines_part.split("-"))
                        else:
                            start_l = end_l = int(lines_part)

                        new_member = Member(desc_part.strip(), file_part.strip(), start_l, end_l)

                        # Filter out false positives (exact same file and lines in the same group)
                        is_duplicate = any(
                            m.filepath == new_member.filepath and
                            m.start_line == new_member.start_line and
                            m.end_line == new_member.end_line
                            for m in current_group["members"]
                        )

                        if not is_duplicate:
                            current_group["members"].append(new_member)
                    except ValueError:
                        pass
                        # Ignore malformed lines gracefully
    # Format Output
    for section in sections:
        # Filter groups with < 2 members after deduplication
        valid_groups = [g for g in section["groups"] if len(g["members"]) >= 2]

        if not valid_groups:
            continue

        print(f"{section['name']}")
        print("=" * len(section['name']) + "\n")

        # Display code if requested
        if args.code:
            if args.compact:
                # Compact mode: collect all members with code across all groups in this section
                code_to_occurrences = defaultdict(list)

                for group in valid_groups:
                    for m in group["members"]:
                        raw_code = read_code(m.filepath, m.start_line, m.end_line)
                        if not raw_code:
                            continue

                        code_content = raw_code
                        if args.no_comments:
                            code_content = strip_rust_comments(code_content)

                        key = (m.filepath, m.start_line, m.end_line, code_content)
                        code_to_occurrences[key].append({
                            "group_id": group["id"],
                            "desc": m.desc,
                            "lines": f"{m.start_line}-{m.end_line}",
                            "filepath": m.filepath,
                        })

                # Now print each unique code block once, with all its occurrences
                for (filepath, start_l, end_l, code_content), occurrences in code_to_occurrences.items():
                    print(f"\n    ```rust")
                    print(code_content)
                    print("    ```")

                    # Group occurrences by group id for nicer output
                    by_group = defaultdict(list)
                    for occ in occurrences:
                        by_group[occ["group_id"]].append(occ)

                    for gid, occs in by_group.items():
                        # Sort by description then lines
                        occs_sorted = sorted(occs, key=lambda x: (x["desc"], x["lines"]))
                        # Print one line per occurrence, but merge same desc+lines if they appear multiple times
                        desc_lines = defaultdict(list)
                        for o in occs_sorted:
                            desc_lines[(o["desc"], o["lines"])].append(o["filepath"])

                        for (desc, lines), fpaths in desc_lines.items():
                            # If same desc+lines appears in multiple files, show all
                            if len(fpaths) == 1:
                                print(f"      Group {gid}: {desc} (lines {lines}) in {fpaths[0]}")
                            else:
                                print(f"      Group {gid}: {desc} (lines {lines}) in {', '.join(fpaths)}")
                    print()
            else:
                for group in valid_groups:
                    # Build group header
                    parts = []
                    if args.with_fingerprint:
                        parts.append(f"fingerprint: {group['fingerprint']}")
                    if group['similarity']:
                        parts.append(f"similarity: {group['similarity']}%")
                    parts.append(f"{len(group['members'])} members")

                    header_info = ", ".join(parts)
                    print(f"Group {group['id']} ({header_info}):")

                    # Condense by file path
                    by_file = defaultdict(list)
                    for m in group["members"]:
                        by_file[m.filepath].append(m)

                    for filepath, members in by_file.items():
                        print(f"  File: {filepath}")

                        # Condense repeated messages inside the same file
                        desc_to_lines = defaultdict(list)
                        for m in members:
                            desc_to_lines[m.desc].append(f"{m.start_line}-{m.end_line}")

                        for desc, lines_list in desc_to_lines.items():
                            lines_str = ", ".join(lines_list)
                            print(f"    - {desc} (lines {lines_str})")

                    code_to_members = defaultdict(list)
                    for m in group["members"]:
                        raw_code = read_code(m.filepath, m.start_line, m.end_line)
                        if not raw_code:
                            continue

                        code_content = raw_code
                        if args.no_comments:
                            code_content = strip_rust_comments(code_content)

                        code_to_members[code_content].append(m)

                    for code_content, m_list in code_to_members.items():
                        print("\n    ```rust")
                        print(code_content)
                        print("    ```")
        print()

if __name__ == "__main__":
    parse_and_process()
