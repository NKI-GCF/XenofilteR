#!/usr/bin/env python3
import sys
import re

def parse_cargo_check():
    current_block = []
    is_error = False
    has_note_or_help = False
    seen_keyword = False

    # Regex to strip ANSI color and formatting codes
    ansi_escape = re.compile(r'\x1b\[[0-9;]*[a-zA-Z]')

    # Matches rustc location pointers like "   --> src/lib.rs:155:36"
    file_line_re = re.compile(r'^(\s*--+>\s*)([^:]+):(\d+)(?::\d+)?')

    # Matches "note: ", "help: ", "= note: ", "= help: "
    keyword_re = re.compile(r'^\s*(= )?(note|help):')

    def flush_block():
        # Only print the block if it's an error AND contains note/help
        if is_error and has_note_or_help:
            for b_line in current_block:
                print(b_line, end="")

    for line in sys.stdin:
        # Create a colorless version of the line for our logic checks
        clean_line = ansi_escape.sub('', line)

        # Detect the start of a new diagnostic block
        if clean_line.startswith("error:") or clean_line.startswith("error["):
            flush_block()
            current_block = [line] # Keep the original line with its colors
            is_error = True
            has_note_or_help = False
            seen_keyword = False

        # Ignore warnings and other cargo outputs
        elif (clean_line.startswith("warning:") or clean_line.startswith("warning[") or
              clean_line.startswith("Finished") or clean_line.startswith("Compiling") or
              clean_line.startswith("Blocking")):
            flush_block()
            current_block = []
            is_error = False

        else:
            if not is_error:
                continue

            # Flag if we've encountered a note/help in this specific block
            if keyword_re.match(clean_line):
                has_note_or_help = True
                seen_keyword = True

            # Check if line contains a file path using the clean string
            match = file_line_re.match(clean_line)
            if match:
                prefix = match.group(1)
                filepath = match.group(2)
                linenum = match.group(3)

                # We construct the gvim command text (it will output in default terminal color)
                if seen_keyword:
                    modified_line = f"{prefix}gvim {filepath} +{linenum}\n"
                else:
                    modified_line = f"gvim {filepath} +{linenum}\n"

                current_block.append(modified_line)
            else:
                # Append the original colored line for all other output
                current_block.append(line)

    # Flush the last block when the stream ends
    flush_block()

if __name__ == "__main__":
    parse_cargo_check()
