import os
import re
from collections import defaultdict

def extract_function_name(description):
    """
    Heuristically extract the function or method name from the cargo-mutants description.
    """
    # 1. Look for "in <function_name>" at the end of the line
    match = re.search(r'\bin\s+(.+)$', description)
    if match:
        return match.group(1).strip()
    
    # 2. Look for "replace <function_name> ->" 
    match = re.search(r'^replace\s+(.*?)\s*->', description)
    if match:
        return match.group(1).strip()
    
    # 3. Look for "replace <function_name> with"
    # (Avoid capturing simple operator replacements if they don't have an "in" clause, 
    # though valid cargo-mutants output usually appends "in func" for operators).
    match = re.search(r'^replace\s+(.*?)\s+with\b', description)
    if match:
        func = match.group(1).strip()
        # Filter out simple operator replacements just in case
        if len(func) > 2 and func not in ['==', '!=', '&&', '||', '<', '>']:
            return func
            
    return "<unknown_function>"

def main(mutants_dir="mutants.out"):
    # Data structures
    # file_line_stats: dict[file][line] = {'missed': int, 'caught': int}
    file_line_stats = defaultdict(lambda: defaultdict(lambda: {'missed': 0, 'caught': 0}))
    
    # func_stats: dict[func_name] = {'missed': int, 'caught': int, 'lines': set()}
    func_stats = defaultdict(lambda: {'missed': 0, 'caught': 0, 'lines': set()})
    
    # Map filenames to outcomes
    # 'caught' and 'timeout' both count as successful tests (the test suite caught the mutant)
    # 'missed' means the mutant survived (the test suite failed to catch it)
    # 'unviable' means it didn't compile, so we typically exclude it from coverage metrics
    categories = {
        "missed.txt": "missed",
        "caught.txt": "caught",
        "timeout.txt": "caught" 
    }
    
    pattern = re.compile(r'^([^:]+):(\d+):(\d+):\s*(.*)$')
    
    for filename, outcome in categories.items():
        filepath = os.path.join(mutants_dir, filename)
        if not os.path.exists(filepath):
            continue
            
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                match = pattern.match(line)
                if match:
                    src_file = match.group(1)
                    line_num = int(match.group(2))
                    description = match.group(4)
                    
                    func_name = extract_function_name(description)
                    
                    # Update Line Stats
                    file_line_stats[src_file][line_num][outcome] += 1
                    
                    # Update Function Stats
                    func_stats[func_name][outcome] += 1
                    func_stats[func_name]['lines'].add(f"{src_file}:{line_num}")

    # --- 1. Summarize per file/line ---
    print(f"{'='*60}")
    print("FILE / LINE SUMMARY (Failures vs Total Mutants)")
    print(f"{'='*60}")
    
    for src_file, lines in sorted(file_line_stats.items()):
        print(f"\nFile: {src_file}")
        for line_num, stats in sorted(lines.items()):
            missed = stats['missed']
            caught = stats['caught']
            total = missed + caught
            if total > 0:
                print(f"  Line {line_num:<4}: {missed} failures out of {total} tests "
                      f"(Survived: {missed/total*100:.1f}%)")
    
    # --- 2. Evaluate per function ---
    print(f"\n{'='*60}")
    print("FUNCTION TEST COVERAGE & DENSITY")
    print(f"{'='*60}")
    
    # Sort functions by lowest coverage, then by highest total mutants
    def sort_key(item):
        name, stats = item
        total = stats['missed'] + stats['caught']
        cov = stats['caught'] / total if total > 0 else 0
        return (cov, -total, name)

    for func_name, stats in sorted(func_stats.items(), key=sort_key):
        missed = stats['missed']
        caught = stats['caught']
        total = missed + caught
        unique_lines = len(stats['lines'])
        
        if total == 0:
            continue
            
        coverage = (caught / total) * 100
        # Density metric: Average number of generated mutants per mutated line in this function
        density = total / unique_lines if unique_lines > 0 else 0
        
        print(f"\nFunction: {func_name}")
        print(f"  Coverage Score : {coverage:.1f}% ({caught} caught / {total} total mutants)")
        print(f"  Failures       : {missed} missed")
        print(f"  Mutant Density : {density:.1f} mutants/line (across {unique_lines} unique lines)")

if __name__ == "__main__":
    # You can change the directory path if you are running this from outside the cargo project root
    main(mutants_dir="mutants.out")
