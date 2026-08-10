import sys
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pass1dir", type=Path)
    args = parser.parse_args()

    brace_parts = []
    for run_dir in sorted(args.pass1dir.iterdir()):
        if not run_dir.is_dir():
            continue
        catgt_subdirs = [d for d in run_dir.iterdir() if d.is_dir() and d.name.startswith("catgt_")]
        for catgt_dir in catgt_subdirs:
            brace_parts.append(f'{{{run_dir},{catgt_dir.name}}}')

    if not brace_parts:
        print(f"No catgt_* output folders found under {args.pass1dir}", file=sys.stderr)
        sys.exit(1)

    print("".join(brace_parts))

if __name__ == "__main__":
    main()