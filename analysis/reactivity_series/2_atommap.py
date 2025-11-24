#!/usr/bin/env python3

from pathlib import Path

IN_DIR = Path("out/1_project")
OUT_DIR = Path("out/2_map")
ONLY_TABLES = []

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Stage 2 stub")

if __name__ == "__main__":
    main()
