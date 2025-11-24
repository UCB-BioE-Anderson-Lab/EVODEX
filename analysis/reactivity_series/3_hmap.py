#!/usr/bin/env python3

from pathlib import Path

IN_DIR = Path("out/2_map")
OUT_DIR = Path("out/3_hmap")
ONLY_TABLES = []

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Stage 3 stub")

if __name__ == "__main__":
    main()
