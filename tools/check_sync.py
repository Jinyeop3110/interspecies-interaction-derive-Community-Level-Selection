#!/usr/bin/env python3
"""Report drift between this release and the working repository.

The release is meant to be *generated*, not maintained as a parallel fork. Any
file listed in ``manifest.tsv`` is a verbatim copy of a file in the working
repository; this script tells you when the two have diverged, and can re-copy.

    python tools/check_sync.py --source ~/path/to/working/repo
    python tools/check_sync.py --source ~/path/to/working/repo --update

Exits non-zero if anything is out of sync, so it can gate a release.
"""

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path(__file__).resolve().parent / "manifest.tsv"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_manifest():
    entries = []
    for line_number, raw in enumerate(MANIFEST.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            sys.exit(f"{MANIFEST}:{line_number}: expected two tab-separated paths")
        entries.append((parts[0].strip(), parts[1].strip()))
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path,
                        help="root of the working repository")
    parser.add_argument("--update", action="store_true",
                        help="re-copy drifted files from the working repository")
    args = parser.parse_args()

    source_root = args.source.expanduser().resolve()
    if not source_root.is_dir():
        sys.exit(f"source repository not found: {source_root}")

    drifted, missing, ok = [], [], 0

    for public_rel, source_rel in read_manifest():
        public_path = REPO_ROOT / public_rel
        source_path = source_root / source_rel

        if not source_path.exists():
            missing.append((public_rel, source_rel, "source missing"))
            continue
        if not public_path.exists():
            missing.append((public_rel, source_rel, "release copy missing"))
            continue

        if digest(public_path) == digest(source_path):
            ok += 1
            continue

        drifted.append((public_rel, source_rel))
        if args.update:
            public_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_path, public_path)

    print(f"{ok} file(s) in sync")

    for public_rel, source_rel in drifted:
        verb = "updated" if args.update else "DRIFT"
        print(f"  {verb}: {public_rel}  <-  {source_rel}")

    for public_rel, source_rel, reason in missing:
        print(f"  MISSING ({reason}): {public_rel}  <-  {source_rel}")

    if missing or (drifted and not args.update):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
