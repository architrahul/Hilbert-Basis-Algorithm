#!/usr/bin/env python3
"""
Build a local SQLite database of covering designs from the
La Jolla Covering Repository (LJCR).

Features:
- Resume-safe (can be stopped and restarted)
"""

import urllib.request
import sqlite3
import time
import re
import sys

# ----------------------------
# Configuration
# ----------------------------
DB_FILE = "coverings.db"
REQUEST_DELAY = 0.5

V_RANGE = range(2, 100)   # v < 100
T_RANGE = range(2, 9)     # t <= 8
K_MAX = 25                # k <= 25
# ----------------------------


def download_covering_design(v, k, t):
    """
    Download C(v,k,t) from LJCR and return list of blocks
    (each block is a list of ints, 1-indexed).
    Returns None if not available.
    """
    url = f"https://ljcr.dmgordon.org/cover/show_cover.php?v={v}&k={k}&t={t}"

    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            html = response.read().decode("utf-8")
    except Exception:
        return None

    # Strip HTML tags
    text = re.sub(r"<[^>]+>", "\n", html)

    blocks = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        clean = (
            line.replace("{", "")
                .replace("}", "")
                .replace(",", " ")
        )
        clean = re.sub(r"\s+", " ", clean)

        if re.fullmatch(r"[0-9 ]+", clean):
            parts = clean.split()
            if len(parts) == k:
                try:
                    block = [int(x) for x in parts]
                    if all(1 <= x <= v for x in block):
                        blocks.append(block)
                except ValueError:
                    pass

    return blocks if blocks else None


def design_exists(cur, v, k, t):
    """
    Check if *any* block of C(v,k,t) is already stored.
    If yes, we assume the design was already ingested.
    """
    cur.execute(
        "SELECT 1 FROM coverings WHERE v=? AND k=? AND t=? LIMIT 1",
        (v, k, t)
    )
    return cur.fetchone() is not None


def main():
    conn = sqlite3.connect(DB_FILE)
    cur = conn.cursor()

    # --- schema (resume-safe) ---
    cur.execute("""
    CREATE TABLE IF NOT EXISTS coverings (
        v INTEGER,
        k INTEGER,
        t INTEGER,
        block TEXT,
        UNIQUE(v, k, t, block)
    )
    """)

    cur.execute("""
    CREATE INDEX IF NOT EXISTS idx_vkt
    ON coverings (v, k, t)
    """)

    conn.commit()

    print("Starting / resuming database build")
    print(f"Database: {DB_FILE}")
    print("=" * 70)

    try:
        for v in V_RANGE:
            for t in T_RANGE:
                if t >= v:
                    continue

                for k in range(t + 1, min(v, K_MAX) + 1):
                    # Resume logic
                    if design_exists(cur, v, k, t):
                        print(f"Skipping C({v},{k},{t}) (already stored)")
                        continue

                    print(f"Fetching C({v},{k},{t})...")

                    blocks = download_covering_design(v, k, t)
                    if not blocks:
                        print("  not available")
                        continue

                    for block in blocks:
                        cur.execute(
                            "INSERT OR IGNORE INTO coverings VALUES (?, ?, ?, ?)",
                            (v, k, t, ",".join(map(str, block)))
                        )

                    conn.commit()
                    print(f"  stored {len(blocks)} blocks")

                    time.sleep(REQUEST_DELAY)

    except KeyboardInterrupt:
        print("\nInterrupted by user. Progress saved.")
    except Exception as e:
        print("\nError occurred:", e)
        print("Progress saved. You can safely rerun the script.")
    finally:
        conn.close()
        print("Database connection closed.")


if __name__ == "__main__":
    main()
