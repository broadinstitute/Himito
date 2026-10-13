#!/usr/bin/env python3
"""Fail if a Rust change touches anything but comments.

Compares each changed `.rs` file between a git revision and the working tree
after removing all comments (`//`, `///`, `//!`, nested `/* */`) and
collapsing whitespace. Strings, raw strings, byte strings and char literals
are kept verbatim, so `"a // b"` is code, not a comment.

Doc comments are removed too, so shortened clap help text passes.

usage: check_comment_only_diff.py [REV]   (default REV = HEAD)
exit 0 = comment-only; 1 = code changed (prints the first difference)
"""
import subprocess, sys

def strip(src: str) -> str:
    out, i, n = [], 0, len(src)
    while i < n:
        c = src[i]
        if src.startswith("//", i):
            j = src.find("\n", i); i = n if j < 0 else j; out.append(" "); continue
        if src.startswith("/*", i):
            depth, i = 1, i + 2
            while i < n and depth:
                if src.startswith("/*", i): depth += 1; i += 2
                elif src.startswith("*/", i): depth -= 1; i += 2
                else: i += 1
            out.append(" "); continue
        # raw strings: r"..", r#".."#, br#".."#
        if c in "rb" and (src.startswith('r"', i) or src.startswith("r#", i) or src.startswith('br"', i) or src.startswith("br#", i)):
            k = i + (2 if src[i] == "b" else 1)
            hashes = 0
            while k < n and src[k] == "#": hashes += 1; k += 1
            if k < n and src[k] == '"':
                end = src.find('"' + "#" * hashes, k + 1)
                end = n if end < 0 else end + 1 + hashes
                out.append(src[i:end]); i = end; continue
        if c == '"' or (c == "b" and src.startswith('b"', i)):
            k = i + (2 if c == "b" else 1)
            while k < n and src[k] != '"':
                k += 2 if src[k] == "\\" else 1
            out.append(src[i:k + 1]); i = k + 1; continue
        if c == "'":
            # char literal ('x', '\n', '\u{1F}', b'x' handled as 'x') vs lifetime ('a)
            if i + 2 < n and src[i + 1] == "\\":
                k = src.find("'", i + 2); out.append(src[i:k + 1]); i = k + 1; continue
            if i + 2 < n and src[i + 2] == "'":
                out.append(src[i:i + 3]); i += 3; continue
        out.append(c); i += 1
    return " ".join("".join(out).split())

def main():
    rev = sys.argv[1] if len(sys.argv) > 1 else "HEAD"
    files = subprocess.run(["git", "diff", "--name-only", rev, "--", "*.rs"], capture_output=True, text=True, check=True).stdout.split()
    bad = 0
    for f in files:
        old = subprocess.run(["git", "show", f"{rev}:{f}"], capture_output=True, text=True).stdout
        try:
            new = open(f).read()
        except FileNotFoundError:
            print(f"{f}: deleted"); bad += 1; continue
        a, b = strip(old), strip(new)
        if a != b:
            k = next((x for x in range(min(len(a), len(b))) if a[x] != b[x]), min(len(a), len(b)))
            print(f"{f}: CODE CHANGED near ...{a[max(0,k-60):k+60]!r}\n   vs ...{b[max(0,k-60):k+60]!r}")
            bad += 1
        else:
            print(f"{f}: comment-only")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
