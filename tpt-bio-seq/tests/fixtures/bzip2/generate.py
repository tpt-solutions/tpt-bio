"""Writes the bzip2 fixtures used by the decoder tests in src/cram/bzip2.rs.

Compression uses Python's bz2 module, which wraps libbzip2 (test-only oracle; the
decoder itself is pure Rust). The plaintext is regenerated in the Rust tests from
the same formulas, so only the compressed bytes are stored.
"""

import bz2
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def pattern(n):
    # Compressible but not trivial: a slowly varying index into ten letters.
    return bytes(b"abcdefghij"[((i * i) // 17 + i // 64) % 10] for i in range(n))


def runs():
    out = bytearray()
    for k in range(20):
        out += b"q" * 300
        out += b"b" * (k % 5)
        out.append(65 + k)
    return bytes(out)


def write(name, data, level=9, parts=None):
    compressed = bz2.compress(data, compresslevel=level)
    with open(os.path.join(HERE, name), "wb") as f:
        f.write(compressed)
    print(name, len(data), "->", len(compressed))


write("empty.bz2", b"")
write("one.bz2", b"x")
write("runs.bz2", runs())
write("pattern.bz2", pattern(3000))
write("blocks.bz2", pattern(250000), level=1)
first = bz2.compress(pattern(1000))
second = bz2.compress(b"second stream")
with open(os.path.join(HERE, "concat.bz2"), "wb") as f:
    f.write(first + second)
print("concat.bz2", 1000 + 13, "->", len(first) + len(second))
