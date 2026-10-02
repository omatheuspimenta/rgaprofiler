#!/usr/bin/env bash
# Tests for bin/fasta_qc.sh's stop-codon handling: a terminal '*' is stripped silently,
# an internal '*' becomes X and is logged as a WARNING (never deleted, which would splice
# its neighbours into a sequence that doesn't exist).
# Needs seqkit on PATH, or docker (falls back to the same seqkit image FASTA_QC uses).
# Run: bash tests/bin/test_fasta_qc.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$HERE/../../bin/fasta_qc.sh"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fails=0

check() { # name expected actual
    if [[ "$2" == "$3" ]]; then echo "ok   - $1"; else echo "FAIL - $1: expected '$2', got '$3'"; fails=$((fails+1)); fi
}
seq_of() { awk -v id="$2" '/^>/ { keep = (substr($1, 2) == id); next } keep' "$1"; }

run_qc() { # input output  (stderr -> $TMP/qc.err)
    cp "$SCRIPT" "$TMP/fasta_qc.sh"
    if command -v seqkit >/dev/null 2>&1; then
        (cd "$TMP" && bash fasta_qc.sh "$1" "$2" size 100) 2> "$TMP/qc.err"
    else
        docker run --rm -u "$(id -u):$(id -g)" -v "$TMP:/w" -w /w \
            quay.io/biocontainers/seqkit:2.13.0--he881be0_0 \
            bash fasta_qc.sh "$1" "$2" size 100 2> "$TMP/qc.err"
    fi
}

printf '>trailing\nMKV*\n>internal\nMK*LV\n>none\nMKLV\n>wrapped_trailing\nMKVL\nAY*\n>space_crlf\r\nmkv* \r\n>only_stop\n*\n>internal_and_trailing\nM*K*V*\n' > "$TMP/in.fasta"
run_qc in.fasta out.fasta
check "fasta_qc.sh exits 0" 0 "$?"
OUT="$TMP/out.fasta"

# Output sequences
check "trailing *: stripped"                   MKV    "$(seq_of "$OUT" trailing)"
check "internal *: X, neighbours not joined"   MKXLV  "$(seq_of "$OUT" internal)"
check "no *: unchanged"                        MKLV   "$(seq_of "$OUT" none)"
check "trailing * on a wrapped last line"      MKVLAY "$(seq_of "$OUT" wrapped_trailing)"
check "trailing * before CR/space, uppercased" MKV    "$(seq_of "$OUT" space_crlf)"
check "internal + trailing: only last stripped" MXKXV "$(seq_of "$OUT" internal_and_trailing)"
check "'*'-only sequence dropped"              0      "$(grep -c '^>only_stop' "$OUT")"

# Warnings
check "internal * logged" 1 \
    "$(grep -c '^WARNING: 1 internal stop codon(s) (\*) in sequence internal replaced with X$' "$TMP/qc.err")"
check "two internal * logged as 2" 1 \
    "$(grep -c '^WARNING: 2 internal stop codon(s) (\*) in sequence internal_and_trailing replaced with X$' "$TMP/qc.err")"
check "trailing-only sequences not warned about" 0 \
    "$(grep -Ec 'sequence (trailing|none|wrapped_trailing|space_crlf) ' "$TMP/qc.err")"
check "internal * not double-counted as invalid" 0 \
    "$(grep -c 'invalid character(s) in sequence internal ' "$TMP/qc.err")"

echo
if (( fails )); then echo "$fails test(s) FAILED"; exit 1; fi
echo "all tests passed"
