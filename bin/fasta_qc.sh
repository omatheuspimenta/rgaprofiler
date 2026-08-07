#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.fasta> <output.fasta>" >&2
    exit 1
fi

INPUT=$1
OUTPUT=$2
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

VALID_AA="ACDEFGHIKLMNPQRSTVWYBXZJUO"

# 1. Remove duplicate IDs (seqkit rmdup dedups by ID by default, keeping the
#    first occurrence). Log the removed ones.
seqkit rmdup -D "$TMPDIR/duplicates.log" -o "$TMPDIR/step1.fasta" "$INPUT"
if [ -s "$TMPDIR/duplicates.log" ]; then
    sed 's/^/WARNING: duplicate ID removed -> /' "$TMPDIR/duplicates.log" >&2
fi

# 2. Strip stray whitespace/CR/asterisks from the sequence only (-s), then uppercase
seqkit replace -s -p '[ \t\r\*]' -r '' "$TMPDIR/step1.fasta" \
    | seqkit seq -u -o "$TMPDIR/step2.fasta"

# 3. Replace any remaining invalid characters with X, logging how many were
#    changed per sequence. seqkit has no built-in "replace + log" combo, so
#    this step drops to tabular format, does the substitution in awk, and
#    converts back to FASTA.
seqkit fx2tab "$TMPDIR/step2.fasta" > "$TMPDIR/step2.tab"

awk -v valid="$VALID_AA" '
BEGIN { FS="\t"; OFS="\t" }
{
    header = $1
    seq    = $2
    split(header, parts, /[ \t]/)
    id = parts[1]

    if (seq == "") {
        print "WARNING: sequence " id " ignored. Empty sequence." > "/dev/stderr"
        next
    }

    count = gsub("[^" valid "]", "X", seq)
    if (count > 0) {
        print "WARNING: " count " invalid character(s) in sequence " id " replaced with X" > "/dev/stderr"
    }

    print header, seq
}
' "$TMPDIR/step2.tab" > "$TMPDIR/step3.tab"

# 4. Convert back to FASTA. -w 0 keeps sequences unwrapped (single line),
#    matching the original script's output. Drop -w 0 if you'd rather have
#    standard 60-column wrapping.
seqkit tab2fx "$TMPDIR/step3.tab" -w 0 -o "$OUTPUT"

if [ ! -s "$OUTPUT" ]; then
    echo "ERROR: no valid sequences were written to $OUTPUT" >&2
    exit 1
fi