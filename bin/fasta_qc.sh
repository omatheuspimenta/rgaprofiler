#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.fasta> <output.fasta> <chunk_size>" >&2
    echo "  chunk_size: number of sequences per chunk (positive integer)" >&2
    exit 1
fi

INPUT=$1
OUTPUT=$2
CHUNK_SIZE=$3
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

if ! [[ "$CHUNK_SIZE" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: chunk_size must be a positive integer, got '$CHUNK_SIZE'" >&2
    exit 1
fi

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

# 5. Split the cleaned FASTA into chunks of $CHUNK_SIZE sequences each, so
#    downstream Nextflow processes can fan out over them. Chunks are written
#    to a dedicated directory next to $OUTPUT, named after it, e.g.
#    "results.fasta" -> "results_chunks/results.part_001.fasta"
OUTBASE=$(basename "$OUTPUT")
OUTDIR=$(dirname "$OUTPUT")
CHUNK_PREFIX="${OUTBASE%.*}"
CHUNK_DIR="$OUTDIR/${CHUNK_PREFIX}_chunks"

rm -rf "$CHUNK_DIR"
seqkit split2 -s "$CHUNK_SIZE" -O "$CHUNK_DIR" -o "$CHUNK_PREFIX" -e .fasta -w 0 "$OUTPUT"

N_CHUNKS=$(find "$CHUNK_DIR" -maxdepth 1 -name '*.fasta' | wc -l)
if [ "$N_CHUNKS" -eq 0 ]; then
    echo "ERROR: chunking produced no files in $CHUNK_DIR" >&2
    exit 1
fi

echo "Split $OUTPUT into $N_CHUNKS chunk(s) of up to $CHUNK_SIZE sequences in $CHUNK_DIR" >&2