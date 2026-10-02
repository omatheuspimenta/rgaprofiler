#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input.fasta> <output.fasta> <split_mode> <split_value>" >&2
    echo "  split_mode:  'size' (split into chunks of split_value sequences each)," >&2
    echo "               'parts' (split into exactly split_value chunks, balanced by seqkit) or" >&2
    echo "               'auto' (choose the chunk count from the sequence count and the host's CPUs)" >&2
    echo "  split_value: positive integer -- sequences-per-chunk for 'size', chunk count for 'parts'," >&2
    echo "               the host's CPU count for 'auto'" >&2
    exit 1
fi

INPUT=$1
OUTPUT=$2
SPLIT_MODE=$3
SPLIT_VALUE=$4
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

if [ "$SPLIT_MODE" != "size" ] && [ "$SPLIT_MODE" != "parts" ] && [ "$SPLIT_MODE" != "auto" ]; then
    echo "ERROR: split_mode must be 'size', 'parts' or 'auto', got '$SPLIT_MODE'" >&2
    exit 1
fi

if ! [[ "$SPLIT_VALUE" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: split_value must be a positive integer, got '$SPLIT_VALUE'" >&2
    exit 1
fi

VALID_AA="ACDEFGHIKLMNPQRSTVWYBXZJUO"

# 1. Remove duplicate IDs (seqkit rmdup dedups by ID by default, keeping the
#    first occurrence). Log the removed ones.
seqkit rmdup -D "$TMPDIR/duplicates.log" -o "$TMPDIR/step1.fasta" "$INPUT"
if [ -s "$TMPDIR/duplicates.log" ]; then
    sed 's/^/WARNING: duplicate ID removed -> /' "$TMPDIR/duplicates.log" >&2
fi

# 2. Strip stray whitespace/CR from the sequence only (-s), then the terminal stop
#    codon(s) -- only at the very end, after the whitespace is gone so 'MKV* ' still
#    counts as terminal -- then uppercase. An internal '*' is NOT deleted here:
#    that would splice the residues on either side into a sequence that doesn't
#    exist. Step 3 turns it into X and logs it instead.
seqkit replace -s -p '[ \t\r]' -r '' "$TMPDIR/step1.fasta" \
    | seqkit replace -s -p '\*+$' -r '' \
    | seqkit seq -u -o "$TMPDIR/step2.fasta"

# 3. Replace internal stop codons and any remaining invalid characters with X,
#    logging how many were changed per sequence. seqkit has no built-in
#    "replace + log" combo, so this step drops to tabular format, does the
#    substitution in awk, and converts back to FASTA.
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

    # Counted apart from other invalid characters: an internal stop usually means a
    # pseudogene or a broken gene model, which is worth flagging on its own.
    stops = gsub(/\*/, "X", seq)
    if (stops > 0) {
        print "WARNING: " stops " internal stop codon(s) (*) in sequence " id " replaced with X" > "/dev/stderr"
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

# 5. Split the cleaned FASTA into chunks (each a complete, self-contained set of
#    FASTA records -- never an arbitrary line split), so downstream Nextflow
#    processes can fan out over them. Chunks are written to a dedicated directory
#    next to $OUTPUT, named after it, e.g. "results.fasta" -> "results_chunks/
#    results.part_001.fasta". Two mutually exclusive modes, both backed by seqkit
#    split2 (which always keeps whole records together):
#      - 'size':  a fixed number of sequences per chunk (split_value each);
#                 the resulting chunk *count* depends on the input's sequence total.
#      - 'parts': a fixed chunk *count* (split_value chunks); seqkit distributes
#                 sequences across them as evenly as possible.
OUTBASE=$(basename "$OUTPUT")
OUTDIR=$(dirname "$OUTPUT")
CHUNK_PREFIX="${OUTBASE%.*}"
CHUNK_DIR="$OUTDIR/${CHUNK_PREFIX}_chunks"

# 'auto' (--num_blocks auto): pick the chunk count here, now that the cleaned sequence
# count is known, then split exactly as an explicit --num_blocks of that count would.
# Chunk size = enough chunks for ~2 waves of InterProScan tasks (4 CPUs each, the
# longest per-chunk step) across the host's CPUs, so no core idles through the tail, but
# never below AUTO_MIN_SEQS sequences per chunk -- every chunk task pays a fixed model
# load (~15-25s per GPU tool, far more on CPU) -- nor above AUTO_MAX_SEQS, which keeps
# DeepCoil2's GPU memory (grows with chunk size) and each task's runtime bounded.
AUTO_MIN_SEQS=100
AUTO_MAX_SEQS=500
if [ "$SPLIT_MODE" = "auto" ]; then
    N_SEQS=$(grep -c '^>' "$OUTPUT")
    SLOTS=$(( SPLIT_VALUE / 4 )); (( SLOTS < 1 )) && SLOTS=1
    SIZE=$(( (N_SEQS + 2 * SLOTS - 1) / (2 * SLOTS) ))
    (( SIZE < AUTO_MIN_SEQS )) && SIZE=$AUTO_MIN_SEQS
    (( SIZE > AUTO_MAX_SEQS )) && SIZE=$AUTO_MAX_SEQS
    PARTS=$(( (N_SEQS + SIZE - 1) / SIZE ))
    echo "num_blocks auto: ${N_SEQS} sequences, ${SPLIT_VALUE} CPUs -> ${PARTS} chunk(s) of ~$(( (N_SEQS + PARTS - 1) / PARTS )) sequences (rerun with --num_blocks ${PARTS} to reproduce exactly)" >&2
    SPLIT_MODE=parts
    SPLIT_VALUE=$PARTS
fi

rm -rf "$CHUNK_DIR"
if [ "$SPLIT_MODE" = "size" ]; then
    seqkit split2 -s "$SPLIT_VALUE" -O "$CHUNK_DIR" -o "$CHUNK_PREFIX" -e .fasta -w 0 "$OUTPUT"
else
    seqkit split2 -p "$SPLIT_VALUE" -O "$CHUNK_DIR" -o "$CHUNK_PREFIX" -e .fasta -w 0 "$OUTPUT"
fi

N_CHUNKS=$(find "$CHUNK_DIR" -maxdepth 1 -name '*.fasta' | wc -l)
if [ "$N_CHUNKS" -eq 0 ]; then
    echo "ERROR: chunking produced no files in $CHUNK_DIR" >&2
    exit 1
fi

echo "Split $OUTPUT into $N_CHUNKS chunk(s) (mode=$SPLIT_MODE, value=$SPLIT_VALUE) in $CHUNK_DIR" >&2