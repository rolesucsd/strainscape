#!/usr/bin/env bash

# --------------------------
# 1) Define your filters:
#    You can hardcode them here or read them from arguments.
# --------------------------
STRING="G000436455"        # The string to match in the 1st column
COL2_MIN=$(echo "4499379 * 0.25" | bc -l)
COL2_MAX=$(echo "4499379 * 0.38" | bc -l)

# Construct output filename dynamically
OUTPUT="${STRING}_${COL2_MIN}-${COL2_MAX}.txt"

# Remove old summary file if it exists
rm -f "$OUTPUT"

echo "Filtering all *.cov files for lines where:"
echo "  1st col == $STRING"
echo "  2nd col in [$COL2_MIN, $COL2_MAX]"
echo "  3rd col in [$COL2_MIN, $COL2_MAX]"
echo "Writing results to: $OUTPUT"

# --------------------------
# 2) Loop through all .cov files
# --------------------------
for FILE in *.cov; do
  if [[ -f "$FILE" ]]; then
    # For each line in $FILE, check conditions with awk
    # If conditions pass, print the entire line + a tab + filename
    awk -v str="$STRING" \
        -v c2min="$COL2_MIN" -v c2max="$COL2_MAX" \
        -v c3min="$COL2_MIN" -v c3max="$COL2_MAX" \
        -v fname="$FILE" '
        ($1 == str && $2 >= c2min && $2 <= c2max && $3 >= c3min && $3 <= c3max) {
            print $0 "\t" fname
        }
    ' "$FILE" >> "$OUTPUT"
  fi
done

echo "Done. Summary written to $OUTPUT."