#!/bin/bash

# Exit on any error
set -e

# Paths
ROUNDINGSAT_EXEC="./build_debug/roundingsat"

# Usage check
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <instance.opb> [optimum_value]"
    exit 1
fi

INSTANCE_PATH="$1"
OPTIMUM_VALUE="$2"

# Check executable
if [ ! -f "$ROUNDINGSAT_EXEC" ]; then
    echo "Error: RoundingSat not found at $ROUNDINGSAT_EXEC"
    exit 1
fi

echo "Starting backbone extraction for instance: $INSTANCE_PATH"

# If optimum value was given, go directly to backbone extraction
if [ -n "$OPTIMUM_VALUE" ]; then
    echo "Using provided optimum value: $OPTIMUM_VALUE"
    "$ROUNDINGSAT_EXEC" "$INSTANCE_PATH" --calculate-backbones=1 --already-solved=1 --optimum-value="$OPTIMUM_VALUE"
    exit 0
fi

# First call to solver to either extract optimum value or detect decisional instance
echo "No optimum value provided. Running RoundingSat to detect instance type..."
OUTPUT=$("$ROUNDINGSAT_EXEC" "$INSTANCE_PATH" --calculate-backbones=1 2>&1)
echo "$OUTPUT"

# Check if it's decisional
if echo "$OUTPUT" | grep -q "s CALCULATED BACKBONE"; then
    echo "Backbone successfully extracted (instance is decisional)."
    exit 0
fi

# Try to extract optimum value (from line before 's OPTIMUM FOUND')
OPTIMUM_LINE=$(echo "$OUTPUT" | grep -B1 "s OPTIMUM FOUND" | head -n 1)
OPTIMUM_VALUE=$(echo "$OPTIMUM_LINE" | awk '/^o / {print $2}')

if [ -z "$OPTIMUM_VALUE" ]; then
    echo "Error: Failed to extract optimum value from solver output."
    exit 1
fi

echo "Extracted optimum value: $OPTIMUM_VALUE"
echo "Running RoundingSat again to extract backbone..."

# Run with optimum value
"$ROUNDINGSAT_EXEC" "$INSTANCE_PATH" --calculate-backbones=1 --already-solved=1 --optimum-value="$OPTIMUM_VALUE"

echo "Backbone extraction complete."
