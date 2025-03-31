#!/bin/bash

set -euo pipefail
ROUNDINGSAT_EXEC="./build_debug/roundingsat"
DEFAULT_TIMEOUT=3600  # 1 hour

print_help() {
    echo "Usage: $0 <instance.opb> [--opt <optimum_value>] [--timeout <seconds>]"
    echo ""
    echo "Arguments:"
    echo "  <instance.opb>            Path to the .opb instance"
    echo "  --opt <value>             (Optional) Provide the optimum value directly"
    echo "  --timeout <seconds>       (Optional) Time limit in seconds (default: 3600s)"
    echo ""
    echo "Examples:"
    echo "  $0 instance.opb"
    echo "  $0 instance.opb --timeout 300"
    echo "  $0 instance.opb --opt 12345"
    echo "  $0 instance.opb --opt 12345 --timeout 1200"
}

if [ "$#" -lt 1 ]; then
    print_help
    exit 1
fi

INSTANCE_PATH=""
OPTIMUM_VALUE=""
TOTAL_LIMIT=""

# Parse parameters
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --opt)
            OPTIMUM_VALUE="$2"
            shift; shift
            ;;
        --timeout)
            TOTAL_LIMIT="$2"
            shift; shift
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        -*|--*)
            echo "Unknown option $1"
            print_help
            exit 1
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

# Get instance path
if [ ${#POSITIONAL[@]} -lt 1 ]; then
    echo "Error: Missing instance file."
    print_help
    exit 1
fi
INSTANCE_PATH="${POSITIONAL[0]}"

# Set default timeout if none was provided
if [ -z "${TOTAL_LIMIT:-}" ]; then
    echo "No timeout specified. Using default timeout: ${DEFAULT_TIMEOUT}s (1 hour)"
    TOTAL_LIMIT="$DEFAULT_TIMEOUT"
fi

# Check if RoundingSat binary exists
if [ ! -f "$ROUNDINGSAT_EXEC" ]; then
    echo "Error: RoundingSat not found at $ROUNDINGSAT_EXEC"
    exit 1
fi

# RoundingSat options
OPTIONS=(
    "--calculate-backbones=1"
    "--cg-encoding=reified"
    "--cg-resprop=1"
    "--lp-cut-gomory=0"
    "--verbosity=0"
)

START_TIME=$(date +%s)
echo "Starting backbone extraction for instance: $INSTANCE_PATH"

# If optimum is provided, use single-call flow
if [ -n "$OPTIMUM_VALUE" ]; then
    echo "Using provided optimum value: $OPTIMUM_VALUE"
    "$ROUNDINGSAT_EXEC" "$INSTANCE_PATH" "${OPTIONS[@]}" --already-solved=1 --optimum-value="$OPTIMUM_VALUE" --timeout="$TOTAL_LIMIT"
    exit 0
fi

# First call: detect type or extract optimum
echo "No optimum value provided. Running RoundingSat to detect instance type..."
OUTPUT=$("$ROUNDINGSAT_EXEC" "$INSTANCE_PATH" "${OPTIONS[@]}" --timeout="$TOTAL_LIMIT" 2>&1)
echo "$OUTPUT"

# Check if it's a decisional instance
if echo "$OUTPUT" | grep -q "s CALCULATED BACKBONE"; then
    echo "Backbone successfully extracted (instance is decisional)."
    exit 0
fi

# Check if optimum was found
if ! echo "$OUTPUT" | grep -q "s OPTIMUM FOUND"; then
    echo "Error: No optimum value found in output — likely due to timeout or infeasibility."
    exit 1
fi

# Try to extract optimum value (line before 's OPTIMUM FOUND')
OPTIMUM_LINE=$(echo "$OUTPUT" | grep -B1 "s OPTIMUM FOUND" | head -n 1)
OPTIMUM_VALUE=$(echo "$OPTIMUM_LINE" | awk '/^o / {print $2}')

if [ -z "$OPTIMUM_VALUE" ]; then
    echo "Error: Failed to extract optimum value from solver output (unexpected format)."
    exit 1
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
REMAINING_TIME=$((TOTAL_LIMIT - ELAPSED))

if [ "$REMAINING_TIME" -le 0 ]; then
    echo "Time budget exhausted before second call. Skipping backbone extraction."
    exit 1
fi

echo "Extracted optimum value: $OPTIMUM_VALUE"
echo "Running RoundingSat again to extract backbone (remaining time: ${REMAINING_TIME}s)..."
"$ROUNDINGSAT_EXEC" "$INSTANCE_PATH" "${OPTIONS[@]}" --already-solved=1 --optimum-value="$OPTIMUM_VALUE" --timeout="$REMAINING_TIME"

echo "Backbone extraction complete."
