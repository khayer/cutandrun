#!/bin/bash
# run_local_report.sh - Local test runner for CUT_and_RUN report generator

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
TEST_STAGING="${PROJECT_DIR}/tests/report_smoke/staging"
TEST_OUTPUT="${PROJECT_DIR}/tests/report_smoke/CUT_and_RUN_report.html"

# Install dependencies if not already installed
echo "Installing dependencies..."
pip install -q -r "${SCRIPT_DIR}/requirements.txt"

# Run the report generator
echo "Generating report from staging directory: $TEST_STAGING"
python3 "${SCRIPT_DIR}/generator.py" \
    --staging-dir "$TEST_STAGING" \
    --out "$TEST_OUTPUT"

echo "Report generated successfully: $TEST_OUTPUT"
echo "To view the report, open it in a web browser:"
echo "  open '$TEST_OUTPUT'  # macOS"
echo "  xdg-open '$TEST_OUTPUT'  # Linux"
