#!/bin/bash
set -euo pipefail

PACKAGE_NAME="pcms"
WHEEL_DIR="dist"

echo "Cleaning previous build artifacts..."
rm -rf build/ "$WHEEL_DIR"/ "${PACKAGE_NAME}"*.egg-info

echo "Uninstalling previous installation of '$PACKAGE_NAME' (if installed)..."
if pip show "$PACKAGE_NAME" &> /dev/null; then
    pip uninstall -y "$PACKAGE_NAME"
else
    echo "Package '$PACKAGE_NAME' not currently installed."
fi

echo "Building wheel distribution using setup.py..."
python setup.py bdist_wheel

echo "Installing the package from wheel..."
WHEEL_FILE=$(ls "$WHEEL_DIR"/"${PACKAGE_NAME}"-*.whl | head -n 1)
pip install "$WHEEL_FILE" --force-reinstall

echo "Done installing '$PACKAGE_NAME' from: $WHEEL_FILE"
