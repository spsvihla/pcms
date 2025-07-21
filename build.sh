#!/bin/bash
set -euo pipefail

PACKAGE_NAME="pcms"
WHEEL_DIR="dist"
DEBUG_BUILD=0

usage() {
    echo "Usage: $0 [-d|--debug] [-r|--rebuild]"
    echo "  -d, --debug    Build with debug flags"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--debug)
            DEBUG_BUILD=1
            shift
            ;;
        -*|--*)
            echo "Unknown option $1"
            usage
            ;;
        *)
            break
            ;;
    esac
done

echo "Cleaning previous build artifacts..."
rm -rf build/ "$WHEEL_DIR"/ "${PACKAGE_NAME}"*.egg-info

echo "Uninstalling previous installation of '$PACKAGE_NAME' (if installed)..."
if pip show "$PACKAGE_NAME" &> /dev/null; then
    pip uninstall -y "$PACKAGE_NAME"
else
    echo "Package '$PACKAGE_NAME' not currently installed."
fi

echo "Building wheel distribution..."
if [[ $DEBUG_BUILD -eq 1 ]]; then
    python setup.py bdist_wheel --debug
else
    python setup.py bdist_wheel
fi

echo "Installing the package from wheel..."
WHEEL_FILE=$(ls "$WHEEL_DIR"/"${PACKAGE_NAME}"-*.whl | head -n 1)
pip install "$WHEEL_FILE" --force-reinstall

echo "Done installing '$PACKAGE_NAME' from: $WHEEL_FILE"
