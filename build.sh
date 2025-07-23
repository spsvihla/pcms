#!/bin/bash
set -euo pipefail

PACKAGE_NAME="pcms"
WHEEL_DIR="dist"
BUILD_TYPE="release"

usage() {
    echo "Usage: $0 [--build-type <debug|profile|release>]"
    echo "  --build-type   Build type: debug, profile, or release (default: release)"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-type)
            if [[ $# -lt 2 ]]; then
                echo "Error: --build-type requires an argument"
                usage
            fi
            BUILD_TYPE="$2"
            shift 2
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
BUILD_CMD="python setup.py bdist_wheel"

if [[ "$BUILD_TYPE" != "release" ]]; then
    BUILD_CMD+=" --build-type=$BUILD_TYPE"
fi

echo "Running: $BUILD_CMD"
eval "$BUILD_CMD"

echo "Installing the package from wheel..."
WHEEL_FILE=$(ls "$WHEEL_DIR"/"${PACKAGE_NAME}"-*.whl | head -n 1)
pip install "$WHEEL_FILE" --force-reinstall

echo "Done installing '$PACKAGE_NAME' from: $WHEEL_FILE"
