#!/bin/bash
set -euo pipefail

PACKAGE_NAME="pcms"
WHEEL_DIR="dist"
DEBUG_BUILD=0
OPTIMIZE=0
OPT_LEVEL=""

usage() {
    echo "Usage: $0 [-d|--debug] [-O <level>] "
    echo "  -d, --debug                Build with debug flags"
    echo "  -o, --optimize             Build with optimization flags"
    echo "  -O, --opt-level <level>    Optimization level to pass to setup.py (e.g., -O2, -O3)"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--debug)
            DEBUG_BUILD=1
            shift
            ;;
        -o|--optimize)
            OPTIMIZE=1
            shift
            ;;
        -O|--opt-level)
            if [[ $# -lt 2 ]]; then
                echo "Error: --opt-level requires an argument"
                usage
            fi
            OPT_LEVEL="$2"
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
if [[ $DEBUG_BUILD -eq 1 ]]; then
    BUILD_CMD+=" --debug"
fi
if [[ -n "$OPT_LEVEL" ]]; then
    BUILD_CMD+=" --opt=$OPT_LEVEL"
fi
if [[ $OPTIMIZE -eq 1 ]]; then
    BUILD_CMD+=" --optimize"
fi

echo "Running: $BUILD_CMD"
eval "$BUILD_CMD"

echo "Installing the package from wheel..."
WHEEL_FILE=$(ls "$WHEEL_DIR"/"${PACKAGE_NAME}"-*.whl | head -n 1)
pip install "$WHEEL_FILE" --force-reinstall

echo "Done installing '$PACKAGE_NAME' from: $WHEEL_FILE"
