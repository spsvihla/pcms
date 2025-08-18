#!/bin/bash
set -euo pipefail

PACKAGE_NAME="pcms"
WHEEL_DIR="dist"
BUILD_TYPE="Release"   # Default: Release

usage() {
    echo "Usage: $0 [--build-type <debug|profile|release>] [--clean]"
    echo "  --build-type   Build type: debug, profile, or release (default: release)"
    echo "  --clean        Clean build artifacts and exit"
    exit 1
}

CLEAN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-type)
            if [[ $# -lt 2 ]]; then
                echo "Error: --build-type requires an argument"
                usage
            fi
            case "$2" in
                debug) BUILD_TYPE="Debug" ;;
                profile) BUILD_TYPE="RelWithDebInfo" ;; # map "profile" to CMake
                release) BUILD_TYPE="Release" ;;
                *) echo "Unknown build type: $2" ; usage ;;
            esac
            shift 2
            ;;
        --clean)
            CLEAN=true
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

if $CLEAN; then
    echo "Cleaning build artifacts..."
    rm -rf build/ "$WHEEL_DIR"/ "${PACKAGE_NAME}"*.egg-info UNKNOWN.egg-info
    echo "Clean complete."
    exit 0
fi

echo "Cleaning previous build artifacts..."
rm -rf build/ "$WHEEL_DIR"/ "${PACKAGE_NAME}"*.egg-info UNKNOWN.egg-info

echo "Uninstalling previous installation of '$PACKAGE_NAME' (if installed)..."
if pip show "$PACKAGE_NAME" &> /dev/null; then
    pip uninstall -y "$PACKAGE_NAME"
else
    echo "Package '$PACKAGE_NAME' not currently installed."
fi

echo "Building wheel distribution with CMake build type: $BUILD_TYPE"

# Use scikit-build / pip wheel with cmake.build-type
pip wheel . -w "$WHEEL_DIR" --config-settings=cmake.build-type="$BUILD_TYPE"

echo "Installing the package from wheel..."
WHEEL_FILE=$(ls "$WHEEL_DIR"/"${PACKAGE_NAME}"-*.whl | head -n 1)
pip install "$WHEEL_FILE" --force-reinstall

echo "Done installing '$PACKAGE_NAME' from: $WHEEL_FILE"
