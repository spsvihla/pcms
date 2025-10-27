#!/usr/bin/env bash
set -euo pipefail

PACKAGE_NAME="pcms"
WHEEL_DIR="dist"
BUILD_TYPE="Release"   # Default build type
CLEAN=false
EXTRAS=""

usage() {
    echo "Usage: $0 [--build-type <debug|profile|release>] [--clean] [--extras <extras>]"
    echo "  --build-type   Build type: debug, profile, or release (default: release)"
    echo "  --clean        Clean build artifacts and exit"
    echo "  --extras       Install extras (comma-separated) after wheel installation"
    exit 1
}

# --- Parse arguments ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-type)
            if [[ $# -lt 2 ]]; then usage; fi
            case "$2" in
                debug) BUILD_TYPE="Debug" ;;
                profile) BUILD_TYPE="RelWithDebInfo" ;;
                release) BUILD_TYPE="Release" ;;
                *) echo "Unknown build type: $2"; usage ;;
            esac
            shift 2
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        --extras)
            if [[ $# -lt 2 ]]; then usage; fi
            EXTRAS="$2"
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

ENV_FILE=".env"
LINE='PYTHONPATH=.venv/lib/python3.10/site-packages:$PYTHONPATH'

# Check if .env exists
if [ ! -f "$ENV_FILE" ]; then
    echo "Creating $ENV_FILE..."
    echo "$LINE" > "$ENV_FILE"
    echo ".env created with PYTHONPATH entry."
else
    # Check if the line already exists
    if ! grep -qxF "$LINE" "$ENV_FILE"; then
        echo "Appending PYTHONPATH to $ENV_FILE..."
        echo "$LINE" >> "$ENV_FILE"
    else
        echo "$ENV_FILE already contains the PYTHONPATH entry."
    fi
fi

# --- Clean build artifacts and exit ---
if $CLEAN; then
    echo "Cleaning build artifacts..."
    rm -rf build/ "$WHEEL_DIR"/ "${PACKAGE_NAME}"*.egg-info UNKNOWN.egg-info
    echo "Clean complete."
    exit 0
fi

# --- Clean previous builds ---
echo "Cleaning previous build artifacts..."
rm -rf build/ "$WHEEL_DIR"/ "${PACKAGE_NAME}"*.egg-info UNKNOWN.egg-info

# --- Uninstall previous installation ---
echo "Uninstalling previous installation of '$PACKAGE_NAME' (if installed)..."
if pip show "$PACKAGE_NAME" &> /dev/null; then
    pip uninstall -y "$PACKAGE_NAME"
else
    echo "Package '$PACKAGE_NAME' not currently installed."
fi

# --- Build wheel ---
echo "Building wheel distribution with CMake build type: $BUILD_TYPE"
pip wheel . -w "$WHEEL_DIR" --config-settings=cmake.build-type="$BUILD_TYPE" --verbose

# --- Install wheel ---
echo "Installing the package from wheel..."
WHEEL_FILE=$(ls "$WHEEL_DIR"/"${PACKAGE_NAME}"-*.whl | head -n 1)
pip install "$WHEEL_FILE" --force-reinstall
echo "Installed '$PACKAGE_NAME' from: $WHEEL_FILE"

# --- Install extras from pyproject.toml ---
if [[ -n "$EXTRAS" ]]; then
    echo "Installing extras: [$EXTRAS]"
    pip install ".[${EXTRAS}]"
fi

echo "Done."
