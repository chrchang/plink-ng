# ---------- System dependencies ----------
echo "Installing system dependencies..."
apt-get update && apt-get install -y \
    liblapack-dev \
    libopenblas-dev \
    liblapacke-dev \
    build-essential

# Verify lapacke.h location
echo "Searching for lapacke.h..."
find /usr -name "lapacke.h" 2>/dev/null || echo "lapacke.h not found!"


# ---------- Python environment ----------
echo "Setting up Python environment..."

# Install python3-venv if not already installed
apt-get update && apt-get install -y python3-venv

