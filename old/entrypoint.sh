#!/bin/bash --login
# The --login ensures the bash configuration is loaded,
# enabling Conda.

# Enable strict mode.
set -euo pipefail
# ... Run whatever commands ...

# Temporarily disable strict mode and activate conda:
set +euo pipefail
conda activate myenv

# Re-enable strict mode:
set -euo pipefail

# exec the final command:
#exec conda run --no-capture-output -n streamlit run test.py
#exec streamlit run test.py
exec streamlit run gui.py --server.port=8501 --server.address=0.0.0.0