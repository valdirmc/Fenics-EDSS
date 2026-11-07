"""
Script para rodar os testes dentro do container.
Uso: docker compose run --rm fenicsx python3 /app/tests/run_tests.py
"""

import os
import subprocess
import sys


def main():
    print("=== Testes FeniCS-EDSS ===")
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(tests_dir)
    result = subprocess.run(
        [sys.executable, "-m", "pytest", tests_dir, "-v", "--tb=short"],
        cwd=repo,
    )
    sys.exit(result.returncode)

if __name__ == "__main__":
    main()
