"""
Script para rodar os testes dentro do container.
Uso: docker compose run --rm fenicsx python3 /app/tests/run_tests.py
"""

import subprocess
import sys

def main():
    print("=== Testes FeniCS-EDSS ===")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "/app/tests/test_templates.py", "-v", "--tb=short"],
        cwd="/app"
    )
    sys.exit(result.returncode)

if __name__ == "__main__":
    main()
