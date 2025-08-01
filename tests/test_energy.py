import os
import math

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "logs")
TOTAL_ENERGY_PATH = os.path.join(DATA_DIR, "total_energy_distrib.txt")

with open(TOTAL_ENERGY_PATH) as f:
    energies = [float(line.strip()) for line in f if line.strip()]

assert len(energies) >= 51, "Expected at least 51 energy records"
start, end = energies[0], energies[-1]
assert math.isfinite(start) and math.isfinite(end)
relative_change = abs(end - start) / (abs(start) + 1e-12)
# Energy should remain roughly stable for a short test run
assert relative_change < 1.0, f"Energy drift too high: {relative_change}"
print("Energy log test passed")
