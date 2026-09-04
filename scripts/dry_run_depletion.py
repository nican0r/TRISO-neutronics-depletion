"""One-step depletion dry run to validate the full pipeline before the real run.

Run from the output directory you intend to use for the full depletion:
    mkdir -p output/depletion && cd output/depletion
    python ../../scripts/dry_run_depletion.py
"""
import sys
from pathlib import Path

import openmc
import openmc.deplete

sys.path.insert(0, str(Path(__file__).parents[1] / 'src'))
from triso.depletion import build_depletion_model, POWER_DENSITY

CHAIN = Path(__file__).parents[1] / 'data' / 'chain_endfb71_thermal.xml'

model = build_depletion_model()
model.settings.batches = 12
model.settings.inactive = 2
model.settings.particles = 100

op = openmc.deplete.CoupledOperator(
    model,
    chain_file=str(CHAIN),
    normalization_mode='fission-q',
)

integrator = openmc.deplete.CECMIntegrator(
    op, [1.0], power_density=POWER_DENSITY, timestep_units='d'
)
integrator.integrate()

r = openmc.deplete.Results('depletion_results.h5')
t, k = r.get_keff()
print(f'\nDry run OK')
print(f'  k-eff at t=0:  {k[0].n:.5f} +/- {k[0].s:.5f}')
print(f'  k-eff at t=1d: {k[1].n:.5f} +/- {k[1].s:.5f}')
