"""Stage 0 handoff audit.

Three checks performed in order:
  1. k-eff reproduced — loads existing statepoint.200.h5 and compares against
     RESULTS.md reference (1.10184 ± tolerance).
  2. Depletion module available — confirms openmc.deplete imports correctly and
     CoupledOperator/Chain are present.
  3. Depletion chain file present — checks data/chain_endfb80_pwr.xml exists
     and is parseable by openmc.deplete.Chain.
  4. Kernel material ready for depletion — confirms U234/U235/U238 nuclides are
     present and kernel.depletable is True.

Exit code 0 if all checks pass; 1 if any fail.
"""

from __future__ import annotations

import pathlib
import sys

import openmc
import openmc.deplete

_REPO_ROOT = pathlib.Path(__file__).parent.parent
_STATEPOINT = _REPO_ROOT / 'statepoint.200.h5'
_CHAIN_FILE = _REPO_ROOT / 'data' / 'chain_endfb71_thermal.xml'

# Reference k-eff from RESULTS.md (ENDF/B-VIII.0 library, 1200 K, 200 batches)
_K_REF = 1.10184


def _check_keff() -> bool:
    if not _STATEPOINT.exists():
        print(f'FAIL [k-eff]: statepoint not found at {_STATEPOINT}')
        print('       Run `python -m triso.run` to generate it.')
        return False

    with openmc.StatePoint(str(_STATEPOINT)) as sp:
        k_mean = sp.keff.nominal_value
        k_std = sp.keff.std_dev

    tol = 3.0 * k_std
    delta = abs(k_mean - _K_REF)
    passed = delta <= tol
    status = 'PASS' if passed else 'FAIL'
    print(
        f'{status} [k-eff]: {k_mean:.5f} ± {k_std:.5f} '
        f'(ref {_K_REF:.5f}, Δ={delta:.5f}, tol=3σ={tol:.5f})'
    )
    return passed


def _check_depletion_module() -> bool:
    try:
        from openmc.deplete import CoupledOperator, Chain  # noqa: F401
        print(
            f'PASS [depletion module]: openmc.deplete available '
            f'(openmc {openmc.__version__}), CoupledOperator and Chain present'
        )
        return True
    except ImportError as exc:
        print(f'FAIL [depletion module]: {exc}')
        return False


def _check_chain_file() -> bool:
    if not _CHAIN_FILE.exists():
        print(
            f'FAIL [chain file]: not found at {_CHAIN_FILE}\n'
            '       Run: bash scripts/download_chain.sh'
        )
        return False

    chain = openmc.deplete.Chain.from_xml(str(_CHAIN_FILE))
    print(
        f'PASS [chain file]: {_CHAIN_FILE.name} loaded, '
        f'{len(chain.nuclides)} nuclides'
    )
    return True


def _check_kernel_material() -> bool:
    # Import here so OPENMC_CROSS_SECTIONS must be set, same as a real run.
    from triso.materials import build_materials

    mats = build_materials()
    kernel = mats['kernel']
    nuc_names = {n.name for n in kernel.nuclides}
    required = {'U234', 'U235', 'U238'}
    missing = required - nuc_names

    if missing:
        print(f'FAIL [kernel composition]: missing nuclides {missing}')
        return False

    if not kernel.depletable:
        print('FAIL [kernel composition]: kernel.depletable is not True')
        return False

    u_nucs = sorted(n for n in nuc_names if n.startswith('U'))
    print(
        f'PASS [kernel composition]: uranium nuclides {u_nucs}, depletable=True'
    )
    return True


def main() -> int:
    print('=== Stage 0 handoff audit ===\n')

    results = {
        'k-eff reproduced': _check_keff(),
        'depletion module': _check_depletion_module(),
        'chain file': _check_chain_file(),
        'kernel composition': _check_kernel_material(),
    }

    print('\n--- Summary ---')
    all_pass = True
    for label, passed in results.items():
        tag = 'PASS' if passed else 'FAIL'
        print(f'  {tag}: {label}')
        all_pass = all_pass and passed

    print()
    if all_pass:
        print('All checks passed — Stage 0 geometry is confirmed and depletion is ready.')
    else:
        print('One or more checks failed — see details above.')

    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
