#!/usr/bin/env python3
"""Remove unwanted temperature points from an OpenMC HDF5 cross-section library.

Each nuclide/thermal-scattering HDF5 file stores temperature-dependent data in
groups named like '250K', '294K', '900K', etc. under 'energy/', 'kTs/',
'reactions/reaction_NNN/', and 'urr/'.  Simple h5py deletion only removes the
link; the data bytes remain and the file does not shrink.  This script rewrites
each file from scratch, copying only the desired temperature groups, which
produces a genuinely smaller library.

'0K' groups (windowed-multipole / zero-Kelvin cross sections used for on-the-fly
Doppler broadening) are always preserved regardless of the keep list.

Usage
-----
  python scripts/trim_temperatures.py <library_dir> <temp_K> [<temp_K> ...]

Example — keep only 900 K and 1200 K:
  python scripts/trim_temperatures.py data/endfb80_hdf5 900 1200
"""

from __future__ import annotations

import argparse
import re
import shutil
import tempfile
from pathlib import Path

import h5py

_TEMP_PAT = re.compile(r'^(\d+)K$')


def _is_temp_label(name: str) -> bool:
    return bool(_TEMP_PAT.match(name))


def _copy_filtered(
    src: h5py.Group,
    dst: h5py.Group,
    keep: set[str],
) -> None:
    """Recursively copy src → dst, skipping unwanted temperature groups/datasets."""
    for k, v in src.attrs.items():
        dst.attrs[k] = v

    for name, item in src.items():
        if _is_temp_label(name) and name != '0K' and name not in keep:
            continue  # drop this temperature

        if isinstance(item, h5py.Group):
            child = dst.require_group(name)
            _copy_filtered(item, child, keep)
        else:
            src.copy(name, dst)


def trim_file(path: Path, keep: set[str]) -> tuple[int, int]:
    """Rewrite path keeping only the requested temperature groups.

    Returns (bytes_before, bytes_after).
    """
    before = path.stat().st_size

    tmp_fd, tmp_name = tempfile.mkstemp(dir=path.parent, suffix='.h5')
    tmp_path = Path(tmp_name)
    try:
        import os
        os.close(tmp_fd)
        with h5py.File(path, 'r') as fin, h5py.File(tmp_path, 'w') as fout:
            _copy_filtered(fin, fout, keep)
        path.unlink()
        shutil.move(tmp_path, path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    after = path.stat().st_size
    return before, after


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('library_dir', type=Path,
                        help='Directory containing the OpenMC HDF5 library.')
    parser.add_argument('temperatures', type=int, nargs='+',
                        help='Temperature(s) in Kelvin to keep (e.g. 900 1200).')
    args = parser.parse_args()

    lib_dir: Path = args.library_dir
    keep = {f'{t}K' for t in args.temperatures}

    if not lib_dir.is_dir():
        raise SystemExit(f'ERROR: {lib_dir} is not a directory.')

    h5_files = sorted(p for p in lib_dir.rglob('*.h5')
                      if p.name != 'cross_sections.h5')

    if not h5_files:
        raise SystemExit(f'No .h5 files found in {lib_dir}.')

    print(f'Library : {lib_dir}')
    print(f'Keeping : {sorted(keep)}  (0K windowed-multipole always kept)')
    print(f'Files   : {len(h5_files)}')
    print()

    total_before = total_after = 0

    for i, fpath in enumerate(h5_files, 1):
        before, after = trim_file(fpath, keep)
        total_before += before
        total_after += after
        saved_pct = 100 * (before - after) / before if before else 0
        print(f'  [{i:4d}/{len(h5_files)}]  {fpath.name:<40}  '
              f'{before / 1e6:7.1f} MB → {after / 1e6:7.1f} MB  '
              f'({saved_pct:.0f}% smaller)')

    saved_gb = (total_before - total_after) / 1e9
    print()
    print(f'Total before : {total_before / 1e9:.2f} GB')
    print(f'Total after  : {total_after  / 1e9:.2f} GB')
    print(f'Space saved  : {saved_gb:.2f} GB')


if __name__ == '__main__':
    main()
