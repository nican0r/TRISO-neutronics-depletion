"""
Validates that OpenMC is importable and cross-section data is accessible.
Run after setup_env.sh and download_data.sh.
Usage: python scripts/validate_env.py
"""
import sys
import os

def check(label, fn):
    try:
        result = fn()
        print(f"  [OK] {label}{': ' + result if result else ''}")
        return True
    except Exception as e:
        print(f"  [FAIL] {label}: {e}")
        return False

print(f"Python: {sys.version}")
print()
print("Checking imports ...")

ok = True
ok &= check("import openmc", lambda: __import__("openmc") and None)
ok &= check("import numpy", lambda: __import__("numpy").__version__)
ok &= check("import h5py", lambda: __import__("h5py").__version__)
ok &= check("import matplotlib", lambda: __import__("matplotlib").__version__)

print()
print("Checking nuclear data ...")

xs = os.environ.get("OPENMC_CROSS_SECTIONS")
if not xs:
    print("  [FAIL] OPENMC_CROSS_SECTIONS not set")
    ok = False
elif not os.path.isfile(xs):
    print(f"  [FAIL] OPENMC_CROSS_SECTIONS points to missing file: {xs}")
    ok = False
else:
    print(f"  [OK] OPENMC_CROSS_SECTIONS: {xs}")
    try:
        import openmc
        lib = openmc.data.DataLibrary.from_xml(xs)
        nuclides = [e["materials"] for e in lib.libraries]
        flat = [n for sub in nuclides for n in (sub if isinstance(sub, list) else [sub])]
        print(f"  [OK] Library loaded — {len(lib.libraries)} datasets")
    except Exception as e:
        print(f"  [FAIL] Could not load library: {e}")
        ok = False

print()
if ok:
    print("Environment OK — ready to build the TRISO model.")
else:
    print("Environment has errors — see above.")
    sys.exit(1)
