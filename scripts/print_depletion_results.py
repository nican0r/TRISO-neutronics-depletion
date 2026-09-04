"""Print k-eff vs. time from a completed depletion run.

Run from the directory containing depletion_results.h5:
    cd output/depletion
    python ../../scripts/print_depletion_results.py
"""
import openmc.deplete

r = openmc.deplete.Results('depletion_results.h5')
t, k = r.get_keff()

print(f'{"Step":>5}  {"Time [d]":>10}  {"k-eff":>10}  {"sigma":>8}')
print('-' * 42)
for i, (ti, ki) in enumerate(zip(t / 86400, k)):
    print(f'{i:>5}  {ti:>10.1f}  {ki[0]:>10.5f}  {ki[1]:>8.5f}')
