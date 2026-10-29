@echo off
cd /d %~dp0
python exact_hamiltonian_check.py 1 > exact_hamiltonian_check.out
python make_notebook.py
jupyter nbconvert --to notebook --execute --inplace analysis.ipynb

