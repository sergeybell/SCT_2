@echo off
cd /d %~dp0
python exact_hamiltonian_check.py > exact_hamiltonian_check.out
python make_notebook.py
jupyter nbconvert --to notebook --execute --inplace analysis.ipynb
