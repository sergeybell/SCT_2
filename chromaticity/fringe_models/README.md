# fringe_models - краевые поля Bell против COSY / MAD-X / OptiM

Решётка: `..\BELL\magnetic_2_md.bell` (та же, что для `chrom_bell.jl`).
Bell подключается как библиотека: `Project.toml` содержит `Bell` с
`[sources] path = "C:\\GIT_REPOS\\BELL"`.

Запуск: `run_fringe.bat` (= `julia --project=. run.jl`). Первый раз:
`julia --project=. -e "using Pkg; Pkg.instantiate()"`.

Результаты: `output\fringe_vs_codes.txt`, `output\xi_y_vs_gap.png`.

Сравнение:
1. Bell `NoEdge` против COSY FR 0; Bell `HardEdge(order=2)` против MAD-X
   (жёсткий край, E1=E2=0, FINT=0) и предела COSY FR 3 при GAP -> 0.
2. Bell `SoftEdge(Enge)` с коэффициентами Enge по умолчанию COSY против
   COSY FR 3 при GAP = 50 / 25 / 12.5 / 6.25 мм. В COSY GAP - полуапертура,
   масштаб Enge - полный зазор, поэтому в Bell `Enge(gap = 2*GAP)`.
3. Жёсткий край с FINT профиля Enge против COSY FR 3.

Опорные числа скопированы из `..\COSY\output\cosy_ov321_fr.txt`,
`..\analysis\exact_hamiltonian_check.out`, `..\MADX`.
Пример в самом репозитории BELL (`examples\19_fringe_models`) показывает
только модели Bell и их внутреннюю сходимость, без сравнения с другими кодами.
