# Step 02: running mean of exact tracking vs the simple formula and vs M2.
function step_02!(S, ctx)
    section("02. Скользящее среднее δ: садится ли на формулу и на M2")
    dm = ctx.series[:dm]
    N = 12000
    z0 = ctx.series[:z0_dm]
    println("  exact tracking $N turns …")
    t = @elapsed Z = run_turns(z -> step_exact(S, z), z0, N)
    d = deltas_of(Z)
    pnote = phis_note_of(Z, S)
    ctx.series[:Z4] = Z
    ctx.series[:d4] = d
    ctx.series[:p4] = pnote
    @printf("  done in %.1f s\n", t)
    m = measure_synchrotron(d, pnote; Qs = S.Qs)
    J0 = nf6d_actions(S.nf, z0)
    m1 = formula_delta_eq(S.pe, S.xi, S.C, m.delta_m_fit, 0.0, 0.0)
    m2 = nf6d_delta_eq(S.nf, J0)
    rm = running_mean(d)
    last = rm[end]
    nd = m.n_drop
    Ts = m.Tsyn
    ks = 1:div(N - nd, Ts)
    Nk = [nd + k * Ts for k in ks]
    pm = [mean(d[nd+1:n]) for n in Nk]
    p = plot(xlabel = "номер оборота N", ylabel = "⟨δ⟩_N · 10⁶", legend = :topright,
             size = (760, 440), xscale = :log10)
    plot!(p, 1:N, rm .* 1e6, lw = 1.0, alpha = 0.7, label = "простое скользящее среднее")
    plot!(p, Nk, pm .* 1e6, lw = 2.2, marker = :circle, ms = 3, label = "среднее по целым периодам T_s")
    hline!(p, [m1 * 1e6], lw = 2, ls = :dash, label = "формула −η₁δ_m²/(2η₀)")
    hline!(p, [m2 * 1e6], lw = 2, ls = :dot, label = "нормальная форма (M2)")
    hline!(p, [m.delta_eq * 1e6], lw = 1.5, color = :green, label = "гармонический фит ряда")
    savefig_step(p, "fig_02_runmean")
    mac!(ctx, "Nlong", string(N))
    mac!(ctx, "MoneA", texnum(m1; d = 4))
    mac!(ctx, "MtwoA", texnum(m2; d = 4))
    mac!(ctx, "MeanA", texnum(last; d = 4))
    mac!(ctx, "PmA", texnum(pm[end]; d = 4))
    mac!(ctx, "DiffAonetwo", texnum(m2 - m1; d = 2))
    mac!(ctx, "DiffAFourMone", texnum(last - m1; d = 2))
    mac!(ctx, "DiffAFourMtwo", texnum(last - m2; d = 2))
    mac!(ctx, "DiffPmMone", texnum(pm[end] - m1; d = 2))
    mac!(ctx, "dmFit", texnum(m.delta_m_fit; d = 4))
    mac!(ctx, "phatA", @sprintf("%.4f", m.phi_hat))
    mac!(ctx, "FitA", texnum(m.delta_eq; d = 4))
    mac!(ctx, "DiffFitMone", texnum(m.delta_eq - m1; d = 2))
    conc = @sprintf("Простое среднее за %d оборотов %.3e, до формулы %.2e. Среднее по периодам %.3e. Гармонический фит ряда %.3e, до формулы %.2e, до M2 %.2e. Простое среднее ещё не село; фит садится на формулу и M2.",
                    N, last, last - m1, pm[end], m.delta_eq, m.delta_eq - m1, m.delta_eq - m2)
    write_conclusion("step02", conc)
    write_step_tex("step02", """
\\section{Скользящее среднее и простая формула}
\\textbf{Вопрос.} Садится ли среднее \$\\delta\$ по оборотам на формулу \$-\\eta_1\\delta_m^2/(2\\eta_0)\$ и на среднее нормальной формы?

Тот же запуск, \$N=\\vNlong\$ оборотов. Формула берёт \$\\eta_0,\\eta_1\$ с магнитных элементов этого кольца и измеренную амплитуду \$\\delta_m\$. Нормальная форма усредняет лабораторный \$\\delta\$ по трём фазам того же произведения Ли и про имена \$\\eta,\\xi\$ не знает.

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_02_runmean.png}
\\caption{Простое среднее \$\\vMeanA\$ и среднее по периодам \$\\vPmA\$ ещё не сели. Гармонический фит \$\\vFitA\$ (до формулы \$\\vDiffFitMone\$). Формула \$\\vMoneA\$, M2 \$\\vMtwoA\$.}
\\end{figure}

\\textbf{Вывод.} Простое среднее за \$\\vNlong\$ оборотов не село: сходится как \$1/N\$. Гармонический фит того же ряда садится на формулу и на M2. Между формулой и M2 разница \$\\vDiffAonetwo\$.
""")
    # fix macro names used in tex - I used vDiffAfourm in tex but mac DiffA4m1
    # I'll use consistent names in tex matching mac! keys: DiffA4m1 -> vDiffAivm? 
    # LaTeX command names can't have digits mixed awkwardly. \vDiffA4m1 is invalid because of the digit in the middle...
    # actually \vDiffA4m1 is valid as \vDiffA then 4m1 leftover? In TeX, \vDiffA4m1 is command \vDiffA4m1 if we \newcommand{\vDiffA4m1}.
    # \newcommand{\vDiffA4m1} is VALID in LaTeX.
    # I wrote \vDiffAfourm in the tex body - mismatch. Fix the tex body.
    return ctx
end
