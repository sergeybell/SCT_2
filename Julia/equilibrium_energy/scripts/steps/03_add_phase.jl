# Step 03: add the φ̂²/32 factor to the formula; same series as step 02.
function step_03!(S, ctx)
    section("03. Добавление фазовой поправки 1 − φ̂²/32")
    d = ctx.series[:d4]
    pnote = ctx.series[:p4]
    z0 = ctx.series[:z0_dm]
    m = measure_synchrotron(d, pnote; Qs = S.Qs)
    m1 = formula_delta_eq(S.pe, S.xi, S.C, m.delta_m_fit, 0.0, 0.0)
    m1p = formula_delta_eq(S.pe, S.xi, S.C, m.delta_m_fit, 0.0, 0.0; phat = m.phi_hat)
    m2 = nf6d_delta_eq(S.nf, nf6d_actions(S.nf, z0))
    rm = running_mean(d)
    last = rm[end]
    p = plot(xlabel = "номер оборота N", ylabel = "⟨δ⟩_N · 10⁶", legend = :topright,
             size = (760, 440), xscale = :log10)
    plot!(p, 1:length(d), rm .* 1e6, lw = 1.6, label = "скользящее среднее, точный интегратор")
    hline!(p, [m1 * 1e6], lw = 2, ls = :dash, label = "формула без фазы")
    hline!(p, [m1p * 1e6], lw = 2, ls = :dot, label = "формула с 1 − φ̂²/32")
    savefig_step(p, "fig_03_phase")
    shift = m1p - m1
    mac!(ctx, "MonePhi", texnum(m1p; d = 4))
    mac!(ctx, "PhaseShift", texnum(shift; d = 2))
    mac!(ctx, "phatA", @sprintf("%.4f", m.phi_hat))
    mac!(ctx, "DiffPhiMean", texnum(last - m1p; d = 2))
    mac!(ctx, "DiffNoPhiMean", texnum(last - m1; d = 2))
    conc = @sprintf("Множитель 1−φ̂²/32 при φ̂=%.4f сдвигает формулу на %.2e. Остаток среднего к формуле с фазой %.2e, без фазы %.2e. Сдвиг меньше оставшихся колебаний среднего на этом N.",
                    m.phi_hat, shift, last - m1p, last - m1)
    write_conclusion("step03", conc)
    write_step_tex("step03", """
\\section{Фазовая поправка}
\\textbf{Вопрос.} Сдвигает ли множитель \$1-\\hat\\varphi^2/32\$ предсказанный уровень так, что среднее садится на новую линию?

Тот же ряд. Амплитуда фазы \$\\hat\\varphi=\\vphatA\$. Этот множитель --- следующий порядок (генератор степени 5); в M2 его нет.

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_03_phase.png}
\\caption{Формула без фазы \$\\vMoneA\$, с фазой \$\\vMonePhi\$, сдвиг \$\\vPhaseShift\$. Среднее интегратора к формуле с фазой: \$\\vDiffPhiMean\$.}
\\end{figure}

\\textbf{Вывод.} Сдвиг уровня \$\\vPhaseShift\$. На этих оборотах он сравни́м с остатком скользящего среднего, отдельной посадки на новую линию не видно.
""")
    return ctx
end
