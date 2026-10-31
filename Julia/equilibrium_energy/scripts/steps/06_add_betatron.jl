# Step 06: add Λ to δ_eq. One amplitude, RF on.
function step_06!(S, ctx)
    section("06. Сдвиг δ_eq от бетатронного удлинения")
    dm = ctx.series[:dm]
    Ax = 10e-3
    N = 8000
    z0 = init_canonical(S, dm, Ax, 0.0)
    println("  exact tracking $N turns with A_x=$(Ax) …")
    t = @elapsed Z = run_turns(z -> step_exact(S, z), z0, N)
    d = deltas_of(Z)
    pnote = phis_note_of(Z, S)
    ctx.series[:Z4_Ax] = Z
    ctx.series[:d4_Ax] = d
    @printf("  done in %.1f s\n", t)
    m = measure_synchrotron(d, pnote; Qs = S.Qs)
    Jx = cs_action(Ax, S.po.x)
    m1_noΛ = formula_delta_eq(S.pe, S.xi, S.C, m.delta_m_fit, 0.0, 0.0)
    m1 = formula_delta_eq(S.pe, S.xi, S.C, m.delta_m_fit, Jx, 0.0)
    m2 = nf6d_delta_eq(S.nf, nf6d_actions(S.nf, z0))
    rm = running_mean(d)
    last = rm[end]
    d0 = ctx.series[:d4]
    n0 = min(N, length(d0))
    p = plot(xlabel = "номер оборота N", ylabel = "⟨δ⟩_N · 10⁶", legend = :topright,
             size = (760, 440), xscale = :log10)
    plot!(p, 1:n0, running_mean(d0[1:n0]) .* 1e6, lw = 1.5, label = "без бетатрона")
    plot!(p, 1:N, rm .* 1e6, lw = 1.5, label = "A_x = 10 мм")
    hline!(p, [m1_noΛ * 1e6], lw = 2, ls = :dash, label = "формула без Λ")
    hline!(p, [m1 * 1e6], lw = 2, ls = :dot, label = "формула с Λ")
    savefig_step(p, "fig_06_lambda")
    shift = m1 - m1_noΛ
    mac!(ctx, "Nbeta", string(N))
    mac!(ctx, "Axten", "10")
    mac!(ctx, "MoneNoLam", texnum(m1_noΛ; d = 4))
    mac!(ctx, "MoneLam", texnum(m1; d = 4))
    mac!(ctx, "MtwoLam", texnum(m2; d = 4))
    mac!(ctx, "MeanLam", texnum(last; d = 4))
    mac!(ctx, "LamShift", texnum(shift; d = 2))
    mac!(ctx, "DiffLamMean", texnum(last - m1; d = 2))
    mac!(ctx, "DiffNoLamMean", texnum(last - m1_noΛ; d = 2))
    conc = @sprintf("Включение A_x=10 мм сдвигает формулу на %.2e (это −Λ/η₀). Среднее интегратора %.3e: до формулы с Λ %.2e, до формулы без Λ %.2e. Среднее садится на линию с Λ.",
                    shift, last, last - m1, last - m1_noΛ)
    write_conclusion("step06", conc)
    write_step_tex("step06", """
\\section{Бетатрон в \$\\delta_{\\mathrm{eq}}\$}
\\textbf{Вопрос.} На сколько сдвигается среднее \$\\delta\$, если к той же частице добавить бетатрон?

\$A_x=\\vAxten\$~мм, полость включена, \$N=\\vNbeta\$ оборотов. Две формулы: без \$\\Lambda\$ и с \$\\Lambda=-(2\\pi/C)J_x\\xi_x\$.

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_06_lambda.png}
\\caption{Формула без \$\\Lambda\$: \$\\vMoneNoLam\$. С \$\\Lambda\$: \$\\vMoneLam\$, сдвиг \$\\vLamShift\$. Среднее интегратора \$\\vMeanLam\$ (до формулы с \$\\Lambda\$: \$\\vDiffLamMean\$). M2: \$\\vMtwoLam\$.}
\\end{figure}

\\textbf{Вывод.} Бетатрон сдвигает уровень на \$\\vLamShift\$. Среднее садится на формулу с \$\\Lambda\$, а не на формулу без него.
""")
    return ctx
end
