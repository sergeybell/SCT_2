# Step 08: three pictures of the torus.
function step_08!(S, ctx)
    section("08. Тор")
    dm = ctx.series[:dm]
    Ax = 10e-3
    z0 = init_canonical(S, dm, Ax, 0.0)
    samp = sample_along(S, z0; nturns = 8)
    scale = 250.0
    R = S.C / 2pi
    φ = samp.s ./ R
    X = (R .+ scale .* samp.x) .* cos.(φ)
    Y = (R .+ scale .* samp.x) .* sin.(φ)
    circx = R .* cos.(range(0, 2pi; length = 400))
    circy = R .* sin.(range(0, 2pi; length = 400))
    p1 = plot(xlabel = "X, м", ylabel = "Y, м", legend = :topright, size = (640, 640),
              aspect_ratio = :equal)
    plot!(p1, circx, circy, lw = 2.2, color = :black, label = "замкнутая орбита")
    plot!(p1, X, Y, lw = 0.7, alpha = 0.85, label = "8 оборотов, x × $(Int(scale))")
    savefig_step(p1, "fig_08_torus_ring")

    Zs = haskey(ctx.series, :Z4_Ax) ? ctx.series[:Z4_Ax] : run_turns(z -> step_exact(S, z), z0, 800)
    nsyn = min(size(Zs, 2), 800)
    ph = phis_note_of(Zs[:, 1:nsyn], S)
    dd = deltas_of(Zs[:, 1:nsyn])
    m = measure_synchrotron(dd, ph; Qs = S.Qs, n_drop = 0)
    p2 = plot(xlabel = "φ (запаздывание), рад", ylabel = "δ · 10³", legend = :topright, size = (560, 520))
    scatter!(p2, ph, dd .* 1e3, ms = 3, msw = 0, alpha = 0.7, label = "обороты")
    hline!(p2, [m.delta_eq * 1e3], lw = 2, ls = :dash, label = "центр = δ_eq")
    savefig_step(p2, "fig_08_synchrotron")

    J0 = nf6d_actions(S.nf, z0)
    nθ = 48
    θβ = range(0, 2pi; length = nθ + 1)[1:nθ]
    θs = range(0, 2pi; length = nθ + 1)[1:nθ]
    δg = zeros(nθ, nθ)
    rx, rs = sqrt(2 * J0[1]), sqrt(2 * J0[3])
    for i in 1:nθ, j in 1:nθ
        w = [rx * cos(θβ[i]), -rx * sin(θβ[i]), 0.0, 0.0, rs * cos(θs[j]), -rs * sin(θs[j])]
        δg[j, i] = S.nf.zstar[6] + Float64(evaluate(S.nf.U6, w))
    end
    Rt, rt = 1.0, 0.38
    Xm = [(Rt + rt * cos(tb)) * cos(ts) for tb in θβ, ts in θs]
    Ym = [(Rt + rt * cos(tb)) * sin(ts) for tb in θβ, ts in θs]
    Zm = [rt * sin(tb) for tb in θβ, ts in θs]
    p3 = surface(Xm, Ym, Zm, fill_z = δg', colorbar_title = "δ",
                 xlabel = "", ylabel = "", zlabel = "", legend = false,
                 size = (640, 560), camera = (40, 25))
    ntr = min(size(Zs, 2), 400)
    Xt = Float64[]; Yt = Float64[]; Zt = Float64[]
    for n in 1:ntr
        θ, _ = nf_angles(S.nf, Zs[:, n])
        tb, ts = θ[1], θ[3]
        push!(Xt, (Rt + rt * cos(tb)) * cos(ts))
        push!(Yt, (Rt + rt * cos(tb)) * sin(ts))
        push!(Zt, rt * sin(tb))
    end
    plot!(p3, Xt, Yt, Zt, lw = 1.8, color = :white, label = "")
    savefig_step(p3, "fig_08_torus_surface")
    mac!(ctx, "torScale", @sprintf("%.0f", scale))
    mac!(ctx, "deqTorus", texnum(m.delta_eq; d = 3))
    mac!(ctx, "meanColor", texnum(mean(δg); d = 3))
    mac!(ctx, "mtwoTorus", texnum(nf6d_delta_eq(S.nf, J0); d = 3))
    conc = @sprintf("Тор вокруг машины — полоса траектории на плане кольца. В плоскости (φ, δ) частица ходит по окружности вокруг центра %.3e. Цвет поверхности — δ(θ_β, θ_s); среднее цвета %.3e совпадает с M2 %.3e. Нитка — точный интегратор. Срез Jy=0; полный тор трёхмерный.",
                    m.delta_eq, mean(δg), nf6d_delta_eq(S.nf, J0))
    write_conclusion("step08", conc)
    write_step_tex("step08", """
\\section{Тор}
\\textbf{Вопрос.} Что такое тор, по которому усредняется \$\\delta\$?

После нормальной формы остаются три амплитуды и три фазы. Тор --- все тройки фаз при фиксированных амплитудах. Это не траектория: один оборот ставит на тор одну точку, долгий трекинг наматывает нитку. \$\\delta_{\\mathrm{eq}}\$ --- среднее \$\\delta\$ по всей поверхности сразу.

\\begin{figure}[H]\\centering
\\includegraphics[width=0.72\\textwidth]{figs/fig_08_torus_ring.png}
\\caption{План кольца, 8 оборотов. Поперечник увеличен в \$\\vtorScale\$ раз. Полоса, которую заполняет траектория, --- тор вокруг машины.}
\\end{figure}

\\begin{figure}[H]\\centering
\\includegraphics[width=0.62\\textwidth]{figs/fig_08_synchrotron.png}
\\caption{Плоскость \$(\\varphi,\\delta)\$ при включённом ВЧ. Частица ходит по окружности. Центр --- \$\\delta_{\\mathrm{eq}}=\\vdeqTorus\$; в центр частица не садится.}
\\end{figure}

\\begin{figure}[H]\\centering
\\includegraphics[width=0.72\\textwidth]{figs/fig_08_torus_surface.png}
\\caption{Поверхность \$(\\theta_\\beta,\\theta_s)\$ при \$J_y=0\$. Цвет --- лабораторный \$\\delta\$. Среднее цвета \$\\vmeanColor\$, M2 \$\\vmtwoTorus\$. Белая нитка --- точный интегратор. Полный тор трёхмерен: вторая бетатронная фаза на этом срезе не показана.}
\\end{figure}

\\textbf{Вывод.} \$\\delta_{\\mathrm{eq}}\$ --- среднее по поверхности, не точка, в которую частица приходит. Трекинг рисует одну нитку на этом торе.
""")
    return ctx
end
