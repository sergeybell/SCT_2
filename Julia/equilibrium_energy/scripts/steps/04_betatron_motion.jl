# Step 04: betatron motion on the same ring, RF on.
function step_04!(S, ctx)
    section("04. Бетатронное движение при включённом ВЧ")
    dm = ctx.series[:dm]
    Ax = 20e-3
    z0 = init_canonical(S, dm, Ax, 0.0)
    ctx.series[:z0_Ax] = z0
    ctx.series[:Ax] = Ax
    samp = sample_along(S, z0; nturns = 2)
    Z = run_turns(z -> step_exact(S, z), z0, 250)
    ux = [canonical_to_user(Z[:, n], S.ref, S.C, S.h)[1] for n in 1:size(Z, 2)]
    uxp = [canonical_to_user(Z[:, n], S.ref, S.C, S.h)[2] for n in 1:size(Z, 2)]
    scale = 200.0
    R = S.C / 2pi
    φ = samp.s ./ R
    X = (R .+ scale .* samp.x) .* cos.(φ)
    Y = (R .+ scale .* samp.x) .* sin.(φ)
    circx = R .* cos.(range(0, 2pi; length = 400))
    circy = R .* sin.(range(0, 2pi; length = 400))
    p1 = plot(xlabel = "X, м", ylabel = "Y, м", legend = :topright, size = (640, 640),
              aspect_ratio = :equal)
    plot!(p1, circx, circy, lw = 2, color = :black, label = "замкнутая орбита")
    plot!(p1, X, Y, lw = 1.2, label = "траектория, x увеличен в $(Int(scale)) раз")
    savefig_step(p1, "fig_04_ring")
    p2 = plot(xlabel = "x, мм", ylabel = "x′, мрад", legend = :topright, size = (560, 520),
              aspect_ratio = :equal)
    scatter!(p2, ux .* 1e3, uxp .* 1e3, ms = 3, msw = 0, label = "оборот за оборотом")
    savefig_step(p2, "fig_04_ellipse")
    Jx = cs_action(Ax, S.po.x)
    mac!(ctx, "Axmm", @sprintf("%.0f", Ax * 1e3))
    mac!(ctx, "xScale", @sprintf("%.0f", scale))
    mac!(ctx, "JxCS", texnum(Jx; d = 3))
    mac!(ctx, "betax", @sprintf("%.4f", S.po.x.beta))
    conc = @sprintf("При включённом ВЧ частица ходит по эллипсу (x, x′) с амплитудой %d мм. На плане кольца поперечник увеличен в %.0f раз, иначе %.0f мм на фоне периметра %.1f м не видны.",
                    Ax * 1e3, scale, Ax * 1e3, S.C)
    write_conclusion("step04", conc)
    write_step_tex("step04", """
\\section{Бетатронное движение}
\\textbf{Вопрос.} Как выглядит поперечное движение на том же кольце, где полость включена?

Амплитуда \$A_x=\\vAxmm\$~мм, \$\\delta_m=\\vdm\$. Действие Куранта--Снайдера \$J_x=\\vJxCS\$~м, \$\\beta_x=\\vbetax\$~м в точке наблюдения.

\\begin{figure}[H]\\centering
\\includegraphics[width=0.72\\textwidth]{figs/fig_04_ring.png}
\\caption{Два оборота на плане кольца. Поперечное смещение увеличено в \$\\vxScale\$ раз: иначе миллиметры на фоне \$C=\\vCirc\$~м не видны. Полоса вокруг орбиты --- след бетатрона.}
\\end{figure}

\\begin{figure}[H]\\centering
\\includegraphics[width=0.62\\textwidth]{figs/fig_04_ellipse.png}
\\caption{Эллипс \$(x,x')\$ оборот за оборотом. ВЧ включён: эллипс медленно дышит из-за синхротрона через дисперсию.}
\\end{figure}

\\textbf{Вывод.} Бетатрон есть поверх синхротрона. Это то движение, которое ниже удлиняет путь.
""")
    return ctx
end
