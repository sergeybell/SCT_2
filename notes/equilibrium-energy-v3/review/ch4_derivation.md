# Вывод формул главы 4

Обозначения совпадают с запиской. Индекс \(s\) у \(\gamma_s\), \(\beta_s\), \(E_s\), \(p_s\), \(v_s\) пишется явно. В этой главе \(y=0\).

Кольцевое среднее \(\langle f\rangle_C=(1/C_s)\oint f\,\mathrm{d}s\) берётся вдоль опорной орбиты. Среднее без индекса \(\langle f\rangle\) --- среднее по оборотам за синхротронный период, как в \eqref{eq:deq-def} записки. В постановке знаменатель среднего записан как \(C\); в записке длина опорной орбиты есть \(C_s\).

Удерживаемый порядок: \(\delta^2\) в \(\tau\) и \(u^2\) в \(v_s/v\). Члены \(O(\delta^3)\) и \(O(u^3)\) отбрасываются.

## 1. Потенциал и уравнение Лапласа

Опорная орбита лежит в плоскости, радиус кривизны \(\rho\), \(x\) --- по главной нормали наружу, \(y\) --- по бинормали. Масштабный множитель вдоль орбиты \(h_s=1+x/\rho\), \(h_x=h_y=1\).

Поле дефлектора не зависит от \(s\). Уравнение Лапласа в этих координатах

\[
\frac{\partial^2\Phi}{\partial x^2}+\frac{\partial^2\Phi}{\partial y^2}+\frac{1}{\rho+x}\frac{\partial\Phi}{\partial x}=0.
\]

Нормировка \(\psi=q\Phi/(p_s v_s)\), \(\xi=x/\rho\), \(\zeta=y/\rho\). Поскольку \(p_s v_s=\beta_s^2 E_s\), та же функция \(\psi\) входит в энергию. После перехода к \(\xi,\zeta\) и умножения на \(\rho^2\)

\[
\psi_{\xi\xi}+\psi_{\zeta\zeta}+\frac{1}{1+\xi}\psi_\xi=0.
\]

Ряд до второго порядка, показатель поля \(n\) свободен:

\[
\psi=\xi+a\xi^2+b\zeta^2+O(3).
\]

Тогда \(\psi_\xi=1+2a\xi+O(2)\), \(\psi_{\xi\xi}=2a\), \(\psi_{\zeta\zeta}=2b\), и \(1/(1+\xi)=1+O(\xi)\). В нулевом порядке

\[
2a+2b+1=0.
\]

Поле на опорной орбите. Принимается \(E_x=\partial\Phi/\partial x\), так что

\[
E_0=\frac{\partial\Phi}{\partial x}\Big|_0=\frac{p_s v_s}{q\rho},\qquad qE_0=\frac{p_s v_s}{\rho}:
\]

линейный член \(\psi\) равен \(\xi\), а не \(-\xi\). Показатель поля

\[
n=-\frac{\rho}{E_0}\frac{\partial E_x}{\partial x}\Big|_0.
\]

\[
\frac{\partial E_x}{\partial x}\Big|_0=\frac{p_s v_s}{q\rho^2}\psi_{\xi\xi}=\frac{p_s v_s}{q\rho^2}\,2a,
\]

\[
\frac{\rho}{E_0}\frac{\partial E_x}{\partial x}=2a,
\]

поэтому \(n=-2a\) и \(a=-n/2\). Из \(2a+2b+1=0\)

\[
2b=n-1,\qquad b=\frac{n-1}{2}.
\]

Итог:

\[
\psi=\xi-\frac{n}{2}\xi^2+\frac{n-1}{2}\zeta^2.
\]

Кубические члены потенциала в \(\eta_0\) и \(\eta_1\) не входят: на орбите \(x=O(\delta)\) они дают \(O(\delta^3)\). Вне дефлекторов \(\Phi=0\), то есть \(\psi=0\), и \(1/\rho=0\).

## 2. Локальный \(\gamma\) и \(v_s/v\)

Полная энергия \(E=\gamma mc^2+q\Phi\) сохраняется в статическом поле и меняется только в зазоре. На опорной орбите \(\psi=0\), поэтому \(E_s=\gamma_s mc^2\). По определению \(w\)

\[
E=E_s+w\beta_s^2 E_s,\qquad q\Phi=\psi\,p_s v_s=\psi\beta_s^2 E_s,
\]

\[
E-q\Phi=E_s\bigl[1+\beta_s^2(w-\psi)\bigr].
\]

Слева \(\gamma mc^2\), значит точно

\[
\gamma=\gamma_s\bigl[1+\beta_s^2 u\bigr],\qquad u=w-\psi.
\]

Дальше индекс \(s\) у \(\beta_s\), \(\gamma_s\) в промежуточных строках опущен: \(\beta^2=1-1/\gamma^2\), \(\gamma\) в правой части --- синхронное значение.

\[
\frac{1}{\gamma_{\mathrm{loc}}^2}=\frac{1}{\gamma^2(1+\beta^2 u)^2}=\frac{1}{\gamma^2}\bigl(1-2\beta^2 u+3\beta^4 u^2+O(u^3)\bigr),
\]

\[
\beta_{\mathrm{loc}}^2=1-\frac{1}{\gamma_{\mathrm{loc}}^2}
=\beta^2+\frac{2\beta^2}{\gamma^2}u-\frac{3\beta^4}{\gamma^2}u^2+O(u^3).
\]

\[
\frac{\beta_{\mathrm{loc}}}{\beta}
=\sqrt{1+\frac{2u}{\gamma^2}-\frac{3\beta^2}{\gamma^2}u^2+O(u^3)}.
\]

Ряд \(\sqrt{1+z}=1+z/2-z^2/8+O(z^3)\) при \(z=2u/\gamma^2+O(u^2)\):

\[
\frac{\beta_{\mathrm{loc}}}{\beta}
=1+\frac{u}{\gamma^2}-\frac{3\beta^2}{2\gamma^2}u^2-\frac{u^2}{2\gamma^4}+O(u^3).
\]

Обращение \(1/(1+\varepsilon)=1-\varepsilon+\varepsilon^2+O(\varepsilon^3)\), \(\varepsilon=u/\gamma^2+O(u^2)\), \(\varepsilon^2=u^2/\gamma^4\):

\[
\frac{v_s}{v}=\frac{\beta}{\beta_{\mathrm{loc}}}
=1-\frac{u}{\gamma^2}+\Biggl(\frac{3\beta^2}{2\gamma^2}+\frac{1}{2\gamma^4}+\frac{1}{\gamma^4}\Biggr)u^2+O(u^3).
\]

Подстановка \(\beta^2=1-1/\gamma^2\):

\[
\frac{3\beta^2}{2\gamma^2}=\frac{3}{2\gamma^2}-\frac{3}{2\gamma^4}.
\]

Члены \(1/\gamma^4\) сокращаются: \(-3/(2\gamma^4)+1/(2\gamma^4)+1/\gamma^4=0\). Остаётся

\[
\frac{v_s}{v}=1-\frac{u}{\gamma^2}+\frac{3u^2}{2\gamma^2}+O(u^3).
\]

## 3. Подынтегральное выражение до \(\delta^2\)

### 3.1. Элемент пути

При \(y=0\)

\[
\frac{\mathrm{d}l}{\mathrm{d}s}=\sqrt{\Bigl(1+\frac{x}{\rho}\Bigr)^2+x'^2}
=\sqrt{1+z},
\qquad
z=2\frac{x}{\rho}+\Bigl(\frac{x}{\rho}\Bigr)^2+x'^2.
\]

\[
\sqrt{1+z}=1+\frac{z}{2}-\frac{z^2}{8}+O(z^3),
\qquad
\frac{z}{2}=\frac{x}{\rho}+\frac12\Bigl(\frac{x}{\rho}\Bigr)^2+\frac{x'^2}{2},
\qquad
\frac{z^2}{8}=\frac12\Bigl(\frac{x}{\rho}\Bigr)^2+O(3).
\]

Квадрат \((x/\rho)^2\) сокращается:

\[
\frac{\mathrm{d}l}{\mathrm{d}s}=1+\frac{x}{\rho}+\frac{x'^2}{2}+O(3).
\]

Опорная частица имеет \(\psi=0\) и \(v=v_s\), поэтому \(T_s=C_s/v_s\) и

\[
1+\tau=\frac{1}{C_s}\oint\Bigl(1+\frac{x}{\rho}+\frac{x'^2}{2}\Bigr)\frac{v_s}{v}\,\mathrm{d}s.
\]

### 3.2. \(u\) на замкнутой орбите

Вне дефлекторов \(\psi=0\), и \(w\) связано с \(\delta\) рядом записки

\[
w=\delta+\frac{\delta^2}{2\gamma^2}+O(\delta^3).
\]

Между зазорами \(w\) не зависит от \(s\). На орбите \(x=D\delta+D_2\delta^2\), \(y=0\),

\[
\xi=\frac{D}{\rho}\delta+\frac{D_2}{\rho}\delta^2,
\qquad
\xi^2=\Bigl(\frac{D}{\rho}\Bigr)^2\delta^2+O(\delta^3),
\]

\[
\psi=\xi-\frac{n}{2}\xi^2
=\frac{D}{\rho}\delta+\Biggl[\frac{D_2}{\rho}-\frac{n}{2}\Bigl(\frac{D}{\rho}\Bigr)^2\Biggr]\delta^2+O(\delta^3).
\]

Обозначения на время выкладки: \(d_1=D/\rho\), \(d_2=D_2/\rho\), \(p_2=D'^2\). Тогда

\[
u=w-\psi
=(1-d_1)\delta+\Biggl[\frac{1}{2\gamma^2}-d_2+\frac{n}{2}d_1^2\Biggr]\delta^2+O(\delta^3),
\]

\[
u^2=(1-d_1)^2\delta^2+O(\delta^3).
\]

Вне дефлектора \(d_1=d_2=0\) и \(\psi=0\), формула для \(u\) сохраняется.

### 3.3. Скорость и путь

\[
\frac{v_s}{v}
=1-\frac{(1-d_1)\delta}{\gamma^2}
+\Biggl[
-\frac{1}{2\gamma^4}+\frac{d_2}{\gamma^2}-\frac{n\,d_1^2}{2\gamma^2}+\frac{3(1-d_1)^2}{2\gamma^2}
\Biggr]\delta^2+O(\delta^3).
\]

\[
1+\frac{x}{\rho}+\frac{x'^2}{2}
=1+d_1\delta+\Bigl(d_2+\frac{p_2}{2}\Bigr)\delta^2+O(\delta^3):
\]

член \(D_2'\) входит в \(x'\) множителем \(\delta^2\) и в \(x'^2\) даёт \(O(\delta^3)\).

Произведение \(PQ\) при \(P=1+A\delta+B\delta^2\), \(Q=1+C\delta+F\delta^2\) имеет коэффициенты \(A+C\) и \(B+F+AC\).

\[
A=d_1,
\qquad
C=-\frac{1-d_1}{\gamma^2},
\qquad
A+C=d_1\Bigl(1+\frac{1}{\gamma^2}\Bigr)-\frac{1}{\gamma^2}.
\]

\[
B=d_2+\frac{p_2}{2},
\]

\[
F=-\frac{1}{2\gamma^4}+\frac{d_2}{\gamma^2}-\frac{n\,d_1^2}{2\gamma^2}+\frac{3}{2\gamma^2}-\frac{3d_1}{\gamma^2}+\frac{3d_1^2}{2\gamma^2},
\]

\[
AC=-\frac{d_1(1-d_1)}{\gamma^2}=-\frac{d_1}{\gamma^2}+\frac{d_1^2}{\gamma^2}.
\]

Сумма \(B+F+AC\).

- \(p_2/2\) остаётся.
- \(d_2\) входит как \(d_2+d_2/\gamma^2=d_2(1+1/\gamma^2)\).
- При \(d_1^2\): \(-n/(2\gamma^2)+3/(2\gamma^2)+1/\gamma^2=(5-n)/(2\gamma^2)\).
- При \(d_1\): \(-3/\gamma^2-1/\gamma^2=-4/\gamma^2\).
- Свободный член: \(-1/(2\gamma^4)+3/(2\gamma^2)=(3\gamma^2-1)/(2\gamma^4)\).

До усреднения подынтегральное выражение равно

\[
\begin{aligned}
&1+\Biggl[\Bigl(1+\frac{1}{\gamma^2}\Bigr)d_1-\frac{1}{\gamma^2}\Biggr]\delta
\\
&+\Biggl[
\frac{p_2}{2}
+\Bigl(1+\frac{1}{\gamma^2}\Bigr)d_2
+\frac{5-n}{2\gamma^2}d_1^2
-\frac{4d_1}{\gamma^2}
+\frac{3\gamma^2-1}{2\gamma^4}
\Biggr]\delta^2.
\end{aligned}
\]

### 3.4. Среднее по кольцу

Постоянная \(1\) даёт \(1\). Величина \(1/\gamma^2\) не зависит от \(s\). Поэтому

\[
\eta_0=\alpha_0\Bigl(1+\frac{1}{\gamma^2}\Bigr)-\frac{1}{\gamma^2},
\qquad
\alpha_0=\Biggl\langle\frac{D}{\rho}\Biggr\rangle_C.
\]

Тот же \(\alpha_0\) стоит в длине орбиты: \((1/C_s)\oint(1+x/\rho+x'^2/2)\,\mathrm{d}s=1+\alpha_0\delta+O(\delta^2)\), линейный член есть \(\langle D/\rho\rangle_C\).

Свободный член \((3\gamma^2-1)/(2\gamma^4)\) постоянен, его среднее равно ему самому. Линейный по \(d_1\) член даёт \(-4\alpha_0/\gamma^2\). Остальные средние не сворачиваются:

\[
\eta_1
=\frac12\langle D'^2\rangle_C
+\Bigl(1+\frac{1}{\gamma^2}\Bigr)\Biggl\langle\frac{D_2}{\rho}\Biggr\rangle_C
+\frac{5-n}{2\gamma^2}\Biggl\langle\Bigl(\frac{D}{\rho}\Bigr)^2\Biggr\rangle_C
-\frac{4\alpha_0}{\gamma^2}
+\frac{3\gamma^2-1}{2\gamma^4}.
\]

Переход к \(w\). Из записки \(\delta=w-w^2/(2\gamma^2)+O(w^3)\),

\[
\tau=\eta_0\Bigl(w-\frac{w^2}{2\gamma^2}\Bigr)+\eta_1 w^2+O(w^3)
=\eta_0 w+\Bigl(\eta_1-\frac{\eta_0}{2\gamma^2}\Bigr)w^2.
\]

Это \eqref{eq:tau-w} записки: \(\tilde\eta_1=\eta_1-\eta_0/(2\gamma^2)\). Подстановка использует только ряд \(w(\delta)\), не вид \(\eta_0\) и \(\eta_1\).

## 4. Кольцо из одного дефлектора

Все сечение --- дефлектор, замкнутая орбита круговая: \(D'=0\), величины не зависят от \(s\), кольцевое среднее совпадает со значением. Радиус опорной орбиты \(\rho=1\), \(m=c=q=1\). Импульс \(\delta\) определён полной энергией через ряд \(w=\delta+\delta^2/(2\gamma^2)\), то есть как импульс при том же \(E\) и \(\psi=0\). Сила на круговой орбите радиуса \(r\):

\[
pv=\gamma_{\mathrm{loc}}-\frac{1}{\gamma_{\mathrm{loc}}}=k\,r^{1-n},
\qquad
k=\beta^2\gamma,
\]

поскольку \(qE_r=k r^{-n}\) и \(qE_r=pv/r\), а на опорной орбите \(k=p_s v_s\).

### 4.1. \(n=1\)

При \(n=1\) произведение \(pv=k\) на всех круговых орбитах, локальный \(\gamma\) равен \(\gamma\), весь избыток энергии сидит в потенциале \(\psi=\ln r\). Скорость равна \(v_s\), поэтому \(\tau=r-1\). Из \(w=\ln r\)

\[
r=e^w=1+w+\frac{w^2}{2}+O(w^3),
\qquad
\tau=w+\frac{w^2}{2}+O(w^3).
\]

\[
w=\delta+\frac{\delta^2}{2\gamma^2},
\qquad
\tau=\delta+\frac{\gamma^2+1}{2\gamma^2}\delta^2+O(\delta^3).
\]

Отсюда \(\eta_0=1\), \(\eta_1=(\gamma^2+1)/(2\gamma^2)\). Смещение

\[
\xi=r-1=w+\frac{w^2}{2}
=\delta+\frac{\gamma^2+1}{2\gamma^2}\delta^2,
\]

то есть \(D/\rho=1\), \(D_2/\rho=(\gamma^2+1)/(2\gamma^2)\). Подстановка в общую формулу, \(g=\gamma^2\):

\[
\Bigl(1+\frac{1}{g}\Bigr)\frac{g+1}{2g}=\frac{(g+1)^2}{2g^2},
\qquad
\frac{5-1}{2g}=\frac{2}{g},
\]

\[
\eta_1=\frac{(g+1)^2}{2g^2}+\frac{2}{g}-\frac{4}{g}+\frac{3g-1}{2g^2}
=\frac{(g+1)^2+(3g-1)-4g}{2g^2}
=\frac{g^2+g}{2g^2}
=\frac{g+1}{2g}.
\]

Общая формула воспроизводит \(\eta_1\).

### 4.2. \(n=2\)

\[
\gamma_{\mathrm{loc}}-\frac{1}{\gamma_{\mathrm{loc}}}=\frac{k}{r},
\qquad
\Phi=k\Bigl(1-\frac{1}{r}\Bigr),
\qquad
k=\beta^2\gamma.
\]

Пусть \(\gamma_{\mathrm{loc}}=\gamma(1+\varepsilon)\), \(r=1+\xi\). Равенство силы до \(\xi^2\):

\[
\varepsilon=\varepsilon_1\xi+\varepsilon_2\xi^2,
\qquad
\varepsilon_1=-\frac{\beta^2 g}{g+1},
\qquad
\varepsilon_2=\frac{\beta^2 g}{g+1}+\frac{\beta^4 g^2}{(g+1)^3},
\]

\(g=\gamma^2\). Полная энергия даёт

\[
w=\frac{\varepsilon}{\beta^2}+\xi-\xi^2+O(\xi^3)
=\frac{\xi}{g+1}+\Biggl[-\frac{1}{g+1}+\frac{\beta^2 g^2}{(g+1)^3}\Biggr]\xi^2+O(\xi^3).
\]

Обращение:

\[
\xi=(g+1)w+(3g+1)w^2+O(w^3).
\]

С \(w=\delta+\delta^2/(2g)\)

\[
\frac{D}{\rho}=g+1,
\qquad
\frac{D_2}{\rho}=\frac{g+1}{2g}+3g+1=\frac{6g^2+3g+1}{2g}.
\]

Линейный член уже достаточен для \(\eta_0\):

\[
\eta_0=(g+1)\Bigl(1+\frac{1}{g}\Bigr)-\frac{1}{g}=g+2.
\]

Подстановка \(D_2\) в общую формулу. Общий знаменатель \(2g^2\):

\[
\Bigl(1+\frac{1}{g}\Bigr)\frac{D_2}{\rho}=\frac{(g+1)(6g^2+3g+1)}{2g^2},
\qquad
\frac{3}{2g}(g+1)^2=\frac{3g(g+1)^2}{2g^2},
\]

\[
-\frac{4(g+1)}{g}=\frac{-8g(g+1)}{2g^2},
\qquad
\frac{3g-1}{2g^2}=\frac{3g-1}{2g^2}.
\]

Числитель:

\[
(g+1)(6g^2+3g+1)=6g^3+9g^2+4g+1,
\]

\[
3g(g+1)^2=3g^3+6g^2+3g,
\qquad
-8g(g+1)=-8g^2-8g.
\]

Сумма с \(3g-1\):

\[
9g^3+7g^2+2g=g(9g^2+7g+2).
\]

\[
\eta_1=\frac{9g^2+7g+2}{2g}=\frac{9\gamma^4+7\gamma^2+2}{2\gamma^2}.
\]

### 4.3. Рациональная точка

При \(\gamma=5/4\) замкнутые формы дают \(n=1\): \(\eta_0=1\), \(\eta_1=41/50\); \(n=2\): \(\eta_0=57/16\), \(\eta_1=8937/800\). Период точной круговой орбиты (\(pv=k r^{1-n}\), ряд до \(\delta^2\)) и подстановка полученных \(D\), \(D_2\) в общую формулу совпадают с этими дробями; та же сверка при \(n=3/2\).
