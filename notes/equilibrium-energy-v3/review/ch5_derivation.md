# Вывод формул главы 5

Обозначения совпадают с запиской. \(\langle\,\rangle_\psi\) --- среднее по бетатронным фазам. Скобки без индекса в \(\langle\psi\rangle\) и \(\langle\gamma\rangle\) --- среднее по синхротронным и бетатронным колебаниям в фиксированной точке \(s\). Удерживается второй порядок по амплитуде: \(J\), \(A^2\) и \(\delta_\eq\). Члены третьего порядка отбрасываются.

## 1. Подынтегральное выражение

При \(w=0\) выполнено \(u=w-\psi=-\psi\). Подстановка в ряд \(v_s/v\) даёт

\[
\frac{v_s}{v}=1+\frac{\psi}{\gamma_s^2}+\frac{3\psi^2}{2\gamma_s^2}.
\]

Элемент пути \(\mathrm{d}l/\mathrm{d}s=1+x/\rho+(x'^2+y'^2)/2\). Произведение

\[
I=\Bigl(1+\frac{x}{\rho}+q\Bigr)
\Bigl(1+\frac{\psi}{\gamma_s^2}+\frac{3\psi^2}{2\gamma_s^2}\Bigr),
\qquad
q=\frac{x'^2+y'^2}{2}.
\]

\(x_\beta,y_\beta,x',y'=O(a)\), \(J\sim a^2\), \(x_{\mathrm{nl}}=O(a^2)\), \(q=O(a^2)\). До порядка \(a^2\)

\[
I=1+\frac{x}{\rho}+q+\frac{\psi}{\gamma_s^2}+\frac{1}{\gamma_s^2}\frac{x}{\rho}\,\psi+\frac{3\psi^2}{2\gamma_s^2}.
\]

Произведения \(q\psi\) и \((x/\rho)\psi^2\) начинаются с \(a^3\).

\[
\psi=\frac{x}{\rho}-\frac{n}{2}\Bigl(\frac{x}{\rho}\Bigr)^2+\frac{n-1}{2}\Bigl(\frac{y}{\rho}\Bigr)^2,
\]

\[
x=x_\beta+x_{\mathrm{nl}},\qquad y=y_\beta.
\]

До порядка \(a^2\)

\[
\psi=\frac{x_\beta}{\rho}+\frac{x_{\mathrm{nl}}}{\rho}-\frac{n}{2}\Bigl(\frac{x_\beta}{\rho}\Bigr)^2+\frac{n-1}{2}\Bigl(\frac{y_\beta}{\rho}\Bigr)^2,
\]

\[
\psi^2=\Bigl(\frac{x_\beta}{\rho}\Bigr)^2,
\qquad
\frac{x}{\rho}\,\psi=\Bigl(\frac{x_\beta}{\rho}\Bigr)^2.
\]

Член \((x/\rho)\cdot\bigl(-(n/2)(x/\rho)^2\bigr)\) есть \(O(a^3)\) и в произведение не входит. Квадрат части \(\psi\), содержащей \(y^2\), есть \(O(a^4)\).

Вклад в \(I-1\), квадратичный по \(x_\beta\), коэффициент при \(x_\beta^2/\rho^2\):

- из \(\psi/\gamma_s^2\): \(-n/(2\gamma_s^2)\);
- из \((x/\rho)\psi/\gamma_s^2\): \(1/\gamma_s^2\);
- из \(3\psi^2/(2\gamma_s^2)\): \(3/(2\gamma_s^2)\).

Сумма

\[
\frac{1}{\gamma_s^2}\Bigl(1-\frac{n}{2}+\frac{3}{2}\Bigr)=\frac{5-n}{2\gamma_s^2}.
\]

Вклад \(y_\beta^2/\rho^2\) входит только из \(\psi/\gamma_s^2\):

\[
\frac{n-1}{2\gamma_s^2}.
\]

Линейный по \(x_{\mathrm{nl}}\) вклад: \(x_{\mathrm{nl}}/\rho\) из элемента пути и \(x_{\mathrm{nl}}/(\rho\gamma_s^2)\) из \(\psi/\gamma_s^2\), вместе \((1+1/\gamma_s^2)x_{\mathrm{nl}}/\rho\).

При \(J=0\) и \(w=0\) траектория есть опорная орбита, \(I=1\). Линейные по \(x_\beta\) члены имеют нулевое фазовое среднее. Поэтому

\[
\Lambda=\frac{1}{C_s}\oint\Biggl[
\Bigl(1+\frac{1}{\gamma_s^2}\Bigr)\frac{\langle x_{\mathrm{nl}}\rangle_\psi}{\rho}
+\frac{\langle x'^2+y'^2\rangle_\psi}{2}
+\frac{5-n}{2\gamma_s^2}\frac{\langle x_\beta^2\rangle_\psi}{\rho^2}
+\frac{n-1}{2\gamma_s^2}\frac{\langle y_\beta^2\rangle_\psi}{\rho^2}
\Biggr]\mathrm{d}s.
\]

\(\langle x_\beta'\rangle_\psi=0\) и \(x_{\mathrm{nl}}'^2=O(a^4)\), так что \(\langle x'^2\rangle_\psi=\langle x_\beta'^2\rangle_\psi\). То же для \(y\). С \(\langle x_\beta'^2\rangle_\psi=J_x\gamma_x\), \(\langle y_\beta'^2\rangle_\psi=J_y\gamma_y\), \(\langle x_\beta^2\rangle_\psi=J_x\beta_x\), \(\langle y_\beta^2\rangle_\psi=J_y\beta_y\)

\[
\frac{\langle x'^2+y'^2\rangle_\psi}{2}=\frac{J_x\gamma_x+J_y\gamma_y}{2},
\]

\[
\frac{5-n}{2\gamma_s^2}\frac{\langle x_\beta^2\rangle_\psi}{\rho^2}
=\frac{J_x}{2}\,\frac{(5-n)\beta_x}{\gamma_s^2\rho^2},
\]

\[
\frac{n-1}{2\gamma_s^2}\frac{\langle y_\beta^2\rangle_\psi}{\rho^2}
=\frac{J_y}{2}\,\frac{(n-1)\beta_y}{\gamma_s^2\rho^2}.
\]

Вне дефлекторов \(1/\rho=0\): слагаемые с \(1/\rho\) и \(1/\rho^2\) там отсутствуют, интегралы \(\gamma_x\) и \(\gamma_y\) берутся по всему кольцу.

## 2. Среднее \(\psi\) в дефлекторе

\[
x=D\delta+D_2\delta^2+x_\beta+x_{\mathrm{nl}},\qquad y=y_\beta.
\]

\[
\Bigl\langle\frac{x}{\rho}\Bigr\rangle
=\frac{D}{\rho}\delta_\eq+\frac{D_2}{\rho}\frac{A^2}{2}+\frac{\langle x_{\mathrm{nl}}\rangle_\psi}{\rho},
\]

поскольку \(\langle x_\beta\rangle=0\), \(\langle\delta\rangle=\delta_\eq\), \(\langle\delta^2\rangle=A^2/2\).

\[
\Bigl\langle\Bigl(\frac{x}{\rho}\Bigr)^2\Bigr\rangle
=\Bigl(\frac{D}{\rho}\Bigr)^2\frac{A^2}{2}+\frac{J_x\beta_x}{\rho^2}.
\]

Перекрёстный член \(2(D/\rho)\delta\,(x_\beta/\rho)\) имеет нулевое среднее: фазы \(\delta\) и \(x_\beta\) независимы. Произведения \(x_{\mathrm{nl}}\) с \(\delta\) и с \(x_\beta\) имеют порядок выше второго. Квадрат \(y\):

\[
\Bigl\langle\Bigl(\frac{y}{\rho}\Bigr)^2\Bigr\rangle=\frac{J_y\beta_y}{\rho^2}.
\]

Подстановка в \(\psi\) даёт

\[
\begin{aligned}
\langle\psi\rangle
&=\frac{D}{\rho}\delta_\eq
+\Biggl[\frac{D_2}{\rho}-\frac{n}{2}\Bigl(\frac{D}{\rho}\Bigr)^2\Biggr]\frac{A^2}{2}
+\frac{\langle x_{\mathrm{nl}}\rangle_\psi}{\rho}
\\
&\quad-\frac{n}{2}\frac{J_x\beta_x}{\rho^2}
+\frac{n-1}{2}\frac{J_y\beta_y}{\rho^2}.
\end{aligned}
\]

Из \(\gamma=\gamma_s[1+\beta_s^2(w-\psi)]\) точно \((\gamma-\gamma_s)/(\beta_s^2\gamma_s)=w-\psi\). Кинетическая энергия есть \(\gamma mc^2\), а \(E_s=\gamma_s mc^2\), поэтому её сдвиг равен \(\beta_s^2 E_s(w-\psi)\). В среднем

\[
\Delta K_\eq(s)=\beta_s^2 E_s\bigl[w_\eq-\langle\psi\rangle(s)\bigr].
\]

Вне дефлекторов \(\psi=0\), и \(\Delta K_\eq=\beta_s^2 E_s w_\eq=\Delta E_\eq\).
