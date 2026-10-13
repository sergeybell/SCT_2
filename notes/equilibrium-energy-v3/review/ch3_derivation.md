# Вывод формул главы 3 и кубического члена приложения A

Обозначения совпадают с запиской. Скобки без индекса --- среднее за синхротронный период. Удерживаются поправки порядка \(\hat\varphi^2\) в форме колебания и в \(\langle w^2\rangle\), а также уже известные члены порядка \(\tilde\eta_1\) и \(\tilde\eta_1^2\) из глав 1--2. Произведения \(\tilde\eta_1\hat\varphi^2\) в гармониках \(w(t)\) отбрасываются. Индекс \(s\) у \(\omega_s\), \(\eta_0\), \(\gamma_s\) пишется явно.

## 1. Система и уравнение для фазы

Второе уравнение сглаженной системы есть \(\mathrm{d}w/\mathrm{d}t=K\cos\varphi_s\sin\varphi\). Ряд

\[
\sin\varphi=\varphi-\frac{\varphi^3}{6}+O(\varphi^5)
\]

обрывается на кубическом члене. Производная этого приближения

\[
\frac{\mathrm{d}}{\mathrm{d}t}\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr)
=\Bigl(1-\frac{\varphi^2}{2}\Bigr)\dot\varphi
\]

совпадает с \(\cos\varphi\cdot\dot\varphi\) при

\[
\cos\varphi=1-\frac{\varphi^2}{2}+O(\varphi^4).
\]

При \(\Lambda=0\) система главы:

\[
\dot\varphi=\omega_{\mathrm{rf}}(\eta_0 w+\tilde\eta_1 w^2),
\qquad
\dot w=K\cos\varphi_s\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr).
\]

При \(\tilde\eta_1=0\) первое уравнение есть \(\dot\varphi=\omega_{\mathrm{rf}}\eta_0 w\). Дифференцирование и подстановка второго:

\[
\ddot\varphi
=\omega_{\mathrm{rf}}\eta_0\cdot K\cos\varphi_s\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr).
\]

Из второго равенства \eqref{eq:ws},

\[
K\omega_{\mathrm{rf}}\cos\varphi_s=-\frac{\omega_s^2}{\eta_0},
\]

умножение на \(\eta_0\) даёт

\[
\omega_{\mathrm{rf}}\eta_0 K\cos\varphi_s=-\omega_s^2.
\]

Знак верен при обоих знаках \(\eta_0\): \(\cos\varphi_s=-\operatorname{sgn}\eta_0\), поэтому

\[
\eta_0\cos\varphi_s=-|\eta_0|,
\qquad
\omega_{\mathrm{rf}}\eta_0 K\cos\varphi_s=-K\omega_{\mathrm{rf}}|\eta_0|=-\omega_s^2.
\]

Отсюда

\[
\ddot\varphi+\omega_s^2\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr)=0,
\]

или

\[
\ddot\varphi+\omega_s^2\varphi=\frac{\omega_s^2}{6}\varphi^3.
\]

Это уравнение Дуффинга со смягчающей силой: точный маятник \(\ddot\varphi+\omega_s^2\sin\varphi=0\) имеет ту же кубическую часть.

Полуразмах фазы в этой главе есть \(\hat\varphi=\max|\varphi|\), а не линейная оценка \(\omega_{\mathrm{rf}}|\eta_0|A/\omega_s\) из главы 1. Эти две величины совпадают в ведущем порядке и различаются на \(O(\hat\varphi^3)\).

## 2. Первый порядок метода Боголюбова--Митропольского

Уравнение записывается как

\[
\ddot x+\omega^2 x=\varepsilon c_3 x^3.
\]

Подстановка анзаца \eqref{eq:ansatz} в порядке \(\varepsilon\) даёт \eqref{eq:ordk}. Правая часть:

\[
\cos^3\psi=\frac{3\cos\psi+\cos 3\psi}{4},
\]

что следует из \(\cos 3\psi=4\cos^3\psi-3\cos\psi\). Поэтому

\[
f_1=c_3 a^3\frac{3\cos\psi+\cos 3\psi}{4}
=\frac{3c_3 a^3}{4}\cos\psi+\frac{c_3 a^3}{4}\cos 3\psi.
\]

Члена \(\sin\psi\) нет, проекция \eqref{eq:ordk} на него даёт \(A_1=0\).

Проекция на \(\cos\psi\). Функция \(u_1\) первой гармоники не содержит, оператор \(\partial_{\psi\psi}+1\) на ней равен нулю, и остаётся

\[
-2a\omega B_1=\frac{3c_3 a^3}{4},
\qquad
B_1=-\frac{3c_3 a^2}{8\omega}.
\]

Проекция на \(\cos 3\psi\). Коэффициент Фурье \(U_3\) в \(u_1=U_3\cos 3\psi\) удовлетворяет

\[
\omega^2(1-9)U_3=\frac{c_3 a^3}{4},
\qquad
-8\omega^2 U_3=\frac{c_3 a^3}{4},
\qquad
U_3=-\frac{c_3 a^3}{32\omega^2}.
\]

При \(\varepsilon=1\)

\[
x=a\cos\psi-\frac{c_3 a^3}{32\omega^2}\cos 3\psi,
\qquad
\Omega=\omega+B_1=\omega-\frac{3c_3 a^2}{8\omega}.
\]

Величина \(a\) есть коэффициент первой гармоники. В точке \(\psi=0\) оба косинуса равны единице:

\[
x(0)=a-\frac{c_3 a^3}{32\omega^2}.
\]

Вторая производная в этой точке отрицательна при \(c_3 a^2\ll\omega^2\), так что это максимум. В точке \(\psi=\pi\) значение равно \(-(x(0))\), и полуразмах модуля равен \(x(0)\).

## 3. Подстановка \(\varphi\), сдвиг на четверть периода

Для \eqref{eq:phi-ode}: \(x=\varphi\), \(\omega=\omega_s\), \(c_3=\omega_s^2/6\).

\[
\frac{3c_3 a^2}{8\omega}
=\frac{3(\omega_s^2/6)a^2}{8\omega_s}
=\frac{a^2\omega_s}{16},
\qquad
\frac{\Omega}{\omega_s}=1-\frac{a^2}{16}.
\]

\[
\frac{c_3 a^3}{32\omega^2}
=\frac{(\omega_s^2/6)a^3}{32\omega_s^2}
=\frac{a^3}{192}.
\]

\[
\varphi=a\cos\psi-\frac{a^3}{192}\cos 3\psi,
\qquad
\varphi_{\max}=a-\frac{a^3}{192}=a\Bigl(1-\frac{a^2}{192}\Bigr).
\]

Сдвиг начала времени \(\psi\mapsto\psi-\pi/2\). Дальше \(\psi\) --- новая фаза, частота \(\Omega\) не меняется.

\[
\cos(\psi-\pi/2)=\sin\psi.
\]

\[
\cos 3(\psi-\pi/2)
=\cos(3\psi-3\pi/2)
=\cos 3\psi\cdot\cos\frac{3\pi}{2}+\sin 3\psi\cdot\sin\frac{3\pi}{2}.
\]

\(\cos(3\pi/2)=0\), \(\sin(3\pi/2)=-1\), поэтому

\[
\cos 3(\psi-\pi/2)=-\sin 3\psi.
\]

Подстановка в косинусное решение:

\[
\varphi
=a\sin\psi-\frac{a^3}{192}(-\sin 3\psi)
=a\sin\psi+\frac{a^3}{192}\sin 3\psi.
\]

Минус перед \(\cos 3\psi\) и минус из сдвига дают плюс перед \(\sin 3\psi\).

В точке \(\psi=\pi/2\): \(\sin(\pi/2)=1\), \(\sin(3\pi/2)=-1\),

\[
\varphi(\pi/2)=a-\frac{a^3}{192}.
\]

В точке \(\psi=3\pi/2\): \(\sin(3\pi/2)=-1\), \(\sin(9\pi/2)=\sin(\pi/2)=1\),

\[
\varphi(3\pi/2)=-\Bigl(a-\frac{a^3}{192}\Bigr).
\]

Полуразмах \(\hat\varphi=\max|\varphi|\) равен прежнему максимуму:

\[
\hat\varphi=a\Bigl(1-\frac{a^2}{192}\Bigr).
\]

Обращение. Пусть \(a=\hat\varphi+\delta\) и \(\delta=O(\hat\varphi^3)\). Тогда \(a^3=\hat\varphi^3+O(\hat\varphi^5)\) и

\[
\hat\varphi=a-\frac{a^3}{192}
=\hat\varphi+\delta-\frac{\hat\varphi^3}{192}+O(\hat\varphi^5),
\]

откуда \(\delta=\hat\varphi^3/192\). До членов \(\hat\varphi^5\)

\[
a=\hat\varphi+\frac{\hat\varphi^3}{192},
\qquad
a^2=\hat\varphi^2+O(\hat\varphi^4),
\qquad
\frac{\Omega}{\omega_s}=1-\frac{\hat\varphi^2}{16}.
\]

Подстановка \(a\) в синусное решение:

\[
\frac{a^3}{192}=\frac{\hat\varphi^3}{192}+O(\hat\varphi^5),
\]

\[
\varphi
=\Bigl(\hat\varphi+\frac{\hat\varphi^3}{192}\Bigr)\sin\psi
+\frac{\hat\varphi^3}{192}\sin 3\psi
=\hat\varphi\sin\psi
+\frac{\hat\varphi^3}{192}(\sin\psi+\sin 3\psi).
\]

Уравнение \(\ddot\varphi+\omega_s^2(\varphi-\varphi^3/6)=0\) чётно по \(\varphi\). Обе ветви \(\varphi\) и \(-\varphi\) допустимы. Ветка, согласованная с \eqref{eq:phi}, умножается на \(\sigma=\operatorname{sgn}\eta_0\):

\[
\varphi
=\sigma\Biggl[
\hat\varphi\sin\psi
+\frac{\hat\varphi^3}{192}(\sin\psi+\sin 3\psi)
\Biggr].
\]

Ведущий член есть \(\sigma\hat\varphi\sin\psi\). В том же порядке \(\Omega=\omega_s\) и \(\hat\varphi=\omega_{\mathrm{rf}}|\eta_0|A/\omega_s\), поэтому

\[
\sigma\hat\varphi
=\operatorname{sgn}\eta_0\cdot\frac{\omega_{\mathrm{rf}}|\eta_0|A}{\omega_s}
=\frac{\omega_{\mathrm{rf}}}{\omega_s}\eta_0 A
=\frac{\omega_{\mathrm{rf}}}{\Omega}\eta_0 A.
\]

Это ведущий член \eqref{eq:phi}. На порядке \(\hat\varphi^3\) коэффициент при \(\sin\psi\) уже содержит вклад третьей гармоники в полуразмах и в первую гармонику \(w\), и равенство \(\varphi=(\omega_{\mathrm{rf}}/\Omega)\eta_0 A\sin\psi\) остаётся только ведущим.

## 4. Переход к \(w\) и полуразмах \(A\)

При \(\tilde\eta_1=0\)

\[
w=\frac{\dot\varphi}{\omega_{\mathrm{rf}}\eta_0}.
\]

Пусть \(\Phi\) --- выражение в скобках выше, \(\varphi=\sigma\Phi\), \(\psi=\Omega t\).

\[
\frac{\mathrm{d}\Phi}{\mathrm{d}\psi}
=\hat\varphi\cos\psi
+\frac{\hat\varphi^3}{192}(\cos\psi+3\cos 3\psi)
=\Bigl(\hat\varphi+\frac{\hat\varphi^3}{192}\Bigr)\cos\psi
+\frac{\hat\varphi^3}{64}\cos 3\psi,
\]

где \(3/192=1/64\). Тогда

\[
\dot\varphi=\sigma\Omega\frac{\mathrm{d}\Phi}{\mathrm{d}\psi},
\qquad
\eta_0=\sigma|\eta_0|,
\]

\[
w
=\frac{\sigma\Omega}{\omega_{\mathrm{rf}}\sigma|\eta_0|}\frac{\mathrm{d}\Phi}{\mathrm{d}\psi}
=\frac{\Omega}{\omega_{\mathrm{rf}}|\eta_0|}
\Biggl[
\Bigl(\hat\varphi+\frac{\hat\varphi^3}{192}\Bigr)\cos\psi
+\frac{\hat\varphi^3}{64}\cos 3\psi
\Biggr].
\]

Знак \(\eta_0\) сокращается, и \(w\) имеет вид \(+(\ldots)\cos\psi\), как в \eqref{eq:w-sol}.

В точках \(\psi=0\) и \(\psi=\pi\) третья гармоника имеет тот же знак, что и первая:

\[
w(0)=\frac{\Omega}{\omega_{\mathrm{rf}}|\eta_0|}
\Bigl(\hat\varphi+\frac{\hat\varphi^3}{192}+\frac{\hat\varphi^3}{64}\Bigr),
\qquad
w(\pi)=-w(0).
\]

Сложение \(\hat\varphi^3/64=3\hat\varphi^3/192\):

\[
\hat\varphi+\frac{4\hat\varphi^3}{192}
=\hat\varphi\Bigl(1+\frac{\hat\varphi^2}{48}\Bigr).
\]

Полуразмах \(A=(w_{\max}-w_{\min})/2=w(0)\):

\[
A=\frac{\Omega\hat\varphi}{\omega_{\mathrm{rf}}|\eta_0|}\Bigl(1+\frac{\hat\varphi^2}{48}\Bigr).
\]

Подстановка \(\Omega=\omega_s(1-\hat\varphi^2/16)\):

\begin{align*}
A
&=\frac{\omega_s\hat\varphi}{\omega_{\mathrm{rf}}|\eta_0|}
\Bigl(1-\frac{\hat\varphi^2}{16}\Bigr)\Bigl(1+\frac{\hat\varphi^2}{48}\Bigr)\\
&=\frac{\omega_s\hat\varphi}{\omega_{\mathrm{rf}}|\eta_0|}
\Bigl(1-\frac{\hat\varphi^2}{16}+\frac{\hat\varphi^2}{48}+O(\hat\varphi^4)\Bigr).
\end{align*}

\[
\frac{1}{16}-\frac{1}{48}=\frac{3}{48}-\frac{1}{48}=\frac{2}{48}=\frac{1}{24},
\]

\[
A=\frac{\omega_s\hat\varphi}{\omega_{\mathrm{rf}}|\eta_0|}\Bigl(1-\frac{\hat\varphi^2}{24}\Bigr).
\]

Деление \(w\) на \(A\). Коэффициент при \(\cos\psi\):

\[
\frac{1+\hat\varphi^2/192}{1+\hat\varphi^2/48}
=\Bigl(1+\frac{\hat\varphi^2}{192}\Bigr)\Bigl(1-\frac{\hat\varphi^2}{48}\Bigr)+O(\hat\varphi^4)
=1-\Bigl(\frac{1}{48}-\frac{1}{192}\Bigr)\hat\varphi^2+O(\hat\varphi^4).
\]

\[
\frac{1}{48}-\frac{1}{192}=\frac{4-1}{192}=\frac{3}{192}=\frac{1}{64}.
\]

Коэффициент при \(\cos 3\psi\):

\[
\frac{\hat\varphi^2/64}{1+\hat\varphi^2/48}=\frac{\hat\varphi^2}{64}+O(\hat\varphi^4).
\]

Итог

\[
w=A\Biggl[
\Bigl(1-\frac{\hat\varphi^2}{64}\Bigr)\cos\psi
+\frac{\hat\varphi^2}{64}\cos 3\psi
\Biggr].
\]

Производная скобки по \(\psi\) обращается в нуль только при \(\sin\psi=0\), пока \(\hat\varphi^2/64\) мало: если \(f=(1-\varepsilon)\cos\psi+\varepsilon\cos 3\psi\) и \(\varepsilon=\hat\varphi^2/64\), то

\[
f'=-\sin\psi\bigl(1+8\varepsilon-12\varepsilon\sin^2\psi\bigr),
\]

и скобка в правой части положительна при \(\varepsilon<1/4\). Экстремумы лежат в \(\psi=0,\pi\), и \(A\) есть полуразмах точно, а не только с точностью до \(\hat\varphi^4\).

## 5. Интеграл энергии

### 5.1. Точный маятник

При \(\tilde\eta_1=0\) и полном \(\sin\varphi\) то же исключение \(w\) даёт

\[
\ddot\varphi+\omega_s^2\sin\varphi=0.
\]

Умножение на \(\dot\varphi\):

\[
\frac{\mathrm{d}}{\mathrm{d}t}\Bigl(\frac12\dot\varphi^2-\omega_s^2\cos\varphi\Bigr)=0.
\]

На точке поворота \(\varphi=\hat\varphi\), \(\dot\varphi=0\) постоянная равна \(-\omega_s^2\cos\hat\varphi\). В точке \(\varphi=0\)

\[
\frac12\dot\varphi^2-\omega_s^2=-\omega_s^2\cos\hat\varphi,
\qquad
\dot\varphi^2=2\omega_s^2(1-\cos\hat\varphi)=4\omega_s^2\sin^2\frac{\hat\varphi}{2},
\]

\[
|\dot\varphi|_{\max}=2\omega_s\sin\frac{\hat\varphi}{2}.
\]

Орбита симметрична, \(w=\dot\varphi/(\omega_{\mathrm{rf}}\eta_0)\) меняется от \(-A\) до \(A\), и максимум \(|w|\) достигается при \(\varphi=0\):

\[
A=\frac{2\omega_s}{\omega_{\mathrm{rf}}|\eta_0|}\sin\frac{\hat\varphi}{2}.
\]

Ряд \(\sin x=x-x^3/6+O(x^5)\) при \(x=\hat\varphi/2\):

\[
\sin\frac{\hat\varphi}{2}
=\frac{\hat\varphi}{2}-\frac{\hat\varphi^3}{48}+O(\hat\varphi^5),
\qquad
2\sin\frac{\hat\varphi}{2}
=\hat\varphi\Bigl(1-\frac{\hat\varphi^2}{24}\Bigr)+O(\hat\varphi^5).
\]

Разложение точной связи до членов \(\hat\varphi^2\) совпадает с результатом Боголюбова--Митропольского.

### 5.2. Кубический осциллятор

Умножение \(\ddot\varphi+\omega_s^2(\varphi-\varphi^3/6)=0\) на \(\dot\varphi\):

\[
\frac{\mathrm{d}}{\mathrm{d}t}\Biggl(
\frac12\dot\varphi^2+\omega_s^2\frac{\varphi^2}{2}-\omega_s^2\frac{\varphi^4}{24}
\Biggr)=0.
\]

Потенциал \(\varphi^2/2-\varphi^4/24\) есть начало ряда \(1-\cos\varphi\). На точке поворота

\[
\frac12\dot\varphi_{\max}^2
=\omega_s^2\frac{\hat\varphi^2}{2}-\omega_s^2\frac{\hat\varphi^4}{24}
=\omega_s^2\frac{\hat\varphi^2}{2}\Bigl(1-\frac{\hat\varphi^2}{12}\Bigr),
\]

\[
|\dot\varphi|_{\max}
=\omega_s\hat\varphi\sqrt{1-\frac{\hat\varphi^2}{12}}
=\omega_s\hat\varphi\Bigl(1-\frac{\hat\varphi^2}{24}\Bigr)+O(\hat\varphi^5).
\]

Деление на \(\omega_{\mathrm{rf}}|\eta_0|\) снова даёт множитель \((1-\hat\varphi^2/24)\). Совпадение с точным маятником на этом порядке следует из совпадения потенциалов до \(\varphi^4\).

### 5.3. Частота точного маятника

Период либрации есть \(T=4\mathrm{K}(k)/\omega_s\), \(k=\sin(\hat\varphi/2)\), поэтому

\[
\frac{\Omega}{\omega_s}=\frac{\pi}{2\mathrm{K}(k)}.
\]

Ряд модуля

\[
\mathrm{K}(k)=\frac{\pi}{2}\Bigl(1+\frac14 k^2+O(k^4)\Bigr),
\qquad
\frac{\pi}{2\mathrm{K}}=1-\frac{k^2}{4}+O(k^4).
\]

\[
k^2=\sin^2\frac{\hat\varphi}{2}=\frac{\hat\varphi^2}{4}+O(\hat\varphi^4),
\qquad
\frac{k^2}{4}=\frac{\hat\varphi^2}{16}+O(\hat\varphi^4),
\]

\[
\frac{\Omega}{\omega_s}=1-\frac{\hat\varphi^2}{16}+O(\hat\varphi^4).
\]

Кубическое приближение воспроизводит частоту точного маятника до членов \(\hat\varphi^2\). Следующий член точного ряда есть \(+\hat\varphi^4/3072\); в главе он не удерживается.

## 6. Первая поправка по \(\tilde\eta_1\)

При \(\tilde\eta_1\neq 0\) и \(\Lambda=0\) дифференцирование второго уравнения и подстановка первого дают

\[
\ddot\varphi
=-\omega_s^2\Bigl(1+\frac{2\tilde\eta_1}{\eta_0}w\Bigr)\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr).
\]

Раскрытие до первого порядка по \(\tilde\eta_1\) и отбрасывание \(\tilde\eta_1\varphi^3\):

\[
\ddot\varphi+\omega_s^2\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr)
=-2\omega_s^2\frac{\tilde\eta_1}{\eta_0}w\varphi.
\]

На орбите \(\tilde\eta_1=0\) в ведущем порядке \(w=A\cos\psi\), \(\varphi=\sigma\hat\varphi\sin\psi\), и произведение \(w\varphi\) пропорционально \(\sin 2\psi\). Первой гармоники в нём нет, поэтому линейный по \(\tilde\eta_1\) сдвиг частоты из этого произведения не возникает.

Постоянная часть \(w\) есть \(w_{\mathrm{eq}}=O(\tilde\eta_1 A^2)\). Её вклад в правую часть имеет порядок \(\tilde\eta_1^2\) и уже содержится в члене \((5/12)(\tilde\eta_1 A/\eta_0)^2\) главы 1. Поправки формы порядка \(\tilde\eta_1\hat\varphi^2\) в задании отброшены.

Квадратичная сила главы 1 порождает вторую гармонику \(\tilde\eta_1 A^2\cos 2\psi/(6\eta_0)\) и сдвиг частоты \(-(5/12)(\tilde\eta_1 A/\eta_0)^2\). Кубическая сила порождает третью гармонику и сдвиг \(-\hat\varphi^2/16\). Перекрёстное искажение второй гармоники есть порядок \(\tilde\eta_1 A^2\hat\varphi^2\) и не удерживается. Произведение двух относительных сдвигов частоты есть порядок \((\tilde\eta_1 A/\eta_0)^2\hat\varphi^2\) и тоже не удерживается. Сложение относительных поправок даёт

\[
w(t)
=w_{\mathrm{eq}}
+A\Biggl[
\Bigl(1-\frac{\hat\varphi^2}{64}\Bigr)\cos\psi
+\frac{\hat\varphi^2}{64}\cos 3\psi
\Biggr]
+\frac{\tilde\eta_1 A^2}{6\eta_0}\cos 2\psi,
\]

\[
\frac{\Omega}{\omega_s}
=1-\frac{\hat\varphi^2}{16}
-\frac{5}{12}\Bigl(\frac{\tilde\eta_1 A}{\eta_0}\Bigr)^2.
\]

Вторая гармоника в точках \(\psi=0\) и \(\psi=\pi\) равна одной и той же величине. Вместе с проверкой экстремумов из п. 4 полуразмах остаётся равным \(A\) на удержанном порядке: сдвиг точки экстремума второй гармоникой начинается с \((\tilde\eta_1 A/\eta_0)^2\).

Постоянная \(w_{\mathrm{eq}}\) этим разложением не фиксируется. Её даёт \(\langle\tau\rangle=0\).

## 7. Среднее \(\langle w^2\rangle\)

Равенства \eqref{eq:deq-avg} следуют только из \(\langle\tau\rangle=0\) и ряда для \(\tau\). Второе уравнение системы в них не входит. Фазовая поправка меняет лишь значение \(\langle w^2\rangle\).

На решении с \(\tilde\eta_1=0\)

\[
w=A\bigl[(1-\varepsilon)\cos\psi+\varepsilon\cos 3\psi\bigr],
\qquad
\varepsilon=\frac{\hat\varphi^2}{64}.
\]

Средние \(\langle\cos^2\psi\rangle=\langle\cos^2 3\psi\rangle=1/2\), \(\langle\cos\psi\cos 3\psi\rangle=0\):

\begin{align*}
\langle w^2\rangle
&=\frac{A^2}{2}\bigl[(1-\varepsilon)^2+\varepsilon^2\bigr]
=\frac{A^2}{2}\bigl[1-2\varepsilon+2\varepsilon^2\bigr]\\
&=\frac{A^2}{2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)+O(\hat\varphi^4),
\end{align*}

поскольку \(2\varepsilon=\hat\varphi^2/32\), а \(\varepsilon^2=O(\hat\varphi^4)\).

Механизм поправки: при \(\psi=0\) вклад гармоник есть \(A(1-\varepsilon)+A\varepsilon=A\). Полуразмах фиксирован, а коэффициент первой гармоники уменьшен до \(A(1-\varepsilon)\). Её вклад в среднее квадрата падает на \(A^2\varepsilon\), вклад третьей гармоники есть \(A^2\varepsilon^2/2=O(\hat\varphi^4)\).

Связь \(w\) и \(\delta\). Из \eqref{eq:w-series}

\[
\delta=w-\frac{w^2}{2\gamma_s^2}+O(w^3),
\qquad
\delta^2=w^2+O(w^3).
\]

Разность \(\langle\delta^2\rangle-\langle w^2\rangle\) есть \(O(A^3)\). Этот порядок отброшен уже в главах 1--2. Поправка \(O(A^2\hat\varphi^2)\) лежит внутри удержанного \(O(A^2)\) и одинакова у обоих средних:

\[
\langle\delta^2\rangle=\langle w^2\rangle
=\frac{A^2}{2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)+O(\hat\varphi^4).
\]

Почему среднее при \(\tilde\eta_1\neq 0\) берётся с орбиты \(\tilde\eta_1=0\). Запись

\[
w=w_{\mathrm{eq}}+u_0+D\cos 2\psi,
\qquad
D=\frac{\tilde\eta_1 A^2}{6\eta_0}=O(\tilde\eta_1),
\qquad
w_{\mathrm{eq}}=O(\tilde\eta_1),
\]

где \(u_0\) --- правая часть \eqref{eq:w-ph}. Тогда

\[
\langle w^2\rangle
=\langle u_0^2\rangle+w_{\mathrm{eq}}^2+\frac{D^2}{2}
+2w_{\mathrm{eq}}\langle u_0\rangle
+2D\langle u_0\cos 2\psi\rangle.
\]

Среднее \(u_0\) равно нулю. Произведения \(\langle\cos\psi\cos 2\psi\rangle\) и \(\langle\cos 3\psi\cos 2\psi\rangle\) равны нулю. Оставшиеся добавки суть \(O(\tilde\eta_1^2)\). В первом порядке по \(\tilde\eta_1\)

\[
\langle w^2\rangle=\frac{A^2}{2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr).
\]

Подстановка в \eqref{eq:deq-avg}:

\[
\delta_{\mathrm{eq}}
=-\frac{\eta_1}{\eta_0}\langle\delta^2\rangle
=-\frac{\eta_1 A^2}{2\eta_0}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr),
\]

\[
w_{\mathrm{eq}}
=-\frac{\tilde\eta_1}{\eta_0}\langle w^2\rangle
=-\frac{\tilde\eta_1 A^2}{2\eta_0}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr).
\]

Проверка связи \eqref{eq:w-series}. При том же \(\langle\delta^2\rangle\)

\begin{align*}
\delta_{\mathrm{eq}}+\frac{\langle\delta^2\rangle}{2\gamma_s^2}
&=-\frac{\eta_1 A^2}{2\eta_0}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)
+\frac{A^2}{4\gamma_s^2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)\\
&=-\frac{A^2}{2\eta_0}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)
\Bigl(\eta_1-\frac{\eta_0}{2\gamma_s^2}\Bigr)\\
&=-\frac{\tilde\eta_1 A^2}{2\eta_0}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr).
\end{align*}

Второе равенство \eqref{eq:deq-ph} воспроизводится. Множитель \((1-\hat\varphi^2/32)\) стоит и при слагаемом \(A^2/(4\gamma_s^2)\):

\[
w_{\mathrm{eq}}
=\delta_{\mathrm{eq}}+\frac{A^2}{4\gamma_s^2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr).
\]

Запись \(w_{\mathrm{eq}}=\delta_{\mathrm{eq}}+A^2/(4\gamma_s^2)\) без этого множителя расходится с \eqref{eq:w-series} на

\[
\frac{A^2}{4\gamma_s^2}\cdot\frac{\hat\varphi^2}{32}
=\frac{A^2\hat\varphi^2}{128\gamma_s^2}.
\]

В текст главы поставлен множитель. При \(\Lambda=0\) он нужен, чтобы формула с бетатронным членом переходила в \eqref{eq:deq-ph}.

## 8. Эллиптическая формула

Для маятника \(\ddot\varphi+\omega_s^2\sin\varphi=0\) либрация внутри сепаратрисы, \(0<\hat\varphi<\pi\), записывается через эллиптические функции Якоби. Модуль \(k=\sin(\hat\varphi/2)\), аргумент \(u=\omega_s t\):

\[
\sin\frac{\varphi}{2}=k\,\mathrm{sn}(u,k).
\]

Проверка скорости. Дифференцирование левой части:

\[
\frac{\mathrm{d}}{\mathrm{d}t}\sin\frac{\varphi}{2}
=\frac12\cos\frac{\varphi}{2}\cdot\dot\varphi.
\]

Правая часть: \(k\,\mathrm{cn}\,\mathrm{dn}\cdot\omega_s\). На этой орбите \(\cos(\varphi/2)=\mathrm{dn}(u,k)\), поэтому

\[
\dot\varphi=2k\omega_s\,\mathrm{cn}(u,k).
\]

При \(u=0\) функция \(\mathrm{cn}=1\), \(\varphi=0\) и \(\dot\varphi=2k\omega_s=2\omega_s\sin(\hat\varphi/2)\), что совпадает с интегралом энергии. При \(u=\mathrm{K}\) функция \(\mathrm{sn}=1\), \(\varphi=\hat\varphi\) и \(\mathrm{cn}=0\).

Квадрат скорости \(\dot\varphi^2=4k^2\omega_s^2\mathrm{cn}^2\). Среднее по периоду \(4\mathrm{K}\) сводится к среднему по четверти периода. Интеграл

\[
\int_0^{\mathrm{K}}\mathrm{dn}^2 u\,\mathrm{d}u=\mathrm{E}(k)
\]

есть определение полного интеграла второго рода в переменной Якоби, эквивалентное

\[
\mathrm{K}(k)=\int_0^{\pi/2}\frac{\mathrm{d}\theta}{\sqrt{1-k^2\sin^2\theta}},
\qquad
\mathrm{E}(k)=\int_0^{\pi/2}\sqrt{1-k^2\sin^2\theta}\,\mathrm{d}\theta.
\]

Аргумент --- модуль \(k\), не параметр \(k^2\). Дальше

\[
\mathrm{dn}^2=1-k^2\mathrm{sn}^2,
\qquad
\int_0^{\mathrm{K}}\mathrm{sn}^2 u\,\mathrm{d}u=\frac{\mathrm{K}-\mathrm{E}}{k^2},
\]

\[
\int_0^{\mathrm{K}}\mathrm{cn}^2 u\,\mathrm{d}u
=\int_0^{\mathrm{K}}(1-\mathrm{sn}^2)\,\mathrm{d}u
=\mathrm{K}-\frac{\mathrm{K}-\mathrm{E}}{k^2}
=\frac{\mathrm{E}-(1-k^2)\mathrm{K}}{k^2}.
\]

\[
\langle\mathrm{cn}^2\rangle
=\frac{\mathrm{E}/\mathrm{K}-1+k^2}{k^2}.
\]

\[
\langle\dot\varphi^2\rangle
=4k^2\omega_s^2\langle\mathrm{cn}^2\rangle
=4\omega_s^2\bigl(\mathrm{E}/\mathrm{K}-1+k^2\bigr).
\]

Максимум \(\dot\varphi^2\) равен \(4k^2\omega_s^2\). Отношение

\[
\frac{\langle\dot\varphi^2\rangle}{\dot\varphi_{\max}^2}
=\frac{\mathrm{E}(k)/\mathrm{K}(k)-1+k^2}{k^2}.
\]

Поскольку \(w\propto\dot\varphi\) и \(A=|\dot\varphi|_{\max}/(\omega_{\mathrm{rf}}|\eta_0|)\),

\[
\langle w^2\rangle
=A^2\frac{\mathrm{E}(k)/\mathrm{K}(k)-1+k^2}{k^2}.
\]

Та же величина получается прямой квадратурой по орбите. При \(\omega_s=1\)

\[
\langle\dot\varphi^2\rangle
=\frac{1}{T}\oint\dot\varphi\,\mathrm{d}\varphi,
\qquad
T=4\mathrm{K}(k),
\qquad
\dot\varphi=\sqrt{2(\cos\varphi-\cos\hat\varphi)},
\]

и деление на \(\dot\varphi_{\max}^2=2(1-\cos\hat\varphi)\) воспроизводит эллиптическое выражение.

### 8.1. Разложение

Ряды до \(k^4\):

\begin{align*}
\mathrm{K}
&=\frac{\pi}{2}\Bigl(1+\frac14 k^2+\frac{9}{64}k^4+O(k^6)\Bigr),\\
\mathrm{E}
&=\frac{\pi}{2}\Bigl(1-\frac14 k^2-\frac{3}{64}k^4+O(k^6)\Bigr).
\end{align*}

Деление:

\[
\frac{\mathrm{E}}{\mathrm{K}}
=1-\frac{k^2}{2}-\frac{k^4}{16}+O(k^6),
\]

\[
\frac{\mathrm{E}}{\mathrm{K}}-1+k^2=\frac{k^2}{2}-\frac{k^4}{16}+O(k^6),
\]

\[
\frac{\mathrm{E}/\mathrm{K}-1+k^2}{k^2}
=\frac12-\frac{k^2}{16}+O(k^4).
\]

\[
k^2=\frac{\hat\varphi^2}{4}+O(\hat\varphi^4),
\qquad
\frac{k^2}{16}=\frac{\hat\varphi^2}{64}+O(\hat\varphi^4),
\]

\[
\frac{\langle w^2\rangle}{A^2}
=\frac12-\frac{\hat\varphi^2}{64}+O(\hat\varphi^4)
=\frac12\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)+O(\hat\varphi^4).
\]

Множитель фурье-разложения подтверждён разложением точного маятника. Следующий член есть \(-\hat\varphi^4/1536\) в \(\langle w^2\rangle/A^2\).

## 9. Бетатронное запаздывание

Система с \(\tau=\eta_0 w+\tilde\eta_1 w^2+\Lambda\) и тем же кубическим рядом:

\[
\dot\varphi=\omega_{\mathrm{rf}}(\eta_0 w+\tilde\eta_1 w^2+\Lambda),
\qquad
\dot w=K\cos\varphi_s\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr).
\]

При \(\tilde\eta_1=0\) неподвижная точка: \(\varphi-\varphi^3/6=0\) даёт \(\varphi=0\) внутри сепаратрисы, а \(\eta_0 w+\Lambda=0\) даёт \(w=-\Lambda/\eta_0\). В переменной \(u=w+\Lambda/\eta_0\)

\[
\dot\varphi=\omega_{\mathrm{rf}}\eta_0 u,
\qquad
\dot u=K\cos\varphi_s\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr),
\]

то есть снова \eqref{eq:phi-ode}. Полуразмах \(A\), связь с \(\hat\varphi\) и среднее \(\langle u^2\rangle\) те же, что в пп. 4--8. Сдвиг \(\Lambda\) частоту маятника не меняет: \(\Omega/\omega_s=1-\hat\varphi^2/16\).

В первом порядке по \(\tilde\eta_1\) и по \(\Lambda\) квадрат среднего, \((\Lambda/\eta_0)^2\), в \(\langle w^2\rangle\) умножается на \(\tilde\eta_1\) и даёт в \(w_{\mathrm{eq}}\) член порядка \(\tilde\eta_1\Lambda^2\). Он отбрасывается, как \(\Lambda^2\) в главе 2. Остаётся

\[
\langle\delta^2\rangle
=\frac{A^2}{2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr).
\]

Из \eqref{eq:deq-raw-L}

\[
\delta_{\mathrm{eq}}
=-\frac{1}{\eta_0}\Biggl[
\frac{\eta_1 A^2}{2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)+\Lambda
\Biggr].
\]

Из \eqref{eq:w-series}

\[
w_{\mathrm{eq}}
=\delta_{\mathrm{eq}}+\frac{\langle\delta^2\rangle}{2\gamma_s^2}
=\delta_{\mathrm{eq}}+\frac{A^2}{4\gamma_s^2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr).
\]

Подстановка \(\tilde\eta_1=\eta_1-\eta_0/(2\gamma_s^2)\) и \eqref{eq:tau-w}:

\begin{align*}
w_{\mathrm{eq}}
&=-\frac{\eta_1 A^2}{2\eta_0}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)
-\frac{\Lambda}{\eta_0}
+\frac{A^2}{4\gamma_s^2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)\\
&=-\frac{1}{\eta_0}\Biggl[
\frac{\tilde\eta_1 A^2}{2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)+\Lambda
\Biggr].
\end{align*}

При \(\Lambda=0\) это \eqref{eq:deq-ph}.

Частота. Линейный сдвиг \(-\tilde\eta_1\Lambda/\eta_0^2\) и квадратичный \(-(5/12)(\tilde\eta_1 A/\eta_0)^2\) взяты из \eqref{eq:Omega-L}: они получены при \(\sin\varphi\approx\varphi\). Кубический сдвиг \(-\hat\varphi^2/16\) вычислен при \(\tilde\eta_1=0\). Относительные множители перемножаются,

\[
\Bigl(1-\frac{\tilde\eta_1\Lambda}{\eta_0^2}-\frac{5}{12}\Bigl(\frac{\tilde\eta_1 A}{\eta_0}\Bigr)^2\Bigr)
\Bigl(1-\frac{\hat\varphi^2}{16}\Bigr),
\]

и произведения малых параметров отбрасываются по условиям

\[
\Bigl|\frac{\eta_1 A}{\eta_0}\Bigr|\ll 1,
\qquad
\Bigl|\frac{\tilde\eta_1\Lambda}{\eta_0^2}\Bigr|\ll 1,
\qquad
\hat\varphi^2\ll 1.
\]

Остаётся

\[
\frac{\Omega}{\omega_s}
=1-\frac{\tilde\eta_1\Lambda}{\eta_0^2}
-\frac{5}{12}\Bigl(\frac{\tilde\eta_1 A}{\eta_0}\Bigr)^2
-\frac{\hat\varphi^2}{16}.
\]

Эллиптическая формула для \(\langle u^2\rangle\) от условия \(\hat\varphi^2\ll 1\) не зависит. Она требует только, чтобы траектория лежала внутри сепаратрисы, где период конечен и \(k<1\).

## 10. Расхождение с черновой записью кинематического члена

В черновой формулировке пункта с \(\Lambda\) было

\[
w_{\mathrm{eq}}=\delta_{\mathrm{eq}}+\frac{A^2}{4\gamma_s^2}
\]

без множителя \((1-\hat\varphi^2/32)\). Формулы \(\delta_{\mathrm{eq}}\) и \(w_{\mathrm{eq}}=-\tilde\eta_1\langle w^2\rangle/\eta_0\) при этом множитель содержат. Из п. 7 следует, что кинематический член обязан нести тот же множитель, иначе при \(\Lambda=0\) выражение не совпадает с \eqref{eq:deq-ph}, а разность равна \(A^2\hat\varphi^2/(128\gamma_s^2)\). В главу поставлен согласованный вариант. Остальные заданные коэффициенты --- \(1/16\), \(1/24\), \(1/64\), \(1/192\), \(5/12\), \(1/6\) и эллиптическое отношение --- воспроизводятся выкладками выше и отдельной символьной проверкой.
