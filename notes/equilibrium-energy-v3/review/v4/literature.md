# Литературный обзор: равновесный сдвиг и slip factor

Проверено по JACoW, INSPIRE, страницам издателей, arXiv и, где указано, по тексту PDF. Файлы `sections/` и `references.bib` не менялись. Если том, страницы, DOI или формула не подтверждены текстом, стоит пометка **не проверено**.

Обозначения заметки: \(\tau=(T-T_s)/T_s=\eta_0\delta+\eta_1\delta^2+\Lambda\), \(\Lambda=-2\pi J\xi/C\), в ЭС-кольце \(\eta_0^{\mathrm{E}}=(1+1/\gamma_s^2)\alpha_0-1/\gamma_s^2\).

Главные расхождения с заметкой вынесены в конце.

---

## 1. Senichev2013

**Статус:** выходные данные проверены (JACoW и шапка PDF). Формулы (6), (9), (10), (18), (25) прочитаны с изображения страниц PDF. Символа \(\gamma_{\mathrm{eff}}\) в статье нет.

DOI у IPAC 2013 нет. В шапке PDF автор «N. Kulabukhova»; в карточке JACoW — «N. V. Kulabukhova».

```bibtex
@inproceedings{Senichev2013,
  author    = {Senichev, Yu. and Maier, R. and Zyuzin, D. and Kulabukhova, N.},
  title     = {Spin Tune Decoherence Effects in Electro- and Magnetostatic Structures},
  booktitle = {Proc. 4th International Particle Accelerator Conference (IPAC'13)},
  year      = {2013},
  pages     = {2579--2581},
  address   = {Shanghai, China},
  publisher = {JACoW Publishing, Geneva, Switzerland},
  isbn      = {978-3-95450-122-9},
  note      = {paper WEPEA036},
  url       = {https://proceedings.jacow.org/IPAC2013/papers/wepea036.pdf}
}
```

Slip factor введён по периоду: \(\Delta T_{\mathrm{rev}}/T_{\mathrm{rev}}=\eta\delta\), \(\eta=\alpha_0-1/\gamma^2\). Разложение до \(\delta^2\), уравнение (6):

\[
\frac{\Delta T}{T}
=\Bigl(\alpha_0-\frac{1}{\gamma^2}\Bigr)\delta
+\Bigl(\alpha_1-\frac{\alpha_0}{\gamma^2}+\frac{1}{\gamma^4}\Bigr)\delta^2
+\Bigl(\frac{\Delta L}{L}\Bigr)_\beta.
\]

Равновесный сдвиг, уравнение (10), прочитан дважды и согласован с ограничением (18):

\[
\Delta\delta_{\mathrm{eq}}
=\frac{\gamma_s^2}{\gamma_s^2\alpha_0-1}
\Biggl[
\frac{\delta_m^2}{2}
\Bigl(\alpha_1-\frac{\alpha_0}{\gamma_s^2}+\frac{1}{\gamma_s^4}\Bigr)
+\Bigl(\frac{\Delta L}{L}\Bigr)_\beta
\Biggr].
\]

Предфактор равен \(1/\eta_0\). В (18) спиновый разброс записан через \(\langle\Delta\gamma_{\mathrm{eq}}/\gamma\rangle\), и обращение (10) содержит множитель \(2/\beta^2\). Это и есть связь \(\Delta\gamma/\gamma=\beta^2\Delta\delta\), то есть \(\gamma+\Delta\gamma_{\mathrm{eq}}=\gamma+\beta^2\gamma\,\Delta\delta_{\mathrm{eq}}\), но символа \(\gamma_{\mathrm{eff}}\) нет. В (16) \((\Delta L/L)_\beta\) заменён членом с \(\varepsilon\nu\); численный множитель \(\pi/(2L)\) прочитан один раз и перед цитированием лучше сверить глазами.

В ЭС-дефлекторе, раздел «Features of electrostatic ring», разброс скорости уменьшается множителем \(F_v\) из-за перехода кинетической энергии в потенциальную. В (25) линейный член имеет вид \((\alpha_0-1/(F_v\gamma^2))\delta\); у квадратичных слагаемых \(\alpha_0/\gamma^2\) и \(1/\gamma^4\) на изображении стоит \(1/F_v^2\). Рядом в тексте сказано только «дополнительный множитель \(1/F_v\)». Степени \(F_v\) при \(\delta^2\) стоит перечитать в PDF. Там же: slip factor может сменить знак, и синхронная фаза сдвигается на \(\pi\). Секступоли правят \(\alpha_1\) и прямое бетатронное удлинение; нужны три семейства.

**Куда ссылаться.** Глава 6 — формула коллег: \(\Delta\delta_{\mathrm{eq}}\) и \(F_v\). Глава 1 — их коэффициент при \(\delta^2\) равен \(1/\gamma^4\), а не \((3\gamma^2-1)/(2\gamma^4)\). Главы 2 и 5 — \((\Delta L/L)_\beta\) у них геометрическое (\(\varepsilon\nu\)), связь \(\Lambda=-2\pi J\xi/C\) в этой статье не выведена; секступоли введены отдельно.

---

## 2. Точное название «Spin decoherence in the frozen spin…»

**Статус:** на INSPIRE и в JACoW это **Aksentyev2019**, авторы A. E. Aksentyev и Y. V. Senichev. E. Valetov в авторах нет: в списке литературы статьи он стоит как автор диссертации. Карточка: https://inspirehep.net/literature/1743479.

```bibtex
@inproceedings{Aksentyev2019,
  author    = {Aksentyev, A. E. and Senichev, Y. V.},
  title     = {Spin Decoherence in the Frozen Spin Storage Ring Method of Search for a Particle {EDM}},
  booktitle = {Proc. 10th International Particle Accelerator Conference (IPAC'19)},
  year      = {2019},
  pages     = {864--866},
  address   = {Melbourne, Australia},
  publisher = {JACoW Publishing, Geneva, Switzerland},
  doi       = {10.18429/JACoW-IPAC2019-MOPTS012},
  note      = {paper MOPTS012},
  url       = {https://accelconf.web.cern.ch/ipac2019/papers/mopts012.pdf}
}
```

В JACoW второй автор записан как «V. Senichev»; в шапке PDF — «Y. V. Senichev». Частота прецессии зависит от equilibrium-level energy, разного у частиц с разной длиной орбиты. Секступоли убирают параболическую зависимость спинового тюна от координаты в идеальном кольце; линейный остаток предлагают давить ВЧ. В неидеальном кольце декогеренция переносится в вертикальную плоскость. Формул \(\eta_0\), \(\eta_1\), \(F_v\) в доступном тексте нет.

Близкая по авторам работа Valetov, где Senichev соавтор, — диссертация, а не эта статья:

```bibtex
@phdthesis{Valetov2017,
  author = {Valetov, Eremey Vladimirovich},
  title  = {Field Modeling, Symplectic Tracking, and Spin Decoherence for {EDM} and Muon {$g$-2} Lattices},
  school = {Michigan State University},
  year   = {2017},
  note   = {FERMILAB-THESIS-2017-21},
  url    = {https://inspirehep.net/literature/1647076}
}
```

**Статус Valetov2017:** название, автор, университет, год и номер Fermilab проверены по JuSER, OSTI и титулу PDF. Число страниц 348 — по карточке JuSER.

**Куда ссылаться.** Глава 6: Aksentyev2019 — «equilibrium-level energy» и секступольная компенсация, без вывода \(\delta_{\mathrm{eq}}\). Valetov2017 — глава 4, сохранение энергии на крае ЭС-элемента (несохранение энергии, если край смоделирован неверно).

---

## 3. Berg2007

**Статус:** проверено по INSPIRE, arXiv и HTML ar5iv. Журнальный год 2007, препринт 2006. Ключ Berg2007, не Berg2006.

```bibtex
@article{Berg2007,
  author  = {Berg, J. Scott},
  title   = {Amplitude dependence of time of flight and its connection to chromaticity},
  journal = {Nucl. Instrum. Methods Phys. Res., Sect. A},
  volume  = {570},
  pages   = {15--21},
  year    = {2007},
  doi     = {10.1016/j.nima.2006.09.099},
  eprint  = {physics/0607275},
  archivePrefix = {arXiv},
  note    = {submitted 29 July 2006; report BNL-76759-2006-JA}
}
```

Из гамильтониана в нормальных координатах, уравнение (5) препринта:

\[
\frac{d\bar t}{ds}
=-\partial_E H_T-\frac{2\pi(\partial_E\vec\nu)\cdot\vec J_n}{L}+O(J_n^{3/2}).
\]

Низший член времени пролёта линеен по действию и пропорционален хроматичности. Для стоячего бакета сдвиг энергии неподвижной точки, уравнение (7):

\[
E_0+\frac{2\pi\vec\xi\cdot\vec J_n}{T_0\eta},
\]

где \(\eta\) — frequency slip factor, а \(\vec\xi\) определена как \(\vec\nu+\vec\xi\,(\Delta p/p)\). Действие \(J_n\) у Берга имеет размерность энергия\(\times\)время: на опорный импульс он не делит. Поэтому численный множитель нельзя переписывать в \(\Lambda=-2\pi J\xi/C\) без перевода действия. Секступоль, правящий хроматичность, даёт средний дипольный толчок и укорачивает орбиту большой амплитуды — тот же механизм, что сдвиг центра в главах 2 и 5.

**Куда ссылаться.** Главы 2 и 5: связь запаздывания с \(\xi\). Глава 6: сдвиг равновесной энергии пропорционален \(\xi J/\eta\). Магнитное кольцо и линак; ЭС-множителя \(1+1/\gamma^2\) нет.

---

## 4. Shoji2004 и Shoji2005

**Статус:** обе карточки проверены по INSPIRE и страницам APS. Shoji2005 уже стоит в `references.bib` и с этой карточкой совпадает.

```bibtex
@article{Shoji2004,
  author  = {Shoji, Yoshihiko},
  title   = {Bunch lengthening by a betatron motion in quasi-isochronous storage rings},
  journal = {Phys. Rev. ST Accel. Beams},
  volume  = {7},
  pages   = {090703},
  year    = {2004},
  doi     = {10.1103/PhysRevSTAB.7.090703}
}

@article{Shoji2005,
  author  = {Shoji, Y.},
  title   = {Dependence of average path length on betatron motion in a storage ring},
  journal = {Phys. Rev. ST Accel. Beams},
  volume  = {8},
  pages   = {094001},
  year    = {2005},
  doi     = {10.1103/PhysRevSTAB.8.094001}
}
```

Shoji2004 (опубликована 28 сентября 2004): локальное удлинение сгустка от линейного горизонтального бетатронного движения в дисперсионной секции; длина зависит от точки кольца. Это не среднее за оборот.

Shoji2005 (опубликована 27 сентября 2005), аннотация APS: изменение средней длины пути за оборот есть сумма произведений инварианта бетатронного движения и хроматичности в обеих плоскостях. Из-за этого в квазиизохронном кольце при большой хроматичности появляется заметный разброс импульса. Коэффициент \(-2\pi/C\) в этом проходе из PDF не переписывался: страница APS отдала 403. **Численный предфактор не проверен.**

**Куда ссылаться.** Главы 2 и 5: Shoji2005 — среднее \(\Lambda\) через \(\xi\); Shoji2004 — локальное удлинение до усреднения по обороту. ЭС-кольца нет.

---

## 5. Metodiev2014

**Статус:** проверено по INSPIRE и PDF APS (аннотация и список авторов).

```bibtex
@article{Metodiev2014,
  author  = {Metodiev, E. M. and Huang, K. L. and Semertzidis, Y. K. and Morse, W. M.},
  title   = {Fringe electric fields of flat and cylindrical deflectors in electrostatic charged particle storage rings},
  journal = {Phys. Rev. ST Accel. Beams},
  volume  = {17},
  pages   = {074002},
  year    = {2014},
  doi     = {10.1103/PhysRevSTAB.17.074002}
}
```

Опубликована 18 июля 2014. Аналитические потенциал и поле плоских и цилиндрических пластин, включая краевое поле; цилиндрический край записан в полярных координатах. Дальше — трекинг спина в геометрии proton EDM. Показатель \(n\) в \(E_r\propto r^{-n}\) и утверждение «кинетическая энергия на смещённой орбите ниже» в аннотации не формулируются: цилиндр отвечает частному случаю логарифмического потенциала.

**Куда ссылаться.** Глава 4: реальное краевое поле вместо жёсткого края. К ступеньке кинетической энергии на жёстком крае ближе ValetovBerz2020.

---

## 6. ValetovBerz2020 и Valetov2020

**Статус:** проверено по OSTI, сайту Elsevier и странице автора. Том **213**, не 215. OSTI датирует главы 30 декабря 2019; том вышел 18 марта 2020. Год в ключе — 2020.

```bibtex
@incollection{ValetovBerz2020,
  author    = {Valetov, Eremey and Berz, Martin},
  title     = {Derivation, cross-validation, and comparison of analytic formulas for electrostatic deflector aberrations},
  booktitle = {Advances in Imaging and Electron Physics},
  publisher = {Elsevier},
  volume    = {213},
  pages     = {145--203},
  year      = {2020},
  doi       = {10.1016/bs.aiep.2019.11.007},
  note      = {OSTI 1661696}
}

@incollection{Valetov2020,
  author    = {Valetov, Eremey},
  title     = {Analysis and fringe field scaling of a legacy set of electrostatic deflector aberration formulas},
  booktitle = {Advances in Imaging and Electron Physics},
  publisher = {Elsevier},
  volume    = {213},
  pages     = {205--224},
  year      = {2020},
  doi       = {10.1016/bs.aiep.2019.11.008},
  note      = {OSTI 1661697}
}
```

Формулы Вольника 1965 описывают главное поле и не выполняют симплектичность, потому что не учитывают изменение кинетической энергии при изменении потенциала. Поправка жёсткого края: на входе и выходе наклоны масштабируются, \(\alpha_O=(1-u)\alpha_I\). Продольная скорость на краю скачет из-за \(\delta\)-образного \(E_z\), поперечная непрерывна. После этого масштаба формулы совпадают с COSY INFINITY. В координатах COSY жёсткий край есть тождественное отображение \((x,a)\), и отдельного масштаба импульса не нужно.

Исходная статья Вольника:

```bibtex
@article{Wollnik1965,
  author  = {Wollnik, H.},
  title   = {Second order approximation of the three-dimensional trajectories of charged particles in deflecting electrostatic and magnetic fields},
  journal = {Nucl. Instrum. Methods},
  volume  = {34},
  pages   = {213--221},
  year    = {1965},
  doi     = {10.1016/0029-554X(65)90297-1}
}
```

**Статус Wollnik1965:** название, том, страницы 213–221, год и DOI проверены по ScienceDirect.

**Куда ссылаться.** Глава 4: кинетическая энергия ниже внутри дефлектора, на жёстком крае продольный импульс пересчитывается. Это не вывод \(\eta_0^{\mathrm{E}}\) кольца.

---

## 7. Wollnik1987

**Статус:** проверено по каталогу GSI, Elsevier и каталогу CERN. ISBN-10 и ISBN-13 относятся к одному изданию.

```bibtex
@book{Wollnik1987,
  author    = {Wollnik, Hermann},
  title     = {Optics of Charged Particles},
  publisher = {Academic Press},
  address   = {Orlando},
  year      = {1987},
  isbn      = {0-12-762130-X},
  pages     = {293},
  note      = {1st ed.; also 978-0-12-762130-2; ix+293 pp.}
}
```

Книга задаёт траектории в главном и краевом поле электрических и магнитных диполей, квадруполей и гексаполей в одной системе обозначений. Вывод slip factor кольца и множителя \(1+1/\gamma^2\) по оглавлению и аннотации не проверялся: **не проверено**, есть ли там \(\eta_1\).

**Куда ссылаться.** Глава 4: справочник по потенциалу секторного дефлектора и краевому полю, рядом с Wollnik1965. Не вместо ValetovBerz2020: книга 1987 года старше разбора ошибки 1965 года.

---

## 8. Berz1999

**Статус:** библиография проверена по INSPIRE и ScienceDirect. Два ISBN в источниках расходятся, оба оставлены. Зависимость времени пролёта от амплитуды в оглавлении отдельным разделом не выделена.

```bibtex
@book{Berz1999,
  author    = {Berz, Martin},
  title     = {Modern Map Methods in Particle Beam Physics},
  publisher = {Academic Press},
  year      = {1999},
  volume    = {108},
  series    = {Advances in Imaging and Electron Physics},
  pages     = {1--318},
  doi       = {10.1016/S1076-5670(08)70227-1},
  isbn      = {978-0-12-014750-2},
  note      = {author page also gives ISBN 0-12-014750-5}
}
```

Глава 7.3 — дифференциально-алгебраическая нормальная форма: движение по кругу, частота зависит от амплитуды, есть хроматическая коррекция (7.2.3). Это сдвиг тюна, не \(\Lambda\). Прямой расчёт compaction любого порядка — соседняя статья, формулу времени пролёта из неё я не выписывал:

```bibtex
@article{Berz1999nim,
  author  = {Berz, M. and Erd{\'e}lyi, B. and Wan, W. and Ng, K. Y.},
  title   = {Differential algebraic determination of high-order off-energy closed orbits, chromaticities, and momentum compactions},
  journal = {Nucl. Instrum. Methods Phys. Res., Sect. A},
  volume  = {427},
  pages   = {310--314},
  year    = {1999},
  doi     = {10.1016/S0168-9002(98)01538-1}
}
```

**Статус Berz1999nim:** карточка INSPIRE проверена; текст формул **не проверен**.

**Куда ссылаться.** Главы 2 и 5: Berz1999 — почему в нормальной форме тюн зависит от \(J\), а время пролёта связано с \(\partial\nu/\partial E\) уже у Berg2007. Berz1999nim — если понадобится ссылка на \(\alpha_1\) как коэффициент карты, без ЭС-множителя.

---

## 9. Lee2019

**Статус:** формула проверена по тексту 4-го издания (PDF, copyright 2019, предисловие 10 июня 2018). INSPIRE даёт и 2018, и 2019. Более ранние издания не сравнивались.

```bibtex
@book{Lee2019,
  author    = {Lee, S. Y.},
  title     = {Accelerator Physics},
  edition   = {4},
  publisher = {World Scientific},
  address   = {Singapore},
  year      = {2019},
  pages     = {568},
  doi       = {10.1142/11111},
  isbn      = {978-981-3274-67-9},
  note      = {hardcover; ebook 978-981-3274-68-6; paperback 978-981-3274-78-5. INSPIRE also lists 2018}
}
```

Определение — по частоте, уравнение (3.11): \(\Delta\omega/\omega_0=-\eta(\delta)\,\delta\), \(\eta(\delta)=\eta_0+\eta_1\delta+\eta_2\delta^2+\cdots\). В линейном порядке, (2.164), \(\Delta f/f_0=-\Delta T/T_0=-\eta\delta\). Орбита \(R=R_0(1+\alpha_0\delta+\alpha_1\delta^2+\cdots)\), то есть \(\alpha_1\) — коэффициент при \(\delta^2\), как в заметке. Уравнение (3.12):

\[
\eta_0=\alpha_0-\frac{1}{\gamma_0^2},
\qquad
\eta_1=\frac{3\beta_0^2}{2\gamma_0^2}+\alpha_1-\alpha_0\eta_0.
\]

Это ровно запрошенная формула, и она относится к частоте. Коэффициент периода \(\tau=\eta_0\delta+\eta_1^{(T)}\delta^2\) больше на \(\eta_0^2\):

\[
\eta_1^{(T)}=\eta_1^{\mathrm{Lee}}+\eta_0^2
=\alpha_1-\frac{\alpha_0}{\gamma^2}+\frac{3\gamma^2-1}{2\gamma^4}.
\]

Второе выражение совпадает с \(\eta_1\) главы 1. Подстановка \(\eta_0=\alpha_0-1/\gamma^2\) в формулу Lee это даёт алгебраически; в книге периодный \(\eta_1\) так не назван.

Рядом, справочник Ng (разд. 2.3.11, текст OSTI 1002002): \(\Delta T/T_0=\eta\delta\) и

\[
\eta_1=\alpha_1-\frac{\eta_0}{\gamma_0^2}+\frac{3\beta_0^2}{2\gamma_0^2},
\]

что алгебраически равно периодному коэффициенту заметки. Привязка этого текста к номеру страницы конкретного издания Handbook **не проверена**.

```bibtex
@techreport{NgSlip,
  author = {Ng, K. Y.},
  title  = {Momentum Compaction and Phase Slip Factor},
  year   = {2011},
  note   = {section 2.3.11; OSTI 1002002. Handbook edition and pages not verified},
  url    = {https://www.osti.gov/servlets/purl/1002002}
}
```

Год 2011 у OSTI-файла в этом проходе по странице издателя Handbook не подтверждался: в библиографию он поставлен по идентификатору OSTI, **год файла не проверен отдельно**.

**Куда ссылаться.** Глава 1: Lee2019, уравнение (3.12), с явной оговоркой, что это \(\Delta\omega/\omega\). Периодный \(\eta_1\) заметки — NgSlip или Lee плюс \(\eta_0^2\). ЭС-множителя \(1+1/\gamma^2\) у Lee нет.

---

## 10. Дополнительно: slip factor и дисперсия в ЭС-кольце

### Grieser2017

**Статус:** проверено по JACoW, шапке PDF и IOP. В карточке JACoW «R. Hahn»; в PDF — «R. von Hahn».

```bibtex
@inproceedings{Grieser2017,
  author    = {Grieser, M. and von Hahn, R. and Vogel, S. and Wolf, A.},
  title     = {The Phase Slip Factor of the Electrostatic Cryogenic Storage Ring {CSR}},
  booktitle = {Proc. 8th International Particle Accelerator Conference (IPAC'17)},
  year      = {2017},
  pages     = {2255--2257},
  address   = {Copenhagen, Denmark},
  publisher = {JACoW Publishing, Geneva, Switzerland},
  doi       = {10.18429/JACoW-IPAC2017-TUPVA077},
  note      = {paper TUPVA077. Same text: J. Phys.: Conf. Ser. 874 (2017) 012049, doi 10.1088/1742-6596/874/1/012049},
  url       = {https://proceedings.jacow.org/ipac2017/papers/tupva077.pdf}
}
```

Нерелятивистское ЭС-кольцо с цилиндрическими дефлекторами. Их \(\eta\) — по частоте, \(\eta=\Delta f/f_0\big/(\Delta p/p_0)\). Итог, уравнение (14): \(\eta=1-2\alpha_p\). Для магнитного кольца в том же пределе \(\eta_m=1-\alpha_p\). Разница — ускорение и замедление в дефлекторе в зависимости от \(x\). При \(\gamma\to 1\) периодный коэффициент заметки \(\eta_0^{\mathrm{E}}=2\alpha_0-1\) равен \(-\eta\) Гризера. Измерение: \(\alpha_p=0.163\), \(\eta=0.674\); MAD8 даёт \(0.679\). Зависимость \(\eta(Q_x)\) разная на островах устойчивости и не такая, как в магнитном кольце. Второго порядка \(\eta_1\) нет: разложение оборвано на линейном члене, и скорость нерелятивистская.

**Куда ссылаться.** Глава 4, предел \(\eta_0^{\mathrm{E}}\). Глава 6: измеренный множитель 2 при \(\alpha_p\).

### Mane2008 и Mane2012

**Статус:** карточки INSPIRE и ScienceDirect проверены. Полный PDF за платным доступом. В открытом фрагменте ScienceDirect у Mane2012 виден множитель \(1+1/\gamma_0^2\) в правой части уравнения дисперсии и вклад потенциала в \(\Delta v/v\). Замкнутая формула \(\eta=(1+1/\gamma^2)\alpha-1/\gamma^2\) в этом фрагменте **не прочитана**.

```bibtex
@article{Mane2008,
  author  = {Mane, S. R.},
  title   = {Orbital dynamics in a storage ring with electrostatic bending},
  journal = {Nucl. Instrum. Methods Phys. Res., Sect. A},
  volume  = {596},
  number  = {3},
  pages   = {288--294},
  year    = {2008},
  doi     = {10.1016/j.nima.2008.08.087}
}

@article{Mane2012,
  author  = {Mane, S. R.},
  title   = {Orbital and spin motion in a storage ring with static electric and magnetic fields},
  journal = {Nucl. Instrum. Methods Phys. Res., Sect. A},
  volume  = {687},
  pages   = {40--50},
  year    = {2012},
  doi     = {10.1016/j.nima.2012.05.098}
}
```

Mane2008, по аннотации: линейный градиент, дисперсия, momentum compaction и натуральная хроматичность для радиального электрического поля; точное решение для слабой фокусировки. Поле модели — \(E_r\propto 1/r\). Mane2012 добавляет сильную фокусировку, синхротронные колебания и спин, в том числе для чисто электрического кольца. В открытом фрагменте дисперсия электрического изгиба содержит \((1+1/\gamma_0^2)\), а \(\Delta v/v_0\) содержит \(\Delta\Phi\).

**Куда ссылаться.** Главы 4 и 5: (а) множитель \(1+1/\gamma^2\) в источнике дисперсии, не в готовом \(\eta_1\); (б) хроматичность ЭС-кольца заявлена в аннотации 2008 года, формула в этом проходе не выписана. Для (в) фрагмент даёт скорость через потенциал, без ступеньки жёсткого края.

### Conte2012

**Статус:** текст arXiv проверен по HTML. Журнальная публикация **не найдена**.

```bibtex
@unpublished{Conte2012,
  author = {Conte, Mario},
  title  = {Electrostatic Storage Ring},
  year   = {2012},
  note   = {arXiv:1202.1697},
  url    = {https://arxiv.org/abs/1202.1697}
}
```

Для цилиндра \(n=1\) горизонтальный тюн сектора \(\nu_H^2=2-\beta_0^2=1+1/\gamma_0^2\), вертикальный \(\nu_V^2=0\). Смещение по радиусу меняет полную энергию: \(\Delta\gamma/\gamma_0=-\beta_0^2 x/\rho\). В матрице изгиба появляется элемент \(M_{31}=-1/\rho\), и её определитель на одном секторе не равен единице, потому что энергия на дуге не сохраняется; за полный период бетатронного колебания баланс восстанавливается.

**Куда ссылаться.** Глава 4: кинетическая энергия на внешнем радиусе ниже. Глава 5: \(K_x\) цилиндра содержит \(1+1/\gamma^2\). Это заметка, не рецензированная статья.

### Zyuzin2011 и Senichev2011

**Статус:** страницы проверены по JACoW. DOI у IPAC 2011 нет.

```bibtex
@inproceedings{Zyuzin2011,
  author    = {Zyuzin, D. and Maier, R. and Senichev, Yu.},
  title     = {High Order Non-linear Motion in Electrostatic Rings},
  booktitle = {Proc. 2nd International Particle Accelerator Conference (IPAC'11)},
  year      = {2011},
  pages     = {2172--2174},
  address   = {San Sebasti{\'a}n, Spain},
  publisher = {JACoW Publishing, Geneva, Switzerland},
  isbn      = {978-92-9083-366-6},
  note      = {paper WEPC066},
  url       = {https://proceedings.jacow.org/IPAC2011/papers/wepc066.pdf}
}

@inproceedings{Senichev2011,
  author    = {Senichev, Yu. and Lehrach, A. and Maier, R. and Zyuzin, D.},
  title     = {The Spin Aberration of Polarized Beam in Electrostatic Rings},
  booktitle = {Proc. 2nd International Particle Accelerator Conference (IPAC'11)},
  year      = {2011},
  pages     = {2175--2177},
  address   = {San Sebasti{\'a}n, Spain},
  publisher = {JACoW Publishing, Geneva, Switzerland},
  note      = {paper WEPC067},
  url       = {https://proceedings.jacow.org/IPAC2011/papers/wepc067.pdf}
}
```

Zyuzin2011: у частицы со смещением свой equilibrium energy level внутри дефлектора; кинетическая энергия переходит в потенциальную и обратно. Край подавляет этот обмен множителем \(k\). Если тюн близок к единице, уровень энергии меняется за время порядка синхротронного периода. Senichev2011: у внеосевой частицы equilibrium energy level сдвинут относительно магического импульса, и ВЧ этот сдвиг не снимает. Оба текста — качественные и по трекингу COSY; готового \(\eta_1^{\mathrm{E}}\) нет.

**Куда ссылаться.** Глава 4 и глава 6: «равновесный уровень энергии» в дефлекторе до формулы (10) работы 2013 года.

### Senichev2014 и Melnikov2022

**Статус:** карточки JACoW проверены. Коэффициенты формулы (1) в Melnikov2022 из текстового слоя PDF надёжно не сняты.

```bibtex
@inproceedings{Senichev2014,
  author    = {Senichev, Y. and Andrianov, S. N. and Ivanov, A. N. and Lehrach, A. and Maier, R. and Zyuzin, D.},
  title     = {Spin Tune Decoherence in Multipole Fields},
  booktitle = {Proc. 5th International Particle Accelerator Conference (IPAC'14)},
  year      = {2014},
  pages     = {3017--3019},
  address   = {Dresden, Germany},
  publisher = {JACoW Publishing, Geneva, Switzerland},
  doi       = {10.18429/JACoW-IPAC2014-THPRO062},
  note      = {paper THPRO062},
  url       = {https://proceedings.jacow.org/ipac2014/papers/thpro062.pdf}
}

@inproceedings{Melnikov2022,
  author    = {Melnikov, A. A. and Aksentyev, A. E. and Senichev, Y. and Syresin, E.},
  title     = {Investigation of Polarized Proton Spin Coherence Time at Storage Rings},
  booktitle = {Proc. 13th International Particle Accelerator Conference (IPAC'22)},
  year      = {2022},
  pages     = {1832--1834},
  address   = {Bangkok, Thailand},
  publisher = {JACoW Publishing, Geneva, Switzerland},
  doi       = {10.18429/JACoW-IPAC2022-WEPOPT005},
  note      = {paper WEPOPT005},
  url       = {https://jacow.org/ipac2022/papers/wepopt005.pdf}
}
```

Senichev2014 продолжает Senichev2013 на высшие мультиполи: мультиполь сдвигает средний радиус орбиты и тем самым средний уровень энергии и спиновый тюн. Линейный член гасится ВЧ, остаётся член с \(\Delta\gamma^2\) и бетатронное удлинение. Melnikov2022 прямо цитирует Senichev2013 как источник \(\Delta\delta_{\mathrm{eq}}\) и пишет \(\Delta\nu_s=G\Delta\gamma_{\mathrm{eq}}\). Оптимизируются \(\xi_x\), \(\xi_y\) и \(\alpha_1\), тремя семействами секступолей. Коэффициенты их уравнения (1) в этом проходе **не проверены**.

**Куда ссылаться.** Глава 6: цепочка от Senichev2013 к спиновому тюну. Главы 2 и 5: секступоль двигает и \(\xi\), и \(\alpha_1\).

### Talman2017

**Статус:** библиография проверена по IOP. Формулы slip factor из текста книги **не извлечены**.

```bibtex
@book{Talman2017,
  author    = {Talman, Richard},
  title     = {The Electric Dipole Moment Challenge},
  publisher = {Morgan \& Claypool Publishers},
  year      = {2017},
  pages     = {183},
  doi       = {10.1088/978-1-6817-4509-1},
  isbn      = {978-1-6817-4508-4},
  note      = {print; online ISBN 978-1-6817-4509-1. Chapter All-electric particle dynamics, pp. 6-1--6-13}
}
```

Глава «All-electric particle dynamics» есть в оглавлении IOP. Что именно там написано про \(\eta\) и хроматичность, в этом проходе не проверялось.

**Куда ссылаться.** Глава 6 только как обзор полностью электрического кольца, после того как формула будет сверена по книге. Сейчас на \(\eta_0^{\mathrm{E}}\) не опираться.

### Abusaif2021

**Статус:** DOI, ISBN, номер отчёта и авторский список проверены по титулу PDF. Число страниц 246 взято из вторичной ссылки (слайды N. Nikolaev) и по PDF не пересчитано: **объём не проверен**.

```bibtex
@techreport{Abusaif2021,
  author      = {Abusaif, F. and others},
  title       = {Storage ring to search for electric dipole moments of charged particles: Feasibility study},
  institution = {CERN},
  year        = {2021},
  number      = {CERN-2021-003},
  doi         = {10.23731/CYRM-2021-003},
  isbn        = {978-92-9083-606-3},
  eprint      = {1912.07881},
  archivePrefix = {arXiv},
  note        = {CERN Yellow Reports: Monographs, vol. 3; CPEDM Collaboration. PDF ISBN 978-92-9083-607-0. Page count not verified}
}
```

Среди авторов — Y. Senichev, E. Valetov, M. Berz, R. Talman. Отчёт задаёт контекст: спин вдоль скорости дольше 15 минут, радиальное электрическое поле, прототип PTR. Алгебры \(\eta_1^{\mathrm{E}}\) на титуле и в аннотации нет; по тексту отчёта формулы slip factor в этом проходе не выписывались.

**Куда ссылаться.** Глава 6: зачем нужен \(\gamma_{\mathrm{eff}}\), не вместо вывода \(\delta_{\mathrm{eq}}\).

---

## Что сверять с заметкой

1. **Глава 1, \(\eta_1\).** У Lee, (3.12), запрошенная формула \(\eta_1=3\beta^2/(2\gamma^2)+\alpha_1-\alpha_0\eta_0\) есть, и она определена через частоту \(\Delta\omega/\omega_0=-(\eta_0+\eta_1\delta)\delta\). Коэффициент периода в заметке, \(\alpha_1-\alpha_0/\gamma^2+(3\gamma^2-1)/(2\gamma^4)\), равен \(\eta_1\) из Lee плюс \(\eta_0^2\). Та же периодная комбинация получается из рекуррентности Ng.

2. **Глава 6, формула коллег.** В Senichev2013 напечатано \(\Delta\delta_{\mathrm{eq}}=(1/\eta_0)[\cdots]\), а в заметке перед скобкой стоит минус. Коэффициент при \(\delta^2\) у них \(\alpha_1-\alpha_0/\gamma^2+1/\gamma^4\): не хватает второго порядка \(\Delta v/v\), разница с \(\eta_1\) заметки равна \(3\beta^2/(2\gamma^2)\). Символа \(\gamma_{\mathrm{eff}}\) нет; связь \(\Delta\gamma/\gamma=\beta^2\Delta\delta\) видна по множителю \(2/\beta^2\) в (18). Множитель \(F_v\) в (25) есть.

3. **Главы 2 и 5, \(\Lambda\) и \(\xi\).** Berg2007 и аннотация Shoji2005 дают ту же связь времени пролёта с хроматичностью. Предфактор \(-2\pi/C\) и геометрическое действие заметки у Берга другие: его \(J_n\) не делён на импульс. У Senichev бетатронный кусок — геометрический \((\Delta L/L)_\beta\sim\varepsilon\nu\), секступоли добавлены отдельным членом, не через \(\xi=0\Rightarrow\Lambda=0\).

4. **Глава 4, энергия в дефлекторе.** ValetovBerz2020 и Valetov2020 (том 213): на жёстком крае кинетическая энергия обязана скакнуть вместе с потенциалом, иначе аберрации Вольника 1965 не симплектичны. Conte2012: \(\Delta\gamma/\gamma_0=-\beta_0^2 x/\rho\). Grieser2017 проверяет нерелятивистский предел \(\eta_0^{\mathrm{E}}\to 2\alpha_0-1\) (с обратным знаком, потому что их \(\eta\) по частоте). Релятивистский множитель \(1+1/\gamma^2\) в источнике дисперсии виден во фрагменте Mane2012; готовое \(\eta_1^{\mathrm{E}}\) оттуда не снято.

5. **Авторство пункта 2.** Статья с названием про frozen spin — Aksentyev и Senichev, IPAC 2019, MOPTS012. Valetov — диссертация 2017 года и главы AIEP 213.
