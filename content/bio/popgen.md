+++
title = 'Population Genetics'
tags = ["genetics"]
[menu.main]
  parent = "bio"
+++

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/0974707759/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41Y1PqrWh5L._SX180_.jpg" alt="Coalescent Theory: An Introduction" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/0763757373/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51IjUPvsSVL._SX180_.jpg" alt="Genetics of Populations" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/0198502311/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41%2Bu1SOjvcL._SX180_.jpg" alt="Evolutionary Genetics" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/1932846123/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41f8rXshBRL._SX150_.jpg" alt="An Introduction to Population Genetics Theory" /></a>
## Model

- **Wright-Fisher model**
  -   ランダム交配
  -   世代重複なし
  -   **集団サイズは有限** のNで一定
      （ここがHardy-Weinbergと違う。そしてこのことによる遺伝的浮動に興味がある）
  -   係数sの変異の固定確率 $\frac {1 - e^{-2s}} {1 - e^{-4Ns}}$

- **Moran model**
  -   世代重複あり(1個体が複製して、死ぬ1個体を置き換える)。
  -   Wright-Fisherに合わせるならNステップを1世代と考える。
  -   繁殖成功の標準偏差(ヘテロ接合頻度の減少速度＝遺伝的浮動の強さ)はWright-Fisherの倍。
      有効集団サイズが半分。
  -   増殖率rの変異の固定確率 $\frac {1 - 1/r} {1 - 1/r^N}$

## Statistics

**The unfolded site-frequency counts** $\xi_i$
:   派生型のアリルが *i* 個、祖先型のアリルが *n-i* 個である変異サイトの数

**The folded site-frequency counts** $\eta_i$
:   どっちが祖先型か不明な状態。
    片方のアリルが *i* 個、もう片方のアリルが *n-i* 個である変異サイトの数

    $\eta_i = \frac{\xi_i + \xi _{n - i}}{1 + \delta _{i, n - i}}$

**The number of segregating (polymorphic) sites** $S$
:   配列セットの中で、多型のあるサイトの数

**Nucleotide diversity / 塩基多様度** $\pi$
:   整列済み配列セットについてペアワイズで塩基の異なるサイト数を数え、
    ペアあたりで平均したもの。
    多型サイト数が同じでも、アリル頻度が均等なほど大きくなり、
    少数のアリルが優占してたりすると小さくなる。

**Population mutation rate / 集団突然変異率** $\theta$
:   二倍体常染色体なら $4N_e\mu$、
    二倍体X染色体なら $3N_e\mu$、
    一倍体なら $2N_e\mu$。
    直接測定することができないためほかの値から推定する。

    **Watterson (1975)**: $S$ から推定

    <div>\[\begin{split}
    a_1  &= 1 + \frac 1 2 + \frac 1 3 + ... + \frac 1{n-1}\\
    \mbox E[S] &= \theta L a_1\\
    \theta_w &= \frac{S}{L a_1}
\end{split}\]</div>

    **Tajima (1983)**: $\pi$ から推定

    <div>\[\begin{split}
    \mbox E[\pi] &= \theta L\\
    \theta_\pi   &= \frac \pi L
\end{split}\]</div>

## Selection

### 表現型の頻度分布に着目

**directional selection / 方向性選択**
:   形質値の頻度分布が一方向的に動くように働く選択。
    その方向に表現型変化をもたらす変異に対しては
    positive selectionがかかり、
    逆方向の表現型変化をもたらす変異に対しては
    purifying selectionがかかる。

**stabilizing selection / 安定化選択**
:   有利で頻度の高い形質値を中心として、
    頻度分布が広がらないように働く選択。
    disrupbtiveの逆。
    結果的に配列に対して purifying selection がかかることは多いと思われるが、
    より安定してその形質値を実現できるようなアリルに対して
    positive selectionがかかることもあるだろう。

**disruptive selection / 分断化選択**
:   中間的な表現型が不利で、形質値の頻度分布に谷ができるように働く選択。
    stabilizingの逆。
    表現型可塑性や表現型多型で対処される場合もあり、
    必ずしもbalancing selectionや種分化をもたらさない。

### 遺伝子型頻度に着目

**positive selection**
:   有益アリルの頻度を上げるように働く選択

**negative selection = purifying selection / 純化選択**
:   有害アリルを集団から取り除くように働く選択

**background selection** ([Charlesworth et al. 1993](http://www.ncbi.nlm.nih.gov/pubmed/8375663)):
:   有害変異に対する purifying selection によって
    近傍配列まで遺伝的多様度が減少する

### どっちでも

**balancing selection / 平衡選択**
:   多型を維持するように働く選択\
    **heterozygote advantage**\
    **temporally varying selection**\
    **spatially varying selection**\
    **frequency-dependent selection / 頻度依存選択**\
    antagonistic pleiotropy, disassortative mating, self-incompatibility
