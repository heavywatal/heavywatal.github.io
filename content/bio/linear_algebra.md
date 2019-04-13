+++
title = '線形代数'
tags = ["math"]
[menu.main]
  parent = "bio"
+++

## 用語

**正方行列 (square)**
:   行と列の数が等しい行列。

**上三角行列 (upper triangular)**, **下三角行列(lower triangular)**
:   対角線の左下、右上が0となる正方行列。
    対角成分の積が行列式になる。

**対角行列 (diagonal)**
:   対角成分以外が0の正方行列。上三角かつ下三角。
    軸方向の伸縮だけで歪まない写像。
    べき乗が対角成分それぞれのべき乗だけで計算できる。

**単位行列 (identity)** $I, E$
:   対角成分が全て1の対角行列。

**正則行列 (invertible, non-singular, non-degenerate, regular)**
:   逆行列を持つ。行列式が0じゃない。固有値0を持たない。

------------------------------------------------------------------------

**行列式 (determinant)** $\text{det}\, A, |A|$
:   正方行列による変換の体積拡大率。
    行列式が0 ⇔ 写像がぺちゃんこになる ⇔ 逆行列が存在しない。
    固有値の積と等しい。

**跡 (trace)** $\text{tr}\, A$
:   正方行列の対角成分の和。
    固有値の和と等しい。

**核 (kernel)** $\text{Ker}\, A$
:   $A\mathbf{x} = \mathbf{o}$ で原点に移るような
    $\mathbf{x}$ の集合。 核が原点だけ(0次元) ⇔ ランク＝元の次元数 ⇔ 写像は **単射**。

**像 (image)** $\text{Im}\, A$
:   $\mathbf{x}$ を目一杯いろいろ動かしたときの
    $\mathbf{y} = A\mathbf{x}$ の集合。
    像が行き先の全空間 ⇔ ランク＝行き先の次元数 ⇔ 写像は **全射**。

**ランク (rank)** $\text{rank}\, A$
:   像の次元数。

## 行列のべき乗

行列は写像。行列のべき乗は写像の繰り返し。

ベクトル $\mathbf{x}$ に正方行列 $A$ を
$t$ 回かけたらどうなるか知りたい。

<div>\[
\mathbf{x}(t) = A\mathbf{x}(t-1) = A^t\mathbf{x}(0)
\]</div>

そのまま行列計算をするのではなく、適当な正則行列で
$\mathbf{x}(t) = P\mathbf{y}(t)$ という変数変換をしてみると

<div>\[\begin{aligned}
\mathbf{y}(t) &= P^{-1}\mathbf{x}(t) \\
              &= P^{-1}A\mathbf{x}(t-1) \\
              &= P^{-1}AP\mathbf{y}(t-1) \\
              &= (P^{-1}AP)^t\mathbf{y}(0) \\
              &= \Lambda^t\mathbf{y}(0) \\
\mathbf{x}(t) &= P\mathbf{y}(t) \\
              &= P\Lambda^t\mathbf{y}(0) \\
              &= P\Lambda^tP^{-1}\mathbf{x}(0)
\end{aligned}\]</div>

このとき $\Lambda = P^{-1}AP$ が対角行列になってくれてると
$t$ 乗する計算がすごく楽チン。

<div>\[\begin{aligned}
\Lambda^t &= \text{diag}(\lambda _1, ..., \lambda _n)^t \\
          &= \text{diag}(\lambda _1^t, ..., \lambda _n^t)
\end{aligned}\]</div>

この **対角化 (diagonalization)** をもたらす変換行列 $P$ とはどういうものか

<div>\[\begin{aligned}
P^{-1}AP &= \text{diag}(\lambda _1, ..., \lambda _n) \\
      AP &= P \text{diag}(\lambda _1, ..., \lambda _n)
\end{aligned}\]</div>

$P = (\mathbf{p}_1, ..., \mathbf{p}_n)$ として列ごとに見ると

<div>\[\begin{aligned}
A\mathbf{p}_1 &= \lambda_1 \mathbf{p}_1 \\
\vdots \\
A\mathbf{p}_n &= \lambda_n \mathbf{p}_n
\end{aligned}\]</div>

$A$ をかけても長さが変わるだけで方向は変わらない。
この伸縮率 $\lambda$ が **固有値 (eigenvalue)** で、それぞれに対応する
$\mathbf{o}$ でない $\mathbf{p}$ が **固有ベクトル (eigenvector)**。

つまり変換行列 $P$ は $A$ の固有ベクトルを並べたもので、
$\Lambda$ は対角成分に $A$ の固有値を並べたもの。
そうするとさっきの $\mathbf{x} = P\mathbf{y}$ は、
$\mathbf{x}$ を $A$ の固有ベクトルの線形結合として表し、
$A$ をかけても方向が変わらないように変数変換しておくということに相当する。

<div>\[\begin{aligned}
A^t\mathbf{x} &= A^t P \mathbf{y} \\
              &= A^t (y_1\mathbf{p_1} + ... + y_n\mathbf{p_n}) \\
              &= y_1 A^t \mathbf{p_1} + ... + y_n A^t \mathbf{p_n} \\
              &= y_1 \lambda _1^t \mathbf{p_1} + ... + y_n \lambda _n^t \mathbf{p_n} \\
              &= \lambda_k ^t (
                   y_1 \left(\frac {\lambda _1} {\lambda _k}\right) ^t \mathbf{p_1}
                   + ... + y_k \mathbf{p_k}
                   + ... + y_n \left(\frac {\lambda _n} {\lambda _k}\right) ^t \mathbf{p_n}) \\
              &\sim y_k \lambda_k ^t \mathbf{p_k}
\end{aligned}\]</div>

$t$ が大きくなるにつれて最大の固有値 $\lambda_k$
に対応する固有ベクトル $\mathbf{p_k}$ の向きに近づいていく。
その極限には行かないにしても、固有値の大きな固有ベクトルの方向に寄っていく傾向があるってこと。

固有値が重解を含んでいてもその重複と同じ分だけ固有ベクトルが取れれば
(代数的重複度＝幾何的重複度ならば) 対角化可能。
対角化できない正方行列でも $P^{-1}AP = J$ となる $P$ を見つけて
**Jordan標準形** まで持っていくことは可能。
$J^t$ の計算は $\Lambda^t$ ほどじゃないにせよそれなりに楽。

## 参考文献

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4274065782/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51QTY7RSFRL._SX160_.jpg" alt="プログラミングのための線形代数" /></a>
