+++
title = 'NumPy, SciPy'
subtitle = "数値配列演算・科学技術計算ライブラリ"
tags = ["python"]
[menu.main]
  parent = "python"
+++

<https://www.scipy.org/>

## Installation

<https://www.scipy.org/install.html>

```sh
pip install numpy scipy
```

BLAS/LAPACK関連の確認:
```py
import numpy as np
np.show_config()
```

## `NumPy`

数値計算・配列演算ライブラリ

<https://numpy.org/doc/stable/reference/>

### [`ndarray`](https://numpy.org/doc/stable/reference/arrays.ndarray.html): N次元配列クラス

2D行列に特殊化したものとして
[`numpy.matrix`](https://numpy.org/doc/stable/reference/generated/numpy.matrix.html)
もあったが2018年ごろから非推奨。

`data.frame`/`tibble` のようなものがほしいときは
[`pandas.DataFrame`]({{< relref "pandas.md" >}}) の出番。

<https://numpy.org/doc/stable/reference/routines.matlib.html>

### [Universal functions](https://numpy.org/doc/stable/reference/ufuncs.html#available-ufuncs)

`ndarray` に対して効率よく element-wise な演算をする関数


### [Routine](https://numpy.org/doc/stable/reference/routines.html)

`ndarray` を作ったり操作したりする関数全般


## `SciPy`

高度な科学技術計算ライブラリ

<https://docs.scipy.org/doc/scipy/reference/>

### [`scipy.stats`](https://docs.scipy.org/doc/scipy/reference/stats.html)

<https://docs.scipy.org/doc/scipy/tutorial/stats.html>


## 書籍

<a href="https://www.amazon.co.jp/dp/487311845X/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=72a416f5d10a9e84aaab4b3ee9613329&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311845X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=487311845X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873118417/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=6b1a04ec880b6c730bd6e80273e30e9c&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873118417&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873118417" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873117488/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=2181a50362009e68f507d44fc38716b4&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
