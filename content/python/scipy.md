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

`wheel` のおかげで `pip` いっぱつ

```sh
% pip install numpy scipy
```

[Anaconda]({{< relref "install.md" >}}) なら最初からMKL付きで入ってる。

### test

`nose` というパッケージが必要:
```py
>>> import numpy
>>> numpy.test()
```

BLAS/LAPACK関連の確認:
```py
>>> import numpy
>>> numpy.show_config()
```

## `NumPy`

数値計算・配列演算ライブラリ

<https://docs.scipy.org/doc/numpy/reference/>

### [`ndarray`](https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html): N次元配列クラス

それを2D行列に特殊化したものが
[`numpy.matrix`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.matrix.html#numpy.matrix)
で、Rでの `matrix` に相当。
`data.frame`/`tibble` のようなものがほしいときは
[`pandas.DataFrame`]({{< relref "pandas.md" >}}) の出番。

<https://docs.scipy.org/doc/numpy/reference/routines.matlib.html>

### [Universal functions](https://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs)

`ndarray` に対して効率よく element-wise な演算をする関数


### [Routine](https://docs.scipy.org/doc/numpy/reference/routines.html)

`ndarray` を作ったり操作したりする関数全般

- <https://docs.scipy.org/doc/numpy/reference/routines.array-creation.html>
- <https://docs.scipy.org/doc/numpy/reference/routines.array-manipulation.html>
- <https://docs.scipy.org/doc/numpy/reference/routines.math.html>
- <https://docs.scipy.org/doc/numpy/reference/routines.random.html>



## `SciPy`

高度な科学技術計算ライブラリ

<https://docs.scipy.org/doc/scipy/reference/>

### [`scipy.stats`](https://docs.scipy.org/doc/scipy/reference/stats.html)

<https://docs.scipy.org/doc/scipy/reference/tutorial/stats.html>



## 書籍

<a href="https://www.amazon.co.jp/dp/479738946X/ref=as_li_ss_il?ie=UTF8&qid=1485612008&sr=8-6&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=5ea5e48ecc83b9439f21406b6f57c062" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=479738946X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=479738946X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/IPython%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9%E3%82%AF%E3%83%83%E3%82%AF%E3%83%96%E3%83%83%E3%82%AF-%E5%AF%BE%E8%A9%B1%E5%9E%8B%E3%82%B3%E3%83%B3%E3%83%94%E3%83%A5%E3%83%BC%E3%83%86%E3%82%A3%E3%83%B3%E3%82%B0%E3%81%A8%E5%8F%AF%E8%A6%96%E5%8C%96%E3%81%AE%E3%81%9F%E3%82%81%E3%81%AE%E3%83%AC%E3%82%B7%E3%83%94%E9%9B%86-Cyrille-Rossant/dp/4873117488/ref=as_li_ss_il?_encoding=UTF8&psc=1&refRID=X16VFSS3W75RMTG7VGCH&linkCode=li3&tag=heavywatal-22&linkId=b79e2290571289b02621392257a4ac1c" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/Python%E3%81%AB%E3%82%88%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E5%88%86%E6%9E%90%E5%85%A5%E9%96%80-_NumPy%E3%80%81pandas%E3%82%92%E4%BD%BF%E3%81%A3%E3%81%9F%E3%83%87%E3%83%BC%E3%82%BF%E5%87%A6%E7%90%86-Wes-McKinney/dp/4873116554/ref=as_li_ss_il?s=books&ie=UTF8&qid=1485612076&sr=1-16&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=f024f5e24f16c38402149f97591c8aab" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873116554&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873116554" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/%E3%82%BC%E3%83%AD%E3%81%8B%E3%82%89%E4%BD%9C%E3%82%8BDeep-Learning-Python%E3%81%A7%E5%AD%A6%E3%81%B6%E3%83%87%E3%82%A3%E3%83%BC%E3%83%97%E3%83%A9%E3%83%BC%E3%83%8B%E3%83%B3%E3%82%B0%E3%81%AE%E7%90%86%E8%AB%96%E3%81%A8%E5%AE%9F%E8%A3%85-%E6%96%8E%E8%97%A4-%E5%BA%B7%E6%AF%85/dp/4873117585/ref=as_li_ss_il?s=books&ie=UTF8&qid=1485612076&sr=1-3&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=cde2e18b945af3ead43dc15e51b00af6" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117585&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117585" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
