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

### `ndarray`: N次元配列クラス

<https://numpy.org/doc/stable/reference/arrays.ndarray.html>

2D行列に特殊化したものとして
[`numpy.matrix`](https://numpy.org/doc/stable/reference/generated/numpy.matrix.html)
もあったが2018年ごろから非推奨。

`data.frame`/`tibble` のようなものがほしいときは
[`pandas.DataFrame`]({{< relref "pandas.md" >}}) の出番。

<https://numpy.org/doc/stable/reference/routines.matlib.html>

### Universal functions

<https://numpy.org/doc/stable/reference/ufuncs.html#available-ufuncs>

`ndarray` に対して効率よく element-wise な演算をする関数


### Routine

<https://numpy.org/doc/stable/reference/routines.html>

`ndarray` を作ったり操作したりする関数全般


## `SciPy`: 高度な科学技術計算ライブラリ

<https://docs.scipy.org/doc/scipy/reference/>

### `scipy.stats`

- <https://docs.scipy.org/doc/scipy/reference/stats.html>
- <https://docs.scipy.org/doc/scipy/tutorial/stats.html>
