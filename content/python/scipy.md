+++
title = 'NumPy, SciPy'
subtitle = "数値配列演算・科学技術計算ライブラリ"
tags = ["python"]
[menu.main]
  parent = "python"
+++

<https://www.scipy.org/>

## Installation

<https://www.scipy.org/install.html>\
<https://www.scipy.org/scipylib/download.html>\
<https://www.scipy.org/scipylib/building/>\
<https://docs.scipy.org/doc/numpy/user/install.html>

1.  [Homebrew]({{< relref "mac/homebrew.md" >}}) などで `gfortran` をインストール:

        % brew install gcc

    Ubuntuなら:

        % sudo apt-get install build-essential gfortran
        % sudo apt-get install python3-dev
        % sudo apt-get install libatlas-base-dev

2.  `NumPy` をインストール:

        % pip install numpy

3.  `SciPy` をインストール:

        % CC=clang CXX=clang++ FFLAGS=-ff2c pip install scipy

### test

`nose` というパッケージが必要:

    >>> import numpy
    >>> numpy.test()

BLAS/LAPACK関連の確認:

    >>> import numpy
    >>> numpy.show_config()

## `NumPy`

数値計算・配列演算ライブラリ

<https://pypi.python.org/pypi/numpy>\
<http://www.numpy.org/>\
<https://docs.scipy.org/doc/numpy/reference/>

### `ndarray` クラス

<https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html>

N次元配列

------------------------------------------------------------------------

`numpy.matrix` クラスはそれを2D行列に特殊化したもの

<https://docs.scipy.org/doc/numpy/reference/routines.matlib.html>

### Universal functions (`ufunc`)

<https://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs>

`ndarray` に対して効率よく element-wise な演算をする関数

### Routine

<https://docs.scipy.org/doc/numpy/reference/routines.html>

`ndarray` を作ったり操作したりする関数全般

<https://docs.scipy.org/doc/numpy/reference/routines.array-creation.html>

<https://docs.scipy.org/doc/numpy/reference/routines.array-manipulation.html>

<https://docs.scipy.org/doc/numpy/reference/routines.math.html>

<https://docs.scipy.org/doc/numpy/reference/routines.random.html>

## `SciPy`

高度な科学技術計算ライブラリ

<https://pypi.python.org/pypi/scipy>\
<https://www.scipy.org/scipylib/index.html>\
<https://docs.scipy.org/doc/scipy/reference/>

### `scipy.stats`

<https://docs.scipy.org/doc/scipy/reference/tutorial/stats.html>

<https://docs.scipy.org/doc/scipy/reference/stats.html>
