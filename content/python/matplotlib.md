+++
title = 'matplotlib + seaborn'
subtitle = "Pythonでグラフ描画"
tags = ["python", "graph"]
[menu.main]
  parent = "python"
+++

<http://matplotlib.org/>

行列を含む数値計算やグラフ描画を行うためのソフトウェアとしては
有償の [MATLAB](http://www.mathworks.com/products/matlab/) が広く利用されているが、
Python にいくつかのモジュールを導入することで
同等かそれ以上のことを無料で実現することができる。

-   `numpy`: 配列演算や基本的な数学関数
-   `scipy`: 高度な科学技術計算
-   `matplotlib`: グラフ描画
-   `ipython`: 高機能なシェル環境
-   `pylab`: これらを組み合わせてMATLABっぽいインターフェイスにする(不要)

直接 `matplotlib` を触るのは大変なので、
`seaborn` というラッパーを介して使う。
ただし `matplotlib` を全く知らずに `seaborn` を使うのは無理っぽい。
[ggplot2]({{< relref "rstats/ggplot2.md" >}}) は `grid` を知らなくても使えるのに。。。

## 基本

-   <http://matplotlib.org/faq/usage_faq.html>
-   <http://matplotlib.org/users/pyplot_tutorial.html>
-   <http://matplotlib.org/faq/howto_faq.html>
-   <http://www.scipy-lectures.org/intro/matplotlib/matplotlib.html>

```py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

iris = sns.load_dataset('iris')

fig, ax = plt.subplots(figsize=(7, 7))
sns.regplot('sepal_width', 'sepal_length', data=iris, ax=ax)

plt.show()  # opens a window
fig.savefig('example.png')
```


### Figure, Axes

<http://matplotlib.org/faq/usage_faq.html#general-concepts>

<http://matplotlib.org/api/figure_api.html>

<http://matplotlib.org/api/axes_api.html>

## Style

<http://stanford.edu/~mwaskom/software/seaborn/tutorial/aesthetics.html>

背景色や補助線などの設定。

`sns.set_style(style, rc=None)`

`style`:
:   {`darkgrid`, `whitegrid`, `dark`, `white`, `ticks`}

`rc`:
:   see below

```python
>>> sns.axes_style(style='darkgrid', rc=None)
{'axes.axisbelow': True,
 'axes.edgecolor': 'white',
 'axes.facecolor': '#EAEAF2',
 'axes.grid': True,
 'axes.labelcolor': '.15',
 'axes.linewidth': 0,
 'figure.facecolor': 'white',
 'font.family': ['sans-serif'],
 'font.sans-serif': ['Arial',
                     'Liberation Sans',
                     'Bitstream Vera Sans',
                     'sans-serif'],
 'grid.color': 'white',
 'grid.linestyle': '-',
 'image.cmap': 'Greys',
 'legend.frameon': False,
 'legend.numpoints': 1,
 'legend.scatterpoints': 1,
 'lines.solid_capstyle': 'round',
 'text.color': '.15',
 'xtick.color': '.15',
 'xtick.direction': 'out',
 'xtick.major.size': 0,
 'xtick.minor.size': 0,
 'ytick.color': '.15',
 'ytick.direction': 'out',
 'ytick.major.size': 0,
 'ytick.minor.size': 0}
```

`with` 文でも使える:

    with sns.axes_style('white'):
        # plot

`sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)`
:   指定した枠線を消す。
    生の `matplotlib` だと `ax.spines['top'].set_visible(False)`

## Context

<http://stanford.edu/~mwaskom/software/seaborn/tutorial/aesthetics.html#scaling-plot-elements-with-plotting-context-and-set-context>

ラベルや点・線などのスケール調整。

`sns.set_context(context, font_scale=1, rc=None)`

`context`
:   {`notebook`: 1.0, `paper`: 0.8, `talk`, 1.3, `poster`, 1.6}

`rc`:
:   see below

```python
>>> sns.plotting_context(context='notebook', font_scale=1, rc=None)
{'axes.labelsize': 11,
 'axes.titlesize': 12,
 'figure.figsize': array([ 8. ,  5.5]),
 'grid.linewidth': 1,
 'legend.fontsize': 10,
 'lines.linewidth': 1.75,
 'lines.markeredgewidth': 0,
 'lines.markersize': 7,
 'patch.linewidth': 0.3,
 'xtick.labelsize': 10,
 'xtick.major.pad': 7,
 'xtick.major.width': 1,
 'xtick.minor.width': 0.5,
 'ytick.labelsize': 10,
 'ytick.major.pad': 7,
 'ytick.major.width': 1,
 'ytick.minor.width': 0.5}
```

`with` 文でも使える:

    with sns.plotting_context('talk', font_scale=1.2):
        # plot

## Color

<http://stanford.edu/~mwaskom/software/seaborn/tutorial/color_palettes.html>

個々の描画関数の `palette` や `cmap` に指定するか、
`with sns.color_palette():` のブロック内で描画するか、
`sns.set_palette()` でデフォルトを指定する。

-   [ColorBrewer](http://colorbrewer2.org/) の名前で指定。
    e.g., `sns.color_palette('RdBu', n_colors=7)`

    `_r` をつけると逆に、`_d` をつけると暗めになる。

-   [xkcd](http://xkcd.com/color/rgb/) の名前リストを
    `sns.xkcd_palette()` に渡す
-   `cubehelix_palette()` はgrayscaleでもいい感じで印刷できる
-   自分で作る:
    `sns.hls_palette()`, `sns.husl_palette()`,
    `sns.light_palette()`, `sns.dark_palette()`,
    `sns.diverging_palette()`

## Text, Annotation, Legend

<http://matplotlib.org/users/text_intro.html>

<http://matplotlib.org/users/annotations_guide.html>

<http://matplotlib.org/users/legend_guide.html>

## Plot

-   <http://stanford.edu/~mwaskom/software/seaborn/tutorial/regression.html>
-   <http://stanford.edu/~mwaskom/software/seaborn/tutorial/distributions.html>
-   <http://stanford.edu/~mwaskom/software/seaborn/tutorial/categorical.html>
-   <http://stanford.edu/~mwaskom/software/seaborn/examples/index.html>
-   <http://matplotlib.org/users/screenshots.html>

### Axes-level plot

`ax` を受け取ってそこに描画して `ax` を返す。

`sns.regplot(x, y, data, ..., fit_reg=True, ci=95, ..., ax)`
:   散布図。`fit_reg=False` しないと勝手に回帰線が引かれる。

`sns.distplot(a, bins, hist=True, kde=True, rug=False, fit=None, ..., ax)`
:   ヒストグラム、カーネル密度推定、ラグプロットなど分布全般

`sns.boxplot(x, y, hue, data, order, ..., ax)`
:   箱ひげ図

`sns.violinplot(x, y, hue, data, order, ..., ax)`
:   バイオリンプロット

`sns.stripplot(x, y, hue, data, order, ..., ax)`
:   片軸がカテゴリカル変数の散布図

`sns.pointplot(x, y, hue, data, order, ..., ax)`
:   平均値の折れ線グラフ + エラーバー

`sns.barplot(x, y, hue, data, order, ..., ax)`
:   平均値の棒グラフ + エラーバー

`sns.countplot(x, y, hue, data, order, ..., ax)`
:   カテゴリカル変数の頻度棒グラフ

`sns.heatmap(data, vmin, vmax, cmap, center, ..., square, mask, ax)`
:   ヒートマップ

### Figure-level plot

Grid に Axes-level plot を乗せて返す高レベル関数。

`sns.lmplot(x, y, data, hue, col, row, ...)`
:   `regplot()` + `FacetGrid()`

`sns.jointplot(x, y, data, kind, stat_func, ...)`
:   散布図 + 周辺ヒストグラム with `JointGrid()`.\
    `kind`: {`scatter`, `reg`, `resid`, `kde`, `hex`}

`sns.pairplot(data, hue, hue_order, palette, vars, x_vars, y_vars, kind, diag_kind, ...)`
:   ペアワイズ散布図 with `PairGrid()`

`sns.clustermap(data, ...)`
:   `heatmap()` + `ClusterGrid()`

`sns.factorplot(x, y, hue, data, row, col, ..., kind, ...)`
:   カテゴリカル変数全般.\
    `kind`: {`point`, `bar`, `count`, `box`, `violin`, `strip`}

## Grid

<http://stanford.edu/~mwaskom/software/seaborn/tutorial/axis_grids.html>

### `sns.FacetGrid(data, row, col, hue, col_wrap, sharex, sharey, ...)`

カテゴリカル変数でプロットを分ける:
```py
grid = sns.FacetGrid(iris, col='species', col_wrap=2)
grid.map(sns.regplot, 'sepal_width', 'sepal_length')
```

色分けもこれの仕事:
```py
grid = sns.FacetGrid(iris, hue='species')
grid.map(sns.regplot, 'sepal_width', 'sepal_length')
```

### `sns.PairGrid(data, ..., vars, x_vars, y_vars, ...)`

ペアワイズ散布図 + ヒストグラム
```py
grid = sns.PairGrid(iris)
grid = grid.map_diag(sns.distplot)
grid = grid.map_offdiag(sns.regplot)
```

### `sns.JointGrid(x, y, data, ...)`

散布図 + 周辺分布:
```py
grid = sns.JointGrid('sepal_width', 'sepal_length', iris)
grid = grid.plot_joint(sns.regplot)
grid = grid.plot_marginals(sns.distplot, kde=False)
```

------------------------------------------------------------------------

<http://matplotlib.org/users/gridspec.html>

<http://matplotlib.org/users/tight_layout_guide.html>

データに関係なく複数の図を並べる。

### `plt.subplots(nrows, ncols, sharex, sharey, ...)`

等サイズに分割:
```py
fig, axes = plt.subplots(2, 2)
sns.regplot('x', 'y', d, ax=axes[0, 0])
fig.tight_layout()
```

### `sns.gridspec.GridSpec(nrows, ncols, ...)`

e.g., 2x2分割して "品" みたいな配置にする:
```py
gs = sns.gridspec.GridSpec(2, 2)
ax_top = plt.subplot(gs[0, :])
ax_bottom_l = plt.subplot(gs[1, 0])
ax_bottom_r = plt.subplot(gs[1, 1])
```

### `sns.gridspec,GridSpecFromSubplotSpec(nrows, ncols, subplot_spec, ...)`

入れ子で分割。
e.g., 左右に分け、それぞれをさらに3段に分ける:
```py
gs = sns.gridspec.GridSpec(1, 2)
gsl = sns.gridspec.GridSpecFromSubplotSpec(3, 1, gs[0])
gsr = sns.gridspec.GridSpecFromSubplotSpec(3, 1, gs[1])
ax_ltop = plt.subplot(gsl[0])
```

## その他

### インストール

<http://matplotlib.org/faq/installing_faq.html>

[pip]({{< relref "pip.md" >}}) 一発でいけるはず:

    % pip install numpy scipy pandas
    % pip install matplotlib seaborn

### 設定

<http://matplotlib.org/users/customizing.html>

<http://matplotlib.org/faq/troubleshooting_faq.html>

<http://matplotlib.org/faq/environment_variables_faq.html>

`$HOME/.matplotlib/matplotlibrc` が読まれる。

`site-packages/matplotlib/mpl-data/matplotlibrc` にテンプレートがある。

### キャッシュ問題

使用するPythonを変更すると以下のようなランタイムエラーが出ることがある:

    RuntimeError: Could not open facefile /Library/Python/2.6/site-packages/matplotlib/mpl-data/fonts/ttf/Vera.ttf; Cannot_Open_Resource

そういうときはキャッシュを削除してみるとよい:

    % rm -rf ~/.matplotlib/fontList.cache

### バックエンド問題

Macで非Frameworkとしてインストールした自前Pythonを使うと怒られる:

    RuntimeError: Python is not installed as a framework. The Mac OS X backend will not be able to function correctly if Python is not installed as a framework. See the Python documentation for more information on installing Python as a framework on Mac OS X. Please either reinstall Python as a framework, or try one of the other backends.

`~/.matplotlib/matplotlibrc`
でバックエンドを `MacOSX` 以外に変更すればとりあえずOK:

    backend: TkAgg

<http://qiita.com/katryo/items/918667f28301fdec89ba>

<http://matplotlib.org/faq/usage_faq.html#what-is-a-backend>

### ラッパー

-   `seaborn`: <http://stanford.edu/~mwaskom/software/seaborn/>
-   easyplot: <https://github.com/HamsterHuey/easyplot>
-   prettyplotlib: <http://blog.olgabotvinnik.com/prettyplotlib/>
-   uglyplotlib: <https://gitlab.com/padawanphysicist/uglyplotlib>

## 書籍

<a rel="nofollow" href="http://www.amazon.co.jp/gp/product/4873116554/ref=as_li_ss_il?ie=UTF8&camp=247&creative=7399&creativeASIN=4873116554&linkCode=as2&tag=heavywatal-22"><img border="0" src="http://ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873116554&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="http://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=as2&o=9&a=4873116554" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a rel="nofollow" href="http://www.amazon.co.jp/gp/product/4873117356/ref=as_li_ss_il?ie=UTF8&camp=247&creative=7399&creativeASIN=4873117356&linkCode=as2&tag=heavywatal-22"><img border="0" src="http://ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117356&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="http://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=as2&o=9&a=4873117356" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
