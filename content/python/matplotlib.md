+++
title = 'matplotlib + seaborn'
subtitle = "Pythonでグラフ描画"
tags = ["python", "graph"]
[menu.main]
  parent = "python"
+++

- <http://matplotlib.org/>
- <https://seaborn.pydata.org/>

行列を含む数値計算やグラフ描画を行うためのソフトウェアとしては
有償のMATLABが広く利用されているが、
Python にいくつかのモジュールを導入することで
同等かそれ以上のことを無料で実現することができる。

- `matplotlib`: グラフ描画
- [`numpy`]({{< relref "scipy.md" >}}): 配列演算や基本的な数学関数
- [`scipy`]({{< relref "scipy.md#scipy" >}}): 高度な科学技術計算
- [`ipython`]({{< relref "ipython.md" >}}): 高機能なシェル環境

直接 `matplotlib` を触るのは大変なので、
[`seaborn`](https://seaborn.pydata.org/) というラッパーを介して使う。
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

gs = plt.GridSpec(1, 1)
fig = plt.figure()
ax = fig.add_subplot(gs[0])

sns.regplot('sepal_width', 'sepal_length', data=iris, ax=ax)

fig.show()
fig.savefig('example.png')
```

`from matplotlib.pylab import *`
は単にMATLABっぽいインターフェイスにするための乱暴な手段で、
名前空間が汚れるので使用しない。

`pyplot` も同様にあまり使いたくない `import` 主体のモジュールだが、
こちらはbackendのお世話などもしてくれているらしいので、
使わずに済ませるのは難しそう。
例えば `pyplot.figure()` から生成したやつじゃないと `fig.show()` できない、とか。


### Figure, Axes

<http://matplotlib.org/faq/usage_faq.html#general-concepts>

<http://matplotlib.org/api/figure_api.html>

<http://matplotlib.org/api/axes_api.html>

## Style

<https://seaborn.pydata.org/tutorial/aesthetics.html>

背景色や補助線などの設定。

`sns.set_style(style, rc=None)`

`style`:
:   {`darkgrid`, `whitegrid`, `dark`, `white`, `ticks`}

`rc`:
:   see below

```py
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

<https://seaborn.pydata.org/tutorial/aesthetics.html#scaling-plot-elements-with-plotting-context-and-set-context>

ラベルや点・線などのスケール調整。

`sns.set_context(context, font_scale=1, rc=None)`

`context`
:   {`notebook`: 1.0, `paper`: 0.8, `talk`, 1.3, `poster`, 1.6}

`rc`:
:   see below

```py
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

<https://seaborn.pydata.org/tutorial/color_palettes.html>

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

- <https://seaborn.pydata.org/tutorial/regression.html>
- <https://seaborn.pydata.org/tutorial/distributions.html>
- <https://seaborn.pydata.org/tutorial/categorical.html>
- <https://seaborn.pydata.org/examples/index.html>
- <http://matplotlib.org/users/screenshots.html>

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

<https://seaborn.pydata.org/tutorial/axis_grids.html>

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

### `mpl.gridspec.GridSpec(nrows, ncols, ...)`

e.g., 2x2分割して "品" みたいな配置にする:
```py
gs = plt.GridSpec(2, 2)
ax_top = plt.subplot(gs[0, :])
ax_bottom_l = plt.subplot(gs[1, 0])
ax_bottom_r = plt.subplot(gs[1, 1])
```

### `mpl.gridspec.GridSpecFromSubplotSpec(nrows, ncols, subplot_spec, ...)`

入れ子で分割。
e.g., 左右に分け、それぞれをさらに3段に分ける:
```py
gs = plt.GridSpec(1, 2)
gsl = mpl.gridspec.GridSpecFromSubplotSpec(3, 1, gs[0])
gsr = mpl.gridspec.GridSpecFromSubplotSpec(3, 1, gs[1])
ax_ltop = plt.subplot(gsl[0])
```

## その他

### [インストール](http://matplotlib.org/faq/installing_faq.html)

[Anaconda]({{< relref "install.md#anaconda" >}})
には最初から含まれているので楽チン。
Minicondaなら `conda install seaborn` で一発。
そうじゃなくても `pip install seaborn` でいけるはず。

### 設定

- <http://matplotlib.org/users/customizing.html>
- <http://matplotlib.org/faq/troubleshooting_faq.html>
- <http://matplotlib.org/faq/environment_variables_faq.html>

`~/.matplotlib/matplotlibrc` が読まれる。

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

### その他のラッパー

-   easyplot: <https://github.com/HamsterHuey/easyplot>
-   prettyplotlib: <http://blog.olgabotvinnik.com/prettyplotlib/>
-   uglyplotlib: <https://gitlab.com/padawanphysicist/uglyplotlib>


## 書籍

<a href="https://www.amazon.co.jp/dp/479738946X/ref=as_li_ss_il?ie=UTF8&qid=1485612008&sr=8-6&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=5ea5e48ecc83b9439f21406b6f57c062" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=479738946X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=479738946X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/IPython%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9%E3%82%AF%E3%83%83%E3%82%AF%E3%83%96%E3%83%83%E3%82%AF-%E5%AF%BE%E8%A9%B1%E5%9E%8B%E3%82%B3%E3%83%B3%E3%83%94%E3%83%A5%E3%83%BC%E3%83%86%E3%82%A3%E3%83%B3%E3%82%B0%E3%81%A8%E5%8F%AF%E8%A6%96%E5%8C%96%E3%81%AE%E3%81%9F%E3%82%81%E3%81%AE%E3%83%AC%E3%82%B7%E3%83%94%E9%9B%86-Cyrille-Rossant/dp/4873117488/ref=as_li_ss_il?_encoding=UTF8&psc=1&refRID=X16VFSS3W75RMTG7VGCH&linkCode=li3&tag=heavywatal-22&linkId=b79e2290571289b02621392257a4ac1c" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/Python%E3%81%AB%E3%82%88%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E5%88%86%E6%9E%90%E5%85%A5%E9%96%80-_NumPy%E3%80%81pandas%E3%82%92%E4%BD%BF%E3%81%A3%E3%81%9F%E3%83%87%E3%83%BC%E3%82%BF%E5%87%A6%E7%90%86-Wes-McKinney/dp/4873116554/ref=as_li_ss_il?s=books&ie=UTF8&qid=1485612076&sr=1-16&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=f024f5e24f16c38402149f97591c8aab" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873116554&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873116554" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/%E3%82%BC%E3%83%AD%E3%81%8B%E3%82%89%E4%BD%9C%E3%82%8BDeep-Learning-Python%E3%81%A7%E5%AD%A6%E3%81%B6%E3%83%87%E3%82%A3%E3%83%BC%E3%83%97%E3%83%A9%E3%83%BC%E3%83%8B%E3%83%B3%E3%82%B0%E3%81%AE%E7%90%86%E8%AB%96%E3%81%A8%E5%AE%9F%E8%A3%85-%E6%96%8E%E8%97%A4-%E5%BA%B7%E6%AF%85/dp/4873117585/ref=as_li_ss_il?s=books&ie=UTF8&qid=1485612076&sr=1-3&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=cde2e18b945af3ead43dc15e51b00af6" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117585&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117585" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
