+++
title = 'matplotlib + seaborn'
subtitle = "Pythonでグラフ描画"
tags = ["python", "graph"]
[menu.main]
  parent = "python"
+++

[`matplotlib`](https://matplotlib.org/) はPythonにおけるデータ可視化のデファクトスタンダード。
とはいえユーザーが直接これを使ってグラフを描ききるのは難しく、
Rでいうgridパッケージに近い階層と見なしたほうがいいかもしれない。
[`seaborn`](https://seaborn.pydata.org/) 越しに使うのが便利。

## 基本

<https://matplotlib.org/stable/users/explain/quick_start.html>

```py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

iris = sns.load_dataset('iris')

# Create an empty Figure
fig = plt.figure()

# Add an Axes to this fig
ax = fig.subplots()

# Plot on this ax
ax.scatter('sepal_width', 'sepal_length', data=iris)

# Show figure (in Jupyter and other inline IPython environments)
display(fig)

# Show figure in a new window (with non-inline backends)
fig.show()
plt.close(fig)

# Write to a file
fig.savefig('example.png')
```

渡すデータは生のlistとかではなくtidyな
[`pandas.DataFrame`]({{< relref "pandas.md" >}})型にしておく。

`from matplotlib.pylab import *`
は単にMATLABっぽいインターフェイスにするための乱暴な手段で、
例としてよく見かけるけど公式に非推奨とされている。

FigureやAxesを意識せず `plt.plot()` などを使うスタイルは分かりにくいので不採用。


### pyplot

<https://matplotlib.org/stable/api/pyplot_summary.html>

最上位のモジュール。
Figure, Axesインスタンスを明示的に操作するスタイルでは、
最初にFigureを作るくらいしか出番が無いはず。

`plt.figure(num=None, figsize=None, dpi=None, facecolor=None, edgecolor=None, ...)`
: Backendなどを設定してFigureインスタンスを作る

`plt.close(*args)`
: Figureウィンドウを閉じる。
  backend関係も切るっぽいので再び `fig.show()` しても開けない。


### Figure

<https://matplotlib.org/stable/api/figure_api.html>

ウィンドウを表示したり画像ファイルを保存したりする単位となるクラス。
コンストラクタを直接呼ぶのではなく
`plt` や `sns` に作らせるのが普通。
複数の子Axesを保持できる。

`fig.subplots(nrows=1, ncols=1, sharex=False, sharey=False, ...)`
: 子Axesをタイル状に並べて作る。
  デフォルトでは単体のAxesを返す。

`fig.add_subplot(*args, **kwargs)`
: 指定した位置に子Axesを作る。
  少し複雑な配置にしたいときはこれに `GridSpec` を渡すのが便利。

`fig.clear()`
: `fig.axes` を空っぽにする。
  backend関係は切れないので子Axesを追加して再描画可能。

`fig.savefig(fname, dpi=None, facecolor='w', edgecolor='w', **kwargs)`
: 拡張子から画像形式を推定してくれるので `format=` は省略可能。

`fig.show()`
: Windowを開いてFigureを表示する。
  JupyterやHydrogenなどの[IPython]({{< relref "ipython.md" >}})環境では
  `display(fig)` を使う。

`fig.axes`
: 子Axesへの参照


### Axes

<https://matplotlib.org/stable/api/axes_api.html>

軸やラベルを持ったひとつのプロットの単位となるクラス。

描画のためのメソッドが
`ax.plot()`, `ax.scatter()`, `ax.hist()`
などたくさんある。
[pandas.DataFrame]({{< relref "pandas.md" >}}) のメソッドとして呼び出すのもあり:
```py
# Axes method with data
ax.scatter('sepal_width', 'sepal_length', data=iris)

# DataFrame method with ax
iris.plot.scatter('sepal_width', 'sepal_length', ax=ax)
```

設定のためのメソッドも
`ax.set_title()`, `ax.set_xlim()` などたくさん。
`ax.set(**kwargs)` でまとめて設定することもできる。

`ax.clear()`
: プロットしたものを消す。親Figureへの参照は残る。

`ax.figure`
: 親Figureへの参照


### 複数のAxesを配置する

- <https://matplotlib.org/stable/api/gridspec_api.html>
- <https://matplotlib.org/stable/tutorials/intermediate/gridspec.html>

`plt.subplots(nrows, ncols, sharex, sharey, ...)`
:   等サイズに分割:

    ```py
    fig = plt.figure()
    axes = fig.subplots(2, 2)
    sns.regplot('x', 'y', d, ax=axes[0, 0])
    fig.tight_layout()
    ```

`mpl.gridspec.GridSpec(nrows, ncols, ...)`
:   e.g., 2x2分割して "品" みたいな配置にする:

    ```py
    fig = plt.figure()
    gs = plt.GridSpec(2, 2)
    ax_top = fig.add_subplot(gs[0, :])
    ax_bottom_l = fig.add_subplot(gs[1, 0])
    ax_bottom_r = fig.add_subplot(gs[1, 1])
    ```

`mpl.gridspec.GridSpecFromSubplotSpec(nrows, ncols, subplot_spec, ...)`
:   入れ子で分割。
    e.g., 左右に分け、それぞれをさらに3段に分ける:

    ```py
    fig = plt.figure()
    gs = plt.GridSpec(1, 2)
    gsl = sns.mpl.gridspec.GridSpecFromSubplotSpec(3, 1, gs[0])
    gsr = sns.mpl.gridspec.GridSpecFromSubplotSpec(3, 1, gs[1])

    ax_ltop = fig.add_subplot(gsl[0])
    ```

### Text, Annotation, Legend

<https://matplotlib.org/stable/users/explain/text/text_intro.html>

<https://matplotlib.org/stable/users/explain/axes/legend_guide.html>



## Seaborn

### `seaborn.objects` interface

<https://seaborn.pydata.org/tutorial/objects_interface.html>

v0.12 で導入された新しいAPI。
[ggplot2]({{< relref "ggplot2.md" >}}) に近い感覚で書くことができる。
従来のAPIで描けるものをどれくらいカバーしているかはわからないけど、
今後はこちらを使っていきたい。

```py
import seaborn as sns
import seaborn.objects as so

penguins = sns.load_dataset("penguins")
p = so.Plot(penguins, x="bill_length_mm", y="bill_depth_mm", color="species")
p.add(so.Dot(alpha = 0.5))
p.show()
```

### Axes-level plot

色分けや推定値の追加など、生のmatplotlibでやるにはちょっと大変なことがseabornの関数で簡単にできる。
Axesを受け取ってそこに描画するという単純な構造なので、何か自分で作ってもいい:
```py
def my_scatter(x, y, data, ax):
    ax.scatter(x, y, data=data)
    return ax
```

<https://seaborn.pydata.org/examples>


#### Visualizing statistical relationships

<https://seaborn.pydata.org/tutorial/relational.html>

`sns.scatterplot(data, *, markers, ..., ax)`
:   散布図。

`sns.lineplot(data, *, style, markers, dashes, ..., ax)`
:   折れ線グラフ。


#### Distribution plots

<https://seaborn.pydata.org/tutorial/distributions.html>

`sns.histplot(data, stat='count', bins='auto', binwidth=None, binrange=None, discrete=None, ..., ax)`
:   ヒストグラム。
:   `discrete=True` とするだけで整数をうまく描いてくれる。
:   半端な立ち位置だった [`distplot()`はdeprecated](https://gist.github.com/mwaskom/de44147ed2974457ad6372750bbe5751).

`sns.kdeplot(data, data2=None, shade=False, ..., ax)`
:   カーネル密度推定。

`sns.ecdfplot(data, ..., ax)`
:   empirical cumulative density function.

`sns.rugplot(data, ..., ax)`
:   軸沿いにtickを描く。

`sns.heatmap(data, vmin, vmax, cmap, center, ..., square, mask, ax)`
:   ヒートマップ。入力データはtidyじゃなくて行列の形。

#### Categorical plots

<https://seaborn.pydata.org/tutorial/categorical.html>

`sns.stripplot(x, y, hue, data, order, ..., ax)`
:   片軸がカテゴリカル変数の `geom_jitter` に相当。
:   これよりやや規則的な `sns.swarmplot()` も良い。

`sns.boxplot(x, y, hue, data, order, ..., ax)`
:   箱ひげ図

`sns.violinplot(x, y, hue, data, order, ..., ax)`
:   バイオリンプロット

`sns.boxenplot(x, y, hue, data, order, ..., ax)`
:   箱ひげ図の変種。別名 letter-value plot 。

`sns.pointplot(x, y, hue, data, order, ..., ax)`
:   点推定値(平均値とか)の折れ線グラフ + エラーバー

`sns.barplot(x, y, hue, data, order, ..., ax)`
:   平均値の棒グラフ + エラーバー

`sns.countplot(x, y, hue, data, order, ..., ax)`
:   カテゴリカル変数の頻度棒グラフ

#### Regression plots

<https://seaborn.pydata.org/tutorial/regression.html>

`sns.regplot(x, y, data, ..., fit_reg=True, ci=95, ..., ax)`
:   散布図 + 回帰線。


### Axis Grid

<https://seaborn.pydata.org/tutorial/axis_grids.html>

FigureとAxisをいい感じに初期化して、関連するデータを縦・横・色の方向に並べる土台。
これにAxis-level plotを乗せるところまでショートカットする高級関数がFigure-level plot。
できあがったGridクラスの`.set()`系メソッドとか`.fig`プロパティを通じていろいろ調整できる。
これをまた別のグリッドに埋め込むというRのgrobのような操作はたぶんできない。

#### `sns.FacetGrid`

<https://seaborn.pydata.org/generated/seaborn.FacetGrid.html>

- `row=None`, `col=None`, `hue=None`,
- `col_wrap=None`,
- `sharex=True`, `sharey=True`,
- `height=3`, `aspect=1`,
- `legend_out=True`, `despine=True`, `margin_titles=False`,
- dropna, palette, xlim, ylim,
  row_order, col_order, hue_order,
  hue_kws, subplot_kws, gridspec_kws

カテゴリカル変数でプロットを分けて並べる:
```py
grid = sns.FacetGrid(iris, col='species', col_wrap=2)
grid.map(sns.scatterplot, 'sepal_width', 'sepal_length')
```

変数によって色分けする:
```py
grid = sns.FacetGrid(iris, hue='species')
grid.map(sns.scatterplot, 'sepal_width', 'sepal_length')
```

`map()` メソッドには `plt.scatter` など生のmatplotlib関数も渡せる。

`sns.relplot(data, *, )`
:   散布図・折れ線グラフ + `FacetGrid()` のショートカット。
:   `kind`: {`scatter`, `line`}

`sns.displot(data, kind='hist', rug=False, ...)`
:   Distribution plot + `FacetGrid()` のショートカット。
:   `kind`: {`hist`, `kde`, `ecdf`}

`sns.catplot(x, y, hue, data, row, col, ..., kind, ...)`
:   Categorical plot + `FacetGrid()` のショートカット。
:   `kind`: {`strip`, `swarm`, `box`, `violin`, `boxen`, `point`, `bar`, `count`}
:   昔は `factorplot` という名前だった。

`sns.lmplot(x, y, data, hue, col, row, ...)`
:   `regplot()` + `FacetGrid()` のショートカット。


#### `sns.PairGrid`

<https://seaborn.pydata.org/generated/seaborn.PairGrid.html>

- `hue=None`
- `vars=None`, `x_vars=None`, `y_vars=None`
- `corner=False`: 下半分だけのcorner plotにするには `True`
- `height=2.5`, `aspect=1`, `layout_pad=0.5`
- hue_order=None, palette=None, hue_kws=None,
  diag_sharey=True, despine=True, dropna=False

ペアワイズ散布図 + 対角線ヒストグラム
```py
grid = sns.PairGrid(iris, corner=True)
grid = grid.map_offdiag(sns.scatterplot)
grid = grid.map_diag(sns.histplot)
```

`sns.pairplot(data, hue, hue_order, palette, vars, x_vars, y_vars, kind, diag_kind, ...)`
:   `PairGrid()` のショートカット。
:   `kind`: {`scatter`, `reg`}
:   `diag_kind`: {`hist`, `kde`}

#### `sns.JointGrid`

<https://seaborn.pydata.org/generated/seaborn.JointGrid.html>

- x, y, hue
- `height=6`, `ratio=5`, `space=0.2`
- palette, hue_order, hue_norm,
  dropna, xlim, ylim
- `marginal_ticks=False`

散布図 + 周辺分布:
```py
grid = sns.JointGrid(iris, x='sepal_width', y='sepal_length')
grid = grid.plot_joint(sns.scatterplot)
grid = grid.plot_marginals(sns.histplot)
```

`sns.jointplot(x, y, data, kind, stat_func, ...)`
:   `JointGrid()` のショートカット。
:   `kind`: {`scatter`, `kde`, `hist`, `hex`, `reg`, `resid`}

#### `sns.ClusterGrid()`

このGridクラスをユーザーが直接インスタンス化することは想定されていない。

<https://seaborn.pydata.org/generated/seaborn.clustermap.html>

`sns.clustermap(data, ...)`
:   `sns.heatmap()` + `ClusterGrid()`


### Style

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

### Context

<https://seaborn.pydata.org/tutorial/aesthetics.html#scaling-plot-elements>

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

- <https://matplotlib.org/stable/users/explain/colors/colormaps.html>
- <https://seaborn.pydata.org/tutorial/color_palettes.html>

いくつかの方法で指定できる:

- 個々の描画関数の `palette` や `cmap` に指定
- `with sns.color_palette():` のブロック内で描画
- `sns.set_palette()` でデフォルトを指定


パレットもいくつかある:

-   [ColorBrewer](https://colorbrewer2.org/) の名前で指定。
    `_r` をつけると逆に、`_d` をつけると暗めになる。
    e.g., `sns.color_palette('RdBu_r', n_colors=7)`
-   [xkcd](https://xkcd.com/color/rgb/) の名前リストを
    `sns.xkcd_palette()` に渡す
-   `cubehelix_palette()` はgrayscaleでもいい感じで印刷できる
-   自分で作る:
    `sns.hls_palette()`, `sns.husl_palette()`,
    `sns.light_palette()`, `sns.dark_palette()`,
    `sns.diverging_palette()`



## 設定

- <https://matplotlib.org/stable/users/explain/customizing.html>
- <https://matplotlib.org/stable/install/environment_variables_faq.html>

`~/.matplotlib/matplotlibrc` が読まれる。

`site-packages/matplotlib/mpl-data/matplotlibrc` にテンプレートがある。
