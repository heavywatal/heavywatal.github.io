+++
title = 'matplotlib + seaborn'
subtitle = "Pythonでグラフ描画"
tags = ["python", "graph"]
[menu.main]
  parent = "python"
+++

[`matplotlib`](http://matplotlib.org/) はPythonにおけるデータ可視化のデファクトスタンダード。
基本的には何でもできるけど、基本的な機能しか提供していないので、
いくらかの便利機能を [`seaborn`](https://seaborn.pydata.org/) で補う。

## 基本

-   <http://matplotlib.org/faq/usage_faq.html>
-   <http://matplotlib.org/users/>
-   <http://matplotlib.org/faq/howto_faq.html>
-   <http://www.scipy-lectures.org/intro/matplotlib/matplotlib.html>

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

# Show figure (in Jupyter, Hydrogen, or other inline IPython environments)
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


### [pyplot](http://matplotlib.org/api/pyplot_summary.html)

最上位のモジュール。
Figure, Axesインスタンスを明示的に操作するスタイルでは、
最初にFigureを作るくらいしか出番が無いはず。

`plt.figure(num=None, figsize=None, dpi=None, facecolor=None, edgecolor=None, ...)`
: Backendなどを設定してFigureインスタンスを作る

`plt.close(*args)`
: Figureウィンドウを閉じる。
  backend関係も切るっぽいので再び `fig.show()` しても開けない。


### [Figure](http://matplotlib.org/api/figure_api.html)

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


### [Axes](http://matplotlib.org/api/axes_api.html)

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

- <http://matplotlib.org/users/gridspec.html>
- <http://matplotlib.org/users/tight_layout_guide.html>

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

<http://matplotlib.org/users/text_intro.html>

<http://matplotlib.org/users/annotations_guide.html>

<http://matplotlib.org/users/legend_guide.html>



## Seaborn

### Axes-level plot

色分けや推定値の追加など、生のmatplotlibでやるにはちょっと大変なことがseabornの関数で簡単にできる。
Axesを受け取ってそこに描画するという単純な構造なので、何か自分で作ってもいい:
```py
def my_scatter(x, y, data, ax):
    ax.scatter(x, y, data=data)
    return ax
```

<https://seaborn.pydata.org/examples>

#### [Categorical plots](https://seaborn.pydata.org/tutorial/categorical.html)

`sns.boxplot(x, y, hue, data, order, ..., ax)`
:   箱ひげ図

`sns.violinplot(x, y, hue, data, order, ..., ax)`
:   バイオリンプロット

`sns.stripplot(x, y, hue, data, order, ..., ax)`
:   片軸がカテゴリカル変数の散布図

`sns.pointplot(x, y, hue, data, order, ..., ax)`
:   点推定値(平均値とか)の折れ線グラフ + エラーバー

`sns.barplot(x, y, hue, data, order, ..., ax)`
:   平均値の棒グラフ + エラーバー

`sns.countplot(x, y, hue, data, order, ..., ax)`
:   カテゴリカル変数の頻度棒グラフ

#### [Distribution plots](https://seaborn.pydata.org/tutorial/distributions.html)

`sns.distplot(a, bins, hist=True, kde=True, rug=False, fit=None, ..., ax)`
:   `ax.hist()` + `sns.kdeplot()` + `sns.rugplot()`

`sns.kdeplot(data, data2=None, shade=False, ..., ax)`
:   カーネル密度推定。

`sns.heatmap(data, vmin, vmax, cmap, center, ..., square, mask, ax)`
:   ヒートマップ。入力データはtidyじゃなくて行列の形。

#### [Regression plots](https://seaborn.pydata.org/tutorial/regression.html)

`sns.regplot(x, y, data, ..., fit_reg=True, ci=95, ..., ax)`
:   散布図 + 回帰線。

### [Axis Grid](https://seaborn.pydata.org/tutorial/axis_grids.html)

FigureとAxisをいい感じに初期化して、関連するデータを縦・横・色の方向に並べる土台。
これにAxis-level plotを乗せるところまでショートカットする高級関数がFigure-level plot。
できあがったGridクラスの`.set()`系メソッドとか`.fig`プロパティを通じていろいろ調整できる。
これをまた別のグリッドに埋め込むというRのgrobのような操作はたぶんできない。

#### [`sns.FacetGrid`](https://seaborn.pydata.org/generated/seaborn.FacetGrid.html)

`(data, row, col, hue, col_wrap, sharex, sharey, ...)`

カテゴリカル変数でプロットを分けて並べる:
```py
grid = sns.FacetGrid(iris, col='species', col_wrap=2)
grid.map(sns.regplot, 'sepal_width', 'sepal_length')
```

変数によって色分けする:
```py
grid = sns.FacetGrid(iris, hue='species')
grid.map(plt.scatter, 'sepal_width', 'sepal_length')
```

`map()` メソッドには `plt.scatter` など生のmatplotlib関数も渡せる。

`sns.lmplot(x, y, data, hue, col, row, ...)`
:   `regplot()` + `FacetGrid()` のショートカット。

`sns.factorplot(x, y, hue, data, row, col, ..., kind, ...)`
:   Categorical plot + `FacetGrid()` のショートカット。
:   `kind`: {`point`, `bar`, `count`, `box`, `violin`, `strip`}

#### [`sns.PairGrid`](https://seaborn.pydata.org/generated/seaborn.PairGrid.html)

`(data, hue, ..., vars, x_vars, y_vars, ...)`

ペアワイズ散布図 + 対角線ヒストグラム
```py
grid = sns.PairGrid(iris)
grid = grid.map_offdiag(sns.regplot)
grid = grid.map_diag(sns.distplot)
```

`sns.pairplot(data, hue, hue_order, palette, vars, x_vars, y_vars, kind, diag_kind, ...)`
:   `PairGrid()` のショートカット。
:   `kind`: {`scatter`, `reg`}
:   `diag_kind`: {`hist`, `kde`}

#### [`sns.JointGrid`](https://seaborn.pydata.org/generated/seaborn.JointGrid.html)

`(x, y, data, size, ratio, space, dropna, xlim, ylim)`

散布図 + 周辺分布:
```py
grid = sns.JointGrid('sepal_width', 'sepal_length', iris)
grid = grid.plot_joint(sns.regplot)
grid = grid.plot_marginals(sns.distplot, kde=False)
```

`sns.jointplot(x, y, data, kind, stat_func, ...)`
:   `JointGrid()` のショートカット。
:   `kind`: {`scatter`, `reg`, `resid`, `kde`, `hex`}

#### `sns.ClusterGrid()`

`sns.clustermap(data, ...)`
:   `heatmap()` + `ClusterGrid()`


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

- <https://matplotlib.org/users/colormaps.html>
- <https://seaborn.pydata.org/tutorial/color_palettes.html>

いくつかの方法で指定できる:

- 個々の描画関数の `palette` や `cmap` に指定
- `with sns.color_palette():` のブロック内で描画
- `sns.set_palette()` でデフォルトを指定


パレットもいくつかある:

-   [ColorBrewer](http://colorbrewer2.org/) の名前で指定。
    `_r` をつけると逆に、`_d` をつけると暗めになる。
    e.g., `sns.color_palette('RdBu_r', n_colors=7)`
-   [xkcd](http://xkcd.com/color/rgb/) の名前リストを
    `sns.xkcd_palette()` に渡す
-   `cubehelix_palette()` はgrayscaleでもいい感じで印刷できる
-   自分で作る:
    `sns.hls_palette()`, `sns.husl_palette()`,
    `sns.light_palette()`, `sns.dark_palette()`,
    `sns.diverging_palette()`



## その他

### インストール

[pyenvか何かで最新のPython3系をインストールして]({{< relref "install.md" >}})、
`pip install seaborn` を実行。

### 設定

- <http://matplotlib.org/users/customizing.html>
- <http://matplotlib.org/faq/troubleshooting_faq.html>
- <http://matplotlib.org/faq/environment_variables_faq.html>

`~/.matplotlib/matplotlibrc` が読まれる。

`site-packages/matplotlib/mpl-data/matplotlibrc` にテンプレートがある。

### [バックエンド問題](http://matplotlib.org/faq/usage_faq.html#what-is-a-backend)

Macで非Frameworkとしてインストールした自前Pythonを使うと怒られる:

    RuntimeError: Python is not installed as a framework. The Mac OS X backend will not be able to function correctly if Python is not installed as a framework. See the Python documentation for more information on installing Python as a framework on Mac OS X. Please either reinstall Python as a framework, or try one of the other backends.

`~/.matplotlib/matplotlibrc` に
`backend: tkagg` などと書いて
`macosx` 以外に変更すればとりあえずOK。
[Frameworkでインストールする]({{< relref "install.md" >}})ほうがいいとは思うけど。

### キャッシュ問題

使用するPythonを変更すると以下のようなランタイムエラーが出ることがある:

    RuntimeError: Could not open facefile /Library/Python/2.6/site-packages/matplotlib/mpl-data/fonts/ttf/Vera.ttf; Cannot_Open_Resource

そういうときはキャッシュを削除してみるとよい:

    rm -rf ~/.matplotlib/fontList.cache


## 書籍

<a href="https://www.amazon.co.jp/dp/487311845X/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=72a416f5d10a9e84aaab4b3ee9613329&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311845X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=487311845X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873118417/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=6b1a04ec880b6c730bd6e80273e30e9c&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873118417&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873118417" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873117488/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=2181a50362009e68f507d44fc38716b4&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
