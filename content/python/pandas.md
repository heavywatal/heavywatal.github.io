+++
title = 'Pandas'
subtitle = "data.frame in Python"
tags = ["python"]
[menu.main]
  parent = "python"
+++

<https://pandas.pydata.org/>

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

iris = sns.load_dataset('iris')
```

## 型

<https://pandas.pydata.org/pandas-docs/stable/user_guide/dsintro.html>

[`Series`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.html)
: 1次元の名前付き `np.array()` みたいなもの。

[`DataFrame`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html)
: Seriesを列とする2次元の表。
: `axis=0` が行(index)で、`axis=1` が列(columns)。

DataFrameを時間軸方向に重ねたような3次元構造として
`Panel` ってのもあったけど `MultiIndex` があれば不要ってことでv0.20からdeprecated、v0.25で削除。
["pan(el)-da(ta)-s" が名前の由来だったのに。](https://pandas.pydata.org/pandas-docs/version/0.24/getting_started/dsintro.html#panel)


## 読み書き

[`pd.read_csv(infile)`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html)
: `sep=','`
: `header='infer'`
: `names=None`
: `index_col=None`
: `usecols=None`

[`df.to_csv(outfile)`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html)
: `sep=','`
: `float_format=None`
: `index=True`


## 基本操作

https://pandas.pydata.org/pandas-docs/stable/10min.html

```python
df.columns
df.index
df.values
df.dtypes
df.head()
df.tail()
```

https://pandas.pydata.org/pandas-docs/stable/indexing.html

```python
df[['species']]  # DataFrame with 1 column
df['species']    # Series
df.species       # Series (not recommended)

# label-based
df.loc[0]
df.loc[:,'sepal_width']

# integer-based
df.iloc[0]
df.iloc[:,1]

# fast scalar (single value) lookup
df.at[0,'sepal_width']
df.iat[0,1]
```

`df.species` のような attribute access
は既存のメソッド名(e.g., `min`)と被るとダメなので基本的には使わないほうが良さそう。

番号でもラベルでもアクセスできる `df.ix[]` もあったが、
曖昧で危険なのでdeprecated.

### method chaining

DataFrameを何回も再帰代入するより、メソッドを繋いだほうが書きやすいし読みやすい。
カッコかバックスラッシュを使えばドット前に改行やスペースを入れて整形可能。

```py
(iris.query('species != "setosa"')
     .filter(regex='^(?!sepal)')
     .assign(petal_area=lambda x: x['petal_length'] * x['petal_width'] * 0.5)
     .groupby('species')
     .aggregate(np.mean)
     .sort_values(['petal_length'], ascending=False)
)
```

See also [dplyr + magrittr on R]({{< relref "dplyr.md" >}})

## 変形

- [Reshaping](https://pandas.pydata.org/pandas-docs/stable/reshaping.html)
- [Hierarchical indexing](https://pandas.pydata.org/pandas-docs/stable/advanced.html)
- [`MultiIndex`](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.MultiIndex.html)

### Index

`df.set_index(keys, drop=True, append=False, inplace=False, verify_integrity=False)`
: 指定したcolumnをindexにする。

`df.reset_index(level=None, drop=False, inplace=False, ...)`
: 現在のindexをcolumn化して、新しいものを振り直す。

`df.swaplevel(0, 1, axis=0)`

`df.reorder_levels([1, 0], axis=0)`

`df.rename_axis(mapper, axis=0, copy=True, inplace=False)`


### melt, pivot

Rと同様、列を列として扱う。

`df.melt(id_vars=None, value_vars=None, var_name=None, value_name='value', ...)`
: 縦長に変形。Rでいう `tidyr::gather()`。
  e.g., `iris.melt('species', var_name='variable', value_name='value')`

`df.pivot(index=None, columns=None, values=None)`
: 横長に変形。Rでいう `tidyr::spread()`。
  動かさない列を指定できないのかな？
  `reshape2::dcast()` のように複数の値をaggregationできる亜種として
  `df.pivot_table()` があるけどまあ使わないのが無難か。

```py
molten = iris.reset_index().melt(['index', 'species'])
molten.pivot('index', 'variable')           # species columns are redundant
molten.pivot('index', 'variable', 'value')  # species column is removed
```

### stack, unstack

MultiIndexを中心に考える。
nested tibble的なイメージ？

`df.stack(level=-1, dropna=True)`
: 縦長に変形。
  `melt()` と違って変数名はcolumnにならず、
  新しいindexとして最内側に追加される。
  そのindexに名前をつけるオプションが欲しかった。
  部分的に変形するオプションも無いので、
  残したい列を予めindexにしておく必要がある。

`df.unstack(level=-1, fill_value=None)`, `s.unstack()`
: 横長に変形。
  展開するindexの階層を指定できる(デフォルトの`-1`は最内側)。

```py
stacked = iris.set_index(['species'], append=True).stack()
stacked.unstack()

# more explicitly
stacked.rename_axis(['id', 'species', 'variable']).unstack('variable')
```

## 設定

https://pandas.pydata.org/pandas-docs/stable/options.html

```python
pd.set_option('display.max_rows', 20)
pd.set_option('display.width', None)
```

ディスプレイ幅は `os.get_terminal_size().columns` で明示的に取得してもいい。


## misc.

英語の発音としては [pǽndəz] のように濁るのが自然だけど
[開発者は pan-duss [pǽndəs] と発音している](https://twitter.com/wesmckinn/status/706661972431892483)
とのこと。



## 書籍

<a href="https://www.amazon.co.jp/dp/487311845X/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=72a416f5d10a9e84aaab4b3ee9613329&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311845X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=487311845X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873118417/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=6b1a04ec880b6c730bd6e80273e30e9c&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873118417&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873118417" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873117488/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=2181a50362009e68f507d44fc38716b4&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
