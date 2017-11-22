+++
title = 'Pandas'
subtitle = "data.frame in Python"
tags = ["python"]
[menu.main]
  parent = "python"
+++

http://pandas.pydata.org/

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

iris = sns.load_dataset('iris')
```

## 型

http://pandas.pydata.org/pandas-docs/stable/dsintro.html

[`Series`](http://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.html)
: 1次元の名前付き `np.array()` みたいなもの。

[`DataFrame`](http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html)
: Seriesを列とする2次元の表。
: `axis=0` が行(index)で、`axis=1` が列(columns)。

DataFrameを時間軸方向に重ねたような3次元構造として
`Panel` ってのもあったけど `MultiIndex` があれば不要ってことでdeprecated。


## 読み書き

[`pd.read_table(infile)`](http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_table.html)
: `sep='\t'`
: `header='infer'`
: `names=None`
: `usecols=None`

[`df.to_csv(outfile)`](http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.to_csv.html)
: `sep=','`
: `float_format=None`
: `index=True`


## 基本操作

http://pandas.pydata.org/pandas-docs/stable/10min.html

```python
df.columns
df.index
df.values
df.dtypes
df.head()
df.tail()
```

http://pandas.pydata.org/pandas-docs/stable/indexing.html

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
- [`MultiIndex`](https://pandas.pydata.org/pandas-docs/stable/generated/pandas.MultiIndex.html)

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

http://pandas.pydata.org/pandas-docs/stable/options.html

```python
pd.set_option('display.max_rows', 20)
pd.set_option('display.width', None)
```

ディスプレイ幅は `os.get_terminal_size().columns` で明示的に取得してもいい。


## 書籍

<a href="https://www.amazon.co.jp/dp/479738946X/ref=as_li_ss_il?ie=UTF8&qid=1485612008&sr=8-6&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=5ea5e48ecc83b9439f21406b6f57c062" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=479738946X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=479738946X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/IPython%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9%E3%82%AF%E3%83%83%E3%82%AF%E3%83%96%E3%83%83%E3%82%AF-%E5%AF%BE%E8%A9%B1%E5%9E%8B%E3%82%B3%E3%83%B3%E3%83%94%E3%83%A5%E3%83%BC%E3%83%86%E3%82%A3%E3%83%B3%E3%82%B0%E3%81%A8%E5%8F%AF%E8%A6%96%E5%8C%96%E3%81%AE%E3%81%9F%E3%82%81%E3%81%AE%E3%83%AC%E3%82%B7%E3%83%94%E9%9B%86-Cyrille-Rossant/dp/4873117488/ref=as_li_ss_il?_encoding=UTF8&psc=1&refRID=X16VFSS3W75RMTG7VGCH&linkCode=li3&tag=heavywatal-22&linkId=b79e2290571289b02621392257a4ac1c" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/Python%E3%81%AB%E3%82%88%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E5%88%86%E6%9E%90%E5%85%A5%E9%96%80-_NumPy%E3%80%81pandas%E3%82%92%E4%BD%BF%E3%81%A3%E3%81%9F%E3%83%87%E3%83%BC%E3%82%BF%E5%87%A6%E7%90%86-Wes-McKinney/dp/4873116554/ref=as_li_ss_il?s=books&ie=UTF8&qid=1485612076&sr=1-16&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=f024f5e24f16c38402149f97591c8aab" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873116554&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873116554" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/%E3%82%BC%E3%83%AD%E3%81%8B%E3%82%89%E4%BD%9C%E3%82%8BDeep-Learning-Python%E3%81%A7%E5%AD%A6%E3%81%B6%E3%83%87%E3%82%A3%E3%83%BC%E3%83%97%E3%83%A9%E3%83%BC%E3%83%8B%E3%83%B3%E3%82%B0%E3%81%AE%E7%90%86%E8%AB%96%E3%81%A8%E5%AE%9F%E8%A3%85-%E6%96%8E%E8%97%A4-%E5%BA%B7%E6%AF%85/dp/4873117585/ref=as_li_ss_il?s=books&ie=UTF8&qid=1485612076&sr=1-3&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=cde2e18b945af3ead43dc15e51b00af6" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117585&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117585" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
