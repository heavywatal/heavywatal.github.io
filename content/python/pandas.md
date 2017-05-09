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
```

## 型

`DataFrame`
: `Series` の集まり。

`Series`
: 名前付き `np.array()` みたいなもの。


## 読み書き

http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_table.html

```python
pd.read_table()
pd.read_csv()
df.to_csv(outfile, sep='\t', float_format='%g', index=False)
```


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
df['Species']
df.Species

# label-based
df.loc[0]
df.loc[0,'Sepal.Width']

# integer-based
df.iloc[0]
df.iloc[0,1]

# fast scalar lookup
df.at[0,'Sepal.Width']
df.iat[0,1]
```

`df.Species` のような attribute access
は既存のメソッド名(e.g., `min`)と被るとダメなので基本的には使わないほうが良さそう。

番号でもラベルでもアクセスできる `df.ix[]` もあったが、
曖昧で危険なのでdeprecated.


## 設定

http://pandas.pydata.org/pandas-docs/stable/options.html

```python
pd.set_option('display.max_rows', 20)
pd.set_option('display.width', None)
```

ディスプレイ幅は `os.get_terminal_size().columns` で明示的に取得してもいい。
