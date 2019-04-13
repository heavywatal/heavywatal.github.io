+++
title = 'copy'
tags = ["python"]
[menu.main]
  parent = "python"
+++

-   <http://docs.python.org/library/copy.html>
-   <http://docs.python.org/reference/datamodel.html>

C++からプログラミングを始めて、Pythonにおけるオブジェクトの扱い、
特に代入・参照・コピー・mutable/immutableらへんの理解に苦しんでる人のメモ。

## Pythonの代入は基本的に参照渡し

```py
>>> l = [1, 2, 3]
>>> m = l
>>> l[1] = 0
>>> m
[1, 0, 3]
```

`m` に渡るのは `l` と同じ実体に対する参照。
`l[1] = 0` は `l` を変更しているのではなく、`l` を通してその実体を変更してる。
その変更は `m` を通して見ても同じ。

```py
>>> l = [1, 2, 3]
>>> l = m
>>> l = [6, 6, 6]
>>> m
[1, 2, 3]
```

`l = [6, 6, 6]` は新しい実体への参照を `l` に与える。
`l` が参照していた実体に対する変更ではないので、
`m` から見える実体にも変化は無い。

```py
>>> x = 0
>>> y = x
>>> x = 1
>>> y
0
```

これも同様。`y` には `x` と同じ `0` への参照が渡る。
`x = 1` は `x` を通した `0` への変更ではなく
`x` の参照先を `1` に変えるだけなので、`y` の指す値に変更は無い。

## mutable / immutable

mutable
: `list`
: `set`
: `dict`

immutable
: `0, 1, 2, ...`
: `"strings"`
: `tuple`

immutableなオブジェクトを変えることはできない。
immutableなオブジェクトを参照していた変数に別のオブジェクトの参照を渡すことはできる。

```py
>>> t = (1, 2, 3)
>>> t[1] = 0
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: 'tuple' object does not support item assignment
>>> t = (0, 0, 0)
>>> t
(0, 0, 0)
```

文字列や数字も考え方は同じ。

```py
>>> s = "abc"
>>> s[1] = "z"
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: 'str' object does not support item assignment
>>> s = "xyz"
>>> s
'xyz'
>>> 0 = 1
  File "<stdin>", line 1
SyntaxError: can't assign to literal
>>> x = 0
>>> x = 1
>>> x
1
```

immutableなオブジェクトで変更できないのは「何を参照しているか」という情報。
参照先の実体がmutableなオブジェクトならそいつは変更できる。

```py
>>> l = [1, 2, 3]
>>> t = (l, l)
>>> t
([1, 2, 3], [1, 2, 3])
>>> l[1] = 0
>>> t
([1, 0, 3], [1, 0, 3])
>>> t[0][2] = 0
>>> l
[1, 0, 0]
```

`t` が参照しているのは、`l` が指す実体へのimmutableな参照を持つ `tuple` 。
「その `tuple` が `[0]` 番目で何を参照しているか」を変更することはできないが、
参照している実体はmutableな `list` なので `l` や `t[0]` を通して変更できる。
つまり、Pythonの `tuple` はC++でいうところの
`int* const` なポインタを格納した固定長配列みたいなもの？

## copy {#copy-1}

```py
>>> import copy
>>>
>>> l = [1, 2, 3]
>>> m = copy.copy(l)
>>> l[1] = 0
>>> l
[1, 0, 3]
>>> m
[1, 2, 3]
```

`copy.copy()` によって `l` と同じ中身の新しい `list` 実体が作られ、
それへの参照が `m` に渡される。`l` を通した古い `list` への変更は `m` から見える `list` には影響しない。
中身がすべてimmutableな一重のコンテナならこれでいいが、
そうじゃない場合「`l` と同じ中身」というのが問題になる。
C++で言えば、メンバ変数としてポインタを持つクラスを
デフォルトのコピーコンストラクタでコピーしたような状態になる。

```py
>>> import copy
>>>
>>> l = [1, 2, 3]
>>> m = [l, l]
>>> x = m
>>> y = copy.copy(m)
>>> z = copy.deepcopy(m)
>>> m
[[1, 2, 3], [1, 2, 3]]
>>> m[1][1] = 0
>>> m[1] = 0
>>> x
[[1, 0, 3], 0]
>>> y
[[1, 0, 3], [1, 0, 3]]
>>> z
[[1, 2, 3], [1, 2, 3]]
```

-   `x` はただの代入によって `m` と同じ実体への参照を受け取ったので、
    `m` を通して行った変更はすべて共有している。
-   `y` は `m` の実体と同じ中身（`l = [1, 2, 3]` で作った `list` への参照2つ）
    を持つ新しい `list` への参照を与えられている。
    `m[1] = 0` は `m` が参照している外側の `list` を変更するものなので、
    `y` や `y[1]` には影響しない。
    `m[1][1] = 0` は `l` や `x[0]` や `y[1]` が参照している
    `list` 実体への変更なので `y` から見えるものも変わる。
-   `z` は `y` と同じように新しい `list` への参照が与えられている。
    `y` とは違ってその `list` 実体が持つのは `m` と同じ中身ではなく、
    同じ値になるようにすべての要素が再起的にコピーされたもの。
    `z` の中にある `[1, 2, 3]` は `l` とは異なる実体であり、
    `l` や `m` などを通して行われた変更は全く影響しない。

<!-- workaround blackfriday bug -->
```py
>>> import copy
>>>
>>> class CopyTest(object):
...     def __init__(self):
...         self.l = [1,2,3]
...         self.s = "xyz"
...     def __repr__(self):
...         return str([self.l, self.s])
...     def self_copy(self):
...         self = copy.copy(self)
...     def self_deepcopy(self):
...         self = copy.deepcopy(self)
...     def copy(self):
...         return copy.copy(self)
...     def deepcopy(self):
...         return copy.deepcopy(self)
...
>>> a = CopyTest()
>>> b = a
>>> b.self_copy()
>>> c = a
>>> c.self_deepcopy()
>>> d = a.copy()
>>> e = a.deepcopy()
>>> a.l[1] = 0
>>> a.s = "---"
>>> a
[[1, 0, 3], '---']
>>> b
[[1, 0, 3], '---']
>>> c
[[1, 0, 3], '---']
>>> d
[[1, 0, 3], 'xyz']
>>> e
[[1, 2, 3], 'xyz']
```

`self = copy.copy(self)` なメソッドは成功しないらしい。
