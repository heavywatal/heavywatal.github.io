+++
title = 'IPython'
subtitle = "強力な対話型Python実行環境"
tags = ["python"]
[menu.main]
  parent = "python"
  weight = -90
+++

- <https://ipython.org/>
- <https://ipython.readthedocs.io/>

## 対話型実行環境

### 起動

ターミナルから `ipython` を実行するか、
普通の対話Pythonから:
```py
import IPython
IPython.start_ipython(argv=[])
```

とりあえずドキュメントを読む:
```py
?
%quickref
```

### 文脈を考慮したタブ補完

定義済み変数などはある程度 `rlcompleter` でも補完できるが、
IPythonはさらに文脈を考慮して賢い補完をしてくれる:
```py
import o[TAB]
import os.[TAB]
open('~/.[TAB]')
```

### 履歴

-   上下キーで単純に遡る
-   途中まで入力して <kbd>control-p</kbd> で前方一致する履歴のみ遡る
-   <kbd>control-r</kbd> から部分一致する履歴を検索
-   `%hist`
-   input cache: `_i`, `_ii`, `_iii`, `_ih[n]`, `_i<n>`
-   output cache: `_`, `__`, `___`, `_oh[n]`, `_<n>`
-   directory: `_dh`

### Object introspection

<https://ipython.readthedocs.io/en/stable/interactive/reference.html#dynamic-object-information>

関数がどんな引数をとるか、
クラスがどんなメンバを持っているか、
などをパッと覗くことができる。
`help()` コマンドの強力版。

先頭(か末尾)にクエスチョンをつけるだけ:
```py
?os
??os
```

-   `%pdoc <object>`: docstring
-   `%pdef <object>`: ?
-   `%psource <object>`: 定義しているソースコード
-   `%pfile <object>`: 定義しているファイル全体

### システムコマンド実行

頭にエクスクラメーションをつけるだけ:
```py
!pwd
files = !ls
```

### マジックコマンド

<https://ipython.readthedocs.io/en/stable/interactive/magics.html>

- 一行単位の line magic は `%` で始める。
- 複数行うけつける cell magic は `%%` で始める。
- デフォルトの `automagic = True` では `%` が省略可能で怖い。

### 環境設定

<https://ipython.readthedocs.io/en/stable/config/intro.html>

```sh
ipython help profile
ipython profile create [name]
ipython profile list
ipython profile locate
ipython --profile=<name>
```

[`~/.ipython/profile_default/ipython_config.py`](https://github.com/heavywatal/dotfiles/blob/master/.ipython/profile_default/ipython_config.py)

起動時に自動的に
[pandas]({{< relref "pandas.md" >}}) や
[matplotlib]({{< relref "matplotlib.md" >}})
を読み込むとか。
automagicを切るとか。


## Jupyter

<https://jupyter.org/>

ウェブブラウザ上で動く対話的実行環境。
元はIPython Notebookだったが、
カーネルを入れ替えることで他言語を扱うことが可能になり、
汎用[Jupyter Notebook](https://jupyter-notebook.readthedocs.io/)として生まれ変わった。
それをさらに開発環境として洗練させたのが[JupyterLab](https://jupyterlab.readthedocs.io/)。

ノートブック形式 `.ipynb` は
Markdown/LaTeX記法による見出し・コメント・数式とともに
ソースコードと実行結果をひとまとめに保存しておけるので、
研究ノートのような使い方に向いている。
Mathematica/Mapleの使い勝手に似ている。
GitHub上でも直接閲覧できるし、[VSCode]({{< relref "vscode.md" >}}) でも扱える。
ターミナルやテキストエディタに馴染みのない非プログラマ向けには便利だろう。

ファイルがインプット・アウトプット混在の複雑なテキストになるので、
**コマンドラインやGit上での取り回しが悪いのは致命的な欠点。**
後述のJupytextがこれを回避する救世主かもしれない。


### 始め方

<https://jupyterlab.readthedocs.io/en/stable/getting_started/starting.html>

1.  適当な作業ディレクトリを作って移動: `mkdir -p ~/jupyter; cd $_`
1.  ターミナルから起動: `jupyter lab [file or directory]`
1.  ウェブブラウザで `http://localhost:8888/lab/` が立ち上がる
1.  Launcher内のNotebookにある適当なカーネル (e.g., Python 3) を選択
1.  `[ ]:` の右の箱に適当なコマンド `print('Hello, world!')` を入れて <kbd>shift-return</kbd>
1.  適当に保存してブラウザを閉じる
1.  ターミナルに戻って <kbd>control-c</kbd> で終了


### キーボードショートカット

key                     | action
----------------------- | ----
<kbd>esc</kbd>          | enter command mode
<kbd>enter</kbd>        | enter edit mode
<kbd>h</kbd>            | show keyboard shortcuts
<kbd>y</kbd>            | to code
<kbd>m</kbd>            | to markdown
<kbd>a</kbd>            | insert cell above
<kbd>b</kbd>            | insert cell bellow
<kbd>dd</kbd>           | delete selected cells
<kbd>ctrl-return</kbd>  | run selected cells


### 出力

例えばある Pandas DataFrame `df` を表示したいとき、
単に `df` で評価するのと明示的に `print(df)` を実行するのとでは結果が異なる。
前者はコードセルの出力として扱われ、HTML+CSSで描画される。
後者は副作用扱いで出力番号も与えられず、普通のテキストで表示される。
テキストで出力するオプションが欲しいけど見つからない。

`nbconvert` による変換時は
`~/.jupyter/jupyter_nbconvert_config.py`
に以下のように書いておけば `print()` 無しでもテキスト表示できる。

```py
c.NbConvertBase.display_data_priority = ['text/plain']
```

[Jupyter’s Common Configuration Approach](https://jupyter.readthedocs.io/en/latest/use/config.html)
によれば代入ではなく `.prepend()` でも行けそうだが、
それだとなぜか `text/html` が優先されてしまう。
また、いくつもある `display_data_priority` の中でなぜ
`NbConvertBase` だけが効くのか、というのもよく分からない。


## Jupytext

<https://jupytext.readthedocs.io/>

`.ipynb`, `.py`, `.md` などの形式を相互に変換してくれる。
例えば、手元のソースコードは `.py` としてバージョン管理・編集し、
配布時に `.ipynb` で出力といった使い方ができる。

Jupytextではソースコードを同期することが主眼なので、
knitrのようにコードセルの実行結果を含むMarkdownを書き出す機能は無い。
それは
`jupyter nbconvert --execute --to markdown`
とかでやるべき仕事っぽい。
[jupytext/issues/220](https://github.com/mwouts/jupytext/issues/220)
で議論はある。

### Format

`md` (Jupytext Markdown)
: YAMLヘッダーにメタデータを保存する。
: コードセルにオプションを付けられる。
: それ以外はほぼ普通のMarkdown。
: 現状のVSCodeではコードセル内での補完が貧弱なのは残念だが、
  リポート目的で非コードが多いならこれが扱いやすそう。

`rmarkdown` ([R Markdown](https://rmarkdown.rstudio.com/))
: コードセルに波括弧を付ける: `{python}`
: Rから[knitr](https://yihui.org/knitr/)を使えば結果を含むMarkdownを出力できる。
  RとPythonの橋渡しは[reticulate](https://rstudio.github.io/reticulate/)が担う。

`md:myst` (MyST Markdown; Markedly Structured Text)
: CommonMarkに準拠しつつreStrucutredTextやSphinxの機能をサポートするリッチなMarkdown。
: コードセルの中にメタデータを埋め込むのが難点。

`md:pandoc` (Pandoc Markdown)
: Pandoc divs `:::` という謎要素でセルを区切るらしい。

`py:light`
: `# コメント` のかたまりをテキストセルとして扱う
: コードと隣接するコメントはコードセル内に含められる。

`py:percent`
: セルの頭に明示的に `# %%` を書く。
: VSCode CodeLens でも認識されるが、むしろ表示が邪魔なのでオフにする。
: 現状のVSCodeでコードを書くにはこれが一番扱いやすいか。

`py:nomarker`
: 情報が落ちすぎて round-trip 不可。

### Config

[`${XDG_CONFIG_HOME}/jupytext/jupytext.toml`](https://github.com/heavywatal/dotfiles/blob/master/.config/jupytext/jupytext.toml)



## 書籍

<a href="https://www.amazon.co.jp/dp/487311845X/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=72a416f5d10a9e84aaab4b3ee9613329&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311845X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=487311845X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873118417/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=6b1a04ec880b6c730bd6e80273e30e9c&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873118417&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873118417" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873117488/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=2181a50362009e68f507d44fc38716b4&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
