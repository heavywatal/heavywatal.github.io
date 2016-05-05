+++
title = 'IPython'
subtitle = "強力な対話型Python実行環境"
tags = ["python"]
[menu.main]
  parent = "python"
  weight = -90
+++

<http://ipython.org/>

<http://ipython.readthedocs.org/>

`pip install ipython jupyter` でインストール

## 対話型実行環境

### 起動

ターミナルから `ipython` を実行するか、
普通の対話Pythonから:

    import IPython
    IPython.start_ipython(argv=[])

とりあえずドキュメントを読む:

    ?
    %quickref

### 文脈を考慮したタブ補完

定義済み変数などはある程度 `rlcompleter` でも補完できるが、
`IPython` はさらに文脈を考慮して賢い補完をしてくれる:

    import o[TAB]
    import os.[TAB]
    open('~/.[TAB]')

### 履歴

-   上下キーで単純に遡る
-   途中まで入力して `control + p` で前方一致する履歴のみ遡る
-   `control + r` から部分一致する履歴を検索
-   `%hist`
-   input cache: `_i`, `_ii`, `_iii`, `_ih[n]`, `_i<n>`
-   output cache: `_`, `__`, `___`, `_oh[n]`, `_<n>`
-   directory: `_dh`

### Object introspection

<http://ipython.readthedocs.org/en/stable/interactive/reference.html#dynamic-object-information>

関数がどんな引数をとるか、
クラスがどんなメンバを持っているか、
などをパッと覗くことができる。
`help()` コマンドの強力版。

先頭(か末尾)にクエスチョンをつけるだけ:

    ?os
    ??os

-   `%pdoc <object>`: docstring
-   `%pdef <object>`: ?
-   `%psource <object>`: 定義しているソースコード
-   `%pfile <object>`: 定義しているファイル全体

### システムコマンド実行

頭にエクスクラメーションをつけるだけ:

    !pwd
    files = !ls

### マジックコマンド

<http://ipython.readthedocs.org/en/stable/interactive/magics.html>

一行単位の line magic は `%` で始める。

複数行うけつける cell magic は `%%` で始める。

しかしautomagic設定により省略可能。

## Jupyter Notebook

<http://jupyter.org/>

<http://jupyter.readthedocs.org/>

<http://jupyter-notebook.readthedocs.org/>

ウェブブラウザ上で動く対話的実行環境。
Markdown/LaTeX記法による見出し・コメント・数式とともに
プログラムと実行結果をひと括りに保存しておけるので、
研究ノートや共同研究者へのレポートとして便利。
Rでいうknitr + Rmarkdownのような位置づけだが、
インプットとアウトプットの近さという点で
Mathematica/Mapleの使い勝手に似ている。

元はIPython Notebookだったが、
カーネルを入れ替えることで他言語サポートが可能になり、
汎用Notebookとして生まれ変わった。

### 始め方

<http://jupyter-notebook.readthedocs.org/en/latest/examples/Notebook/rstversions/Notebook%20Basics.html>

1.  `IPython` とともにインストール
2.  ターミナルから起動:

        % jupyter-notebook [file or directory]

3.  ウェブブラウザで `http://localhost:8888/tree` が立ち上がる
4.  右上の New から適当なNotebookカーネル
    (e.g., Python 3) を選択
5.  In [ ]: の右の箱に適当なPythonコマンドを入れて
    `shift + return`
6.  左上の File から適当に保存してブラウザを閉じる
7.  ターミナルから `control + c` で終了

### 設定

<http://jupyter-notebook.readthedocs.org/en/latest/config.html>

<http://jupyter-notebook.readthedocs.org/en/latest/frontend_config.html>
