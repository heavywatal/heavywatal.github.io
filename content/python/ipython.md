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

開発中の[JupyterLab](https://github.com/jupyter/jupyterlab)がさらに良いらしい。

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


## 書籍

<a href="https://www.amazon.co.jp/dp/479738946X/ref=as_li_ss_il?ie=UTF8&qid=1485612008&sr=8-6&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=5ea5e48ecc83b9439f21406b6f57c062" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=479738946X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=479738946X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/IPython%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9%E3%82%AF%E3%83%83%E3%82%AF%E3%83%96%E3%83%83%E3%82%AF-%E5%AF%BE%E8%A9%B1%E5%9E%8B%E3%82%B3%E3%83%B3%E3%83%94%E3%83%A5%E3%83%BC%E3%83%86%E3%82%A3%E3%83%B3%E3%82%B0%E3%81%A8%E5%8F%AF%E8%A6%96%E5%8C%96%E3%81%AE%E3%81%9F%E3%82%81%E3%81%AE%E3%83%AC%E3%82%B7%E3%83%94%E9%9B%86-Cyrille-Rossant/dp/4873117488/ref=as_li_ss_il?_encoding=UTF8&psc=1&refRID=X16VFSS3W75RMTG7VGCH&linkCode=li3&tag=heavywatal-22&linkId=b79e2290571289b02621392257a4ac1c" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
