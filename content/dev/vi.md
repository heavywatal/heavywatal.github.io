+++
date = "2016-05-20T00:14:03+09:00"
tags = ["editor"]
title = "vi"

[menu.main]
  parent = "dev"
+++

## 基礎知識

オリジナルがvi、改良版がvim (Vi IMproved)。
ただし、shの正体がbashであるように、LinuxやMacに入ってるviの正体はvim。
vimのデフォルトは拡張を無効にしたvi互換モード。
つまり設定をいじらなければ起動コマンドはviでもvimでも同じ？

文字入力はインサートモード、それ以外はすべてノーマルモードで行う。
起動時はノーマルモード。

## ノーマルモード

### とりあえず覚えてないと死ぬやつ

key             | action
--------------- | ------
<kbd>hjkl</kbd> | <kbd>←↓↑→</kbd>
<kbd>i</kbd>    | インサートモードに移行
<kbd>:w</kbd>   | 保存
<kbd>:q</kbd>   | 終了
<kbd>:q!</kbd>  | 保存しないで強制終了
<kbd>u</kbd>    | 元に戻す
<kbd>x</kbd>    | 1文字削除 (行は消せない)
<kbd>dd</kbd>   | 行カット
<kbd>ctrl-c</kbd> | 中断 (Emacs <kbd>ctrl-g</kbd>)


### 移動

key           | action
------------- | ------
<kbd>b</kbd>  | 単語頭
<kbd>e</kbd>  | 単語末
<kbd>w</kbd>  | 単語末の空白
<kbd>0</kbd>  | 行頭
<kbd>^</kbd>  | 行頭 (非空白)
<kbd>$</kbd>  | 行末
<kbd>H</kbd>  | 画面先頭
<kbd>L</kbd>  | 画面末尾
<kbd>gg</kbd> | ファイル先頭
<kbd>G</kbd>  | ファイル末尾
<kbd>ctrl-b</kbd> | 前ページ
<kbd>ctrl-u</kbd> | 半ページ上
<kbd>ctrl-d</kbd> | 半ページ下
<kbd>ctrl-f</kbd> | 次ページ

数字と組み合わせられる。
e.g., <kbd>3j</kbd>で3行下

<kbd>shift-b</kbd>, <kbd>shift-e</kbd>, <kbd>shift-w</kbd>
はドットなどの区切り文字をすっ飛ばして大きめ移動。

カーソルの移動だけでなく、コピーなどの動作対象の指定にも使う。

### コピペ

key           | action
------------- | ------
<kbd>d</kbd>  | カット
<kbd>y</kbd>  | コピー
<kbd>shift-p</kbd> | ペースト

カットとコピーは対象を指定する必要がある。

```
d3l  # 右に3文字カット; 3dl でもいい
yy   # 1行まるまるコピー
```

ただの<kbd>p</kbd>でもペーストだが、カーソルの右に入るのが気持ち悪い。

範囲選択を見ながら操作したい場合は下記のヴィジュアルモードを使う。


## インサートモード `-- INSERT --`

キー入力しかしない。
カーソルキーで移動はできるが、基本的にはしないつもりで。

ノーマルモードへの戻り方
: <kbd>esc</kbd>: 標準だが遠すぎるので却下。
: <kbd>ctrl-[</kbd>: それなりに押しやすく、<kbd>esc</kbd>と同じ挙動。
: <kbd>ctrl-c</kbd>: 最も押しやすいし覚えやすい。
  インサートモードでのあらゆる動作を中断して戻るので注意。


## ヴィジュアルモード `-- VISUAL --`

<kbd>v</kbd> で始まる範囲選択モード。
Emacsでいう<kbd>ctrl-space</kbd>。

<kbd>shift-v</kbd> で行単位選択。
<kbd>ctrl-v</kbd> で矩形選択。

