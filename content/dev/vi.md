+++
date = 2016-05-20T00:14:03+09:00
tags = ["editor"]
aliases = ["/dev/nano.html"]
title = "vi"

[menu.main]
  parent = "dev"
+++

## 概要

ターミナル上で軽快に動作するテキストエディタ。
文字入力はインサートモード、それ以外はすべてノーマルモードで行う。
起動時はノーマルモード。

オリジナルがvi、改良版がvim (Vi IMproved)。
ただし、shの正体がbashであるように、LinuxやMacに入ってるviの正体はvim。
vimのデフォルトは拡張を無効にしたvi互換モード。
つまり設定をいじらなければ起動コマンドはviでもvimでも同じ？

ターミナル上で設定ファイルを修正したり、
[git]({{< relref "git.md" >}}) commit でメッセージ入力したり、
といった軽い用途でしか私は使わない。
その用途でも人に紹介するなら後述の[nano](#gnu-nano)から。
それ以上のテキスト編集は[VSCode]({{< relref "vscode.md" >}})で。
[Neovim](https://neovim.io/)も気になってはいる。


## ノーマルモード

### とりあえず覚えてないと死ぬやつ

key             | action
--------------- | ------
<kbd>i</kbd>    | インサートモードに移行
<kbd>:w</kbd>   | 保存
<kbd>:q</kbd>   | 終了
<kbd>:q!</kbd>  | 保存しないで強制終了
<kbd>u</kbd>    | 元に戻す
<kbd>x</kbd>    | 1文字削除 (行は消せない)
<kbd>dd</kbd>   | 行カット
<kbd>ctrl-c</kbd> | 中断


### 移動

key               | action
----------------- | ------
<kbd>h</kbd><kbd>j</kbd><kbd>k</kbd><kbd>l</kbd> | <kbd>←</kbd><kbd>↓</kbd><kbd>↑</kbd><kbd>→</kbd>
<kbd>b</kbd> <kbd>B</kbd> | 単語頭
<kbd>w</kbd> <kbd>W</kbd> | 次の単語頭
<kbd>0</kbd>      | 行頭
<kbd>^</kbd>      | 行頭 (非空白)
<kbd>$</kbd>      | 行末
<kbd>H</kbd>      | 画面先頭
<kbd>L</kbd>      | 画面末尾
<kbd>gg</kbd>     | ファイル先頭
<kbd>G</kbd>      | ファイル末尾
<kbd>ctrl-b</kbd> | 前ページ
<kbd>ctrl-u</kbd> | 半ページ上
<kbd>ctrl-d</kbd> | 半ページ下
<kbd>ctrl-f</kbd> | 次ページ

数字と組み合わせられる。
e.g., <kbd>3j</kbd>で3行下

<kbd>e</kbd>は単語末じゃなくてそのひとつ前に移動するので気持ち悪い。

カーソルの移動だけでなく、コピーなどの動作対象の指定にも使う。

### コピペ

key          | action
------------ | ------
<kbd>d</kbd> | カット
<kbd>y</kbd> | コピー
<kbd>P</kbd> | ペースト

カットとコピーは対象を指定する必要がある。

```nohighlight
d3l  # 右に3文字カット; 3dl でもいい
diw  # カーソル位置の単語をカット
yy   # 1行まるまるコピー
```

<kbd>p</kbd>はカーソルの右に挿入される謎仕様なので却下。
<kbd>P</kbd>も挿入後のカーソル位置が気持ち悪いけど仕方ない。

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




## GNU nano

<https://www.nano-editor.org/>

機能は控えめだが操作が平易なテキストエディタ。
ターミナル初心者に勧めやすいし、軽微なサーバー作業でも頼りになる。

- 画面の下の方に保存や終了のコマンドが書いてあり、
[emacs]({{< relref "emacs.md" >}})や[vi]({{< relref "vi.md" >}})と比べて迷いにくい。
- vimと同じかそれ以上に、OS標準装備として利用可能な場合が多い。

### Installation

Linux にも Mac にも最初からインストールされている。
それが古くてどうしても気に入らない場合は
[Homebrew]({{< relref "homebrew.md" >}})
とかで新しいのを入れる:

```sh
brew install nano
nano -V
```

### Configuration

各種設定は `~/.nanorc` ファイルで。
とはいえ、フルカスタムしたところで機能は高が知れているし、
動かなくなったら困るので最低限の設定にとどめる。
