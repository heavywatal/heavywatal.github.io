+++
title = 'Emacs'
tags = ["editor"]
+++

## 単語

**buffer**
:   ファイル編集領域

**minibuffer**
:   画面の下端でコマンドを受け付けるところ

**kill-ring**
:   paste board

**window**
:   tmuxでいうpane

**frame**
:   tmuxでいうwindow

## Keybind

<kbd>C-</kbd> は左手小指の <kbd>control</kbd>

<kbd>M-</kbd> は <kbd>esc</kbd> または <kbd>C-[</kbd>

key              | command           | action
---------------- | ----------------- | ------------------
<kbd>C-g</kbd>   | `keyboard-quit`   | とにかくキャンセル
<kbd>C-z</kbd>   | `suspend-emacs`   | とりあえずemacsを抜ける
<kbd>C-x</kbd><kbd>u</kbd> | `advertised-undo` | 元に戻す
<kbd>C-_</kbd>   | `advertised-undo` | 元に戻す
<kbd>C-/</kbd>   | `advertised-undo` | 元に戻す
 | |
<kbd>M-h</kbd>   | `help` | ヘルプ
<kbd>M-x</kbd>   |        | ミニバッファをコマンド受付状態にする
<kbd>M-!</kbd>   |        | ミニバッファをシェルコマンド受付状態にする
 | |
<kbd>C-x</kbd><kbd>C-f</kbd> | `find-file`               | 開く
<kbd>C-x</kbd><kbd>C-s</kbd> | `save-file`               | 上書き保存
<kbd>C-x</kbd><kbd>C-w</kbd> | `write-file`              | 別名で保存
<kbd>C-x</kbd><kbd>C-c</kbd> | `save-buffers-kill-emacs` | 終了
<kbd>C-x</kbd><kbd>d</kbd>   | `dired`                   | directory edit に突入 (下記)
 | |
<kbd>C-s</kbd>   | `isearch-forward`     | 下に検索
<kbd>C-r</kbd>   | `isearch-backward`    | 上に検索
<kbd>M-%</kbd>   |                       | 置換
 | |
<kbd>C-@</kbd>   | `set-mark-command`    | 開始位置をマーク
<kbd>C-SPC</kbd> | `set-mark-command`    | (IME切り替え)
<kbd>C-w</kbd>   | `kill-region`         | マークから現在地までカット
<kbd>M-w</kbd>   | `copy-region-as-kill` | マークから現在地までコピー
<kbd>C-y</kbd>   | `yank`                | 貼り付け
 | |
<kbd>C-x</kbd><kbd>r k</kbd> | `kill-rectangle` | 矩形にカット
<kbd>C-x</kbd><kbd>r y</kbd> | `yank-rectangle` | 矩形にペースト
<kbd>C-x</kbd><kbd>r o</kbd> | `open-rectangle` | 矩形にスペース
 | |
<kbd>C-q</kbd><kbd>TAB</kbd> |                  | タブコード \t 入力
<kbd>C-x</kbd><kbd>TAB</kbd> | `indent-rigidly` | 選択された領域を左右キーで手動インデント
 | |
<kbd>C-x</kbd><kbd>k</kbd>   | `kill-buffer`      | バッファを消す = ファイルを閉じる
<kbd>C-x</kbd><kbd>b</kbd>   | `switch-to-buffer` | バッファを切り替える
<kbd>C-x</kbd><kbd>C-b</kbd> | `list-buffers`     | バッファリストを開く(下記)
 | |
<kbd>C-x</kbd><kbd>0</kbd>   | | このwindowを閉じて分割解除
<kbd>C-x</kbd><kbd>1</kbd>   | | 分割解除してこのwindowを最大化
<kbd>C-x</kbd><kbd>2</kbd>   | | 上下分割
<kbd>C-x</kbd><kbd>3</kbd>   | | 左右分割
<kbd>C-x</kbd><kbd>o</kbd>   | | other (next) windowにフォーカスを移す
 | |
<kbd>C-x</kbd><kbd>5 0</kbd> | | delete-frame
<kbd>C-x</kbd><kbd>5 1</kbd> | | delete-other-frame
<kbd>C-x</kbd><kbd>5 2</kbd> | | make-frame-command
<kbd>C-x</kbd><kbd>5 o</kbd> | | other-frame
<kbd>C-x</kbd><kbd>5 b</kbd> | | switch-tobuffer-other-frame
<kbd>C-x</kbd><kbd>5 f</kbd> | | find-file-other-frame
<kbd>C-x</kbd><kbd>5 d</kbd> | | dired-other-frame

### `list-buffers` <kbd>C-x</kbd><kbd>C-b</kbd>

バッファリストを隣のウィンドウで開く。

<kbd>C-u</kbd> を頭につけるとファイルを開いてるbufferに限定して。

バッファメニューをそのペインで開くコマンドはずばり `buffer-menu`

<kbd>?</kbd> ヘルプ\
<kbd>o</kbd> 隣のwindowで開く\
<kbd>f</kbd> そのwindowで開く (`RET`)\
<kbd>k</kbd> 閉じる候補にマーク\
<kbd>s</kbd> 保存候補にマーク\
<kbd>x</kbd> マーク内容の実行\
<kbd>u</kbd> マーク取り消し\
<kbd>q</kbd> list終了

### `dired` <kbd>C-x</kbd><kbd>d</kbd>

ファイラのようなバッファ。
コピーや移動もできるけど、まあそれはシェルからやればいい。

<kbd>o</kbd> 別バッファで開く\
<kbd>f</kbd> そのバッファで開く <kbd>return</kbd>\
<kbd>v</kbd> 見てみる\
<kbd>q</kbd> 終了

### 繰り返し入力

<kbd>C-u</kbd><kbd>N</kbd><kbd>X</kbd> あるいは <kbd>M-N</kbd><kbd>X</kbd> で、XをN回入力する。

e.g., <kbd>C-u</kbd><kbd>79</kbd><kbd>-</kbd> と打てば水平線を入力できる。

e.g., <kbd>C-u</kbd><kbd>3</kbd><kbd>C-_</kbd> とすれば3操作分だけ元に戻せる。

## モード

<kbd>C-c</kbd> モード特有コマンドのprefix\
<kbd>M-;</kbd> モードに従ってコメント記号の自動挿入\
<kbd>M-/</kbd> コード補完

`list-faces-display`

`customize-face`

### Markdown

<kbd>C-c</kbd><kbd>C-c</kbd><kbd>p</kbd> preview on browser\
<kbd>C-c</kbd><kbd>C-c</kbd><kbd>m</kbd> preview on buffer\
<kbd>C-c</kbd><kbd>C-c</kbd><kbd>v</kbd> write preview and open in browser\
<kbd>C-c</kbd><kbd>C-c</kbd><kbd>e</kbd> write preview

### R (ESS)

<kbd>M-x</kbd><kbd>R</kbd> R起動\
<kbd>C-c</kbd><kbd>C-r</kbd> `ess-eval-region`\
<kbd>C-c</kbd><kbd>C-q</kbd> R終了

## 設定

`~/.emacs.d/init.el`

<https://github.com/heavywatal/dotfiles/blob/master/.config/emacs/init.el>

### package.el でパッケージ管理

Emacs 24で標準入りしたので、基本的にこれを使うのが良さそう。
まず`init.el`の先頭でこのように宣言

```el
(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/"))

;; Added by Package.el.  This must come before configurations of
;; installed packages.  Don't delete this line.  If you don't want it,
;; just comment it out by adding a semicolon to the start of the line.
;; You may delete these explanatory comments.
(package-initialize)
```

あとは使いたいパッケージを以下のように列挙。

```el
(package-install 'flycheck)
(package-install 'auto-complete)
(package-install 'markdown-mode)
```

initialize後に
`(package-refresh-contents)` を入れておくとパッケージ情報を最新に保てるが、
結構遅いので気が向いたときに手動でやったほうがいい:
<kbd>M-x</kbd> `package-list-packages` <kbd>U</kbd> <kbd>x</kbd>
