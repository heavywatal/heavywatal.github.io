+++
title = 'Emacs'
tags = ["editor"]
[menu.main]
  parent = "dev"
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

    C-g        keyboard-quit    とにかくキャンセル
    C-z        suspend-emacs    とりあえずemacsを抜ける
    C-x u      advertised-undo  元に戻す
    C-_        advertised-undo  元に戻す
    C-/        advertised-undo  元に戻す

    M-h        help    ヘルプ
    M-x                ミニバッファをコマンド受付状態にする
    M-!                ミニバッファをシェルコマンド受付状態にする

    C-x C-f    find-file                 開く
    C-x C-s    save-file                 上書き保存
    C-x C-w    write-file                別名で保存
    C-x C-c    save-buffers-kill-emacs   終了
    C-x d      dired                     directory edit に突入 (下記)

    C-s        isearch-forward     下に検索
    C-r        isearch-backward    上に検索
    M-%                            置換

    C-Space    set-mark-command    開始位置をマーク
    C-w        kill-region         マークから現在地までカット
    M-w        copy-region-as-kill マークから現在地までコピー
    C-y        yank                貼り付け

    C-x r k    kill-rectangle    矩形にカット
    C-x r y    yank-rectangle    矩形にペースト
    C-x r o    open-rectangle    矩形にスペース

    C-q TAB                      タブコード \t 入力
    C-x TAB    indent-rigidly    選択された領域を左右キーで手動インデント

    C-x k      kill-buffer       バッファを消す = ファイルを閉じる
    C-x b      switch-to-buffer  バッファを切り替える
    C-x C-b    list-buffers      バッファリストを開く(下記)

    C-x 0      このwindowを閉じて分割解除
    C-x 1      分割解除してこのwindowを最大化
    C-x 2      上下分割
    C-x 3      左右分割
    C-x o      other (next) windowにフォーカスを移す

    C-x 5 0    delete-frame
    C-x 5 1    delete-other-frame
    C-x 5 2    make-frame-command
    C-x 5 o    other-frame
    C-x 5 b    switch-tobuffer-other-frame
    C-x 5 f    find-file-other-frame
    C-x 5 d    dired-other-frame

### `list-buffers` <kbd>C-x C-b</kbd>

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

### `dired` <kbd>C-x d</kbd>

ファイラのようなバッファ。
コピーや移動もできるけど、まあそれはシェルからやればいい。

<kbd>o</kbd> 別バッファで開く\
<kbd>f</kbd> そのバッファで開く <kbd>return</kbd>\
<kbd>v</kbd> 見てみる\
<kbd>q</kbd> 終了

### 繰り返し入力

<kbd>C-u N X</kbd> あるいは <kbd>M-N X</kbd> で、XをN回入力する。

e.g., <kbd>C-u 79 -</kbd> と打てば水平線を入力できる。

e.g., <kbd>C-u 3 C-_</kbd> とすれば3操作分だけ元に戻せる。

## モード

<kbd>C-c</kbd> モード特有コマンドのprefix\
<kbd>M-;</kbd> モードに従ってコメント記号の自動挿入\
<kbd>M-/</kbd> コード補完

`list-faces-display`

`customize-face`

### Markdown

<kbd>C-c C-c p</kbd> preview on browser\
<kbd>C-c C-c m</kbd> preview on buffer\
<kbd>C-c C-c v</kbd> write preview and open in browser\
<kbd>C-c C-c e</kbd> write preview

### R (ESS)

<kbd>M-x R</kbd> R起動\
<kbd>C-c C-r</kbd> `ess-eval-region`\
<kbd>C-c C-q</kbd> R終了

## 設定

`~/.emacs.d/init.el`

<https://github.com/heavywatal/dotfiles/blob/master/.emacs.d/init.el>

## Caskでパッケージ管理

-   <https://github.com/cask/cask>
-   <http://cask.readthedocs.org/>

Emacs 24以上じゃないと動かない (=CentOS 6.5ではダメ)

`cask` コマンドは基本的に `Cask`
ファイルが置いてあるディレクトリ (e.g., `~/.emacs.d/`)
で実行することが想定されている。
それ以外の場所で実行するときは `--path ~/.emacs.d/`
のようなオプションを付ける。

### セットアップ

<http://cask.readthedocs.org/en/latest/guide/installation.html>

1.  Cask本体を `~/.cask` にインストール:

        % curl -fsSL https://raw.githubusercontent.com/cask/cask/master/go | python

    実行ファイルは `~/.cask/bin/` に置かれるので、
    `.zshrc` などで `PATH` を設定するか、
    既に通ってるところにシムリンクを張るか、フルパスで使う。

2.  `~/.emacs.d/Cask` を書く。 cf. <https://github.com/heavywatal/dotfiles/blob/master/.emacs.d/Cask>
3.  `~/.emacs.d/Cask` に書かれたパッケージをインストール:

        % cask install --path ~/.emacs.d/

4.  `~/.emacs.d/init.el` に設定を書く:

        (require 'cask "~/.cask/cask.el")
        (cask-initialize)

### メンテナンス

パッケージをアップデート:

    % cask update --path ~/.emacs.d/

Cask本体をアップデート:

    % cask upgrade
