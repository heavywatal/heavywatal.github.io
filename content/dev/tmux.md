+++
title = 'tmux'
subtitle = "仮想端末でリモート仕事を安全に"
tags = ["job", "shell"]
[menu.main]
  parent = "dev"
+++

<http://tmux.github.io/>

[GNU screen](http://www.gnu.org/software/screen/)
の後を継ぐ端末多重化ソフト(terminal multiplexer)。

1つの画面の中でウインドウを追加・分割して複数の端末を開く
:   GUIのタブが不要になる
:    1つの `ssh` セッションで複数の端末を持てる

`ssh` 切断後も端末丸ごと継続され、後でまた繋ぎ直せる
:   不意の `ssh` 切断でも作業が失われない
:   別の端末から接続しても同じ作業を継続できる
:   `nohup` とかバックグラウンド化とか考えるより楽チン cf. [nohup]({{< relref "nohup.md" >}})

[Homebrew]({{< relref "mac/homebrew.md" >}}) あるいはLinuxbrewで一発インストール:
`brew install tmux`

## キーバインド

tmux 内で **prefix key** に続けて特定のキーを送信すると、
そのキーに応じたさまざまなコマンドを実行できる。
prefix keyはデフォルトで `C-b` だが
後述の設定で `C-t` に変更することにする。
`control + b` を `C-b` のように表記する
(e.g. `C-t ?` でキーバインドを列挙)。

```nohighlight
? 　list-keys
d 　detach-client
[ 　copy-mode
] 　paste-buffer
c 　new-window
n 　next-window
p 　previous-window
l 　last-window
, 　rename-window
" 　split-window
% 　split-window -h
↑ 　select-pane -U
↓ 　select-pane -D
← 　select-pane -L
→ 　select-pane -R
o 　select-pane -t:.+
: 　command-prompt
x 　confirm-before kill-pane
& 　confirm-before kill-window
```

### コピーモード

上に戻ってスクロールしたり、その内容をコピーしたいときはコピーモードを使う。
コピーモード中のキー操作はデフォルトでは `emacs` 風になっている。

1.  `<prefix> [` でコピーモードに入る
2.  `C-space` でコピー開始点をマーク
3.  `C-w` で終点をマークし、コピーモードを出る
4.  `<prefix> ]` でペースト

設定ファイルに
`bind-key -t emacs-copy C-w copy-pipe "pbcopy"`
と書いておけばコピー内容がMacのクリップボードにも送られる。

## 設定

設定ファイル： `$HOME/.tmux.conf`

<https://github.com/heavywatal/dotfiles/blob/master/.tmux.conf>

prefix 変更
:   `C-b` はキャレット左移動に使われるべきなので、
    `zsh` や `emacs` で使わない `C-t` に変更する。
    tmux の頭文字で覚えやすいし、`b` より若干近い。

起動時ウィンドウサイズ変更 `aggressive-resize`
:   サイズの異なる端末からアクセスしたときに随時ウィンドウサイズ変更

Mac `open` 問題
:   `tmux` 内だと `open` がうまく働かないのでそれを回避するために
    [reattach-to-user-namespace](https://github.com/ChrisJohnsen/tmux-MacOSX-pasteboard)
    をインストールして挟む:

        % brew install reattach-to-user-namespace


## 利用例

1.  リモートサーバーに ssh ログインし、
    tmux の新しいセッションを開始:

        % ssh charles
        % tmux -2u

2.  ウィンドウを縦に分割し、右ペインでPythonインタプリタを起動:

        <prefix> %

        % python

3.  左ペインにフォーカスを戻し、ファイルを閲覧したり何だり:

        <prefix> o

        % less ~/.ssh/config

4.  新しいウィンドウを作って `root` 仕事をしたり何だり:

        <prefix> c

        % su -
        Password:

5.  ウィンドウを切り替える:

        <prefix> l
        <prefix> n
        <prefix> p

6.  このセッションをデタッチし、ログアウトして家に帰る:

        <prefix> d

        % logout

7.  家からサーバーに再び ssh ログインして、
    さっきの tmux セッションをアタッチして作業を再開:

        % ssh charles
        % tmux attach -d

### 備忘

デタッチ後しばらくしてシェルを起動すると
既に存在しているセッションを忘れがちなので、
以下のようなものを `.zshrc` とかに書いておく。

```sh
tmux has-session >/dev/null 2>&1 && if [ -z "${TMUX}" ]; then
    echo '% tmux list-sessions'
    tmux list-sessions
    echo '% tmux list-windows -a'
    tmux list-windows -a
fi
```

### `ssh` 先で即 `tmux`

以下の様なzsh関数を `sshmux` とかいう名前で定義する

```sh
if [ -n "$1" ]; then
    ssh -t $* "tmux -2u attach -d || tmux -2u"
else
    echo "$0: missing hostname"
    return 1
fi
```

さらに `.zshrc` に `compdef sshmux=ssh`
と書いておけば補完もいい感じになる。

## 書籍

<a href="http://www.amazon.co.jp/gp/product/B00A4I3ZVY/ref=as_li_ss_il?ie=UTF8&camp=247&creative=7399&creativeASIN=B00A4I3ZVY&linkCode=as2&tag=heavywatal-22"><img border="0" src="http://ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=B00A4I3ZVY&Format=_SL160_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="http://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=as2&o=9&a=B00A4I3ZVY" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="http://www.amazon.co.jp/gp/product/178398516X/ref=as_li_ss_il?ie=UTF8&camp=247&creative=7399&creativeASIN=178398516X&linkCode=as2&tag=heavywatal-22"><img border="0" src="http://ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=178398516X&Format=_SL160_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="http://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=as2&o=9&a=178398516X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
