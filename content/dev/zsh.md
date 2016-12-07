+++
title = 'zsh'
tags = ["shell"]
[menu.main]
  parent = "dev"
+++

## The most powerful shell

-   <http://www.zsh.org/>
-   <http://zsh.sourceforge.net/Doc/>

## Environmental variables

`ZDOTDIR`
:   設定ファイルを読み込むディレクトリ。デフォルトは `${HOME}`

`ZSH_VERSION`
:   `.zshrc` とかで条件分岐するのに使える

`fpath`
:   zsh関数や補完関数のパス

`HISTFILE`
:   ヒストリーを保存するファイル

## Configuration files

[L] ログインシェルとして実行時\
[Z] 非ログインで `zsh` 起動時\
[S] シェルスクリプト実行時

`.zshenv` [LZS]
:   スクリプトの実行にも必要な環境変数(`PATH` とか)の指定

`.zprofile` [L]
:   ログインシェルとして使うのに必要な設定

`.zshrc` [LZ]
:   インタラクティブシェルとして使うのに必要な設定

`.zlogin` [L]
:   `.zshrc` より後に読まれる以外は `.zprofile` と同じ。

`.zlogout`
:   ログアウト時にしてほしいことが万が一あれば

{{%div class="note"%}}
`$ZDOTDIR/` 以下の個人設定ファイルの前に、
システム全体の設定ファイルとして `/etc/z*` が読み込まれる。
これがしばしば問題を起こすので、とりあえず `unsetopt GLOBAL_RCS` で切っとく。

例えばMacでは `/usr/libexec/path_helper` が
`/usr/bin` などの基本的なPATHを設定してくれる。
Yosemiteまでは `/etc/zshenv` で実行されていたが、
El Capitanからは `/etc/zprofile` に変更されたせいで
`~/.zshenv` での設定がうまく反映されなくなってしまった。
これを防ぐには以下のようにする:
```sh
# ~/.zshenv
unsetopt GLOBAL_RCS
if [ $(uname) = Darwin ]; then
    PATH=''
    eval $(/usr/libexec/path_helper -s)
fi
PATH=${HOME}/local/bin:${PATH}
export PATH
```
{{%/div%}}

## `~/.zsh/` 以下にまとめる

1.  ディレクトリを作ってその中に設定ファイルを入れる:

        % mkdir ~/.zsh

2.  `~/.zsh/.zshenv` に以下の内容を記述して `ZDOTDIR` を設定:

        export ZDOTDIR=${HOME}/.zsh

3.  `~/.zshenv` にシンボリックリンクを張る:

        % cd
        % ln -s .zsh/.zshenv

4.  `zsh` 起動
    1.  `~/.zshenv -> ~/.zsh/.zshenv` が読み込まれる
    2.  `ZDOTDIR=${HOME}/.zsh` が設定される
    3.  `$ZDOTDIR` 以下の設定ファイルが読み込まれる

## Installation

Macでは [Homebrew]({{< relref "mac/homebrew.md" >}})` (あるいは `[MacPorts]({{< relref "mac/macports.md" >}})) を使うと良い:

    % brew install zsh --without-etcdir
    % brew install zsh-completions
