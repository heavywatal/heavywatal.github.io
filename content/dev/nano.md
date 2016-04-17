+++
title = 'nano'
subtitle = "ターミナルで動く軽量テキストエディタ"
tags = ["editor"]
[menu.main]
  parent = "dev"
+++

<http://www.nano-editor.org/>

emacs や vi ほど高機能ではないけど敷居が低い。

## Usage

ターミナルから起動:

    % nano

下の方に保存や終了のコマンドが書いてあるので
emacs や vi のように迷うことは無い。

## Installation

Linux にも Mac にも最初からインストールされている。

ただMac付属の `/usr/bin/nano` は古いし
色設定ファイルも含まれていないので
パッケージマネージャーで新しいのを入れる。
cf. [/mac/homebrew]({{< relref "mac/homebrew.md" >}}) :

    % brew install nano

コンパイルオプションを確認してみる:

    (Mac)% /usr/bin/nano -V
     GNU nano version 2.0.6 (compiled 20:02:10, Jun 24 2011)
     Email: nano@nano-editor.org    Web: http://www.nano-editor.org/
     Compiled options: --disable-nls --enable-color --enable-extra --enable-multibuffer --enable-nanorc --enable-utf8

    (Mac)% nano -V
     GNU nano version 2.2.6 (compiled 19:50:42, Jun 14 2012)
     (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,
     2008, 2009 Free Software Foundation, Inc.
     Email: nano@nano-editor.org    Web: http://www.nano-editor.org/
     Compiled options: --disable-nls --enable-color --enable-extra --enable-multibuffer --enable-nanorc --enable-utf8

## Configuration

各種設定は `~/.nanorc` ファイルで。

Linuxでは `/usr/share/doc/nano/examples/nanors.samples` や
`/etc/nanorc` を自分のホームにコピーして編集する:

    % cp /etc/nanorc ~/.dotfiles/.nanorc
    % ln -s .dotfiles/.nanorc ~/
    % nano ~/.nanorc

Macの homebrew によるインストールでは生成されない？ :

    brew list nano

最新のソースコードからビルドして生成される `nanorc.sample`
を参考にするのがよさそう。

### Syntax highlight

<https://github.com/scopatz/nanorc>

有志によってメンテされている定義ファイルを利用する:

    % git clone https://github.com/scopatz/nanorc.git ~/.nano

`~/.nanorc` でワイルドカード指定:

    include "~/.nano/*.nanorc"
