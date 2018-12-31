+++
title = 'nano'
subtitle = "ターミナルで動く軽量テキストエディタ"
tags = ["editor"]
[menu.main]
  parent = "dev"
+++

<http://www.nano-editor.org/>

[emacs]({{< relref "emacs.md" >}})や[vi]({{< relref "vi.md" >}})ほど高機能ではないけど敷居が低い。

## Usage

ターミナルから起動:

    nano

下の方に保存や終了のコマンドが書いてあるので
emacs や vi のように迷うことは無い。

## Installation

Linux にも Mac にも最初からインストールされている。
それが古くて気に入らない場合は
[Homebrew]({{< relref "homebrew.md" >}})
とかで新しいのを入れる:

```sh
brew install nano
nano -V
```

## Configuration

各種設定は `~/.nanorc` ファイルで。

Linuxでは `/usr/share/doc/nano/examples/nanors.samples` や
`/etc/nanorc` を自分のホームにコピーして編集する:

    cp /etc/nanorc ~/.dotfiles/.nanorc
    ln -s .dotfiles/.nanorc ~/
    nano ~/.nanorc

Macの homebrew によるインストールでは生成されない？ :

    brew list nano

最新のソースコードからビルドして生成される `nanorc.sample`
を参考にするのがよさそう。

### Syntax highlight

<https://github.com/scopatz/nanorc>

有志によってメンテされている定義ファイルを利用する:

    git clone https://github.com/scopatz/nanorc.git ~/.nano

`~/.nanorc` でワイルドカード指定:

    include "~/.nano/*.nanorc"
