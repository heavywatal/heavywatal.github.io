+++
title = '開発環境'
tags = ["package"]
[menu.main]
  parent = "dev"
+++

## Linux

### ユーザ追加

ユーザーを作って管理権限を与える:

    # adduser USERNAME
    # gpassed -a USERNAME sudo

### すぐ入れるパッケージ

[/linux/apt]({{< relref "linux/apt.md" >}})

[/linux/centos]({{< relref "linux/centos.md" >}})

`sudo apt-get install`:

    build-essential
    zsh
    tmux
    git
    mercurial

これらシステム標準のものが古すぎたり、
管理者権限がなくて自由にインストールできない場合は
次の Linuxbrew を利用してユーザのホームに入れる。

### Linuxbrew

<https://github.com/Linuxbrew/linuxbrew>

Macの [Homebrew]({{< relref "mac/homebrew.md" >}}) をLinuxに移植したパッケージマネージャ。
使い勝手はほぼ一緒だけど依存関係の処理などがやや甘い。

1.  `git --version` を確認して 1.7.12 未満だったら
    [最新版](https://github.com/git/git/releases)をソースコードからインストール:

        % wget -O- https://github.com/git/git/archive/v2.8.3.tar.gz | tar xz
        % cd git-2.8.3/
        % autoreconf -i
        % ./configure --prefix=${HOME}/local
        % make
        % make install

2.  Linuxbrewをインストール:

        % git clone https://github.com/Linuxbrew/linuxbrew.git ~/.homebrew

3.  必要なもろもろをインストール:

        % brew install tmux

        git
        mercurial
        zsh --without-etcdir
        tmux
        gcc
        python3
        llvm --with-clang --with-libcxx --with-compiler-rt --with-rtti
        boost --with-c++11
        emacs

## Mac

### Command Line Tools

コンパイラや `make` などはOSに付いてこないので別途インストールが必要。
ターミナルから以下のコマンドを実行:

    % xcode-select --install

あるいは <https://developer.apple.com/downloads/> からダウンロード。

{{%div class="note"%}}
総合開発環境 Xcode をインストールしたければ

1.  App Store から
    [Xcode](https://itunes.apple.com/jp/app/xcode/id497799835)
    を選択
2.  Xcode 上での補完のためのindexingが重いらしいので切っとく:

        defaults write com.apple.dt.Xcode IDEIndexDisable 1

{{%/div%}}
### パッケージ管理ツール

-   [/mac/homebrew]({{< relref "mac/homebrew.md" >}})
-   [/mac/macports]({{< relref "mac/macports.md" >}})

### その他のプログラム

-   Activity Monitor --- OS標準
-   OS X Server
    --- [App Store](https://itunes.apple.com/jp/app/os-x-server/id537441259)
-   MenuMeters --- <http://www.ragingmenace.com/software/menumeters/>

### リソースフォークを無視

`tar` などでリソースフォークを無視:

    if [ $(uname) = Darwin ]; then
            export COPYFILE_DISABLE=true
            export COPY_EXTENDED_ATTRIBUTES_DISABLE=true
    fi

## 共通

[/dev/etc]({{< relref "dev/etc.md" >}})

### Python

[/python/install]({{< relref "python/install.md" >}})

[/python/pip]({{< relref "python/pip.md" >}})

`Pillow` をインストールする前に:

    % sudo apt-get install libtiff5-dev libwebp-dev libfreetype6-dev liblcms2-dev libopenjpeg-dev

### C++

[/cxx/gcc]({{< relref "cxx/gcc.md" >}})

[/cxx/boost]({{< relref "cxx/boost.md" >}})

[SFMT]({{< relref "cxx/random.md" >}})

### R

[/rstats/config]({{< relref "rstats/config.md" >}})

### エディタ

[/dev/emacs]({{< relref "dev/atom.md" >}})

[/dev/emacs]({{< relref "dev/emacs.md" >}})

[/dev/nano]({{< relref "dev/nano.md" >}})

### Trash

`rm` はゴミ箱を経由せず削除してしまうので、
間違って消してしまっても基本的には元に戻せない。
以下のような対策によりその危険が少しは減るかも。

`alias rm='rm -i'`
:   `rm` コマンドそのものを少しだけ安全にしておくため、
    `.zshrc` でエイリアス設定しておく。
    するとファイルごとにホントに消していいかどうか確認してくれるようになる:

        % rm .DS_Store
        rm: remove regular file '.DS_Store'?

    でも結局ろくすっぽ確認せず `y` を押す癖がついてしまうという...

`trash-cli`
:   Python製なので [pip]({{< relref "python/pip.md" >}}) で
    `pip install trash-cli` して入れる。
    すると以下のようなコマンドがインストールされる:

        restore-trash
        trash
        trash-empty
        trash-list
        trash-put
        trash-rm

    ちなみにLinuxのゴミ箱は `~/.local/share/Trash`

`rmtrash`
:   Macでも `trash-cli` を使えないことはないが、
    ゴミ箱のパスがMac標準の `~/.Trash` ではなく
    Linuxのものになってしまうので
    [Homebrew]({{< relref "mac/homebrew.md" >}}) で `rmtrash` を入れる:

        % brew install rmtrash
