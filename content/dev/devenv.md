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

`sudo apt-get install`:

    build-essential
    zsh
    tmux
    git

これらシステム標準のものが古すぎたり、
管理者権限がなくて自由にインストールできない場合は
次の Linuxbrew を利用してユーザのホームに入れる。

### Linuxbrew

<http://linuxbrew.sh/>

Macの [Homebrew]({{< relref "homebrew.md" >}}) をLinuxに移植したパッケージマネージャ。

1.  RHEL/CentOS 6系の場合まずlibcurlが古すぎるので、
    [最新のcurl](https://curl.haxx.se/download.html)をソースコードからインストール:
    ```sh
    % wget -O- https://curl.haxx.se/download/curl-7.59.0.tar.gz | tar xz
    % cd curl-7.59.0/
    % ./configure --prefix=${HOME}/opt/local
    % make -j4
    ```

1.  `git --version` を確認して 1.7.12 未満だったら
    [最新のgit](https://github.com/git/git/releases)をソースコードからインストール:
    ```sh
    % wget -O- https://github.com/git/git/archive/v2.16.3.tar.gz | tar xz
    % cd git-2.16.3/
    % autoreconf -i
    % ./configure --prefix=${HOME}/opt/local --with-curl=${HOME}/opt/local
    % make -j4
    % make install
    ```

1.  上記の自前curl/gitを利用するために環境変数をセット:
    ```sh
    export PATH=${HOME}/opt/local/bin:$PATH
    export HOMEBREW_NO_ENV_FILTERING=1
    ```

1.  Linuxbrewを `~/.linuxbrew` にインストール:
    ```sh
    sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
    ```

1.  gccをインストール: `brew install gcc`

    OS標準のgccやglibcが古すぎて
    `glibc cannot be built with any available compilers` と怒られる場合は、
    まずglibc抜きでgccをインストールし、そのgccでglibcをビルドして、
    gccを入れ直す、という手間が必要:
    ```sh
    HOMEBREW_BUILD_FROM_SOURCE=1 brew install gcc --without-glibc
    brew install glibc
    brew remove gcc
    brew install gcc
    ```
    それでもうまくいかないときは
    `HOMEBREW_NO_ENV_FILTERING=1 brew install --force-bottle glibc`
    とかでいいのか...

1.  あとは欲しいものを `brew install ___`
    ```
    zsh --without-etcdir
    tmux
    pyenv
    boost
    ```

## Mac

### Command Line Tools

コンパイラや `make` などはOSに付いてこないので別途インストールが必要。
<https://developer.apple.com/downloads/> からダウンロードするか、
ターミナルから以下のコマンドを実行:

```sh
% xcode-select --install
```

インストールされているバージョンなどを確認するには:
```sh
% pkgutil --pkg-info=com.apple.pkg.CLTools_Executables
% clang -v
```

総合開発環境 Xcode をインストールしたければ、App Store から [Xcode](https://itunes.apple.com/jp/app/xcode/id497799835) を選択。

### パッケージ管理ツール

-   [/mac/homebrew]({{< relref "homebrew.md" >}})
-   [/mac/macports]({{< relref "macports.md" >}})

### その他のプログラム

-   [MenuMeters](https://member.ipmu.jp/yuji.tachikawa/MenuMetersElCapitan/)
-   [QuickLook plugins]({{< relref "quicklook.md" >}})
-   [`defaults`コマンドで各種設定]({{< relref "command.md#defaults" >}})

### リソースフォークを無視

`tar` などでリソースフォークを無視:

```sh
if [ $(uname) = Darwin ]; then
    export COPYFILE_DISABLE=true
    export COPY_EXTENDED_ATTRIBUTES_DISABLE=true
fi
```

## 共通

### Python

- [install]({{< relref "/python/install.md" >}})
- [pip]({{< relref "pip.md" >}})
- `Pillow` をインストールする前に:

      % sudo apt-get install libtiff5-dev libwebp-dev libfreetype6-dev liblcms2-dev libopenjpeg-dev

### C++

- [boost]({{< relref "boost.md" >}})
- [SFMT]({{< relref "random.md" >}})

### R

[/rstats/config]({{< relref "/rstats/config.md" >}})

### エディタ

- [atom]({{< relref "atom.md" >}})
- [emacs]({{< relref "emacs.md" >}})
- [nano]({{< relref "nano.md" >}})

### Trash

`rm` はゴミ箱を経由せず削除してしまうので、
間違って消してしまっても基本的には元に戻せない。
以下のような対策によりその危険が少しは減るかも。

`alias rmi='rm -i'`
:   ホントに消していいかどうか確認してくれるようなオプションつきのエイリアスを
    `.zshrc` に設定しておく。

        % rmi .DS_Store
        rm: remove regular file '.DS_Store'?

    エイリアス名を `rm` そのものにしてしまうと、
    結局ろくすっぽ確認せず `y` を押す癖や、
    いちいち確認されないように `rm -rf` する癖がつくので逆に危険。
    普段は `rmi` を使う癖をつけ、必要なときたまに `rm` を使い、
    `-f` はよほどのことが無い限り使わないようにする。

`trash-cli`
:   Python製なので [pip]({{< relref "pip.md" >}}) で
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
    [Homebrew]({{< relref "homebrew.md" >}}) で `rmtrash` を入れる:

        % brew install rmtrash
