+++
title = 'Homebrew'
tags = ["mac", "package"]
[menu.main]
  parent = "mac"
+++

Unixツールをパッケージとして手軽にインストールできるMac用パッケージ管理ソフト。

<http://brew.sh/>

## Installation

https://github.com/Homebrew/brew/blob/master/share/doc/homebrew/Installation.md

1.  Command Line Tools をインストールする。 cf. [/dev/devenv]({{< relref "dev/devenv.md" >}})
2.  公式では `/usr/local/` へのインストールが推奨されているが、
    個人的にあまり好ましくないのでホームディレクトリに
    `~/.homebrew/` を作ってインストールする:

        % cd
        % mkdir .homebrew
        % curl -L https://github.com/mxcl/homebrew/tarball/master | tar xz --strip 1 -C .homebrew

3.  `.zshrc` を書き換えてシェル環境を整える:

        # パスを通す
        export PATH=$HOME/.homebrew/bin:$PATH
        export MANPATH=$HOME/.homebrew/share/man:$MANPATH
        if [ -d $HOME/.homebrew/share/zsh-completions ]; then
            fpath=($HOME/.homebrew/share/zsh-completions $fpath)
        fi

        # お好みでデフォルトコンパイラを変更
        export HOMEBREW_CC=clang

## Usage

https://github.com/Homebrew/brew/blob/master/share/doc/homebrew/FAQ.md

- Homebrew本体とカタログをアップデートし、アップグレード可能なパッケージを表示:

        % brew update && brew outdated

- `outdated` なものを全てアップグレード:

        % brew upgrade --all

- パッケージのバージョンを固定し、`brew upgrade --all` の適用外にする。
  頻繁に更新され、やたらCPUを使うやつらに。

        % brew pin imagemagick

- パッケージ検索:

        % brew search text

- パッケージ情報の表示:

        % brew info formula

- パッケージのインストール・アンインストール:

        % brew install formula
        % brew uninstall formula

- インストール済みパッケージ、またはパッケージ内ファイルの一覧:

        % brew list [formula]


## brew install

公式リポジトリから明示的にインストールしたものメモ:

    aspell
    autoconf
    binutils
    boost
    colordiff
    coreutils
    doxygen
    emacs
    fftw
    findutils
    gcc
    gibo
    git --without-completions
    gnu-sed
    gnu-tar
    gnuplot
    graphviz --with-bindings
    gsl
    jags
    lesspipe
    mercurial
    nkf
    rmtrash
    tmux
    tree
    wakeonlan
    wget
    xz
    zsh --without-etcdir
    zsh-completions

{{%div class="note"%}}
コンパイル済みパッケージ(bottle)が提供されてる場合はダウンロードだけで済むが、
オプションを付けるとソースからビルドする羽目になるので注意。
gccやboostなどのデカいやつはとりあえずデフォルトで入れたほうが早い。

`--force-bottle` と `--build-from-source` で明示的に切り替えられる。
{{%/div%}}

{{%div class="note"%}}
`coreutils`, `gnu-sed`, `gnu-tar`
などは既存のコマンドとごっちゃにならないよう頭に `g`
を付けてインストールしてくれる。
元の名前でアクセスする方法はいくつかあるが、
`$(brew --prefix)/opt/{{coreutils,gnu-sed,gnu-tar}}/libexec/gnubin` に
`PATH` を通すのがよい。
`brew unlink coreutils gnu-sed gnu-tar` してもそれらのディレクトリは残る。
{{%/div%}}

## brew tap

https://github.com/Homebrew/brew/blob/master/share/doc/homebrew/brew-tap.md

https://github.com/Homebrew/brew/blob/master/share/doc/homebrew/Interesting-Taps-&-Branches.md

非公式フォーミュラを公開しているリポジトリを追加する:

    % brew tap homebrew/science
    % brew install samtools

`brew tap` せずに直接インストールも可能:

    % brew install homebrew/science/samtools

-   `brew tap homebrew/dupes` - <https://github.com/Homebrew/homebrew-dupes>\
    システムに既にあるものと重複してでも新しいのを持っておきたいツール群。例えば

        diffutils
        grep --with-default-names
        make
        nano
        rsync

-   `brew tap homebrew/versions` - <https://github.com/Homebrew/homebrew-versions>\
    既に存在しているパッケージのバージョン違いを提供してくれている。例えば

        gcc5

-   `brew tap homebrew/science` - <https://github.com/Homebrew/homebrew-science>\
    バイオインフォマティクスなど科学計算のツール群。例えば

        bcftools
        blast
        bowtie2
        bwa
        cd-hit
        clustal-w
        cufflinks
        emboss
        fastqc
        fwdpp
        igv
        igvtools
        libsequence
        mafft
        paml
        phylip
        repeatmasker
        samtools
        snpeff
        tophat
        varscan

-   `brew tap homebrew/python` - <https://github.com/Homebrew/homebrew-python>\
    ライブラリ依存性などにより [pip]({{< relref "python/pip.md" >}}) からインストールしにくいPythonライブラリ。例えば

        matplotlib
        numpy
        pillow
        scipy

## Cask

<http://caskroom.io/>

GUIアプリケーションもHomebrewで管理してしまおうという野心的な拡張機能。

インストールは1行:

    % brew tap caskroom/cask

使うときは普通の `brew` コマンドに `cask` を挟むだけ:

    % brew cask install libreoffice
    % brew cask list

アプリ側でアップデートを実行するとCask内でのバージョンと食い違っちゃうけど使用上は問題ないらしい。

alfred, amazon-cloud-drive, atom, audacity,
basictex, bibdesk, caffeine,
dia, dropbox, evernote, firefox, gephi,
google-chrome, google-drive, google-earth, google-japanese-ime,
inkscape, iterm2, karabiner, kindle, libreoffice,
macfusion, megasync, mendeley-desktop,
onyx, openoffice, osxfuse, pandoc, quicksilver,
r, rstudio, seashore, skim, skitch, skype, spideroakone, sshfs,
the-unarchiver, vlc, yujitach-menumeters
