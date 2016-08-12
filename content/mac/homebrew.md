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

1.  Command Line Tools をインストールする。
    cf. [/dev/devenv]({{< relref "dev/devenv.md" >}})

2.  ターミナルから下記のコマンドを実行し、指示に従ってパスワードを入力する:
    ```sh
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    ```

    デフォルトの `/usr/local/` にインストールするのが嫌なら、
    例えばホーム以下の `~/.homebrew/` にインストールすることもできる:
    ```sh
    % cd
    % mkdir .homebrew
    % curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1 -C .homebrew
    ```
    が、`/usr/local/` 以外にインストールすると、
    bottle機能を封じられて毎回自前ビルドする仕様になってしまったので、
    非力なラップトップとかでは結構厳しい。

3.  `.zshenv` (もしくは`.zshrc`) でパスを通す:
    ```sh
    PATH=${HOME}/.homebrew/bin:/usr/local/bin:${PATH}
    brew_prefix=$(brew --prefix 2>/dev/null)
    if [ -n "${brew_prefix}" ]; then
        MANPATH=${brew_prefix}/share/man:${MANPATH}
        fpath=(${brew_prefix}/share/zsh-completions ${fpath})
    fi
    unset brew_prefix
    ```

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
    findutils
    gcc
    gibo
    git
    gnu-sed
    gnu-tar
    gnuplot
    graphviz
    gsl
    lesspipe
    mercurial
    nkf
    pandoc
    rmtrash
    tmux
    tree
    wakeonlan
    wget
    xz
    zsh
    zsh-completions

{{%div class="note"%}}
コンパイル済みパッケージ(bottle)が提供されてる場合はダウンロードだけで済むが、
オプションを付けるとソースからビルドする羽目になるので注意。
gccやboostなどのデカいやつはとりあえずデフォルトで入れたほうが早い。

`--force-bottle` と `--build-from-source` で明示的に切り替えられる。
しかしいつの間にか、Homebrewが`/usr/local/`
以外にインストールしてあるとbottleを使えない仕様になってしまった。
{{%/div%}}

{{%div class="note"%}}
`coreutils`, `gnu-sed`, `gnu-tar`
などは既存のコマンドとごっちゃにならないよう頭に `g`
を付けてインストールしてくれる。
元の名前でアクセスする方法はいくつかあるが、
`$(brew --prefix)/opt/{coreutils,gnu-sed,gnu-tar}/libexec/gnubin` に
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
        r
        repeatmasker
        samtools
        snpeff
        tophat
        varscan

-   `brew tap homebrew/python` - <https://github.com/Homebrew/homebrew-python>\
    ライブラリ依存性などにより [pip]({{< relref "python/pip.md" >}})
    からインストールしにくいPythonライブラリ。例えば

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
onyx, openoffice, osxfuse, quicksilver,
rstudio, seashore, skim, skitch, skype, spideroakone, sshfs,
the-unarchiver, vlc, yujitach-menumeters
