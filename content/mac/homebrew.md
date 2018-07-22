+++
title = 'Homebrew'
tags = ["mac", "package"]
[menu.main]
  parent = "mac"
+++

Unixツールをパッケージとして手軽にインストールできるMac用パッケージ管理ソフト。

<http://brew.sh/>

## Installation

<http://docs.brew.sh/Installation.html>

1.  Command Line Tools をインストールする。
    cf. [/dev/devenv]({{< relref "devenv.md" >}})

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

http://docs.brew.sh/FAQ.html

- Homebrew本体とカタログをアップデートし、アップグレード可能なパッケージを表示:

        % brew update && brew outdated

- `outdated` なものを全てアップグレード:

        % brew upgrade

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
    automake
    binutils
    boost
    cmake
    coreutils
    diffutils
    doxygen
    eigen
    emacs
    findutils
    git
    gnu-sed
    gnu-tar
    go
    grep --with-default-names
    imagemagick
    less
    make --with-default-names
    mercurial
    nano
    nkf
    pandoc
    pyenv
    rmtrash
    rsync
    tmux
    tree
    wget
    xz
    zopfli
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
Rをここからインストールするとバイナリ版のパッケージが利用できず、
毎回ソースからビルドすることになるので、
後述のcaskで `r-app` を入れるほうが簡単。
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

- http://docs.brew.sh/brew-tap.html
- http://docs.brew.sh/Interesting-Taps-and-Forks.html

明示的にリポジトリを追加する:

    % brew tap brewsci/bio
    % brew install libsequence

暗黙に `brew tap` しつつ直接インストールも可能:

    % brew install brewsci/bio/libsequence

バイオインフォマティクスなど科学計算のツール群はHomebrew公式タップ
[`homebrew/science`](https://github.com/Homebrew/homebrew-science)
に収録されていたがdeprecatedになった。
[`brewsci/science`](https://brewsci.github.io/homebrew-bio/)
がフォーミュラを一旦引き継いで、 `homebrew/core` と
[`brewsci/bio`](https://brewsci.github.io/homebrew-bio/)
に振り分けて移行を進めている。


### Tapを作る

https://docs.brew.sh/How-to-Create-and-Maintain-a-Tap.html

GitHubに `homebrew-nicetap` のような名前のリポジトリを作成し、
ルート直下の `Formula/` に `goodtool.rb` のようなファイルを置くだけ。
使うときは `brew install <username>/nicetap/goodtool`
のようにリポジトリ名から `homebrew-` を削った形で指定する。


### Formulaを作る

https://docs.brew.sh/Formula-Cookbook.html

新規作成するには `brew create <URL>` コマンドが便利。

`url` には [`git tag`]({{< relref "git.md#tag" >}})
でGitHubに作られるバージョン付きアーカイブを指定してやるのが楽ちん。

`head` にリポジトリを登録してバージョン無しで運用することも可能。
ただし `brew upgrade` ではHEADの更新をチェックしてくれないので要注意。
`brew reinstall` するしかないのかな？


## Cask

<http://caskroom.io/>

GUIアプリケーションもHomebrewで管理してしまおうという野心的な拡張機能。

インストールは1行:

    % brew tap caskroom/cask

使うときは普通の `brew` コマンドに `cask` を挟むだけ:

    % brew cask install libreoffice
    % brew cask list

アプリ側でアップデートを実行するとCask内でのバージョンと食い違っちゃうけど使用上は問題ないらしい。

alfred, amazon-drive, amazon-music, atom,
basictex, betterzipql, bibdesk, boostnote, caffeine, cmd-eikana,
dia, dropbox, firefox, gephi, gitter,
google-backup-and-sync, google-chrome, google-earth, google-japanese-ime,
inkscape, iterm2, kindle, libreoffice, megasync,
openoffice, osxfuxe, qlstephen, quicklook-csv, r-app, rstudio,
skim, skyfonts, skype, slack, spideroakone, sshfs,
the-unarchiver, virtualbox, vlc,
whatsapp, xquartz, yujitach-menumeters
