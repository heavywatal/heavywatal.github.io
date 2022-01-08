+++
title = 'Homebrew'
tags = ["mac", "package"]
[menu.main]
  parent = "mac"
+++

Unixツールをパッケージとして手軽にインストールできるパッケージ管理ソフト。

<https://brew.sh/>

## Installation

1.  Command Line Tools をインストールする。
    cf. [/dev/devenv]({{< relref "devenv.md" >}}):
    ```sh
    xcode-select --install
    ```

    [Xcodeを丸ごとインストールしてある場合でも独立CLTが必要](https://github.com/Homebrew/brew/issues/11250)らしい。


1.  ターミナルから下記のコマンドを実行し、指示に従う:
    ```sh
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
    ```

1.  ちゃんと入ったか確認:
    ```sh
    brew doctor
    ```

---

<https://docs.brew.sh/Installation.html>

デフォルトの `/usr/local/` にインストールするのが嫌なら、
例えばホーム以下の `~/.homebrew/` にインストールすることもできる:
```sh
mkdir ~/.homebrew
curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1 -C ~/.homebrew
```
が、`/usr/local/` 以外にインストールすると、
bottle機能を封じられて毎回自前ビルドすることになるみたいなので、
非力なラップトップとかでは結構厳しい。


## Usage

<https://docs.brew.sh/FAQ.html>

-   Homebrew本体とカタログをアップデートし、アップグレード可能なパッケージを表示:

        brew update && brew outdated

-   `outdated` なものを全てアップグレード:

        brew upgrade

-   パッケージのバージョンを固定し、一括 `brew upgrade` の適用外にする。
    頻繁に更新され、やたらCPUを使うやつらに。

        brew pin imagemagick

-   パッケージ検索:

        brew search text

-   パッケージ情報の表示:

        brew info formula

-   パッケージのインストール・アンインストール:

        brew install formula
        brew uninstall formula

-   インストール済みパッケージ、またはパッケージ内ファイルの一覧:

        brew list [formula]


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
    grep
    imagemagick
    less
    lftp
    make
    nano
    nkf
    pandoc
    parallel
    rsync
    sshfs
    tmux
    vim
    wget
    xz
    zopfli
    zsh
    zsh-completions

Rをここからインストールするとバイナリ版のパッケージが利用できず、
毎回ソースからビルドすることになるので、
後述のcaskのほうの `r` を入れるほうが簡単。

`coreutils`, `gnu-sed`, `gnu-tar`, `grep`
などは既存のコマンドとごっちゃにならないよう頭に `g`
を付けてインストールしてくれる。
元の名前でアクセスする方法はいくつかあるが、
`$(brew --prefix)/opt/{coreutils,gnu-sed,gnu-tar,grep}/libexec/gnubin` に
`PATH` を通すのが楽ちん。
`brew unlink coreutils gnu-sed gnu-tar grep` してもそれらのディレクトリは残る。


## brew tap

- https://docs.brew.sh/brew-tap.html
- https://docs.brew.sh/Interesting-Taps-and-Forks.html

明示的にリポジトリを追加する:

```sh
brew tap brewsci/bio
brew install libsequence
```

暗黙に `brew tap` しつつ直接インストールも可能:
`brew install brewsci/bio/libsequence`

バイオインフォマティクスなど科学計算のツール群はHomebrew公式タップ
[`homebrew/science`](https://github.com/Homebrew/homebrew-science)
に収録されていたがdeprecatedになった。
[`brewsci/science`](https://github.com/brewsci/homebrew-science)
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

GUIアプリケーションもHomebrewで管理してしまおうという野心的な拡張機能。
昔は `brew cask install` のような形で使っていたが、
今は `brew` 本体に統合されている。
同名のformulaがある場合などは `--cask` で限定できる:

```sh
brew install --cask r rstudio
brew list --cask
```

アプリ側でアップデートを実行するとCask内でのバージョンと食い違っちゃうけど使用上は問題ないらしい。

adobe-acrobat-reader amazon-photos atom
basictex bibdesk
discord docker drawio dropbox firefox gimp gitter
google-chrome google-drive google-japanese-ime
inkscape iterm2 julia kindle libreoffice macfuse
marshallofsound-google-play-music-player
megasync menumeters monitorcontrol
r rstudio
skim skype slack spideroakone
the-unarchiver virtualbox visual-studio-code vlc
whatsapp xquartz zoom

### [Quicklook]({{< relref "quicklook.md" >}})

### Fonts

```sh
brew tap homebrew/cask-fonts
```

font-ubuntu font-ubuntu-mono font-ubuntu-mono-nerd-font
font-noto-sans font-noto-serif font-noto-sans-mono
font-source-sans-3 font-source-serif-4
font-open-sans font-roboto font-dejavu
font-libertinus font-lora font-merriweather
