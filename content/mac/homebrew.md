+++
title = 'Homebrew'
tags = ["mac", "package"]
[menu.main]
  parent = "mac"
+++

Unixツールを手軽にインストールできるパッケージ管理ソフト。

<https://brew.sh/>

## Installation

<https://docs.brew.sh/Installation>

1.  Command Line Tools をインストールする。
    cf. [/dev/devenv]({{< relref "devenv.md" >}}):
    ```sh
    xcode-select --install
    ```

    [Xcodeを丸ごとインストールしてある場合でも独立CLTが必要](https://github.com/Homebrew/brew/issues/11250)らしい。


1.  ターミナルから下記のコマンドを実行し、指示に従う:
    ```sh
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    ```

1.  ちゃんと入ったか確認:
    ```sh
    brew doctor
    ```


## Usage

<https://docs.brew.sh/FAQ>

-   パッケージのインストール・アンインストール:

        brew install cmake
        brew uninstall cmake

-   Homebrew本体とカタログをアップデートし、アップグレード可能なパッケージを表示:

        brew update && brew outdated

-   `outdated` なものを全てアップグレード:

        brew upgrade

-   パッケージ検索:

        brew search cmake

-   パッケージ情報の表示:

        brew info cmake

-   インストール済みパッケージ、またはパッケージ内ファイルの一覧:

        brew list [formula]


## brew install

公式リポジトリから明示的にインストールしたものメモ:\
boost
cmake
doxygen
eigen
exiftool
fzf
go
make
miller
nkf
pandoc
parallel
qpdf
rbenv
rsync
tmux
webp
wget
zsh-completions
zstd

Rをここからインストールするとバイナリ版のパッケージが利用できず、
毎回ソースからビルドすることになるので、
後述のように `--cask r` で入れるほうが簡単。

`coreutils`, `gnu-tar`
などは既存のシステムコマンドとごっちゃにならないよう頭に `g`
を付けてインストールしてくれる。
元の名前でアクセスする方法はいくつかあるが、
`$(brew --prefix)/opt/{coreutils,gnu-sed,gnu-tar,grep}/libexec/gnubin` に
`PATH` を通すのが楽ちん。
`brew unlink coreutils gnu-tar` してもそれらのディレクトリは残る。

rust製ツールもcargoで自前ビルドするよりこちらで入れてしまった方が楽ちん:\
as-tree
bat
diskus
dust
eza
fd
git-delta
hck
hexyl
hyperfine
lsd
monolith
oxipng
procs
qsv
ripgrep
sd



## brew tap

- <https://docs.brew.sh/Taps>
- <https://docs.brew.sh/Interesting-Taps-and-Forks>

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

<https://docs.brew.sh/How-to-Create-and-Maintain-a-Tap>

GitHubに `homebrew-nicetap` のような名前のリポジトリを作成し、
ルート直下の `Formula/` に `goodtool.rb` のようなファイルを置くだけ。
使うときは `brew install <username>/nicetap/goodtool`
のようにリポジトリ名から `homebrew-` を削った形で指定する。


### Formulaを作る

<https://docs.brew.sh/Formula-Cookbook>

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

アプリ側にアップデート機能があって `auto_updates true` とされているものは普通の
`brew outdated` や `brew upgrade` には出てこない。
`brew outdated --greedy-auto-updates`
で一覧して明示的に `brew upgrade <cask>` することも可能だが、
そのままCaskのバージョンから離れていっても問題ない。
`$(brew --cache)` 以下にインストーラーが保持されるので、
ストレージ不足で気になる場合は確認して消す。

aldente amazon-photos
basictex bibdesk
discord drawio dropbox equinox firefox
google-chrome google-drive
inkscape joplin julia
macfuse megasync menumeters monitorcontrol
orbstack quarto r rstudio
skim skype slack
the-unarchiver visual-studio-code vlc
wezterm xquartz zoom

### Quicklook

[See Quicklook]({{< relref "quicklook.md" >}}).

### Fonts

```sh
brew tap homebrew/cask-fonts
```

font-sf-mono font-sf-pro
font-ubuntu font-ubuntu-mono font-ubuntu-mono-nerd-font
font-ubuntu-sans font-ubuntu-sans-mono font-ubuntu-sans-nerd-font
font-noto-sans font-noto-serif font-noto-sans-mono
font-source-sans-3 font-source-serif-4
font-open-sans font-roboto font-dejavu
font-lora font-merriweather
