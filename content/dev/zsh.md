+++
title = 'zsh'
tags = ["shell"]
[menu.main]
  parent = "dev"
+++

## The most powerful shell

-   <http://www.zsh.org/>
-   <http://zsh.sourceforge.net/Doc/>

## Installation

基本的にはOSに入ってる `/bin/zsh` を使う。
新しいのを入れるなら[Homebrew]({{< relref "homebrew.md" >}})を使うのが楽:

```sh
brew install zsh
brew install zsh-completions
```

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

`$ZDOTDIR` 以下の個人設定ファイルが場合に応じて下記の順で読まれる。
システム全体の設定ファイルとして `/etc/z*` が個人設定ファイルの前に読み込まれる。
`unsetopt GLOBAL_RCS` で切れるが `/etc/zshenv` だけは絶対最初。

`.zshenv`
:   スクリプトの実行時も含めてあらゆる場合に読み込まれる。
    インタラクティブ用の設定などはせず、最低限の記述に留める。
    例えば `ZDOTDIR`, `unsetopt NOMATCH` など。

`.zprofile`
:   ログインシェルとして立ち上げるときのみ読まれる。
    `export` する環境変数(`PATH` とか)を設定するのに適している。
    `.bash_profile` に対応するので共通設定を
    `.profile` に書いておいて `source` するとか。
:   例えばローカル環境Mac + リモート環境Linux CUIで開発する場合、
    ターミナルも[tmux]({{< relref "tmux.md" >}})もデフォルトでログインシェルを立ち上げるので、
    `.zshrc` に一本化してしまっても構わない。
:   使い分けるのはLinux GUIを使う場合とか、
    よほど重い初期化をログインシェル1回で済ませたい場合とか。

`.zshrc`
:   ログイン・非ログイン問わず、インタラクティブシェルとして立ち上げるときに読まれる。
    だいたいどの設定もこれに書いておけば問題ない。
:   `.zprofile` と使い分けるなら
    `setopt` や `autoload` など、親シェルから引き継がれないものはこちら。
    `alias` などは別ファイルを読み込む形にして `.bashrc` と共有。

`.zlogin`
:   `.zshrc` より後に読まれる以外は `.zprofile` と同じ。使わない。

`.zlogout`
:   ログアウト時にしてほしいことが万が一あれば。

Macでは `/usr/libexec/path_helper` が
`/usr/bin` などの基本的なPATHを設定してくれる。
Yosemiteまでは `/etc/zshenv` で実行されていたが、
El Capitanからは `/etc/zprofile` で実行されるため、
`.zshenv` で `PATH` を設定しようとするとうまく反映されない。
`unsetopt GLOBAL_RCS` で `/etc/zprofile` をスキップして `path_helper` を手動実行するか、
素直に `.zprofile` 以降のファイルで設定する。


### 起動時間短縮

まずはプロファイリングしてボトルネックを知る:

```sh
# head of .zshenv
zmodload zsh/zprof

# tail of .zshrc
zprof
```

`compinit` とかが遅かったり複数回呼ばれていたりするので順番やオプションを変えてみる。
