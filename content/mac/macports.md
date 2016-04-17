+++
title = 'MacPorts'
[menu.main]
  parent = "mac"
+++

{{%div class="note"%}}
[homebrew]({{< relref "homebrew.md" >}})

最近は [Homebrew]({{< relref "mac/homebrew.md" >}}) ばっかり使っててこっちの情報は更新してない
{{%/div%}}

## Package Manager for Mac OS X

<http://www.macports.org/>

MacPorts はLinuxで言うところの
`apt` や `yum` のようなもので、
さまざまなソフトウェアをコマンドラインから簡単に管理できるパッケージ管理システム。
ソースコードからコンパイルしてインストールという難しい手順も自動でやってくれるし、
アップデートやアンインストールもコマンドひとつで行える。
それぞれのパッケージのことを port と呼ぶ。
port は `/opt/local/` 以下にインストールされる。

## Installation

1.  Command Line Tools をインストールする。 cf. [/dev/devenv]({{< relref "dev/devenv.md" >}})
2.  [公式サイト](http://www.macports.org/) から
    OSに合ったディスクイメージ `MacPorts-***.dmg` をダウンロード
3.  それをマウントしてインストーラを実行
4.  コマンドがインストールされる `/opt/local/bin` にパスを通す。
    例えば `.zshenv` や `.bashrc` に以下のように記述。:

        export PATH=/opt/local/bin:/opt/local/sbin:${PATH}

5.  Terminal を開いて
    MacPorts 本体とportカタログをアップデート:

        % sudo port selfupdate

## Usage

ほかにもいろんなコマンドがあるけど代表的なものだけ。詳しくは `man port` で。

-   欲しいportを検索:

        % port search bio

-   気になったportの詳細を表示:

        % port info [port]

-   portをインストールする時に選択できるオプションを表示。:

        % port variants [port]

-   そのportを入れるために必要な（依存している）portを表示。
    下はrecursive版で、`--full` オプションでフル表示できる。:

        % port deps [port]
        % port rdeps [port]

-   そのportに依存しているやつらを表示。:

        % port dependents [port]
        % port rdependents [port]

-   portのインストール/アンインストール。
    そいつが依存しているほかのportも自動的にインストールされる。:

        % sudo port install [port] [+variant]
        % sudo port uninstall [port]

-   インストール済みのportを一覧表示:

        % port installed

-   MacPorts 本体とportカタログをアップデートし、
    アップデート可能なものを一覧表示:

        % sudo port selfupdate && port outdated

-   アップデート可能なものをすべてアップデート。
    古いものを自動でアンインストールするには `-u` オプションを付ける。:

        % sudo port upgrade outdated

-   インストール済みのportを再インストール:

        % sudo port -n upgrade --force [port]

-   

    過去のバージョンのportをインストール
    :   1.  <http://trac.macports.org/browser/trunk/dports> から目的のportのページを開く
        2.  右上の"Revision Log"からお目当てのバージョンのリビジョン番号を確認
        3.  以下のような svn コマンドで
            `Portfile` をダウンロードし、インストール:

                % svn checkout -r 74577 http://svn.macports.org/repository/macports/trunk/dports/shells/zsh-devel zsh-devel-4.3.11
                % cd zsh-devel-4.3.11
                % sudo port install +mp_completion +doc +examples

### pseudo-portnames

各コマンドの対象となるportを、実名だけでなく状態によってまとめて指定できる。
e.g.:

    % port list leaves
    % sudo port uninstall $(port echo inactive)

`all`
:   all the ports in each ports tree listed in sources.conf

`current`
:   the port in the current working directory.

`active`
:   set of installed and active ports.

`inactive`
:   set of installed but inactive ports.

`actinact`
:   set of installed ports that have both an active version and one or more inactive versions.

`installed`
:   set of all installed ports.

`uninstalled`
:   ports in the ports tree(s) that aren't installed.

`outdated`
:   installed ports that are out of date with respect to their current version/revision in the ports tree(s)

`obsolete`
:   set of ports that are installed but no longer exist in any port tree

`requested`
:   installed ports that were explicitly asked for.

`unrequested`
:   installed ports that were installed only to satisfy dependencies.

`leaves`
:   installed ports that are unrequested and have no dependents.

## port installed

### coreutils

MacはFreeBSDの上にできてるので各種コマンドラインツールもBSD製のものが多い。
一方LinuxのコマンドはGNU製が多いので、
同じ名前でも微妙にオプションや挙動が異なったりして戸惑うことがある。
そこで、MacPorts を使って
GNU製のコマンドラインツール群であるcoreutilsをMacに入れる。

各プログラムはデフォルトで接頭辞 `g` のついた状態で
`/opt/local/bin/` にインストールされる。
variantとして `+with_default_names` を指定してインストールすれば元の名前で入るが、
そうすべきではないと思う。
必要なものを `/usr/local/bin/` にシムリンク張るか、
`.zshrc` などにエイリアスを定義して使う。:

    % sudo port install coreutils
    % sudo ln -s /opt/local/bin/gls /usr/local/bin/ls

    alias ls="gls -vF --show-control-chars --color=auto"

### GNU tar & xz

GNU `tar` 1.22からxz圧縮をサポートするようになった。
Macにプリインストールされてる `bsdtar` や
`gnutar` (1.17) では使えない。
gzipは圧縮も展開も高速で（CPU負荷が小さくて）圧縮率が低い。
bzipは圧縮も展開も遅くて（CPU負荷が大きくて）圧縮率が高い。
xzは最も圧縮率が高く、圧縮にはbzip2よりも時間がかかる一方で展開はgzip並みに速い。
ということで、多くのユーザーに配布・展開されるようなファイルの圧縮に効果的。
あるいは、書き換える予定は無いが長期保存しておかなければいけないデカいファイルとか。:

    % sudo port install xz gnutar
    % tar cJf archive.tar.xz archive/

### zsh

Macには元から入ってるし、ソースからインストールするのも簡単だけど、
port コマンドを補完できるvariantがあるのでこれを使う。
ログインシェルにするにはひと手間必要。
その後の設定は [こちらのページ参照]({{< relref "dev/zsh.md" >}}) 。:

    % sudo port install zsh-devel +mp_completion +doc +examples
    % sudo emacs -nw /etc/shells # 末尾に/opt/local/bin/zshを追加
    % chsh -s /opt/local/bin/zsh
    % exit

### misc.

そのほか意識的に入れるもの。依存パッケージとして勝手に入るものではなく。:

    clang-3.1
    emacs
    gcc47
    gnuplot +no_x11
    graphviz +no_x11
    grep
    gsed
    llvm-3.1
    nkf
    rmtrash
    rsync
    tmux
    wakeonlan
    wget

意識的には入れないけど、ぼーっとしてるとX11関連のportと共に勝手にインストールされてしまうもの。
明示的に `+no_x11` をつけてインストールしておき、それを防ぐ。:

    cairo +no_x11 +quartz
    gd2 +no_x11
    ghostscript +no_x11

あるいは `/opt/local/etc/macports/variants.conf`
に以下のように書いておくと自動的にそうしてくれる。:

    +no_x11 +quartz
