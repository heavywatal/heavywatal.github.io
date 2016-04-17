+++
title = 'Installation'
[menu.main]
  parent = "python"
+++

<https://www.python.org/downloads/>

## Ubuntu

`/usr/bin/python` が既にインストールされている。
他のバージョンも Synaptic や apt-get で簡単に入れられる。

新しいや古いのが欲しい場合はソースからビルドするか
[ppa:fkrull/deadsnakes](https://launchpad.net/~fkrull/+archive/deadsnakes)
などのリポジトリを登録するとよい。:

    % sudo add-apt-repository ppa:fkrull/deadsnakes
    % sudo apt-get update
    % sudo apt-get install python2.7

## Mac

`/usr/bin/python` が既にインストールされている。
他のバージョンも [Homebrew]({{< relref "mac/homebrew.md" >}})` や `[MacPorts]({{< relref "mac/macports.md" >}}) で簡単に入れられる。
昔は tk や quartz 周りで面倒があったが、
tk が別パッケージに分離されたので楽になった:

    % brew install python3
    % sudo port install python27

いろんなライブラリも提供されてるけどそれらは利用せず
[pip]({{< relref "pip.md" >}}) とかを使ったほうが良い:

    % port search py27

## Source

1.  必要なパッケージをインストールしておく

    Ubuntuなら:

        % sudo apt-get install build-essential
        % sudo apt-get install libreadline6-dev
        % sudo apt-get install libsqlite3-dev
        % sudo apt-get install libgdbm-dev
        % sudo apt-get install zlib1g-dev
        % sudo apt-get install libbz2-dev
        % sudo apt-get install liblzma-dev

    CentOSなら:

        % sudo yum groupinstall "Development Tools"
        % sudo yum install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel gdbm-devel xz-devel

    Macなら:

        % brew install gdbm openssl readline sqlite xz

2.  ダウンロードして展開:

        % wget -O- http://www.python.org/ftp/python/3.4.1/Python-3.4.1.tar.xz | tar xJ

3.  configure してビルド:

        % cd Python-3.4.1/
        % ./configure --help
        % ./configure --with-threads
        % make

    {{%div class="note"%}}
ユニコードにはバイト幅の異なる UCS-4 と UCS-2 という2種類があり、
Python 2の configure のデフォルトは UCS-2。
`sys.maxunicode` で確認できる。
Python 3.3以降ではUCS-4のみ。
Python 2をucs4でビルドするには
`./configure --with-threads --enable-unicode=ucs4`
    {{%/div%}}

    {{%div class="note"%}}
モジュールをビルドするのに必要なヘッダファイルが見つからなかったとかで
警告メッセージが表示されるが、だいたい問題ない。
使いそうなモジュールが含まれている場合は、
必要なヘッダファイルを持ってそうなパッケージ (`libXXX-dev` のようなもの) を
Synaptic などからインストールして `make` し直すとよい。
Homebrew で入れたライブラリを利用する場合はオプション付きで `configure`
(特に readline や sqlite は keg-only なので注意):

```sh
./configure --with-threads --prefix=${HOME}/.virtualenv/python CPPFLAGS="-I$(brew --prefix)/include -I$(brew --prefix)/opt/readline/include -I$(brew --prefix)/opt/sqlite/include -I$(brew --prefix)/opt/openssl/include" LDFLAGS="-L$(brew --prefix)/lib -L$(brew --prefix)/opt/readline/lib -L$(brew --prefix)/opt/sqlite/lib -L$(brew --prefix)/opt/openssl/lib"
```
    {{%/div%}}

4.  インストール
    (古いバージョンに上書きせず共存させるため `altinstall`):

        % sudo make altinstall

## 環境設定

See [pip]({{< relref "pip.md" >}})

### `PYTHONPATH`

自分で書いたプログラムをいつでも `import` できるようにする。

-   ファイルはまとめて `$HOME/local/lib/python/` 以下に置く
-   環境変数 `PYTHONPATH` にそのディレクトリを指定する。
    例えば `.zshenv` に以下のように記述するとか:

        export PYTHONPATH=$HOME/local/lib/python

### `PYTHONSTARTUP`

インタラクティブモードで起動するときに読み込むファイルを指定する環境変数。
例えば `.zshrc` に:

    export PYTHONSTARTUP=$HOME/local/lib/python/pythonstartup.py

以下のようなものを書いておくと、`tab` とか `^i` で補完できるようになる。:

    try:
        import readline
    except ImportError:
        print("Module readline not available.")
    else:
        import rlcompleter
        if 'libedit' in readline.__doc__:
            readline.parse_and_bind("bind ^I rl_complete")
        else:
            readline.parse_and_bind("tab: complete")
