+++
title = 'Installation'
tags = ["python"]
[menu.main]
  parent = "python"
  weight = -99
+++

<https://www.python.org/downloads/>

## Installer

MacやLinuxならシステムの一部として
`/usr/bin/python` が既にインストールされているので基本的にはそれを使えばよい。
Scientificな用途で使いたい場合は
[Numpy/Scipy]({{< relref "scipy.md" >}})
などの主要パッケージもまとめて面倒みてくれる
[Anaconda](https://docs.continuum.io/anaconda/)
で最新版を入れると良い。
GUIもあってかなり親切だが、
`PATH`上でシステムコマンド(`curl`など)を上書きしちゃうヤンチャな面もあるので、
それが気になる人は
[pyenv](https://github.com/yyuu/pyenv)
越しに入れることで汚染を防げる。


## Source

1.  必要なパッケージをインストールしておく:

    Ubuntuなら
    ```sh
    % sudo apt-get install build-essential libreadline6-dev libsqlite3-dev libgdbm-dev zlib1g-dev libbz2-dev liblzma-dev
    ```

    CentOSなら
    ```sh
    % sudo yum groupinstall "Development Tools"
    % sudo yum install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel gdbm-devel xz-devel
    ```

    Macなら
    ```sh
    % brew install gdbm libressl readline sqlite xz
    ```

2.  ダウンロードして展開:

        % wget -O- https://www.python.org/ftp/python/3.5.1/Python-3.5.1.tar.xz | tar xJ

3.  configure してビルド:
    ```sh
    % cd Python-3.5.1/
    % ./configure --help
    % ./configure --prefix=${HOME}/.virtualenv/Python
    % make
    ```

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
パッケージマネージャからインストールして `make` し直すとよい。

Macの場合は `--enable-framework`
を付けてビルドしておかないと使えないモジュールが出てくるので注意。
Homebrewで入れたライブラリを利用する場合は明示的に位置指定が必要。
(特に readline, sqlite, openssl/libressl は keg-only なので注意):

```sh
DST=${HOME}/.virtualenv
./configure --enable-framework=${DST} --prefix=${DST} CPPFLAGS="-I$(brew --prefix)/include -I$(brew --prefix)/opt/readline/include -I$(brew --prefix)/opt/sqlite/include -I$(brew --prefix)/opt/libressl/include" LDFLAGS="-L$(brew --prefix)/lib -L$(brew --prefix)/opt/readline/lib -L$(brew --prefix)/opt/sqlite/lib -L$(brew --prefix)/opt/libressl/lib"
```
    {{%/div%}}

4.  インストール
    (古いバージョンに上書きせず共存させるため `altinstall`):

        % make altinstall


## 環境設定

### パッケージ管理

See [pip]({{< relref "pip.md" >}})

### `PYTHONPATH`

自分で書いたプログラムをいつでも `import` できるようにする。

-   ファイルはまとめて `~/local/lib/python/` 以下に置く
-   環境変数 `PYTHONPATH` にそのディレクトリを指定する。
    例えば `.zshenv` に以下のように記述するとか:

        export PYTHONPATH=${HOME}/local/lib/python

### `PYTHONSTARTUP`

インタラクティブモードで起動するときに読み込むファイルを指定する環境変数。
例えば `.zshrc` に:

    export PYTHONSTARTUP=${HOME}/local/lib/python/pythonstartup.py

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

対話モードをさらに便利にするには [IPython]({{< relref "ipython.md" >}}) を使う。
