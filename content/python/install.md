+++
title = 'Pythonインストール'
tags = ["python"]
[menu.main]
  parent = "python"
  weight = -99
+++

## 標準・公式

MacやLinuxならシステムの一部として
`/usr/bin/python` が既にインストールされているので、
何もしなくても使える。
違うバージョンを使いたければ
[python.org公式のインストーラ](https://www.python.org/downloads/)
で入れるのも悪くない。


## Anaconda

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
全部入りに抵抗がある場合はMinicondaで小さくスタートすることも可能。

例えば Homebrew + pyenv + Miniconda で基本的な環境を整える手順:

```sh
% brew install pyenv
% exec ${SHELL} -l
% pyenv install -l | less
% pyenv install miniconda3-latest
% pyenv global minoconda3-latest
% exec ${SHELL} -l
% python --version
% conda install seaborn biopython flake8
```

その後のパッケージ管理も `conda` で。


## Source

万がいちソースコードからビルドしたい場合の手順

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

        % wget -O- https://www.python.org/ftp/python/3.5.3/Python-3.5.3.tar.xz | tar xJ

3.  configure してビルド:
    ```sh
    % cd Python-3.5.3/
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


## 関連書籍

<a href="https://www.amazon.co.jp/dp/479738946X/ref=as_li_ss_il?ie=UTF8&qid=1485612008&sr=8-6&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=5ea5e48ecc83b9439f21406b6f57c062" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=479738946X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=479738946X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/IPython%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9%E3%82%AF%E3%83%83%E3%82%AF%E3%83%96%E3%83%83%E3%82%AF-%E5%AF%BE%E8%A9%B1%E5%9E%8B%E3%82%B3%E3%83%B3%E3%83%94%E3%83%A5%E3%83%BC%E3%83%86%E3%82%A3%E3%83%B3%E3%82%B0%E3%81%A8%E5%8F%AF%E8%A6%96%E5%8C%96%E3%81%AE%E3%81%9F%E3%82%81%E3%81%AE%E3%83%AC%E3%82%B7%E3%83%94%E9%9B%86-Cyrille-Rossant/dp/4873117488/ref=as_li_ss_il?_encoding=UTF8&psc=1&refRID=X16VFSS3W75RMTG7VGCH&linkCode=li3&tag=heavywatal-22&linkId=b79e2290571289b02621392257a4ac1c" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
