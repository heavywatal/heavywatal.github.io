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


## pyenv

管理者権限なしでホーム以下にインストールするには
[pyenv](https://github.com/pyenv/pyenv)
が便利。

[matplotlib]({{< relref "matplotlib.md" >}}) で
`macosx` backend を使いたい場合などは環境変数
`PYTHON_CONFIGURE_OPTS="--enable-framework"`
をセットしてFramework型でビルドする。

[Mojaveで "zlib not available" と怒られる問題](https://github.com/pyenv/pyenv/issues/1219)は
`CFLAGS="-I$(xcrun --show-sdk-path)/usr/include"` を定義して回避。

```sh
git clone https://github.com/pyenv/pyenv.git ~/.pyenv
~/.pyenv/bin/pyenv install -l | less
~/.pyenv/bin/pyenv install 3.7.2
~/.pyenv/bin/pyenv global 3.7.2
```

シェルのPATHを設定する。
Pythonのバージョンを頻繁に切り替える場合は、
[公式の説明](https://github.com/pyenv/pyenv#installation)どおりに
`eval "$(pyenv init -)"`
を設定してshimsを使う方法のほうがいいかもしれないけど、
そうでなければ以下のようにPATHだけ設定するほうが単純:

```sh
export PYENV_ROOT=${HOME}/.pyenv
if [ -d $PYENV_ROOT ]; then
  PATH=${PYENV_ROOT}/bin:$PATH
  PATH=${PYENV_ROOT}/versions/$(pyenv global)/bin:$PATH
fi
export PATH
```

シェルを再起動して [pip]({{< relref "pip.md" >}}) でパッケージを入れる:

```sh
exec $SHELL -l
pip install -U setuptools pip wheel
pip install -U flake8 psutil requests
pip install -U seaborn ipython biopython
```


## Anaconda

Scientificな用途で使いたい場合は
[Numpy/Scipy]({{< relref "scipy.md" >}})
などの主要パッケージもまとめて面倒みてくれる
[Anaconda](https://docs.continuum.io/anaconda/)
で最新版を入れるという選択肢もある。
GUIのインストーラでもいいし、Homebrewでもいける:

```sh
brew cask install anaconda
export PATH=/usr/local/anaconda3/bin:"$PATH"
```

ただし`PATH`上でシステムコマンド(`curl`など)を上書きしちゃうヤンチャな面もあるので、
それが気になる人はpyenv越しに入れることで汚染を防げる。
全部入りに抵抗がある場合は
`pyenv install miniconda3-latest`
から小さくスタートすることも可能。
パッケージ管理では `pip` の代わりに非公式の `conda` を使うことになる。


## Source

万がいちソースコードからビルドしたい場合の手順

1.  必要なパッケージをインストールしておく:

    Ubuntuなら
    ```sh
    sudo apt install build-essential libreadline6-dev libsqlite3-dev libgdbm-dev zlib1g-dev libbz2-dev liblzma-dev
    ```

    CentOSなら
    ```sh
    sudo yum groupinstall "Development Tools"
    sudo yum install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel gdbm-devel xz-devel
    ```

    Macなら
    ```sh
    brew install gdbm libressl readline sqlite xz
    ```

2.  ダウンロードして展開:

        wget -O- https://www.python.org/ftp/python/3.5.3/Python-3.5.3.tar.xz | tar xJ

3.  configure してビルド:
    ```sh
    cd Python-3.5.3/
    ./configure --help
    ./configure --prefix=${HOME}/Python
    make -j2
    ```

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
DST=${HOME}/Python
./configure --enable-framework=${DST} --prefix=${DST} CPPFLAGS="-I$(brew --prefix)/include -I$(brew --prefix)/opt/readline/include -I$(brew --prefix)/opt/sqlite/include -I$(brew --prefix)/opt/libressl/include" LDFLAGS="-L$(brew --prefix)/lib -L$(brew --prefix)/opt/readline/lib -L$(brew --prefix)/opt/sqlite/lib -L$(brew --prefix)/opt/libressl/lib"
```
    {{%/div%}}


    {{%div class="note"%}}
ユニコードにはバイト幅の異なる UCS-4 と UCS-2 という2種類があり、
Python 2の configure のデフォルトは UCS-2。
`sys.maxunicode` で確認できる。
Python 3.3以降ではUCS-4のみ。
Python 2をucs4でビルドするには
`./configure --with-threads --enable-unicode=ucs4`
    {{%/div%}}

4.  インストール
    (古いバージョンに上書きせず共存させるため `altinstall`):

        make altinstall


## 環境変数

https://docs.python.org/3/using/cmdline.html#environment-variables

シェルの設定ファイル(`~/.zshrc` とか)で `export` しておく。
一時的に無効したいときは `python -E` で起動。

### `PYTHONPATH`

`import` の探索パス (`sys.path`) の先頭付近に場所を追加できる。
例えば自分で書いたモジュールやパッケージの置き場所を指定しておけば、
いつでも優先的に `import` できるようになる。


### `PYTHONUSERBASE`

[`pip install`]({{< relref "pip.md" >}}) や `setup.py install` における
`--user` オプションの目的地を指定できる。
デフォルトでは `${HOME}/.local` 。
MacのFramework buildでは `${HOME}/Library/Python/2.7` とかになる。

現在の設定は
[`site`](https://docs.python.org/3/library/site.html)
の `USER_BASE` で確認できる (`python -m site`)。
`USER_SITE` はその下の `{BASE}/lib/python*.*/site-packages` に配置され、
`sys.path` に含まれる。
除外したいときは `PYTHONNOUSERSITE` をセットするか `python -s` で起動。


### `PYTHONSTARTUP`

インタラクティブモードで起動するときに読み込むファイルを指定できる。
例えば以下のようなものを書いておくと、
3.4未満の古いPythonでも `tab` とか `^i` で補完できるようになる:

```py
import sys

if sys.version_info < (3, 4):
    import rlcompleter
    import readline
    rlcompleter.__name__  # suppress F401
    if 'libedit' in readline.__doc__:
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")
    del readline, rlcompleter
else:
    del sys
```

対話モードをさらに便利にするには [IPython]({{< relref "ipython.md" >}}) を使う。


## 関連書籍

<a href="https://www.amazon.co.jp/dp/479738946X/ref=as_li_ss_il?ie=UTF8&qid=1485612008&sr=8-6&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=5ea5e48ecc83b9439f21406b6f57c062" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=479738946X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=479738946X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/487311845X/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=72a416f5d10a9e84aaab4b3ee9613329&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311845X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=487311845X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873118417/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=6b1a04ec880b6c730bd6e80273e30e9c&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873118417&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873118417" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873117488/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=2181a50362009e68f507d44fc38716b4&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
