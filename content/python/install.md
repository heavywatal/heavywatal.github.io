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
何もしなくても使えっちゃ使える。
でも大概そこに入ってるのは古い2.7とかなので、
ちゃんと使える3.xを使いたければ
[python.jp/環境構築ガイド](https://www.python.jp/install/install.html)
に従って最新版を入れるのが簡単。


## pyenv

管理者権限なしでホーム以下にインストールするには
[pyenv](https://github.com/pyenv/pyenv)
が便利。

1.  [Homebrew]({{< relref "homebrew.md" >}}) か
    [Git]({{< relref "git.md" >}}) を使ってpyenvをインストール:

    ```sh
    brew install pyenv
    # or
    git clone https://github.com/pyenv/pyenv.git ~/.pyenv
    mkdir -p ~/.pyenv/cache
    ```

1.  Pythonのインストール先を決める環境変数
    `PYENV_ROOT` を公式推奨の `~/.pyenv` に設定し、
    ついでにPATHも追加しておく。
    シェルの設定ファイル (e.g., `~/.bashrc`) に次のように追記:

    ```sh
    export PYENV_ROOT=${HOME}/.pyenv
    PATH=${PYENV_ROOT}/bin:$PATH
    PATH=${PYENV_ROOT}/versions/$(pyenv global)/bin:$PATH
    export PATH
    ```

    `pyenv shell` や `pyenv local`
    を使ってPythonのバージョンを頻繁に切り替える場合は、
    [公式の説明](https://github.com/pyenv/pyenv#installation)どおりに
    `eval "$(pyenv init -)"`
    を設定してshimsを使う方法のほうがいいかもしれないけど、
    そうでなければこのようにPATHだけ設定するほうが単純で、
    起動時間も短くなる。

1.  シェルを再起動して設定を反映し、
    目当てのバージョンを探してインストール:

    ```sh
    exec $SHELL -l
    pyenv install -l | less
    pyenv install 3.7.3
    ```

1.  インストールしたものを常に使うように設定:

    ```sh
    pyenv global 3.7.3
    exec $SHELL -l
    ```

1.  [pip]({{< relref "pip.md" >}}) のパスを確認し、パッケージを入れる:

    ```sh
    which pip3
    pip3 install -U setuptools pip wheel
    pip3 install -r /path/to/requirements.txt
    ```

    よく使うパッケージは
    [`requirements.txt`](https://github.com/heavywatal/dotfiles/blob/master/.config/python/requirements.txt)
    の形でまとめておくと楽。


### [既知の問題](https://github.com/pyenv/pyenv/wiki/Common-build-problems)

-   3.1.0より古い[matplotlib]({{< relref "matplotlib.md" >}}) で
    `macosx` backend を使いたい場合などは環境変数
    `PYTHON_CONFIGURE_OPTS="--enable-framework"`
    をセットしてFramework型でビルドする。

-   [Mojaveで "zlib not available" と怒られる問題](https://github.com/pyenv/pyenv/issues/1219)は
    `CFLAGS="-I$(xcrun --show-sdk-path)/usr/include"` を定義して回避。

-   [PEP 394](https://www.python.org/dev/peps/pep-0394/)
    になるべく沿うように、
    `python` では `/usr/bin/python` が呼び出される状態を維持しつつ、
    `python3`, `pip3` を明示的に使うようにしたい。
    が、いまのところできない？


## Anaconda

Scientificな用途で使いたい場合は
[Numpy/Scipy]({{< relref "scipy.md" >}})
などの主要パッケージもまとめて面倒みてくれる
[Anaconda](https://docs.continuum.io/anaconda/)
を使うという選択肢もある。
GUIのインストーラでもいいし、Homebrewでもいける:

```sh
brew install anaconda
export PATH=/usr/local/anaconda3/bin:"$PATH"
```

ただし`PATH`上でシステムコマンド(`curl`など)を上書きしちゃうヤンチャな面もあるので、
それが気になる人はpyenv越しに入れることで汚染を防げる。
全部入りに抵抗がある場合は
`pyenv install miniconda3-latest`
から小さくスタートすることも可能。
パッケージ管理では `pip` の代わりに非公式の `conda` を使うことになる。


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
