+++
title = 'pip'
subtitle = "パッケージ管理"
tags = ["python"]
[menu.main]
  parent = "python"
  weight = -95
+++

Python 2 (3.2以下) ではパッケージ管理のために
外部ツールを別途インストールする必要がある(結構多くて複雑...)。
Python 3.4 以降では `venv` と `ensurepip`
が標準ライブラリに入るため楽チン。

## Installation (Python 2)

1.  ベースのPythonインタープリタをひとつ決め
    (例えばOSに標準装備されてる `/usr/bin/python`)、
    それを使って `setuptools` をインストール:

        % curl -O https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
        % /usr/bin/python ez_setup.py --user

2.  続いて `pip` をインストール:

        % curl -O https://bootstrap.pypa.io/get-pip.py
        % /usr/bin/python get-pip.py --user

3.  Macなら `~/Library/Python/2.7/` 以下に、
    Linuxなら `~/.local/` 以下に上記のライブラリがインストールされるので、
    その中にある `bin` を `PATH` に追加する。
4.  `pip` を使って `virtualenv` をインストール:

        % pip install --user virtualenv

## Usage

### `venv` / `virtualenv`

Doc
:   <http://www.virtualenv.org/>

PyPI
:   <http://pypi.python.org/pypi/virtualenv>

Python3
:   <http://docs.python.org/3/library/venv>

Python実行環境を仮想化するパッケージ。
これで作った仮想環境内で `pip` を使ってパッケージ管理する。
Python 3.3 以降では `venv` が標準ライブラリ入りしたので
`virtualenv` は不要になった。

-   仮想環境を作る:

        % python3.4 -m venv [OPTIONS] /path/to/projectx
        % pyvenv3.4         [OPTIONS] /path/to/projectx
        % virtualenv        [OPTIONS] /path/to/projectx

-   仮想環境に入る、から出る:

        % source /path/to/projectx/bin/activate
        (projectx) % deactivate

-   Pythonを更新したりするとシムリンクが壊れるので一旦消して張り直す:

        % find /path/to/projectx/ -type l
        % find /path/to/projectx/ -type l -delete
        % find /path/to/projectx/ -type l
        % virtualenv /path/to/projectx

### `ensurepip`

Python3
:   <http://docs.python.org/3/library/ensurepip>

`pip` を自動インストールするための公式ライブラリ。
手動でも使える:

    % python3.4 -m ensurepip -h
    % python3.4 -m ensurepip --upgrade

### `pip`

Doc
:   <https://pip.pypa.io/>

PyPI
:   <http://pypi.python.org/pypi/pip>

[PyPI](http://pypi.python.org) からのパッケージのダウンロード、ビルド、インストールを
簡単にできるようにするツール。
アンインストール機能の無い easy\_install に取って代わり、
現在では公式推奨となっている。
Python 3.4 以降では標準ライブラリの `ensurepip`
によって自動的にインストールされる。
`setuptools` に依存している。

-   `pip` コマンド
    -   `install`: Install packages.
    -   `uninstall`: Uninstall packages.
    -   `freeze`: Output installed packages in requirements format.
    -   `list`: List installed packages.
    -   `show`: Show information about installed packages.
    -   `search`: Search PyPI for packages.
    -   `help`: Show help for commands.
-   全体のヘルプ、コマンド毎の詳細ヘルプ:

        % pip help
        % pip install --help

-   `setup.py` にオプションを渡す:

        % pip install --install-option="--prefix=/usr/local" numpy

-   よく使うコマンド:

        % pip list && pip list --outdated
        % pip search numpy
        % pip install --upgrade numpy

-   設定ファイルは `~/.pip/pip.conf` :

        [global]
        download-cache = ~/.pip/cache

### `setuptools`

パッケージ管理の基本となるライブラリ。
コマンドラインツール easy\_install
はこれの一部として含まれているが、直接使うことはない。
(`pip` を使う)

Website
:   <http://pythonhosted.org/setuptools/>

PyPI
:   <https://pypi.python.org/pypi/setuptools/>

### `distribute`

`setuptools` の改良版としてしばらく利用されていたが、
その成果が `setuptools` にマージされたのでもう使われない。

Website
:   <http://packages.python.org/distribute>

PyPI
:   <http://pypi.python.org/pypi/distribute>

## 書籍

<a href="https://www.amazon.co.jp/dp/479738946X/ref=as_li_ss_il?ie=UTF8&qid=1485612008&sr=8-6&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=5ea5e48ecc83b9439f21406b6f57c062" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=479738946X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=479738946X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/IPython%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9%E3%82%AF%E3%83%83%E3%82%AF%E3%83%96%E3%83%83%E3%82%AF-%E5%AF%BE%E8%A9%B1%E5%9E%8B%E3%82%B3%E3%83%B3%E3%83%94%E3%83%A5%E3%83%BC%E3%83%86%E3%82%A3%E3%83%B3%E3%82%B0%E3%81%A8%E5%8F%AF%E8%A6%96%E5%8C%96%E3%81%AE%E3%81%9F%E3%82%81%E3%81%AE%E3%83%AC%E3%82%B7%E3%83%94%E9%9B%86-Cyrille-Rossant/dp/4873117488/ref=as_li_ss_il?_encoding=UTF8&psc=1&refRID=X16VFSS3W75RMTG7VGCH&linkCode=li3&tag=heavywatal-22&linkId=b79e2290571289b02621392257a4ac1c" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
