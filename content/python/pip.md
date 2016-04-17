+++
title = 'pip'
subtitle = "パッケージ管理"
tags = ["python"]
[menu.main]
  parent = "python"
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

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4797371595/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51uQpDtF%2BdL._SX160_.jpg" alt="みんなのPython 第3版" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4798032948/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41bZDMzzeTL._SX160_.jpg" alt="Pythonプロフェッショナルプログラミング" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4048686291/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51wSTTIQtgL._SX160_.jpg" alt="エキスパートPythonプログラミング" /></a>
