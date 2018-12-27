+++
title = 'pip'
subtitle = "パッケージ管理"
tags = ["python", "package"]
[menu.main]
  parent = "python"
  weight = -95
+++

古いPythonではパッケージ管理のためにツールを別途インストールする必要があった。
Python 3.4 以降では `venv` と `ensurepip` が標準ライブラリに入って少しマシに。

科学技術系の利用だけなら、Python本体のインストールからパッケージ管理までぜーんぶ
[Anaconda]({{< relref "install.md#anaconda" >}}) に任せるのが楽ちん。
その場合、本記事はほぼ無用。


## [`pip`](https://pip.pypa.io/)

[PyPI](http://pypi.python.org)
からの簡単にパッケージをインストールできるようにするツール。
アンインストール機能の無い `easy_install` に取って代わり、
現在では公式に推奨されている。
Python 2.7.9以降、3.4以降では標準ライブラリの
[`ensurepip`](https://docs.python.org/3/library/ensurepip.html)
によって自動的にインストールされる。

-   手動インストール:
    ```sh
    % python -m ensurepip --user
    または
    % curl -O https://bootstrap.pypa.io/get-pip.py
    % python get-pip.py --user
    ```
    `--user` を付けた場合のインストール先は環境変数
    [`PYTHONUSERBASE`]({{< relref "install.md#pythonuserbase" >}})
    で指定できる。
    その中の `bin/` を `PATH` に追加するか、
    絶対パスで `pip` コマンドを使う。

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
        % pip search psutil
        % pip install --upgrade setuptools

-   設定ファイルは `~/.pip/pip.conf` :

        [global]
        download-cache = ~/.pip/cache

-   環境内のパッケージを引き継ぐには:

        & pip freeze --local >requirements.txt
        % pip install -r requirements.txt


## [`venv`](https://docs.python.org/3/library/venv.html) / `virtualenv`

Python実行環境を仮想化するパッケージ。
これで作った仮想環境内で `pip` を使ってパッケージ管理する。
Python 3.3 以降では `venv` が標準ライブラリ入りしたので
[`virtualenv`](https://virtualenv.pypa.io/)
の個別インストールは不要になった。

-   仮想環境を作る:

        % python3.6 -m venv [OPTIONS] /path/to/projectx
        % virtualenv        [OPTIONS] /path/to/projectx

-   仮想環境に入る、から出る:

        % source /path/to/projectx/bin/activate
        (projectx) % deactivate

-   Pythonを更新したりするとシムリンクが壊れるので一旦消して張り直す:

        % find /path/to/projectx/ -type l
        % find /path/to/projectx/ -type l -delete
        % find /path/to/projectx/ -type l
        % virtualenv /path/to/projectx


## [`setuptools`](https://github.com/pypa/setuptools)

パッケージ管理・作成の基本となるライブラリ。
コマンドラインツール `easy_install`
はこれの一部として含まれているが、直接使うことはない。
(`pip` を使う)

See also ["setuptools --- Pythonパッケージ作成"]({{< relref "setuptools.md" >}})

`setuptools` の改良版としてしばらく `distribute` も利用されていたが、
その成果が `setuptools` にマージされたので忘れていい。


## 関連書籍

<a href="https://www.amazon.co.jp/dp/479738946X/ref=as_li_ss_il?ie=UTF8&qid=1485612008&sr=8-6&keywords=python&linkCode=li3&tag=heavywatal-22&linkId=5ea5e48ecc83b9439f21406b6f57c062" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=479738946X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=479738946X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/487311845X/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=72a416f5d10a9e84aaab4b3ee9613329&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311845X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=487311845X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873118417/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=6b1a04ec880b6c730bd6e80273e30e9c&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873118417&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873118417" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873117488/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=2181a50362009e68f507d44fc38716b4&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117488&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873117488" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
