+++
title = 'Pythonパッケージ作成'
tags = ["python", "package"]
aliases = ["setuptools.html"]
[menu.main]
  parent = "python"
+++

- <https://docs.python.org/3/tutorial/modules.html>
- <https://docs.python.org/3/reference/import.html>
- <https://packaging.python.org/>
- <https://setuptools.pypa.io/>

## ファイル構成

GitHubやローカルの開発環境から `pip` で直接インストールできる形。

```sh
pywtl/
├── LICENSE
├── README.md
├── pyproject.toml
└── wtl/
    ├── __init__.py
    └── hello.py
```

リポジトリ名(`pywtl`)とパッケージ名(`wtl`)は必ずしも一致してなくてもよい。
開発向けの `-e,--editable` オプションをつけたローカルインストールではコピーが起こらず、
編集後に再インストールしなくてもそのまま反映される。

```sh
pip install -v --user -e ~/git/pywtl/
python -m wtl.hello
```

### `pyproject.toml`

`setuptools` に依存しない形式として
[PEP 517](https://www.python.org/dev/peps/pep-0517),
[PEP 621](https://www.python.org/dev/peps/pep-0621)
で決められた。
[PyPA/Flit](https://flit.readthedocs.io/) (setuptools後継？),
[PDM](https://pdm.fming.dev/),
[Poetry](https://python-poetry.org/),
など後発のツールは大概このファイルだけ書けば済む。

`setuptools` は未対応なので書けるのはこれだけ:
```toml
[build-system]
requires = ["setuptools", "wheel"]
build-backend = "setuptools.build_meta"
```
残りの項目は `setup.cfg` へ。


### `setup.py`

`setup.cfg` の設定を読んでくれるようになってからは
`setup()` の引数にいろいろ渡す必要は無くなった。

```py
from setuptools import setup
setup()
```

さらに新しいバージョンではこれすらも不要になった。
ただしこれが無いと `pip install --editable` がサポートされないらしい。


### `setup.cfg`

<https://setuptools.pypa.io/en/latest/userguide/declarative_config.html>

```ini
[metadata]
name = wtl
version = attr: wtl.__version__
url = https://github.com/heavywatal/pywtl
author = Watal M. Iwasaki
author_email = heavywatal@gmail.com
license_files = LICENSE
description = wtl: Personal Python package
long_description = file: README.md
long_description_content_type = text/markdown

[options]
install_requires =
  psutil
  requests
python_requires = >=3.9
entry_points = file: entry_points.cfg
packages = find:
```

`version` は `__init__.py` に `__version__ = "0.1.0"`
などと書いてあるものを参照できる。
比較したいときは
[`packaging.version.parse()`](https://packaging.pypa.io/en/latest/version.html)
を利用する。

`install_requires` に列挙された依存パッケージは
`pip install` で自動的にインストールされる。
一方 `requirements.txt` はsetuptoolsではなくpipの機能で、
能動的に `pip install -r requirements.txt`
を打たなきゃインストールされないし、
そこに列挙されたパッケージが本当に必要になるまで怒られない。


### `entry_points`

- <https://setuptools.pypa.io/en/latest/userguide/quickstart.html#automatic-package-discovery>
- <https://setuptools.pypa.io/en/latest/userguide/entry_point.html>
- <https://setuptools.pypa.io/en/latest/pkg_resources.html#entry-points>

用途はいろいろあるけど
`${prefix}/bin/` に実行可能ファイルを配置するのによく使われる。
設定は下記のように別ファイル `entry_points.cfg` として
`setup.cfg` から読ませるのが楽チン。

```ini
[console_scripts]
hello.py = wtl.hello:main
```

引数を取らない関数のみ利用可能。
コマンドライン引数を受け取りたい場合はその関数の中で標準の
[`argparse`](https://docs.python.org/3/library/argparse.html)
を使って処理する。


### ソースコード

`wtl/__init__.py`
: このディレクトリがひとつのパッケージであることを示すファイル。
  空でもいいし、初期化処理やオブジェクトを記述してもよい。
  文字列変数 `__version__ = "0.1.2"` を定義して
  `wtl.__verison__` のように参照できるようにしておくのが慣例。

`wtl/hello.py`
```py
"""Sample module
"""


def main():
    import getpass
    print("Hello, " + getpass.getuser() + "!")


if __name__ == "__main__":
    main()
```


### `MANIFEST.in`

Pythonモジュール以外で配布物に含めたいもの、除外したいものがあれば
`include`, `recursive-include`, `exclude`
で指定する。


## コマンド

以前は
[`setup.py`](https://setuptools.pypa.io/en/latest/userguide/commands.html)
を使って操作してたけどそれは廃れつつある。
最小限のビルド機能は
[`build`](https://pypa-build.readthedocs.io/) パッケージに引き継がれた。

```sh
pip install build
python3 -m build --help
python3 -m build
```
