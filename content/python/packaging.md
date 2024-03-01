+++
title = 'Pythonパッケージ作成'
tags = ["python", "package"]
[menu.main]
  parent = "python"
+++

- <https://docs.python.org/3/tutorial/modules.html>
- <https://docs.python.org/3/reference/import.html>
- <https://packaging.python.org/>

## ファイル構成

GitHubやローカルの開発環境から `pip` で直接インストールできる形。

```sh
pywtl/
├── LICENSE
├── README.md
├── pyproject.toml
├── src/wtl/
│   ├── __init__.py
│   └── hello.py
└── tests/
```

リポジトリ名(`pywtl`)とパッケージ名(`wtl`)は必ずしも一致してなくてもよい。

ソースコードは `src` の中に入れる流派と、ルート直下に置く流派がある。
[落とし穴が少なくて推奨されているのは前者](https://docs.pytest.org/explanation/goodpractices.html)。

開発向けの `-e,--editable` オプションをつけたローカルインストールではコピーが起こらず、
編集後に再インストールしなくてもそのまま反映される。

```sh
pip3 install -v --user -e ~/git/pywtl/
python3 -m wtl.hello
python3 -m site
```

### `pyproject.toml`

パッケージ作成に関わる全てのメタ情報を書いておくファイル。
`setuptools` に依存しない形式として
[PEP 517](https://www.python.org/dev/peps/pep-0517),
[PEP 621](https://www.python.org/dev/peps/pep-0621)
で決められた。
過去によく使われていた `setup.py`,
[`setup.cfg`](https://setuptools.pypa.io/en/latest/userguide/declarative_config.html),
[`MANIFEST.in`](https://setuptools.pypa.io/en/latest/userguide/miscellaneous.html)
などは非推奨になった。

[PyPA/Flit](https://flit.readthedocs.io/) (setuptools後継？),
[PDM](https://pdm.fming.dev/),
[Poetry](https://python-poetry.org/),
など後発のツールは早くから対応していて、
`setuptools` も[ようやく61.0から使えるようになった](https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html)。
`[project]` テーブルは [PEP 621](https://www.python.org/dev/peps/pep-0621)
で項目が決められているためツールによらず共通。
それ以外の `[build-system]` などは使うツールによって異なる。

```toml
[build-system]
requires = ["flit_core >=3.6,<4"]
build-backend = "flit_core.buildapi"

[project]
name = "wtl"
authors = [
  {name = "Watal M. Iwasaki", email = "heavywatal@gmail.com"}
]
license = {file = "LICENSE"}
readme = "README.md"
dynamic = ["description", "version"]
classifiers = [
  "Development Status :: 2 - Pre-Alpha",
  "Environment :: Console",
  "Intended Audience :: Science/Research",
  "License :: OSI Approved :: MIT License",
  "Topic :: Scientific/Engineering :: Bio-Informatics",
]
requires-python = ">=3.12"
dependencies = [
  "tomli-w",
]

[project.optional-dependencies]
dev = [
  "pytest",
  "pytest-cov",
  "ruff",
]

[project.urls]
Source = "https://github.com/heavywatal/pywtl"

[project.scripts]
"hello.py" = "wtl.hello:main"

[tool.pyright]
typeCheckingMode = "strict"

[tool.ruff.lint]
select = ["F", "E", "W", "I", "N", "UP", "S", "B", "A", "PL"]
ignore = ["S101"]

[tool.pytest.ini_options]
pythonpath = ["src"]
testpaths = ["tests"]

[tool.coverage.run]
source = ["src"]

[tool.coverage.report]
exclude_also = [
  "if __name__ == .__main__.:",
]
```

`dynamic` に指定したものは `__init__.py` に `__version__ = "0.1.0"`
などと書いてあるものを参照できる。
比較したいときは
[`packaging.version.parse()`](https://packaging.pypa.io/en/latest/version.html)
を利用する。

`dependencies` に列挙された依存パッケージは
`pip3 install` で自動的にインストールされる。
一方 `requirements.txt` はインストール過程には関与せず、
能動的に `pip3 install -r requirements.txt`
を打たなきゃインストールされない。

optional dependencies もインストールしたい場合は
`pip3 install -v -e .[dev]` のように `[key]` を使ってパッケージを指定する。

`project.scripts` で設定したものは
`${prefix}/bin/` に実行可能ファイルが配置される。
以前は `console_scripts` で設定していた。

コード整形やテストのような各種開発ツールの設定も `[tool.***]` に記述できる。
linterとしては
[`pyproject.toml` 対応拒否のflake8](https://github.com/PyCQA/flake8/issues/234)
を捨てて超高速Rust製[ruff](https://docs.astral.sh/ruff/)を使う。
0.2からはformatterとしても使えるようになり、
[https://black.readthedocs.io](black)も不要になった。


### ソースコード

`wtl/__init__.py`
: このディレクトリがひとつのパッケージであることを示すファイル。
  空でもいいし、初期化処理やオブジェクトを記述してもよい。
  文字列変数 `__version__ = "0.1.2"` を定義して
  `wtl.__version__` のように参照できるようにしておくのが慣例。

`wtl/hello.py`
```py
"""Simple module to say hello
"""
import getpass


def main():
    print("Hello, " + getpass.getuser() + "!")


if __name__ == "__main__":
    main()
```

ソースツリーの中にあるファイルを参照するには
[`importlib.resources`](https://docs.python.org/library/importlib.html#module-importlib.resources)
が使える。
Pythonスクリプトではない設定ファイルなどを同梱して読み込むのに便利。
