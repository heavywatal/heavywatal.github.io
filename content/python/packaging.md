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
uv pip install -e .
python3 -m wtl.hello
python3 -m site
```

手動で作ってもいいけど [`uv`]({{< relref "install.md#uv" >}}) に任せるのが楽。
いくつかの形式があるけどとりあえず `--lib`:
- `--app`: パッケージとして扱われることを想定しないスクリプトやウェブアプリなど。
- `--package`: `src/` レイアウトで `[build-system]` も設定される。
- `--lib`: 上記に加えて `py.typed` も作成される。
```sh
uv init --lib example-lib
cd example-lib/
uv run python -c 'import example_lib; print(example_lib.hello())'
uv version
```
See <https://docs.astral.sh/uv/concepts/projects/init/>.


### `pyproject.toml`

- <https://packaging.python.org/en/latest/specifications/pyproject-toml/>
- <https://packaging.python.org/en/latest/guides/writing-pyproject-toml/>

パッケージ作成に関わる全てのメタ情報を書いておくファイル。
`setuptools` に依存しない形式として
[PEP 517](https://www.python.org/dev/peps/pep-0517),
[PEP 621](https://www.python.org/dev/peps/pep-0621)
で決められた。
過去によく使われていた `setup.py`,
[`setup.cfg`](https://setuptools.pypa.io/en/latest/userguide/declarative_config.html),
[`MANIFEST.in`](https://setuptools.pypa.io/en/latest/userguide/miscellaneous.html)
などは非推奨になった。

`[build-system]`, `[project]`, `[tool]`
という3つのテーブルから成る。後に
[PEP 735](https://www.python.org/dev/peps/pep-0735) で
`[dependency-groups]` が追加された。


#### `build-system`

必須ではないけど推奨。
[PyPA/Flit](https://flit.readthedocs.io/) (setuptools後継？),
[PDM](https://pdm.fming.dev/),
[Poetry](https://python-poetry.org/),
など後発のツールは早くから対応していて、
`setuptools` も[ようやく61.0から使えるようになった](https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html)。
とりあえず `uv init --lib` 初期設定の `uv_build` を使い、もし不満を感じたら考える。
```toml
[build-system]
requires = ["uv_build>=0.10.0,<1.0.0"]
build-backend = "uv_build"
```

#### `project`

```toml
[project]
name = "wtl"
version = "0.1.0"
description = "Personal Python Package"
authors = [{ name = "Watal M. Iwasaki", email = "heavywatal@gmail.com" }]
license = "MIT"
license-files = ["LICENSE"]
readme = "README.md"
classifiers = [
  "Development Status :: 2 - Pre-Alpha",
  "Environment :: Console",
  "Intended Audience :: Science/Research",
  "License :: OSI Approved :: MIT License",
  "Topic :: Scientific/Engineering :: Bio-Informatics",
]
requires-python = ">=3.14"
dependencies = [
  "tomli-w",
]

[project.urls]
source = "https://github.com/heavywatal/pywtl"

[project.scripts]
"hello.py" = "wtl.hello:main"

[dependency-groups]
dev = [
  "pytest",
  "pytest-cov",
  "ruff",
]
```

`project.dynamic` に `["description", "version"]` と指定して
`__init__.py` のdocstringや `__version__`
を参照できるかはbackend次第。
`uv_build` は今のところサポートしていないので
([uv#8714](https://github.com/astral-sh/uv/issues/8714))
`pyproject.toml` に書いたバージョンを
[`importlib.metadata`](https://docs.python.org/3/library/importlib.metadata.html)
で `__init__.py` に取り込む。
```py
import importlib.metadata

assert __package__
__version__ = importlib.metadata.version(__package__)
__doc__ = importlib.metadata.metadata(__package__)["Summary"]
```
バージョンを比較したいときは
[`packaging.version.parse()`](https://packaging.pypa.io/en/latest/version.html)
を利用する。

依存関係を書けるところはいくつかある。
See <https://docs.astral.sh/uv/concepts/projects/dependencies/>:
- `project.dependencies`:
  普通の依存関係。
  `uv pip install` で自動的にインストールされる。
- `project.optional-dependencies`:
  通称"extras"。
  ユーザー向けに公開されるけどデフォルトではインストールされない。
  `uv add altair --optional plot` のように追加し、
  `uv pip install polars[plot]` のようにインストールする。
- [`dependency-groups`](https://packaging.python.org/en/latest/specifications/dependency-groups/):
  開発者向けで `[project]` の外にある。
  パッケージ化しないプロジェクトの依存関係を記述するのにも使える。
  `uv add --group dev ruff` のようにして追加。
  名前は何でもいいけど `dev` は `uv` で特別扱いされていて、
  `--dev` オプションがあったり、デフォルトで `uv sync` 対象だったりする。
- `requirements.txt`:
  インストール過程には関与せず、能動的に
  `pip install -r requirements.txt` のように参照するためのもの。

`project.scripts` で設定したものは
`${prefix}/bin/` に実行可能ファイルが配置される。
以前は `console_scripts` で設定していた。


#### `tool`

コード整形やテストのような各種開発ツールの設定を記述する。

```toml
[tool.pyright]
typeCheckingMode = "strict"

[tool.ruff.lint]
select = ["ALL"]
ignore = [
  "D1",   # missing docstring
  "D203", # incompatible
  "D213", # incompatible
  "ANN401", # Any
  "T201", # print
  "S101", # assert
  "DTZ",  # timezone
  "COM812", # trailing comma
  "TD",   # todo
  "FIX",  # todo
]

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

linterとしては
[`pyproject.toml` 対応拒否のflake8](https://github.com/PyCQA/flake8/issues/234)
を捨てて超高速Rust製[ruff](https://docs.astral.sh/ruff/)を使う。
0.2からはformatterとしても使えるようになり、
[black](https://black.readthedocs.io)も不要になった。
各プロジェクトの `pyproject.toml` で設定するのが基本だが、次のようなルールで読み込まれる。
See <https://docs.astral.sh/ruff/configuration/>.
- 明示的にCLIで指定するのが最優先。
- 最初に見つかったファイルを読み込んで探索終了。
  マージしたい場合は `extends` で明示的に指定する。
- ファイル名は `.ruff.toml`, `ruff.toml`, `pyproject.toml` で、この順に優先。
- カレントディレクトリから上に向かってファイルを探索。
- ユーザー設定を探索: `${XDG_CONFIG_HOME}/ruff/`, `${HOME}/.config/ruff/`, `${HOME}/.ruff/`, etc.


### ソースコード

`wtl/__init__.py`
: このディレクトリがひとつのパッケージであることを示すファイル。
  空でもいいし、初期化処理やオブジェクトを記述してもよい。
  文字列変数 `__version__ = "0.1.2"` を定義して
  `wtl.__version__` のように参照できるようにしておくのが慣例。

`wtl/hello.py`
```py
"""Simple module to say hello."""

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
