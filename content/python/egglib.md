+++
title = 'EggLib'
subtitle = "Tools for Evolutionary Genetics and Genomics"
tags = ["python"]
[menu.main]
  parent = "python"
+++

<http://egglib.sourceforge.net/>

前は `seqlib` という名前だったけどいつの間にか生まれ変わった。

## Usage

FASTA形式の配列セットを読み込んで多型に関する統計量を計算

```python
import egglib
with open("sequences.fasta", 'rU') as fin:
    align = egglib.Align(string=fin.read())
stats = align.polymorphism()
print(stats)
```

## Installation

<http://egglib.sourceforge.net/installation.html>

1.  C++で書かれたライブラリをインストール:

        % wget -O- http://sourceforge.net/projects/egglib/files/current/egglib-cpp-2.1.7.tar.gz | tar xz
        % cd egglib-cpp-2.1.7/
        % ./configure --help
        % make
        % sudo make install

2.  Python用のライブラリをインストール:

        % wget -O- http://sourceforge.net/projects/egglib/files/current/egglib-py-2.1.7.tar.gz | tar xz
        % cd egglib-py-2.1.7/
        % python setup.py build
        % python setup.py install
