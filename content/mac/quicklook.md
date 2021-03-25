+++
title = 'QuickLook'
tags = ["mac"]
[menu.main]
  parent = "mac"
+++

Finder上で <kbd>Space</kbd> を押すだけで、
いろんなファイルを開かずに覗き見ることができる。
`/Library/QuickLook/` または `~/Library/QuickLook/`
以下にプラグインを配置することで、機能を拡張できる。
[Homebrew]({{< relref "homebrew.md" >}}) で入れるのが楽ちん。

## Plugins

-   [BetterZipQL.qlgenerator](https://macitbetter.com/BetterZip-Quick-Look-Generator/):
    圧縮ファイルを展開せずに一覧表示
-   [Suspicious Package.qlgenerator](https://mothersruin.com/software/SuspiciousPackage/):
    pkgからインストールする前にの中身をチェックできるように
-   [QuickLookCSV.qlgenerator](https://github.com/p2/quicklook-csv):
    CSVファイルをセルで表示
-   [QLStephen.qlgenerator](https://whomwah.github.io/qlstephen/):
    READMEやMakefileのような拡張子無しのファイルを表示
-   [QLColorCode.qlgenerator](https://github.com/anthonygelibert/QLColorCode):
    さまざまなソースコードを色分け表示。
-   [FigTree](http://tree.bio.ed.ac.uk/software/figtree/):
    NEXUSあるいはNEWICK形式のtreeファイルから系統樹を描く

## Commands

```sh
brew install qlcommonmark qlcolorcode qlstephen quicklook-csv webpquicklook betterzip suspicious-package

# コマンドラインから利用
qlmanage -p some.csv

# キャッシュ削除してプラグイン再読込:
qlmanage -r cache && qlmanage -r

# 認識されているプラグインを一覧表示
qlmanage -m
```
