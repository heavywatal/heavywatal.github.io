+++
title = 'QuickLook'
tags = ["mac"]
[menu.main]
  parent = "mac"
+++

Leopardからの新機能。
Finder 上で <kbd>Space</kbd> を押すだけで、
いろんなファイルを開くことなくチョイ見することができる。
`/Library/QuickLook/` または `~/Library/QuickLook/`
以下にプラグインを配置することで、機能を拡張できる。
[Homebrew]({{< relref "homebrew.md" >}}) の Cask で入れるのが楽ちん。

## Plugins

-   [BetterZipQL.qlgenerator](http://macitbetter.com/BetterZip-Quick-Look-Generator/):
    圧縮ファイルを展開せずに一覧表示
-   [Suspicious Package.qlgenerator](http://www.mothersruin.com/software/SuspiciousPackage/):
    pkgからインストールする前にの中身をチェックできるように
-   [QuickLookCSV.qlgenerator](http://code.google.com/p/quicklook-csv/):
    CSVファイルをセルで表示
-   [QLStephen.qlgenerator](http://whomwah.github.com/qlstephen/):
    READMEやMakefileのような拡張子無しのファイルを表示
-   [QLColorCode.qlgenerator](http://code.google.com/p/qlcolorcode/):
    さまざまなソースコードを色分け表示。
-   [ScriptQL.qlgenerator](http://kainjow.com/):
    AppleScript
-   [qlImageSize.qlgenerator](https://github.com/Nyx0uf/qlImageSize):
    画像表示時のタイトルにピクセル数を表示
-   [FigTree](http://tree.bio.ed.ac.uk/software/figtree/):
    NEXUSあるいはNEWICK形式のtreeファイルから系統樹を描く

## Commands

-   キャッシュ削除してプラグイン再読込:

        qlmanage -r cache && qlmanage -r

-   認識されているプラグインを一覧表示:

        qlmanage -m
