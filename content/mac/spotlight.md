+++
title = 'Spotlight'
tags = ["mac"]
[menu.main]
  parent = "mac"
+++

## 除外設定

余計なインデックスはメモリとCPUの無駄遣いになるので、
検索が必要ないフォルダは積極的に指定しとく。
(e.g. `~/Music`, `~/Pictures`, `~/Movies`)

-   `System Preferences... --> Spotlight --> Privacy --> "+"`
    or Finder からD&D

[Homebrew]({{< relref "homebrew.md" >}})` や `[MacPorts]({{< relref "macports.md" >}})
などのパッケージ管理ツールを入れている人は
それらのルートディレクトリ(`/opt` など)を除外しといたほうがいい。
普通は表示されてないので、以下のコマンドで
Finder に表示させてから、D&Dで放り込む:

    open /opt

## `mds`, `mdworker`

外付けのボリュームを除外したい場合、
上記の方法では外す度に設定が戻ってしまうのでダメで、
下記のように `mdutil` を使う必要がある:

    # インデックス情報を表示
    sudo mdutil -s /Volumes/Macintosh\ HD

    # インデックスサービスを切る
    sudo mdutil -i off /Volumes/Macintosh\ HD

    # インデックスを作る・更新する
    sudo mdiutil -p /Volumes/Macintosh\ HD

    # インデックスを一旦削除して作り直し
    sudo mdutil -E /Volumes/Macintosh\ HD

日本語のせいで暴走してそうな場合は以下を消去:

    ~/Library/Preferences/com.apple.japaneseAnalysys/AppleContextualKKC.index/AdaptiveMap
    ~/Library/Preferences/com.apple.japaneseAnalysys/AppleContextualKKC.index/InputHistory.plist
