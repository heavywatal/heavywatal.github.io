+++
title = 'Mac固有コマンド'
tags = ["mac"]
[menu.main]
  parent = "mac"
+++

## 簡単・便利

### `open`

Finder でダブルクリックするのと同じように、
Terminal から一発でファイルを開くことができる。

```sh
# 関連付けられたデフォルトのアプリケーションで開く
open hudson1992g.pdf

# アプリケーションを指定して開く
open -a Skim hudson1992g.pdf

# カレントディレクトリをFinderで開く
open .
```

### `pbcopy`, `pbpaste`

ターミナルからクリップボード(pasteboard)を操作する。

```sh
pwd | pbcopy
pbpaste | wc
```


### `say` 音声読み上げ

```sh
say hello world
pbpaste | say
```

声は環境設定の Accessibility > Speech で変更可能。

カーソル行あるいは選択中のテキストを読み上げるショートカットも設定できる。
<kbd>option</kbd><kbd>esc</kbd>


### `killall`

Finder や Dock など、GUIから終了させにくいアプリケーションを再起動する:

```sh
killall Finder
killall Dock
killall Kotoeri
```

動作が不安定になったとき、設定変更を反映させたいとき、メモリを開放したいときなどに。


## 設定関連

### `defaults`

`~/Library/Preferences/` 以下にある各種設定ファイル.plistを編集する。
`true` に設定した項目を元に戻すには、
項目自体を `delete` するか、`false` に設定する。

```sh
defaults [write/delete] DOMAIN KEY -TYPE VALUE

# Finderのタイトルバーにフルパスを表示
defaults write com.apple.finder _FXShowPosixPathInTitle -boolean true

# Launchpadの並び順をリセット
defaults write com.apple.dock ResetLaunchPad -bool true; killall Dock

# 使わないXcodeのための重いindexingを切っておく
defaults write com.apple.dt.Xcode IDEIndexDisable 1

# Quicklook上でコピペできるようにする
defaults write com.apple.finder QLEnableTextSelection -bool true

# 特定アプリをダークモードから除外する
# osascript -e 'id of app "Google Chrome"'
defaults write com.google.Chrome NSRequiresAquaSystemAppearance -bool true

# スクリーンショットの保存先を変更
defaults write com.apple.screencapture location "~/Pictures/Screenshot"
# 影なしで撮る
defaults write com.apple.screencapture disable-shadow -bool true
# 撮影後のサムネイル表示(=遅延)を無効にする
defaults write com.apple.screencapture show-thumbnail -bool false
```

[OnyX](https://www.titanium-software.fr/en/onyx.html) や
[TinkerTool](https://www.bresink.com/osx/TinkerTool.html)
などのGUIアプリを使うほうが簡単で安心かも


### `launchctl`

サービスの起動、終了

MacBook Pro の Touch Bar は完全なる失敗作であり欠陥品である。
目視を要求してtouch typingを阻害するばかりでなく、
意図せぬ接触による誤動作を誘発する。
例えば、再生ボタンをうっかり触るとMusic.appが起動してしまう。
次のコマンドでそれを防げる:
```sh
launchctl unload -w /System/Library/LaunchAgents/com.apple.rcd.plist
```
ただし副作用としてそれ以外の用途でも再生ボタンを使えなくなる。
Music.appの起動だけを防ぐオプションは存在しないのかな？


### `lsregister`

Open with で表示されるアプリケーションが重複しまくったときなど、
ファイルとアプリケーションの関連付けに関する古い情報を消して再構築。
まっさらに戻るわけではない。:

    /System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister -kill -r -domain local -domain system -domain user

設定ファイルは `~/Library/Preferences/com.apple.LaunchServices.plist`


### `mds`, `mdworker`

Spotlightの設定やインデックス作成に関わる。
外付けのボリュームを除外したい場合、
GUIから環境設定してもejectするたびに戻ってしまうので、
下記のように `mdutil` を使う必要がある:

```sh
# インデックス情報を表示
sudo mdutil -s /Volumes/Macintosh\ HD

# インデックスサービスを切る
sudo mdutil -i off /Volumes/Macintosh\ HD

# インデックスを作る・更新する
sudo mdiutil -p /Volumes/Macintosh\ HD

# インデックスを一旦削除して作り直し
sudo mdutil -E /Volumes/Macintosh\ HD
```


### インストールなど

`.dmg` のマウント、`.pkg` からのインストール、
システムのソフトウェア・アップデートなどを
`ssh` 越しにやらねばならぬときもある:

```sh
hdiutil mount SomeDiskImage.dmg

sudo installer -pkg SomePackage.pkg -target /

softwareupdate -l
softwareupdate -i -a
```
