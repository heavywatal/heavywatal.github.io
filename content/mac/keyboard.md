+++
title = 'Keyboard'
tags = ["mac", "writing"]
[menu.main]
  parent = "mac"
+++

## Shortcuts

https://support.apple.com/en-us/HT201236

スクリーンショット（全画面）
:   <kbd>shift</kbd><kbd>command</kbd><kbd>3</kbd>

スクリーンショット（範囲指定）
:   <kbd>shift</kbd><kbd>command</kbd><kbd>4</kbd>

スクリーンショット（ウィンドウ選択）
:   <kbd>shift</kbd><kbd>command</kbd><kbd>4</kbd> then <kbd>space</kbd>

スクリーンショット（メニューから選択。動画も撮れる）
:   <kbd>shift</kbd><kbd>command</kbd><kbd>5</kbd>

Dialog (Restart/Sleep/Shut Down)
:   <kbd>control</kbd><kbd>eject</kbd>

Sleep Display
:   <kbd>shift</kbd><kbd>control</kbd><kbd>eject</kbd>

Sleep
:   <kbd>option</kbd><kbd>command</kbd><kbd>eject</kbd>

Restart
:   <kbd>control</kbd><kbd>command</kbd><kbd>eject</kbd>

Shut Down
:   <kbd>control</kbd><kbd>option</kbd><kbd>command</kbd><kbd>eject</kbd>

Cocoaアプリのキーバインドを変更
:   `~/Library/KeyBindings/DefaultKeyBinding.dict`


## US Keyboard

左手小指Aの隣、<kbd>caps lock</kbd>を<kbd>control</kbd>として使う
:   System Preferences > Keyboard > Modifier Keys...

入力言語切り替え
:   <kbd>control</kbd><kbd>space</kbd>

半角英数入力
:   <kbd>control</kbd><kbd>shift</kbd><kbd>;</kbd>

ひらがな入力
:   <kbd>control</kbd><kbd>shift</kbd><kbd>j</kbd>

[⌘英かな](https://github.com/iMasanari/cmd-eikana)
もしくは
[Karabiner-Elements](https://pqrs.org/osx/karabiner/)
を使うことで左右の<kbd>command</kbd>で入力言語切り替えができるようになる。
ただし、前者は頻繁に落ちるためアプリ再起動を要し、
後者は他の設定に干渉する上に反応が遅い。
また、画面共有やVirtualBoxでの<kbd>command</kbd>の挙動が厄介になる。


## ウムラウト、アクセント符号

<kbd>option</kbd> キーを使って修飾つきアルファベットを入力できる

Umlaut (ä ï ü ë ö ÿ)
:   <kbd>option</kbd><kbd>u</kbd> + 母音

Acute accent (á í ú é ó)
:   <kbd>option</kbd><kbd>e</kbd> + 母音

Grave accent (à ì ù è ò)
:   <kbd>option</kbd><kbd>`</kbd> + 母音

Circumflex (â î û ê ô)
:   <kbd>option</kbd><kbd>i</kbd> + 母音

Tilde (ã õ ñ)
:   <kbd>option</kbd><kbd>n</kbd> + 母音

## Entering Unicode Text and Symbols

1.  System Preferences > Keyboard > Input Sources で Unicode Hex Input にチェック
1.  右上メニューバーから Unicode Hex Input を選択
1.  <kbd>option</kbd> を押しつつ番号を入力
    (e.g. <kbd>option</kbd><kbd>2318</kbd> で ⌘ が入力される)
