+++
title = 'Keyboard'
tags = ["mac", "writing"]
[menu.main]
  parent = "mac"
+++

## Shortcuts

https://support.apple.com/en-us/HT201236

スクリーンショット（全画面）
:   <kbd>⇧shift</kbd><kbd>⌘command</kbd><kbd>3</kbd>

スクリーンショット（範囲指定）
:   <kbd>⇧shift</kbd><kbd>⌘command</kbd><kbd>4</kbd>

スクリーンショット（ウィンドウ選択）
:   <kbd>⇧shift</kbd><kbd>⌘command</kbd><kbd>4</kbd> then <kbd>space</kbd>

スクリーンショット（メニューから選択。動画も撮れる）
:   <kbd>⇧shift</kbd><kbd>⌘command</kbd><kbd>5</kbd>

Dialog (Restart/Sleep/Shut Down)
:   <kbd>⌃control</kbd><kbd>⏏︎eject</kbd>

Sleep Display
:   <kbd>⇧shift</kbd><kbd>⌃control</kbd><kbd>⏏︎eject</kbd>

Sleep
:   <kbd>⌥option</kbd><kbd>⌘command</kbd><kbd>⏏︎eject</kbd>

Restart
:   <kbd>⌃control</kbd><kbd>⌘command</kbd><kbd>⏏︎eject</kbd>

Shut Down
:   <kbd>⌃control</kbd><kbd>⌥option</kbd><kbd>⌘command</kbd><kbd>⏏︎eject</kbd>

Cocoaアプリのキーバインドを変更
:   `~/Library/KeyBindings/DefaultKeyBinding.dict`

### Display Brightness

明るさ調整の<kbd>🔅<sub>F1</sub></kbd><kbd>🔆<sub>F2</sub></kbd>キーは通常iMacやMacBookにのみ有効で、外部ディスプレイには効かない。
[ExternalDisplayBrightness](https://github.com/fnesveda/ExternalDisplayBrightness/)
というアプリを用いることで外部ディスプレイもキーボードで調節できるようになる。


## System Preferences > Keyboard

Key Repeat & Delay Until Repeat を最速にする
:   Keyboard > スライダーを両方とも右端に

自動修正をオフにする
:   Text > 右のやつ全部チェック外す

ショートカットはデフォルト設定を使う
:   Shortcuts > Restore Defaults (左の項目ごとに、全部)

### US Keyboard

左手小指<kbd>A</kbd>の隣、<kbd>⇪caps lock</kbd>を<kbd>⌃control</kbd>として使う
:   Keyboard > Modifier Keys...

入力言語切り替え
:   <kbd>⌃control</kbd><kbd>space</kbd>

半角英数入力
:   <kbd>⌃control</kbd><kbd>⇧shift</kbd><kbd>;</kbd>

ひらがな入力
:   <kbd>⌃control</kbd><kbd>⇧shift</kbd><kbd>j</kbd>

[⌘英かな](https://github.com/iMasanari/cmd-eikana)
もしくは
[Karabiner-Elements](https://pqrs.org/osx/karabiner/)
を使うことで左右の<kbd>⌘command</kbd>で入力言語切り替えができるようになる。
ただし、前者は頻繁に落ちるためアプリ再起動を要し、
後者は他の設定に干渉する上に反応が遅い。
また、画面共有やVirtualBoxでの<kbd>⌘command</kbd>の挙動が厄介になる。

## 特殊な文字・記号の入力

### ウムラウト、アクセント符号

Umlaut (ä ï ü ë ö ÿ)
:   <kbd>⌥option</kbd><kbd>u</kbd> + 母音

Acute accent (á í ú é ó)
:   <kbd>⌥option</kbd><kbd>e</kbd> + 母音

Grave accent (à ì ù è ò)
:   <kbd>⌥option</kbd><kbd>`</kbd> + 母音

Circumflex (â î û ê ô)
:   <kbd>⌥option</kbd><kbd>i</kbd> + 母音

Tilde (ã õ ñ)
:   <kbd>⌥option</kbd><kbd>n</kbd> + 母音

### Chinese Pinyin

System Preferences > Keyboard > Input Sources > Add "ABC - Extended"

1st tone (ā ē ī ō ū ǖ)
:   <kbd>⌥option</kbd><kbd>a</kbd> + vowel

2nd tone (á é í ó ú ǘ)
:   <kbd>⌥option</kbd><kbd>e</kbd> + vowel

3rd tone (ǎ ě ǐ ǒ ǔ ǚ)
:   <kbd>⌥option</kbd><kbd>v</kbd> + vowel

4th tone (à è ì ò ù ǜ)
:   <kbd>⌥option</kbd><kbd>`</kbd> + vowel

u with umlaut and tone (ǖ ǘ ǚ ǜ)
:   <kbd>⌥option</kbd><kbd>`</kbd> + <kbd>v</kbd>

### 絵文字・記号

<kbd>⌃control</kbd><kbd>⌘command</kbd><kbd>space</kbd>
でポップアップした窓から検索・入力すればキーボードから手を離さずに済む。

e.g.,
beer🍺🍻, metal🤘, muscle💪, thumb👍, smile😁🤣, option ⌥, schwa ə

## Entering Unicode Text and Symbols

1.  System Preferences > Keyboard > Input Sources で Unicode Hex Input にチェック
1.  右上メニューバーから Unicode Hex Input を選択
1.  <kbd>⌥option</kbd> を押しつつ番号を入力
    (e.g. <kbd>⌥option</kbd><kbd>2318</kbd> で ⌘ が入力される)
