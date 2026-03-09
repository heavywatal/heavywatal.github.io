+++
title = 'Keyboard'
tags = ["mac", "writing"]
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
:   <kbd>⌃control</kbd><kbd>⏻power</kbd>

Sleep Display
:   <kbd>⇧shift</kbd><kbd>⌃control</kbd><kbd>⏏︎eject</kbd>
:   <kbd>⇧shift</kbd><kbd>⌃control</kbd><kbd>⏻power</kbd>

Sleep
:   <kbd>⌥option</kbd><kbd>⌘command</kbd><kbd>⏏︎eject</kbd>
:   <kbd>⌥option</kbd><kbd>⌘command</kbd><kbd>⏻power</kbd>

Restart
:   <kbd>⌃control</kbd><kbd>⌘command</kbd><kbd>⏏︎eject</kbd>

Shut Down
:   <kbd>⌃control</kbd><kbd>⌥option</kbd><kbd>⌘command</kbd><kbd>⏏︎eject</kbd>

[Startup key combinations](https://support.apple.com/en-us/HT201255)
:   機種によって違うので要確認。
:   Safe mode: キャッシュ削除、諸々のチェックと修復など
:   Reset NVRAM/PRAM: 画面や音の設定
:   [Reset SMC](https://support.apple.com/en-us/HT201295): 電源、ファン、内蔵カメラなどハードウェア関係


### Display Brightness

明るさ調整の<kbd>🔅<sub>F1</sub></kbd><kbd>🔆<sub>F2</sub></kbd>キーは通常iMacやMacBookにのみ有効で、外部ディスプレイには効かない。
[MonitorControl](https://github.com/MonitorControl/MonitorControl)
というアプリを用いることで外部ディスプレイもキーボードで調節できるようになる。


## System Preferences > Keyboard

Key Repeat & Delay Until Repeat を最速にする
:   Keyboard > スライダーを両方とも右端に

自動修正をオフにする
:   Text > 右のやつ全部チェック外す

ショートカットはデフォルト設定を使う
:   Shortcuts > Restore Defaults (左の項目ごとに、全部)
:   まっさらセットアップ直後に Restore Defaults を押して変わるところがあるのは謎。

### US Keyboard

左手小指<kbd>A</kbd>の隣、<kbd>⇪caps lock</kbd>を<kbd>⌃control</kbd>として使う
:   Keyboard > Modifier Keys...

入力言語切り替え
:   <kbd>⌃control</kbd><kbd>space</kbd>

半角英数入力
:   <kbd>⌃control</kbd><kbd>⇧shift</kbd><kbd>;</kbd>

ひらがな入力
:   <kbd>⌃control</kbd><kbd>⇧shift</kbd><kbd>j</kbd>

[Karabiner-Elements](https://github.com/pqrs-org/Karabiner-Elements)
を使うことで左右の<kbd>⌘command</kbd>単体押下を
<kbd>英数</kbd><kbd>かな</kbd>として扱えるようになる。
が、<kbd>⌃control</kbd><kbd>space</kbd>も慣れればそれほど面倒じゃなくなる。


### Input Sources

[Google日本語入力](https://www.google.co.jp/ime/)を使う。
`brew install google-japanese-ime`

macOS標準の日本語入力がどうしても使いにくいポイント:

- 変換候補ウィンドウを開いた場合 <kbd>return</kbd> を2回押さないと確定されない。
- "Windows-like shortcuts" をオンにすれば
  <kbd>return</kbd>1回で確定できるようになるが、
  そうすると
  <kbd>⌃n</kbd>, <kbd>⌃p</kbd>, <kbd>⌃k</kbd>
  などのEmacs/Cocoaキーバインドが崩れる。


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


### ⌥ option

|     | <kbd>⌥option</kbd> | <kbd>⌥option</kbd><kbd>⇧shift</kbd> |
| --- | --------- | ----------------- |
| <kbd>8</kbd> | • bullet | ° degree |
| <kbd>-</kbd> | – en dash | — em dash |
| <kbd>=</kbd> | ≠ not equal | ± plus minus |
| <kbd>[</kbd> | “ left double quotation | ” right double quotation |
| <kbd>]</kbd> | ‘ left single quotation | ’ right single quotation |
| <kbd>\\</kbd> | « left double angle quotation | » right double angle quotation |
| <kbd>&lt;</kbd> | ≤ less than or equal to | |
| <kbd>&gt;</kbd> | ≥ greater than or equal to | |

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


## Magic Trackpad

### 応答しなくなった場合に試すこと

1.  スイッチ off → on (触覚フィードバックの有無)
1.  System Preferences → Bluetooth ("Not Connected" or "Connection Rejected"?)
1.  有線接続
1.  一旦デバイスを削除して再登録。
    System Preferences はキーボード操作で完結できないので、
    有線マウスが無い場合は[blueutil](https://github.com/toy/blueutil)で操作。
    ```sh
    brew install blueutil
    blueutil --paired
    blueutil --unpair 00-00-00-00-00-00
    ```
1.  セーフモード、SMCリセット、PRAM/NVRAMリセット
