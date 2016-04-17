+++
title = 'ことえり'
tags = ["mac", "writing"]
[menu.main]
  parent = "mac"
+++

{{%div class="note"%}}
ことえり は使えないにもほどがあるので
[Google日本語入力](http://www.google.co.jp/ime/)
をインストールして使うべし。
以下は、まだ ことえり を使ってたときのメモ。
{{%/div%}}

学習辞書の本体:
:   `~/Library/Preferences/com.apple.JapaneseAnalysis/LearningDictionray.dict`

ユーザー辞書の所在:
:   `~/Library/Dictionary`

ライフサイエンス辞書
:   [ライフサイエンス辞書プロジェクト](http://lsd.pharm.kyoto-u.ac.jp/ja/index.html) から
    ことえり用の辞書をダウンロード。インストール方法は同梱のpdfに書いてある。「卵母細胞」と「扁桃体」が一発変換できるようになればインストール完了。

## 学習辞書の修正

1.  右上のことえり「単語登録／辞書編集」
2.  上バーの「辞書」→「新規ユーザー辞書作成」で適当な名前の辞書作成
3.  上バーの「辞書」→「テキストや辞書から取り込む」で `LearningDictionary.dict` を開く
4.  上バーの「辞書」→「テキストに書き出す」
5.  拡張子をcsvに変え、エクセルなどで編集し、拡張子をtxtに戻す
6.  上バーの「辞書」→「新規ユーザー辞書作成」で `LearningDictionary.dict` を作成
7.  上バーの「辞書」→「テキストや辞書から取り込む」で編集後のtxtを取り込む
8.  `~/library/Dictionaries/` にできた `LearningDictionary.dict` を本物に重ねて再起動

## DictionaryTrainer

たまにCPUを食うプロセス。
噂によると、メールにあるカタカナ語をせっせと辞書に登録してくださってるらしい。
もちろん、余計なお世話。実行権限を取り除いて無効化しておく:

    % sudo chmod 600 /System/Library/Input\ Methods/Kotoeri.app/Contents/Support/DictionaryTrainer.app

## 日本語モードでも半角スペース

普通、日本語モードのときは全角のスペースが入力される。
でも、ほとんど使わないし、むしろ半角スペースを入力したい場面のほうが多い。
ということで、いつでも半角スペースが入力されるように設定する。
全角スペースを入力したいときは `option + space`:

    % defaults write com.apple.inputmethod.Kotoeri zhsy -dict-add " " -bool false
    % killall Kotoeri

と思ったけど、やっぱり日本語モードのときは全角スペースのほうが自然だな。
でもメモとして残しておく。
あと、スペースだけじゃなくていろんな記号で設定できるらしい:

    !"#$%&'()*+,-./:;<=>?@[\]^_`{|}~¥

## 半角英数はRomajiではなくU.S.で

半角英数はことえりのRomajiでも打てるが、U.S.でやったほうがいい。
ことえりはアプリケーション（iTunesとか）との相性が悪くて問題を起こす場合があるので。
半角英数モードは最低1つオンにしておかなきゃいけないようなので、以下の順で設定する。

1.  System Preferences --&gt; Language & Text --&gt; Input Sources
2.  U.S.にチェック
3.  KotoeriのRomajiを外す
