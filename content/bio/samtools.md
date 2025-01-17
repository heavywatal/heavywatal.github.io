+++
title = 'SAMtools'
subtitle = 'Utilities for the Sequence Alignment/Map (SAM)'
tags = ["genetics"]
[menu.main]
  parent = "bio"
+++

<https://www.htslib.org/>

## 操作

<https://www.htslib.org/doc/samtools.html>

### 閲覧・要約

[`view`](https://www.htslib.org/doc/samtools-view.html)
:   BAM/CRAMをプレーンテキストのSAMとして閲覧:
    ```sh
    samtools view -h --no-PG aln.bam | less
    ```
    ここで `-h` をつけないとヘッダーが削れてしまうし、
    `--no-PG` をつけないと元ファイルに含まれていなかった情報
    (`@PG`) が付加されてしまうので要注意。
:   ファイルに含まれるヘッダー情報を改変せずに閲覧する
    `samtools view --header-only --no-PG` のショートカットとして
    `samtools head` サブコマンドが追加された。
:   BAM/CRAMに変換するのも、リードをフィルターするのもこのコマンド。名前が悪い。

[`tview`](https://www.htslib.org/doc/samtools-tview.html)
:   マッピングされた形でインタラクティブに閲覧:
    ```sh
    samtools tview aln.bam [ref.fa]
    ```
    参照配列を与えると表示方法に選択肢が増える。
    予めインデックスを作っておかなくても初回実行時に
    `faidx` が勝手に走って `ref.fa.fai` が作成される。
:   ヘルプ: <kbd>shift</kbd><kbd>/</kbd>
:   ジャンプ: <kbd>g</kbd> または <kbd>/</kbd> して `chr1:10000` のような形式

[`idxstats`](https://www.htslib.org/doc/samtools-idxstats.html)
:   インデックス作成済みBAMの統計量を表示:
    ```
    chr1    249250621       6343976 0
    chr10   135534747       2407204 0
    chr11   135006516       3773511 0
    chr12   133851895       3696141 0
    ...
    ```
    seqnames, seqlengths, mapped reads, unmapped reads

[`flagstat`](https://www.htslib.org/doc/samtools-flagstat.html)
:   フラグの要約統計:
    ```
    65182282 + 0 in total (QC-passed reads + QC-failed reads)
    6895859 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    65182282 + 0 mapped (100.00%:nan%)
    58286423 + 0 paired in sequencing
    29382202 + 0 read1
    28904221 + 0 read2
    52470794 + 0 properly paired (90.02%:nan%)
    55907370 + 0 with itself and mate mapped
    2379053 + 0 singletons (4.08%:nan%)
    456098 + 0 with mate mapped to a different chr
    169500 + 0 with mate mapped to a different chr (mapQ>=5)
    ```


### 下ごしらえ

-   PCR duplicatesを除去
    1. [`samtools collate`](https://www.htslib.org/doc/samtools-collate.html)
       でリード名ごとに並べる。
       順番は関係ないので `sort -n` より `collate` のほうが効率的。
       アラインメント直後は大概こうなっていて省略可能。
    1. [`samtools fixmate -m`](https://www.htslib.org/doc/samtools-fixmate.html)
       で `MC`, `ms` タグを付加。
    1. [`samtools sort`](https://www.htslib.org/doc/samtools-sort.html)
       で位置順にソート。
    1. [`samtools markdup -r`](https://www.htslib.org/doc/samtools-markdup.html)
       で重複を除去。

-   適当な条件でフィルタリング:
    ```
    samtools view -hb -f3 -q2 aln.bam -o filtered.bam
    ```

-   [`bgzip`](https://www.htslib.org/doc/bgzip.html)でFASTAやGFFを圧縮。
    インデックス(`.gzi`)を利用して部分的に展開して高速アクセスすることが可能。
    普通の `gzip` としても展開可能。
    拡張子はデフォルトで `.gz` だけど `.bgz` にすることもある。

-   インデックス作成:
    - [`samtools index`](https://www.htslib.org/doc/samtools-index.html)
      → BAMインデックス (`.bam.bai`)
    - [`samtools faidx`](https://www.htslib.org/doc/samtools-faidx.html)
      → 参照配列インデックス (`.fa.fai`)
    - [`tabix`](https://www.htslib.org/doc/tabix.html)
      → タブ区切りゲノムポジションインデックス (`.bgz.tbi`)\
      いろんな形式を扱える(`-p gff|bed|sam|vcf`)。
      位置順ソート且つbgzip圧縮されている必要がある。
    - [`bgzip -r`](https://www.htslib.org/doc/bgzip.html)
      → BGZFインデックス (`.gz.gzi`)\
      bgzip済みfastaを `faidx` するとついでに作ってもらえるし、
      `tabix` にはおそらく込み込みなので、明示的に作ることは少ない。

### variant calling

- https://samtools.github.io/bcftools/howtos/variant-calling.html
- https://samtools.github.io/bcftools/bcftools.html

VCF/BCFを書き出す:

```sh
bcftools mpileup -f ref.fa aln.bam | bcftools call -mv -Ob -o calls.bcf
```

以前はsamtoolsの機能だったが、bcftoolsが担うことになった。


## SAM形式

- <https://www.htslib.org/doc/sam.html>
- <https://samtools.github.io/hts-specs/>

1.  `QNAME`: リード名
1.  `FLAG`: マッピングの状況をビット表現:

        0x001     1  paired in sequencing
        0x002     2  mapped in proper pair
        0x004     4  unmapped
        0x008     8  mate unmapped
        0x010    16  reverse strand
        0x020    32  mate reverse strand
        0x040    64  first read in pair
        0x080   128  second read in pair
        0x100   256  not primary alignment
        0x200   512  not passing platform/vendor quality checks
        0x400  1024  PCR or optical duplicate
        0x800  2048  supplementary alignment

    paired-end でだいたいちゃんと張り付いたものをとるには
    `-f 3` (= 1 + 2)

    マップされたリード数を数えたいときなどは
    `-F 256` で各リードが1度だけ登場するようにフィルタできるが、
    `0x100` が立っててもmultiple hitで同点優勝してる可能性があり、
    このときどれがprimaryになるかはマッパー依存。
    重複遺伝子などを考慮する場合はむやみに捨ててはいけない。

    数字とフラグを変換してくれる便利 web app:
    <https://broadinstitute.github.io/picard/explain-flags.html>

    strandも場合分けして詳細に数え上げた例:
    <https://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/>

1.  `RNAME`: 参照配列・染色体の名前
1.  `POS`: 位置
1.  `MAPQ`: マッピングクオリティ
    = round($-10\log_{10}\Pr[\text{wrong}]$)

    マッパーによって微妙に違うらしい。
    例えばTopHatでは:

        50:        unique
         3:   2 locations
         2:   3 locations
         1: 4–9 locations
         0: >10 locations

1.  `CIGAR`: マッピング状況とその長さ
    e.g., `101M`, `18M200I83M`, `65M36S`:
    - `M`: alignment match
    - `I`: insertion to the reference
    - `D`: deletion from the reference
    - `N`: skipped (= intron for mRNA-to-genome; undefined otherwise)
    - `S`: soft clipping = `SEQ` の中で無視された部分
    - `H`: hard clipping = `SEQ` の外で無視された部分。先頭と末尾にのみ存在。
    - `P`: padding
    - `=`: sequence match
    - `X`: sequence mismatch

    `M|I|S|=|X` の長さを足せば `SEQ` の長さになる。

1.  `RNEXT`: paired-endの他方が張り付いた染色体。
    同じときは `=` で不明なときは `*`
1.  `PNEXT`: paired-endの他方が張り付いた位置
1.  `TLEN`: inferred Template LENgth。
    左腕の左端から右腕の右端までの距離。
    左腕なら正、右腕なら負、single-endなら0。
1.  `SEQ`: 塩基配列
1.  `QUAL`: 塩基クオリティ
1.  それ以降はマッパー依存。
    形式は `TAG:VTYPE:VALUE`

### CRAM

<https://www.htslib.org/workflow/cram.html>

参照配列からの差分だけを保持することで、BAMよりもコンパクトになりやすい。
裏を返せば、このCRAM単体では完結できない操作も出てくるので扱いに注意が必要。
BAMを置き換えて一般ユーザーの主流になるにはキャッシュの設計がイマイチな気がするけど、
種数あたりのサンプル数・リード数が多くなるほど恩恵も大きくなるからオッケー、なのかなぁ。

#### 参照配列を探しに行く優先順位

<https://www.htslib.org/doc/samtools.html#REFERENCE_SEQUENCES>

1.  samtoolsを呼ぶときの明示的なオプション, e.g., `--reference`.
1.  `M5` タグのハッシュ値 → 環境変数 `$REF_CACHE`.
    - 次の `$REF_PATH` 参照でダウンロードが生じた場合の保存先。
      中身は無圧縮の生FASTAなので容量に注意。
    - ハッシュ値のプレースホルダ `%s` を含める。
      `%2s/%s` のようにするとMD5先頭2文字を消費してディレクトリ構造を作れる。
      これはファイル数が増えすぎて起こる問題を回避するため。
    - デフォルトは `${XDG_CACHE_HOME}/hts-ref/%2s/%2s/%s` or `${HOME}/.cache/hts-ref/%2s/%2s/%s`.
      見えない場所で勝手に容量が膨らんでいくのも恐ろしいし、
      すぐ消してしまいそうな場所にあるのも恐ろしいので変えたほうがいい。
1.  `M5` タグのハッシュ値 → 環境変数 `$REF_PATH`.
    - 普通の `PATH` と同じようにコロン区切りで `$REF_CACHE` 同様 `%s` を含む。
    - デフォルトは `http://www.ebi.ac.uk/ena/cram/md5/%s`.
      つまり `REF_PATH` が空だと `UR` タグより先にインターネットに読みに行った挙げ句、
      無圧縮の生FASTAを勝手に保存する凶悪仕様。
1.  `UR` タグに書かれたファイル。ローカルのみ可、リモートは不可。
    無理やり相対パスにすることも可能だが、
    CRAMファイルからではなくコマンド実行時の `$PWD` からの相対なので実質使えない。

#### 設定例

-   <https://www.ebi.ac.uk/ena>
    に参照配列があって扱う種数も多くなければ、特に設定せずともそれなりに使える。
    新しい参照配列が必要になるたびに大きめのダウンロードが発生することと、
    それがホーム以下の見えないところに溜まっていくことだけ我慢。

-   参照配列を自分で管理するなら `M5` を無効化して
    `UR` のみで運用すればトラフィックやキャッシュの心配が無くなり単純。
    ```sh
    export REF_CACHE=/dev/null
    export REF_PATH=$REF_CACHE
    ```
    パスの変更などでタグを編集したい場合は
    [`samtools reheader`](https://www.htslib.org/doc/samtools-reheader.html) が使える。

-   勝手にEBIを見に行くのは止めたいけど `M5` は使いたい場合、
    何らかの文字を `REF_PATH` に入れておけばいい。
    公式例に習って `$REF_CACHE` を入れておく:
    ```sh
    export REF_CACHE="${HOME}/db/hts-ref-cache/%2s/%2s/%s"
    export REF_PATH="$REF_CACHE"
    ```
    手元のファイルを `REF_CACHE` に配置するには付属のスクリプトが使える。
    FASTAに複数の配列が入っていてもMD5は1本ずつ計算される。
    ```sh
    seq_cache_populate.pl -root ${REF_CACHE%%/\%*} <(gunzip -c genome.fa.gz)
    ```
    大元の参照配列置き場はテキトーに決めてもよくなるけど、
    キャッシュを管理するのもそれはそれで難しそう。
