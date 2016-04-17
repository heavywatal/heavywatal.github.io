+++
title = 'SAMtools'
tags = ["genetics"]
[menu.main]
  parent = "bio"
+++

<http://www.htslib.org/>

## 操作

<http://www.htslib.org/doc/samtools.html>

### 下ごしらえ

-   適当な条件でフィルタリング:

        % samtools view -hb -f3 -q2 aln.bam -o filtered.bam

-   ソート:

        % samtools sort -T tmpsam -@2 -o sorted.bam aln.bam

    一時ファイルの名前なんか適当にユニークに決めてくれたらいいのに、
    `-T {PREFIX}` を指定しないと使えない。

    デフォルトでは位置でソートされるが、
    オプション `-n` で名前順(=ペアが隣り合う)にもできる。

-   標準入力stdinからパイプする場合は `-` を引数に:

        % samtools view -hb -f3 aln.bam | samtools sort -T tmpsam -@2 -o piped.bam

-   BAMインデックス作成:

        % samtools index aln.bam

    `aln.bam.bai` が書き出される

-   参照配列インデックス作成:

        % samtools faidx ref.fa

    `tview` や `pileup` などで必要になる
    `ref.fa.fai` が書き出される

-   PCR duplicatesを除去:

        % samtools rmdup aln.bam unique.bam

fixmate

### 閲覧

-   生のSAMを `less` で閲覧:

        % samtools view -h aln.bam | less

    `-o out.sam` とすればファイルに書き出せる。

-   マッピングされた形で閲覧:

        % samtools tview aln.bam [ref.fa]

    参照配列を与えると表示方法に選択肢が増える。
    予めインデックスを作っておかなくても初回実行時に
    `faidx` が勝手に走って `ref.fa.fai` が作成される。

    ヘルプ
    :   `shift+/`

    ジャンプ
    :   `g` または `/`
        `chr1:10000` のように染色体を付けて

-   インデックス作成済みBAMの統計量を表示:

        % samtools idxstats aln.bam
        chr1    249250621       6343976 0
        chr10   135534747       2407204 0
        chr11   135006516       3773511 0
        chr12   133851895       3696141 0
        ...

    seqnames, seqlengths, mapped reads, unmapped reads

-   フラグの要約統計:

        % samtools flagstat aln.bam
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

-   サイトごとのdepthをタブ区切りで:

        % samtools depth aln.bam
        chr1 12345 1
        chr1 12346 1
        ...

### SNP calling

VCF/BCFを書き出す:

    % samtools mpileup -uv aln.bam
    chr1 12345

`-f ref.fa`

calmd

------------------------------------------------------------------------

merge

reheader

cat

targetcut

phase

## SAM形式

<http://www.htslib.org/doc/sam.html>

<https://samtools.github.io/hts-specs/>

1.  `QNAME`: リード名
2.  `FLAG`: マッピングの状況をビット表現:

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
    <http://broadinstitute.github.io/picard/explain-flags.html>

    strandも場合分けして詳細に数え上げた例:
    <https://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/>

3.  `RNAME`: 参照配列・染色体の名前
4.  `POS`: 位置
5.  `MAPQ`: マッピングクオリティ
    = round($-10\log_{10}\mathrm{Pr[mapping~is~wrong]}$)

    マッパーによって微妙に違うらしい。
    例えばTopHatでは:

        50:        unique
         3:   2 locations
         2:   3 locations
         1: 4–9 locations
         0: >10 locations

6.  `CIGAR`: マッピング状況とその長さ
    e.g., `101M`, `18M200I83M`, `65M36S`
    :

        M: alignment match
        I: insertion
        D: deletion
        N: skipped = intron
        S: soft clipping = 部分マッチで無視された部分
        H: hard clipping = 部分マッチで無視された部分
        X: sequence mismatch

7.  `RNEXT`: paired-endの他方が張り付いた染色体。
    同じときは `=` で不明なときは `*`
8.  `PNEXT`: paired-endの他方が張り付いた位置
9.  `TLEN`: inferred Template LENgth。
    左腕の左端から右腕の右端までの距離。
    左腕なら正、右腕なら負、single-endなら0。
10. `SEQ`: 塩基配列
11. `QUAL`: 塩基クオリティ
12. それ以降はマッパー依存。
    形式は `TAG:VTYPE:VALUE`

## Binding

### pysam

<https://pysam.readthedocs.org/>

<http://bi.biopapyrus.net/python/modules/pysam.html>

### Rsamtools

<http://bioconductor.org/packages/release/bioc/html/Rsamtools.html>
