+++
title = '遺伝研スパコン'
tags = ["job"]
[menu.main]
  parent = "bio"
+++

<https://sc.ddbj.nig.ac.jp/>

各ページのURLがコロコロ変わってリンク切れになりやすい。

## 利用開始

1.  利用規定等を熟読。
1.  手元のコンピュータで[SSH鍵ペアを生成]({{< relref "ssh.md" >}})しておく。
    [公式ドキュメント](https://sc.ddbj.nig.ac.jp/application/ssh_keys)
    に従って特殊なファイル名の鍵を新規作成してもいいけど、
    他の用途で作ったものが既にある場合は作り直さなくてもいい。
1.  [利用登録申請](https://sc.ddbj.nig.ac.jp/application/registration)のフォームを埋める。
    - 申請者
    - 所属機関
    - アカウント: 作っておいたSSH公開鍵をここでコピペ。
    - 責任者
1.  申請者と責任者にメールが届くので、それに従って誓約書PDFを管理者に送信。
1.  アカウント登録証が手元に届く。
1.  手元の `~/.ssh/config` に設定を追加:
    ```
    Host *.ddbj.nig.ac.jp
      User heavywatal
      IdentityFile ~/.ssh/id_ed25519
    ```
    ユーザ名と鍵ファイル名は適宜調整。
    これでsshコマンドを短く済ませられる。
1.  ゲートウェイノードにSSH接続:
    ```sh
    ssh gw.ddbj.nig.ac.jp
    ```
1.  それから `qlogin` コマンドで計算ノードにログインしていたのは過去の話。
    今後の運用はわからないが、システム移行中の2025年3月現在、
    ひとまずインタラクティブノード `a001`, `a002`, `a003`
    のどれかにSSH接続してから作業せよとのこと。
    ```sh
    ssh a001
    ```


## ファイルの送受信

[公式「システムへのファイル転送方法」](https://sc.ddbj.nig.ac.jp/guides/using_general_analysis_division/ga_data_transfer/)
にはscp, sftp, Asperaを使えと書かれてるけど、
[rsync]({{< relref "rsync.md" >}})
を使うのが簡単。

```sh
# send
rsync -auvC ~/input/ gw.ddbj.nig.ac.jp:~/input/

# receive
rsync -auvC gw.ddbj.nig.ac.jp:~/output/ ~/output/
```

ソースコードは当然[Git]({{< relref "git.md" >}})で管理。


## 環境整備

- [ハードウェア構成](https://sc.ddbj.nig.ac.jp/guides/hardware)
- [ソフトウェア構成](https://sc.ddbj.nig.ac.jp/guides/software/)
  (Ubuntu 24.04)
    - [Apptainer](https://sc.ddbj.nig.ac.jp/guides/software/Container/Apptainer/)


## ジョブ投入

SSH接続して、そのまま直に重いコマンドを実行してはいけない。
ジョブスケジューラを使って計算ノードに仕事を投げる必要がある。

PBS系の Sun/Univa/Altair Grid Engine が長らく使われていたが、2025年からSlurmに移行。

- <https://www.schedmd.com/slurm/>
- <https://sc.ddbj.nig.ac.jp/guides/software/JobScheduler/Slurm/>





### コマンド抜粋

<https://slurm.schedmd.com/man_index.html>

[`sbatch`](https://slurm.schedmd.com/sbatch.html)
:   [バッチスクリプト](#バッチスクリプト)の形でジョブを投入する。
    PBS系の `qsub` に相当。
:   [詳細は後述](#sbatch)

[`srun`](https://slurm.schedmd.com/srun.html)
:   実行可能ファイルを直に指定する形でジョブを投入する。
:   e.g., `srun --pty bash` でインタラクティブなシェルを起動。

[`sinfo`](https://slurm.schedmd.com/sinfo.html)
:   システム全体の使用状態を表示。
    `sinfo -alN` でノード毎の詳細表示。

[`squeue`](https://slurm.schedmd.com/squeue.html)
:   現在実行中のジョブ一覧。
    `squeue --me` で自分のジョブだけ表示。
:   [JOB-STATE-CODES](https://slurm.schedmd.com/squeue.html>):
    `CA:CANCELLED`, `CD:COMPLETED`, `F:FAILED`, `PD:PENDING`, `R:RUNNING`

[`sacct`](https://slurm.schedmd.com/sacct.html)
:   ジョブ実行履歴の表示。

[`scontrol`](https://slurm.schedmd.com/scontrol.html)
:   ジョブの詳細表示と条件変更など。
:   `scontrol show job <JOBID>` でジョブの詳細を表示。

[`scancel`](https://slurm.schedmd.com/scancel.html)
:   ジョブ削除。


### `sbatch`

<https://slurm.schedmd.com/sbatch.html>

[バッチスクリプト](#バッチスクリプト)の形でジョブを投げる。

```sh
man sbatch
sbatch --help
sbatch test.sh
```

オプションは `sbatch` コマンドに渡してもいいし、
後述のようにスクリプトの中に書いてもいい。
以下はよく使いそうなオプション抜粋。

`-h, --help`

`-D, --chdir=<directory>`
:   ジョブを実行するディレクトリ。
    デフォルトではカレント `${PWD}`。

`-a, --array=<indexes>`
:   アレイジョブとして投入する。
    タスクごとにパラメータを変えてプログラムを走らせたい場合は、
    スクリプトの中で環境変数 `SLURM_ARRAY_TASK_ID` を拾ってどうにかする。
:   e.g., `-a 1-4` で4個のタスクを持つアレイジョブ。
:   `-a 1-4%2` で同時実行数の上限を2にする。

`-c, --cpus-per-task=<ncpus>`
:   タスクあたりに使用するCPUコア数。
    デフォルトは1コア。

`--export=[ALL,]<variables>`
:   ジョブに引き継ぐ環境変数。
    デフォルトは `ALL` で、`sbatch` を実行したシェルにあるものすべて。
:   追加したり上書きしたりしたい場合は
    `--export=ALL,VAR1=value1`
    のようにカンマ区切りで指定する。

`-J, --job-name=<jobname>`
:   ジョブに名前をつける。
    デフォルトではスクリプト名が採用される。

`--mem=<size>[units]`
:   計算ノードに要求するRAM容量。
:   単位には K, M, G, T が使えて、省略するとメガバイト。
:   `--mem-per-cpu` でCPUコアあたりのRAMを指定することもできる。

`-o, --output=<filename_pattern>`
:   標準出力の書き出し先。
    デフォルトはジョブIDを使って `slurm-%j.out` に書き出される。
    ワーキングディレクトリからの相対パスだと思うけど明記されておらず不明。
:   `-e, --error=<filename_pattern>`
    を指定しなければ標準エラーも同じところに出力される。

`-p, --partition=<partition_names>`
:   どのグループの計算ノードに投げるか。
    PBS系でいうqueueに相当。
    何が選べるかはシステムの設定次第で、例えば遺伝研の場合は
    - `-p short`: 1時間以内に終わる軽いもの
    - `-p epyc`: 124日以内に終わる長いもの
    - `-p medium`: 124日以内に終わる長さで、メモリを多めに使うもの

`-t, --time=<time>`
:   時間を制限する。
    ここでの宣言が上記partitionの設定より長いと一生pending。
:   形式は `MM`, `MM:SS`, `HH:MM:SS`, `DD-HH`, `DD-HH:MM`, `DD-HH:MM:SS` のいずれか。


### バッチスクリプト

`sbatch` コマンドに渡すスクリプト。
`#SBATCH` で始まる行は `sbatch` へのオプションと見なされる。

スクリプト内で参照可能な環境変数をプリントしてみるジョブの例:

```sh
#!/bin/bash
#SBATCH -t 00-00:01:00
#SBATCH --mem-per-cpu 1g
#SBATCH -J print

date -Iseconds

echo SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID-}
echo SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID-}
echo SLURM_ARRAY_TASK_COUNT: ${SLURM_ARRAY_TASK_COUNT-}
echo SLURM_ARRAY_TASK_MIN: ${SLURM_ARRAY_TASK_MIN-}
echo SLURM_ARRAY_TASK_MAX: ${SLURM_ARRAY_TASK_MAX-}
echo SLURM_ARRAY_TASK_STEP: ${SLURM_ARRAY_TASK_STEP-}
echo SLURM_JOB_ID: ${SLURM_JOB_ID-}
echo SLURM_JOB_NAME: ${SLURM_JOB_NAME-}
echo SLURM_JOB_NODELIST: ${SLURM_JOB_NODELIST-}
echo SLURM_JOB_PARTITION: ${SLURM_JOB_PARTITION-}
echo SLURM_JOB_START_TIME: ${SLURM_JOB_START_TIME-}
echo SLURM_MEM_PER_CPU: ${SLURM_MEM_PER_CPU-}
echo SLURM_MEM_PER_NODE: ${SLURM_MEM_PER_NODE-}
echo SLURM_SUBMIT_DIR: ${SLURM_SUBMIT_DIR-}
echo SLURM_SUBMIT_HOST: ${SLURM_SUBMIT_HOST-}
echo SLURM_TASK_PID: ${SLURM_TASK_PID-}
echo SLURMD_NODENAME: ${SLURMD_NODENAME-}

echo HOME: ${HOME-}
echo USER: ${USER-}
echo PWD: ${PWD-}
echo PATH: ${PATH-}

date -Iseconds
```


## Apptainer (Singularity)

- <https://sc.ddbj.nig.ac.jp/guides/software/Container/Apptainer/>
- <https://sc.ddbj.nig.ac.jp/guides/software/Container/BioContainers/>
- [Apptainer]({{< relref "apptainer.md" >}})

`/usr/local/biotools/` 以下に各種ソフトウェアが用意されている。
[BioContainers](https://biocontainers.pro) のものをほぼそのまま置いているらしい。

利用可能なソフトウェアとバージョンを探す:
```sh
find /usr/local/biotools/ -name 'blast*' | sort
```

イメージとプログラム名を指定して実行:
```sh
apptainer exec -e /usr/local/biotools/f/fastp:0.23.4--hadf994f_2 fastp --help
```

そこらに落ちてるイメージを拾ってきて使うこともできる。
[例えばTrinityの公式最新版を使いたい場合](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker#trinity_singularity):
```sh
find /usr/local/biotools/ -name 'trinity*' | sort
mkdir -p ~/image
wget -P ~/image/ https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/trinityrnaseq.v2.11.0.simg
apptainer exec -e ~/image/trinityrnaseq.v2.11.0.simg Trinity --help
```


## R

<https://sc.ddbj.nig.ac.jp/guides/software/DevelopmentEnvironment/R/>

### Apptainer R

存在するけどエラーで動かない:
```sh
find /usr/local/biotools/ -name 'r-base:*' | sort
apptainer exec /usr/local/biotools/r/r-base:4.4.1 R --no-save --no-restore-data
```
```
 *** caught segfault ***
address (nil), cause 'memory not mapped'
```
