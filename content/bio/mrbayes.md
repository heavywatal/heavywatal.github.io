+++
title = 'MrBayes'
tags = ["genetics"]
[menu.main]
  parent = "bio"
+++

<http://mrbayes.sourceforge.net/>

## Installation

1.  OpenMPI をインストール:

        % sudo apt-get install libopenmpi1.5-dev openmpi1.5-bin

2.  ソースコードをダウンロードして展開:

        % wget -O- http://sourceforge.net/projects/mrbayes/files/mrbayes/3.2.2/mrbayes-3.2.2.tar.gz/download | tar xz

3.  コンパイルの仕方を調べて実行:

        % cd mrbayes-3.2.2/src/
        % less CompileInstructions.txt
        % autoconf
        % ./configure --with-beagle=no --enable-mpi=yes

4.  ビルド:

        % make

    このようなエラーが出たら:

        `undefined reference to __exp_finite`
        `undefined reference to __log_finite`

    `LDFLAGS` がオブジェクトのあとに来るように、
    `Makefile` のターゲット mb を書き換える:

        mb: $(OBJECTS)
                $(CC) $(CFLAGS) $^ $(LDFLAGS) $(OUTPUT_OPTION)

5.  インストール:

        % sudo cp mb /usr/local/bin/
        % sudo cp mb /nfsroot/usr/local/bin/
