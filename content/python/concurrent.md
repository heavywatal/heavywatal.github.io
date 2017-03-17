+++
date = "2016-12-14T18:22:07+09:00"
title = "concurrent.futures"
subtitle = "並行処理 in Python"
tags = ["python", "concurrent"]
[menu.main]
  parent = "python"
+++

## 概要

**Global Interpreter Lock (GIL)** の制約により、
1つのPythonインタープリタでは同時に1つのスレッドしかコードを実行できない。
したがってCPUバウンドなピュアPythonコードを
[`threading`](https://docs.python.org/3/library/threading.html)
でマルチスレッド化しても速くならない。
[`subprocess`](https://docs.python.org/3/library/subprocess.html)
による外部プログラム実行やI/OなどGIL外の処理を待つ場合には有効。

一方 [`multiprocessing`](https://docs.python.org/3/library/multiprocessing.html)
は新しいインタプリタを `os.fork()` で立ち上げるので、
CPUバウンドなPythonコードもGILに邪魔されず並列処理できる。
ただし通信のため関数や返り値がpicklableでなければならない。

それらの低級ライブラリを使いやすくまとめたのが
[`concurrent.futures`](https://docs.python.org/3/library/concurrent.futures.html)
(since 3.2) なので、とりあえずこれを使えばよい。
新しい [`asyncio`](https://docs.python.org/3/library/asyncio.html)
(since 3.4) は勝手が違いすぎてとっつきにくい。

並列化対象の関数の例:

```py
import time
import random

def target_func(x):
    time.sleep(random.uniform(0, 1))
    return x + 1
```

## `concurrent.futures`

```py
import os
import concurrent.futures as confu

# 呼び出し順に拾う
with confu.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
    futures = [executor.submit(target_func, x) for x in range(8)]
    (done, notdone) = confu.wait(futures)
    for future in futures:
        print(future.result())

# 終わったやつから拾う
with confu.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
    futures = [executor.submit(target_func, x) for x in range(8)]
    for future in confu.as_completed(futures):
        print(future.result())
```

デフォルトの `max_workers=None` では `5 * os.cpu_count()` になるらしい。

`multiprocessing` 由来の `ProcessPoolExecutor` も使い方は同じ。
こちらは `os.cpu_count()` がデフォルト。


## `threading`

`target` を指定して作った `Thread` インスタンスで `start()` するのが基本。

```py
import threading

threads = []
for i in range(8):
    th = threading.Thread(target=target_func, args=[i])
    th.start()
    threads.append(th)
for th in threads:
    th.join()
```

返り値を得たい場合などはクラスを継承していじる必要がある。
その場合は必ず `run()` メソッドをoverrideする。

```py
class Worker(threading.Thread):
    def __init__(self, target, name=None, args=(), kwargs={}):
        threading.Thread.__init__(self, None, target, name, args, kwargs)
        self._return = None
    def run(self):
        self._return = self._target(*self._args, **self._kwargs)
    def get(self, timeout=None):
        self.join(timeout)
        return self._return

threads = []
for i in range(8):
    th = Worker(target=target_func, args=[i])
    th.start()
    threads.append(th)
for th in threads:
    print(th.get())
```

スレッド数の上限値を設けたい場合は
`threading.Semaphore` でうまくロックしてやる必要がある。


## `multiprocessing`

`threading.Thread` とほぼ同じインターフェイスの
`Process` クラスも用意されているが、
`Pool` を使ったほうが楽チン。

```py
import multiprocessing as mp

with mp.Pool(processes=mp.cpu_count()) as pool:
    results = [pool.apply_async(target_func, [x]) for x in range(8)]
    for res in results:
        print(res.get())
```

### `mp.cpu_count()`

このためだけに `multiprocessing` をimportするのは億劫だったが、
3.4で `os.cpu_count()` が追加された。

Hyper-Threading (HT)が有効な場合は論理コア数が返ってくることに注意。
CPUを100%使い続ける数値計算とかだとそんなに並列化しても早くならない。
物理コア数を取得したい場合は
[`psutil`](https://github.com/giampaolo/psutil)
の `psutil.cpu_count(logical=False)` を使う。
標準ライブラリではないが、広く使われてるらしい。
