+++
date = "2016-12-14T18:22:07+09:00"
title = "concurrent.futures"
subtitle = "並列処理 in Python"
tags = ["python", "concurrent"]
[menu.main]
  parent = "python"
+++

標準ライブラリの `concurrent.futures` で簡単に並列化できる。

```py
import time
import random

def target_func(x):
    time.sleep(random.uniform(0, 1))
    return x + 1
```

## `threading`

https://docs.python.org/3/library/threading.html

**Global Interpreter Lock (GIL)** の制約により、
1つのインタープリタでは同時に1つのスレッドしかコードを実行できない。
したがってCPUバウンドなピュアPythonコードをマルチスレッド化しても速くならない。
外部プログラムやI/OなどGIL外の処理を待つ場合には有効。

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

このモジュールには `multiprocessing.Pool` のような機能が含まれていないので、
スレッド数の上限値を設けたい場合は
`threading.Semaphore` でうまくロックしてやる必要がある。


## `multiprocessing`

https://docs.python.org/3/library/multiprocessing.html

`os.fork()` を使って新しいインタプリタを立ち上げるので、
CPUバウンドなPythonコードもGILに邪魔されず並列処理できる。
ただし通信のため関数や返り値がpicklableでなければならない。

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


## `concurrent.futures` (since 3.2)

https://docs.python.org/3/library/concurrent.futures.html

上記のライブラリを直感的に使いやすくラップした高級インターフェイス。
立ち上げ方も簡単だし、 `as_completed()` で待てるのが特に便利。

```py
import concurrent.futures as confu

with confu.ThreadPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(target_func, x) for x in range(8)]
    (done, notdone) = confu.wait(futures)
    for future in futures:
        print(future.result())

with confu.ThreadPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(target_func, x) for x in range(8)]
    for future in confu.as_completed(futures):
        print(future.result())
```

`ProcessPoolExecutor` もある。

`concurrent.futures.process.multiprocessing.cpu_count()` もある。
