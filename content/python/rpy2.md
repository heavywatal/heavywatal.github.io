+++
title = 'rpy2'
subtitle = "PythonからRを呼び出す"
[menu.main]
  parent = "python"
+++

<http://rpy.sourceforge.net/rpy2.html>

## Example

```python
import rpy2.robjects as robjects
R = robjects.r

R('''
png("plot1.png")
plot(rnorm(1000), rnorm(1000))
dev.off()
''')

R.png("plot2.png")
R.plot(robjects.FloatVector(R.runif(10)), ann=False)
R["dev.off"]()
```

## Installation

1.  [pip]({{< relref "pip.md" >}}) を使ってインストール:

        % pip install rpy2
