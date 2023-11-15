

## 开始

### blockinfo

* 软件支持 Fix，Variable，Index 等几种模式的 block 类型。

* Fix 指的是给定的碱基序列，需要提前提供，legoseq 会使用比对算法确定 Fix 序列在 read 中的位置；

* Variable 指的是可变序列，需要根据上下游的 Fix 序列来确定；

* Index 指的barcode 序列，需要提前提供，可以据此进行序列 demultiplex。