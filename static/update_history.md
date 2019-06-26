# 16s pipelines implemented with qiime2 API
用以记录更新历史.

2019.04.04 上午

第一次提交,并完善了API

目前:
 1. 只能处理demuxed的pairend的数据
 2. 通过填写配置文件可以一键启动deblur和dada2的流程.基于qiime2官方推荐流程.
 3. 最终输出多个qza和qza,并且以纯文本形式输出otu table和代表序列(其中以sOTU进行了重命名)
 
2019.04.04 下午

1. 发现不能处理未demux的数据实在太麻烦,由于dada2和deblur需要的数据是unjoined的,所以最好先分库(并且处理好mix-orientation的问题)
2. qiime2 内部的demux太不完善,使用的是cutadapt, 但无法接收双向的barcodes.
3. 决定使用自己写的demux方法,并且提供fix orientation的方法,并嵌入该流程