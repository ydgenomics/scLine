[Markdown+plantuml:最强流程图绘制](https://mp.weixin.qq.com/s/YbXcBRphGjPEaDgEvfXblA)
```plantuml
start
:10x matrix;
if (SoupX?) then (yes)
  :Run SoupX;
else (no)
endif

:Run scrublet;
if (have annotated scRNA data(.rds) as ref) then (yes)
  :Run singleR;
(no) elseif (have specific genes of specific cell-type as markers) then (yes)
  :Run scType;
else
  :nothing;
endif
:Run Integration;

stop
```