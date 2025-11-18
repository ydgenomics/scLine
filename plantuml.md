[Markdown+plantuml:最强流程图绘制](https://mp.weixin.qq.com/s/YbXcBRphGjPEaDgEvfXblA)
```plantuml
:10x matrix;
if (01_qc) then (doSoupX)
  :SoupX;
  note left
    rawPaths
    filterPaths
    sampleValues
  end note
  :scrublet;
else (not doSoupX)
  :scrublet;
  note right
    test
  end note
endif


if (02_anno) then (auto)
  if (haveRef?) then (yes)
    :singleR;
  (no) elseif (haveMarker?) then (yes)
    :scType;
  else (no)
    stop
  endif
  :summary;
  note left
    weight.txt
  end note
else (manual)
  if (03_enrich) then (haveGo)
    :clusterprofiler;
  else (no)
    stop
  endif

endif

:Run Integration;

stop
```

```mermaid
flowchart LR
A[01_qc] --> A0{doSoupX?}
A0 -->|yes| A1.1[SoupX] --> A1.2[scrublet]
A0 -->|no| A2.1[scrublet]
A --> B[02_anno]
B --> B0{haveRef?}
B0 -->|yes| B1[singleR]
B --> B01{haveMarker?} -->|yes| B2[scType]
B1 --> B3[summary]
B2 --> B3
A --> C[03_enrich]
C --> C0{haveProtein?} -->|yes| C1[eggnog-mapper] --> C2[clusterprofiler/gofigure] --> B3
B --> D[04_integration] --> D1[scVI,harmony,CCA,BBKNN,RLIGER] --> D2[scib-metrics]
B3 --> D
D --> E[05_metaneighbor]
D --> F[06_dea]
D --> H[07_pseudotime]
classDef mainNode fill:#ffcccc,stroke:#ff0000,stroke-width:3px
class A,B,C,D,E,F,H mainNode
linkStyle 4 stroke:#ff0000,stroke-width:3px
linkStyle 11 stroke:#ff0000,stroke-width:3px
linkStyle 16 stroke:#ff0000,stroke-width:3px
linkStyle 20 stroke:#ff0000,stroke-width:3px
linkStyle 21 stroke:#ff0000,stroke-width:3px
linkStyle 22 stroke:#ff0000,stroke-width:3px
```