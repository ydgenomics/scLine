digraph "pipelinedag" {
rankdir=TB;
v0 [shape=point,label="",fixedsize=true,width=0.1];
v4 [label="CONVERT:CONVERT_0"];
v0 -> v4 [label="inputFile"];

v1 [shape=point,label="",fixedsize=true,width=0.1];
v4 [label="CONVERT:CONVERT_0"];
v1 -> v4 [label="assay"];

v2 [shape=point,label="",fixedsize=true,width=0.1];
v4 [label="CONVERT:CONVERT_0"];
v2 -> v4 [label="python_env"];

v3 [shape=point,label="",fixedsize=true,width=0.1];
v4 [label="CONVERT:CONVERT_0"];
v3 -> v4 [label="outdir"];

v4 [label="CONVERT:CONVERT_0"];
v5 [shape=point];
v4 -> v5 [label="results"];

}