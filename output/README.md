# example1

<details> 

```shell
$ tree /data/work/scline/output/example1
/data/work/scline/output/example1
├── anno
│   ├── DotPlot_group1.hr.rds.pdf
│   ├── group1.hr_anno.csv
│   ├── group1.hr_nodes.csv
│   ├── group1.hr_sctype.rds
│   ├── group1.hr_sctype_scores_sorted.csv
│   ├── group1.hr_sctype_umap.pdf
│   └── VlnPlot_group1.hr.rds.pdf
├── cytotrace
│   ├── example1_ct.h5ad
│   ├── example1_violin_pseudotime.pdf
│   └── figures
│       ├── scvelo_example1_ct_pseudotime.pdf
│       ├── umapexample1_ct_all.pdf
│       └── umapexample1_umap.pdf
├── dea
│   └── markers_example1.hr.rds_celltype2_vs_celltype3.csv
├── enrich
│   ├── group1_leiden_res_0.50.csv_0_enrich.txt
│   ├── group1_leiden_res_0.50.csv_0_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_1_enrich.txt
│   ├── group1_leiden_res_0.50.csv_1_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_2_enrich.txt
│   ├── group1_leiden_res_0.50.csv_2_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_3_enrich.txt
│   ├── group1_leiden_res_0.50.csv_3_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_4_enrich.txt
│   ├── group1_leiden_res_0.50.csv_4_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_5_enrich.txt
│   ├── group1_leiden_res_0.50.csv_5_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_6_enrich.txt
│   ├── group1_leiden_res_0.50.csv_6_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_7_enrich.txt
│   ├── group1_leiden_res_0.50.csv_7_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_8_enrich.txt
│   ├── group1_leiden_res_0.50.csv_8_plot1.pdf
│   ├── kegg_info.RData
│   ├── markers_example1.hr.rds_celltype2_vs_celltype3.csv_celltype2_vs_celltype3_enrich.txt
│   ├── markers_example1.hr.rds_celltype2_vs_celltype3.csv_celltype2_vs_celltype3_plot1.pdf
│   └── org.Gspecies.eg.db
│       ├── DESCRIPTION
│       ├── inst
│       │   └── extdata
│       │       └── org.Gspecies.eg.sqlite
│       ├── man
│       │   ├── org.Gspecies.egBASE.Rd
│       │   ├── org.Gspecies.eg_dbconn.Rd
│       │   └── org.Gspecies.egORGANISM.Rd
│       ├── NAMESPACE
│       └── R
│           └── zzz.R
├── example1.hr.rds
├── execution-report.html
├── execution-timeline.html
├── group1.hr.rds
├── pipeline-dag.html
└── qc
    └── group1
        ├── figures
        │   ├── dotplot_leiden_res_0.50_marker.pdf
        │   ├── pca_potentially_undesired_features.pdf
        │   ├── umap_batch.pdf
        │   ├── umap_leiden_clus.pdf
        │   └── umap_quality.pdf
        ├── group1.h5ad
        ├── markers.csv
        │   └── group1_leiden_res_0.50.csv
        ├── qc.pdf
        └── summary.txt

15 directories, 55 files
```

</details>


# example2

<details>

```shell
$ tree /data/work/scline/output/example2
/data/work/scline/output/example2
├── anno
│   ├── DotPlot_group1.hr.rds.pdf
│   ├── DotPlot_group2.hr.rds.pdf
│   ├── group1.hr_anno.csv
│   ├── group1.hr_nodes.csv
│   ├── group1.hr_sctype.rds
│   ├── group1.hr_sctype_scores_sorted.csv
│   ├── group1.hr_sctype_umap.pdf
│   ├── group2.hr_anno.csv
│   ├── group2.hr_nodes.csv
│   ├── group2.hr_sctype.rds
│   ├── group2.hr_sctype_scores_sorted.csv
│   ├── group2.hr_sctype_umap.pdf
│   ├── VlnPlot_group1.hr.rds.pdf
│   └── VlnPlot_group2.hr.rds.pdf
├── cytotrace
│   ├── example2_rliger.INMF_integrated.rh_ct.h5ad
│   ├── example2_rliger.INMF_integrated.rh_violin_pseudotime.pdf
│   └── figures
│       ├── scvelo_example2_rliger.INMF_integrated.rh_ct_pseudotime.pdf
│       ├── umapexample2_rliger.INMF_integrated.rh_ct_all.pdf
│       └── umapexample2_rliger.INMF_integrated.rh_umap.pdf
├── dea
│   └── markers_example2_rliger.INMF_integrated.rds_celltype2_vs_celltype3.csv
├── enrich
│   ├── group1_leiden_res_0.50.csv_0_enrich.txt
│   ├── group1_leiden_res_0.50.csv_0_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_1_enrich.txt
│   ├── group1_leiden_res_0.50.csv_1_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_2_enrich.txt
│   ├── group1_leiden_res_0.50.csv_2_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_3_enrich.txt
│   ├── group1_leiden_res_0.50.csv_3_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_4_enrich.txt
│   ├── group1_leiden_res_0.50.csv_4_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_5_enrich.txt
│   ├── group1_leiden_res_0.50.csv_5_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_6_enrich.txt
│   ├── group1_leiden_res_0.50.csv_6_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_7_enrich.txt
│   ├── group1_leiden_res_0.50.csv_7_plot1.pdf
│   ├── group1_leiden_res_0.50.csv_8_enrich.txt
│   ├── group1_leiden_res_0.50.csv_8_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_0_enrich.txt
│   ├── group2_leiden_res_0.50.csv_0_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_10_enrich.txt
│   ├── group2_leiden_res_0.50.csv_10_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_1_enrich.txt
│   ├── group2_leiden_res_0.50.csv_1_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_2_enrich.txt
│   ├── group2_leiden_res_0.50.csv_2_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_3_enrich.txt
│   ├── group2_leiden_res_0.50.csv_3_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_4_enrich.txt
│   ├── group2_leiden_res_0.50.csv_4_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_5_enrich.txt
│   ├── group2_leiden_res_0.50.csv_5_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_6_enrich.txt
│   ├── group2_leiden_res_0.50.csv_6_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_7_enrich.txt
│   ├── group2_leiden_res_0.50.csv_7_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_8_enrich.txt
│   ├── group2_leiden_res_0.50.csv_8_plot1.pdf
│   ├── group2_leiden_res_0.50.csv_9_enrich.txt
│   ├── group2_leiden_res_0.50.csv_9_plot1.pdf
│   ├── kegg_info.RData
│   ├── markers_example2_rliger.INMF_integrated.rds_celltype2_vs_celltype3.csv_celltype2_vs_celltype3_enrich.txt
│   ├── markers_example2_rliger.INMF_integrated.rds_celltype2_vs_celltype3.csv_celltype2_vs_celltype3_plot1.pdf
│   └── org.Gspecies.eg.db
│       ├── DESCRIPTION
│       ├── inst
│       │   └── extdata
│       │       └── org.Gspecies.eg.sqlite
│       ├── man
│       │   ├── org.Gspecies.egBASE.Rd
│       │   ├── org.Gspecies.eg_dbconn.Rd
│       │   └── org.Gspecies.egORGANISM.Rd
│       ├── NAMESPACE
│       └── R
│           └── zzz.R
├── execution-report.html
├── execution-timeline.html
├── group1.hr.rds
├── group2.hr.rds
├── integrate
│   ├── example2_harmony_integrated_dealplus.h5ad
│   ├── example2_harmony_integrated_plus.pdf
│   ├── example2_rliger.INMF_integrated_dealplus.rds
│   ├── example2_rliger.INMF_integrated_plus.pdf
│   ├── example2_SCTransform.CCA_integrated_dealplus.rds
│   ├── example2_SCTransform.CCA_integrated_plus.pdf
│   ├── example2_SCTransform.harmony_integrated_dealplus.rds
│   ├── example2_SCTransform.harmony_integrated_plus.pdf
│   ├── example2_unintegration_integrated_dealplus.h5ad
│   ├── example2_unintegration_integrated_plus.pdf
│   ├── otherexample2_harmony_integrated_plus.pdf
│   ├── otherexample2_unintegration_integrated_plus.pdf
│   └── scib
│       ├── example2_rliger.INMF_integrated.rh.h5ad
│       ├── example2_scIB.csv
│       ├── example2_scIB.h5ad
│       └── example2_scIB.pdf
├── metaneighbor
│   ├── example2_metaNeighbor.csv
│   ├── example2_metaNeighbor.pdf
│   ├── example2_metaNeighbor_tophits.csv
│   └── example2_sanky.html
├── pipeline-dag.html
└── qc
    ├── group1
    │   ├── figures
    │   │   ├── dotplot_leiden_res_0.50_marker.pdf
    │   │   ├── pca_potentially_undesired_features.pdf
    │   │   ├── umap_batch.pdf
    │   │   ├── umap_leiden_clus.pdf
    │   │   └── umap_quality.pdf
    │   ├── group1.h5ad
    │   ├── markers.csv
    │   │   └── group1_leiden_res_0.50.csv
    │   ├── qc.pdf
    │   └── summary.txt
    └── group2
        ├── figures
        │   ├── dotplot_leiden_res_0.50_marker.pdf
        │   ├── pca_potentially_undesired_features.pdf
        │   ├── umap_batch.pdf
        │   ├── umap_leiden_clus.pdf
        │   └── umap_quality.pdf
        ├── group2.h5ad
        ├── markers.csv
        │   └── group2_leiden_res_0.50.csv
        ├── qc.pdf
        └── summary.txt

21 directories, 113 files
```

</details>


# example3

<details>

```shell
$ tree /data/work/scline/output/example3
/data/work/scline/output/example3
├── anno
│   ├── DotPlot_time1.hr.rds.pdf
│   ├── DotPlot_time2.hr.rds.pdf
│   ├── DotPlot_time3.hr.rds.pdf
│   ├── time1.hr_anno.csv
│   ├── time1.hr_nodes.csv
│   ├── time1.hr_sctype.rds
│   ├── time1.hr_sctype_scores_sorted.csv
│   ├── time1.hr_sctype_umap.pdf
│   ├── time2.hr_anno.csv
│   ├── time2.hr_nodes.csv
│   ├── time2.hr_sctype.rds
│   ├── time2.hr_sctype_scores_sorted.csv
│   ├── time2.hr_sctype_umap.pdf
│   ├── time3.hr_anno.csv
│   ├── time3.hr_nodes.csv
│   ├── time3.hr_sctype.rds
│   ├── time3.hr_sctype_scores_sorted.csv
│   ├── time3.hr_sctype_umap.pdf
│   ├── VlnPlot_time1.hr.rds.pdf
│   ├── VlnPlot_time2.hr.rds.pdf
│   └── VlnPlot_time3.hr.rds.pdf
├── cytotrace
│   ├── example3_harmony_integrated_ct.h5ad
│   ├── example3_harmony_integrated_violin_pseudotime.pdf
│   └── figures
│       ├── scvelo_example3_harmony_integrated_ct_pseudotime.pdf
│       ├── umapexample3_harmony_integrated_ct_all.pdf
│       └── umapexample3_harmony_integrated_umap.pdf
├── dea
│   └── markers_example3_rliger.INMF_integrated.rds_celltype2_vs_celltype3.csv
├── enrich
│   ├── kegg_info.RData
│   ├── markers_example3_rliger.INMF_integrated.rds_celltype2_vs_celltype3.csv_celltype2_vs_celltype3_enrich.txt
│   ├── markers_example3_rliger.INMF_integrated.rds_celltype2_vs_celltype3.csv_celltype2_vs_celltype3_plot1.pdf
│   ├── org.Gspecies.eg.db
│   │   ├── DESCRIPTION
│   │   ├── inst
│   │   │   └── extdata
│   │   │       └── org.Gspecies.eg.sqlite
│   │   ├── man
│   │   │   ├── org.Gspecies.egBASE.Rd
│   │   │   ├── org.Gspecies.eg_dbconn.Rd
│   │   │   └── org.Gspecies.egORGANISM.Rd
│   │   ├── NAMESPACE
│   │   └── R
│   │       └── zzz.R
│   ├── time1_leiden_res_0.50.csv_0_enrich.txt
│   ├── time1_leiden_res_0.50.csv_0_plot1.pdf
│   ├── time1_leiden_res_0.50.csv_1_enrich.txt
│   ├── time1_leiden_res_0.50.csv_1_plot1.pdf
│   ├── time1_leiden_res_0.50.csv_2_enrich.txt
│   ├── time1_leiden_res_0.50.csv_2_plot1.pdf
│   ├── time1_leiden_res_0.50.csv_3_enrich.txt
│   ├── time1_leiden_res_0.50.csv_3_plot1.pdf
│   ├── time1_leiden_res_0.50.csv_4_enrich.txt
│   ├── time1_leiden_res_0.50.csv_4_plot1.pdf
│   ├── time1_leiden_res_0.50.csv_5_enrich.txt
│   ├── time1_leiden_res_0.50.csv_5_plot1.pdf
│   ├── time1_leiden_res_0.50.csv_6_enrich.txt
│   ├── time1_leiden_res_0.50.csv_6_plot1.pdf
│   ├── time1_leiden_res_0.50.csv_7_enrich.txt
│   ├── time1_leiden_res_0.50.csv_7_plot1.pdf
│   ├── time1_leiden_res_0.50.csv_8_enrich.txt
│   ├── time1_leiden_res_0.50.csv_8_plot1.pdf
│   ├── time2_leiden_res_0.50.csv_0_enrich.txt
│   ├── time2_leiden_res_0.50.csv_0_plot1.pdf
│   ├── time2_leiden_res_0.50.csv_1_enrich.txt
│   ├── time2_leiden_res_0.50.csv_1_plot1.pdf
│   ├── time2_leiden_res_0.50.csv_2_enrich.txt
│   ├── time2_leiden_res_0.50.csv_2_plot1.pdf
│   ├── time2_leiden_res_0.50.csv_3_enrich.txt
│   ├── time2_leiden_res_0.50.csv_3_plot1.pdf
│   ├── time2_leiden_res_0.50.csv_4_enrich.txt
│   ├── time2_leiden_res_0.50.csv_4_plot1.pdf
│   ├── time2_leiden_res_0.50.csv_5_enrich.txt
│   ├── time2_leiden_res_0.50.csv_5_plot1.pdf
│   ├── time2_leiden_res_0.50.csv_6_enrich.txt
│   ├── time2_leiden_res_0.50.csv_6_plot1.pdf
│   ├── time2_leiden_res_0.50.csv_7_enrich.txt
│   ├── time2_leiden_res_0.50.csv_7_plot1.pdf
│   ├── time3_leiden_res_0.50.csv_0_enrich.txt
│   ├── time3_leiden_res_0.50.csv_0_plot1.pdf
│   ├── time3_leiden_res_0.50.csv_1_enrich.txt
│   ├── time3_leiden_res_0.50.csv_1_plot1.pdf
│   ├── time3_leiden_res_0.50.csv_2_enrich.txt
│   ├── time3_leiden_res_0.50.csv_2_plot1.pdf
│   ├── time3_leiden_res_0.50.csv_3_enrich.txt
│   ├── time3_leiden_res_0.50.csv_3_plot1.pdf
│   ├── time3_leiden_res_0.50.csv_4_enrich.txt
│   ├── time3_leiden_res_0.50.csv_4_plot1.pdf
│   ├── time3_leiden_res_0.50.csv_5_enrich.txt
│   ├── time3_leiden_res_0.50.csv_5_plot1.pdf
│   ├── time3_leiden_res_0.50.csv_6_enrich.txt
│   ├── time3_leiden_res_0.50.csv_6_plot1.pdf
│   ├── time3_leiden_res_0.50.csv_7_enrich.txt
│   └── time3_leiden_res_0.50.csv_7_plot1.pdf
├── execution-report.html
├── execution-timeline.html
├── integrate
│   ├── example3_harmony_integrated_dealplus.h5ad
│   ├── example3_harmony_integrated_plus.pdf
│   ├── example3_rliger.INMF_integrated_dealplus.rds
│   ├── example3_rliger.INMF_integrated_plus.pdf
│   ├── example3_SCTransform.CCA_integrated_dealplus.rds
│   ├── example3_SCTransform.CCA_integrated_plus.pdf
│   ├── example3_SCTransform.harmony_integrated_dealplus.rds
│   ├── example3_SCTransform.harmony_integrated_plus.pdf
│   ├── example3_unintegration_integrated_dealplus.h5ad
│   ├── example3_unintegration_integrated_plus.pdf
│   ├── otherexample3_harmony_integrated_plus.pdf
│   ├── otherexample3_unintegration_integrated_plus.pdf
│   └── scib
│       ├── example3_harmony_integrated.h5ad
│       ├── example3_scIB.csv
│       ├── example3_scIB.h5ad
│       └── example3_scIB.pdf
├── metaneighbor
│   ├── example3_metaNeighbor.csv
│   ├── example3_metaNeighbor.pdf
│   ├── example3_metaNeighbor_tophits.csv
│   └── example3_sanky.html
├── pipeline-dag.html
├── qc
│   ├── time1
│   │   ├── figures
│   │   │   ├── dotplot_leiden_res_0.50_marker.pdf
│   │   │   ├── pca_potentially_undesired_features.pdf
│   │   │   ├── umap_batch.pdf
│   │   │   ├── umap_leiden_clus.pdf
│   │   │   └── umap_quality.pdf
│   │   ├── markers.csv
│   │   │   └── time1_leiden_res_0.50.csv
│   │   ├── qc.pdf
│   │   ├── summary.txt
│   │   └── time1.h5ad
│   ├── time2
│   │   ├── figures
│   │   │   ├── dotplot_leiden_res_0.50_marker.pdf
│   │   │   ├── pca_potentially_undesired_features.pdf
│   │   │   ├── umap_batch.pdf
│   │   │   ├── umap_leiden_clus.pdf
│   │   │   └── umap_quality.pdf
│   │   ├── markers.csv
│   │   │   └── time2_leiden_res_0.50.csv
│   │   ├── qc.pdf
│   │   ├── summary.txt
│   │   └── time2.h5ad
│   └── time3
│       ├── figures
│       │   ├── dotplot_leiden_res_0.50_marker.pdf
│       │   ├── pca_potentially_undesired_features.pdf
│       │   ├── umap_batch.pdf
│       │   ├── umap_leiden_clus.pdf
│       │   └── umap_quality.pdf
│       ├── markers.csv
│       │   └── time3_leiden_res_0.50.csv
│       ├── qc.pdf
│       ├── summary.txt
│       └── time3.h5ad
├── time1.hr.rds
├── time2.hr.rds
└── time3.hr.rds

24 directories, 140 files
```

</details>