# Output files

All output files are stored in a folder with the same name as the sampleID in the specified output directory. If no output directory is specified, the files are stored in the current working directory.

### Folder structure

```
sampleID/
├── qc/
│   ├── shortread_qc/       (hybrid only)
│   ├── longread_raw/
│   ├── longread_filterd/
│   ├── assembly_qc/
│   └── qc_sampleID.json
├── assembly/
│   ├── graph_plot/
│   ├── miniasm/
│   ├── spades/
│   ├── unicycler/
│   ├── canu/
│   └── flye/
├── assembly_processed/
│   ├── racon/
│   ├── pilon/
│   └── link/
├── genomes/
├── plasmids/
└── plasmIDent
```




