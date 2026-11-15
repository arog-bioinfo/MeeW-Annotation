# DIAMOND index fixtures

These indexes were built from the adjacent synthetic FASTA files with the
pinned runtime executable:

```text
diamond v2.2.6.186
diamond makedb --in small.faa --db small
diamond makedb --in typical.faa --db typical
```

The binary indexes are stored as base64 so the test fixtures remain portable
through text-only patch tooling. Decoded checksums and sizes are:

```text
ced6b1eea34766b32e1eea20ade2f1868f31fed8403d800171d46f595a68e848  small.dmnd    139 bytes
1073bc45a66e51ee799b0abb8e0a95cfcbc19a1d7dfe7708bb1f5480af2b2f2c  typical.dmnd  444 bytes
```
