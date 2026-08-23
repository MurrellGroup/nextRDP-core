# Pinned parity oracle

The target used by the existing locally compiled reference runs is:

```text
RDP5CL.exe  7917d3675f3976b331e8ba8a55cfcb28641bce93844faeb35dc5c836aa5f70b6
DNA5.dll    3e9e4b0b3c4a94eab4a7e8654a0803bc27de138fb0ad58e3481a81a1f199ee70
```

Canonical local location:

```text
../sandbox/source-build/rdp-dll-smoke/rebuilt/
```

The vendored DNA5 source is byte-identical to `originalSource/dna5DLLSource.zip` and to the source used for the locally rebuilt DLL.

RDP-only command-line settings are inherited from the validated local runtime: circular sequences, p=0.05, Bonferroni correction, RDP window 30, RDP method enabled, all other methods disabled, polishing disabled, one worker.

