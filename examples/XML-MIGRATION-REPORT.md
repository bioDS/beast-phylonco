# BEAST3 XML Migration Report — `examples/`

Migration of the top-level `examples/` example-analysis XMLs from BEAST2 to BEAST3,
using `beast3-migration/skills/xml-migration/convert_b2_to_b3.py`. Scope: `examples/`
only (not `phylonco-lphy/examples/`, not FxTemplates).

## Directory layout after this pass

```
examples/
├── legacy/                      ← original BEAST2-format sources (git mv'd here)
│   ├── test_binary_error.xml
│   ├── test_GT16.xml
│   ├── test_GT16_error.xml
│   ├── test_GT16_error_beagle.xml
│   ├── test_SiFit3.xml
│   └── L86/
│       ├── L86_rskycoal_beagle.xml
│       ├── L86_rskycoal_beagle_outgroup.xml
│       ├── L86_rskycoal_em_beagle.xml
│       └── L86_rskycoal_em_beagle_outgroup.xml
├── readcounts/                  ← read-count model XMLs (moved, not yet migrated)
│   ├── test_MA_read_count_model.xml
│   └── test_read_count_model.xml
├── test_binary_error_b3.xml     ← migrated BEAST3 output
├── test_GT16_b3.xml
├── test_GT16_error_b3.xml
├── test_GT16_error_beagle_b3.xml
├── test_SiFit3_b3.xml
├── L86/
│   ├── L86_rskycoal_beagle_b3.xml
│   ├── L86_rskycoal_beagle_outgroup_b3.xml
│   ├── L86_rskycoal_em_beagle_b3.xml
│   ├── L86_rskycoal_em_beagle_outgroup_b3.xml
│   └── reports/                 ← per-file converter reports (auto-generated)
├── reports/
└── data/                        ← unchanged (.nex/.fasta, not XML)
```

## Converted files (9)

| File | Renames | Warnings | TODOs |
|---|---|---|---|
| `test_binary_error.xml` | 17 | 12 | 2 |
| `test_GT16.xml` | 16 | 6 | 1 |
| `test_GT16_error.xml` | 23 | 9 | 0 |
| `test_GT16_error_beagle.xml` | 23 | 9 | 0 |
| `test_SiFit3.xml` | 12 | 8 | 1 |
| `L86/L86_rskycoal_beagle.xml` | 63 | 12 | 1 |
| `L86/L86_rskycoal_beagle_outgroup.xml` | 64 | 12 | 1 |
| `L86/L86_rskycoal_em_beagle.xml` | 70 | 16 | 1 |
| `L86/L86_rskycoal_em_beagle_outgroup.xml` | 71 | 16 | 1 |

Full per-file rename/warning/todo detail is in `tmp/b3migration/log/*-examples.md`
(local, not committed) and the per-file `reports/*.md` next to each converted XML.

**TODOs are expected, not gaps**: every `[todo]` is a class with no `.spec.` twin
(`CompoundDistribution`, `Binary` datatype) — both confirmed non-deprecated core
classes in `beast3`, correctly left unrenamed per the migration skill's U4 rule.

## Validation (`mvn -pl phylonco-beast exec:exec -Dbeast.args="-validate <file>"`)

| Result | Files |
|---|---|
| ✅ PASS | `test_binary_error_b3.xml`, `test_GT16_b3.xml`, `test_GT16_error_b3.xml`, `test_GT16_error_beagle_b3.xml`, `test_SiFit3_b3.xml` |
| ❌ FAIL | all 4 `L86/*_b3.xml` |

### Open issue: `L86` failures — converter gap, not a one-off

All four fail identically:

```
Error 130 parsing the xml input file
Input 102b: type mismatch for input groupSizes
```

**Root cause**: the original `groupSizes` (`spec="parameter.IntegerParameter"`, feeding
`BayesianSkyline`) was converted to `IntVectorParam` — a plain int vector. But BEAST3's
`BayesianSkyline.groupSizeParamInput` requires `IntSimplex<? extends PositiveInt>`: the
group sizes must sum to a declared total (the tree's interval count), not just be a free
vector. This is the int/`groupSizes` analogue of the `RealParameter` → `SimplexParam`
rule already documented in `parameters.md` R3 for `frequencies`-style inputs — the
converter has no equivalent disambiguation rule for `groupSizes`-shaped int inputs, so it
falls through to the generic vector mapping.

**Not hand-patched**: `IntSimplexParam` needs an explicit `sum` (normally taxa count − 1)
that isn't safely inferable from a quick read of the XML — guessing wrong trades a parse
error for a silent validation-logic bug. Needs either:
1. A fix to `convert_b2_to_b3.py`'s disambiguation rules (recognize `groupSizes` inputs
   to `BayesianSkyline`/similar coalescent models and emit `IntSimplexParam` with the
   correct `sum`), or
2. Manual correction of the 4 `L86/*_b3.xml` files once the correct `sum` is known.
