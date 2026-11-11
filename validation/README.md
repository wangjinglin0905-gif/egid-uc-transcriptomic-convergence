# Release checks, 2026-11-11

The portable public copy passed 14 base-R checks recorded in derived_statistics_checks.csv. All six R scripts parsed successfully. Frozen result copies and the matched-null protocol/configuration were hash-checked against the source project; no statistical algorithms were changed.

Text-based release scanning found no local absolute paths or common credential patterns. This bounded scan is not a general guarantee that arbitrary sensitive content can be detected.

The checks were performed on the development computer (Windows, R4.6.1), not an independent machine. No full expression-model or fresh permutation rerun was performed for this public release. Figure binaries are the previously checked 600-dpi PNG/TIFF exports.

MANIFEST_SHA256.csv records the supplied files. The historical freeze hashes apply to the preserved protocol/configuration; public script adaptations have their own provenance hashes.
