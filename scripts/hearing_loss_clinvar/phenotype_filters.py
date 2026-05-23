from __future__ import annotations

import re


def compile_pattern(mode: str) -> re.Pattern[str]:
    """
    Return a compiled regex for identifying hearing-loss-related records in ClinVar PhenotypeList.

    Modes are intentionally conservative and are used for sensitivity analysis rather than claiming
    any single "correct" phenotype definition.
    """
    mode = (mode or "").strip().lower()

    if mode in ("hearing_loss_only", "hl_only", "default"):
        # Baseline definition: explicit phrase match.
        return re.compile(r"\bhearing\s+loss\b", re.IGNORECASE)

    if mode in ("expanded_synonyms", "expanded", "synonyms"):
        # Expanded text definition (still not ontology-grounded).
        # Intentionally avoids very short acronyms (e.g., "HL") to reduce false positives.
        return re.compile(
            r"(\bhearing\s+loss\b|\bsensorineural\s+hearing\s+loss\b|\bhearing\s+impairment\b|\bdeaf(ness)?\b|\bhard\s*of\s*hearing\b|\bauditory\s+neuropathy\b|\bsnhl\b)",
            re.IGNORECASE,
        )

    raise ValueError(f"Unknown phenotype filter mode: {mode}")

