DROP VIEW IF EXISTS unique_candidates;

CREATE VIEW unique_candidates AS
SELECT
    peptide,
    allele,
    MIN(mhcflurry_rank) AS best_mhcflurry_rank,
    MAX(mhcflurry_presentation_score) AS best_presentation_score,
    MIN(consensus_rank) AS best_consensus_rank,
    MAX(high_confidence) AS any_high_confidence,
    COUNT(*) AS isoform_count
FROM final_candidates
GROUP BY peptide, allele;
