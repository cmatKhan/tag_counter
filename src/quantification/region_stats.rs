use statrs::distribution::{Discrete, DiscreteCDF, Poisson};

/// Errors that can occur when computing the Poisson p-value.
#[derive(Debug)]
pub enum PoissonPvalError {
    DivisionByZero,
    InvalidPoissonParameter,
}

/// Errors that can occur when computing enrichment.
#[derive(Debug)]
pub enum EnrichmentError {
    DivisionByZero,
    InvalidComputation,
}

/// Represents total tag counts for control and treatment conditions.
pub struct RegionStats {
    // the total control and treatment tags are reasonably stored as u32. They do
    // get cast to f64 for the enrichment and pvalue calculations. This is also why
    // pseudocount is stored as f64, though it does not need that level of precision.
    total_control_tags: u32,
    total_treatment_tags: u32,
    pub pseudocount: f64,
}

impl RegionStats {
    pub fn new(total_control_tags: u32, total_treatment_tags: u32, pseudocount: f64) -> Self {
        RegionStats {
            total_control_tags,
            total_treatment_tags,
            pseudocount,
        }
    }

    /// Provides read-only access to `total_control_tags`.
    pub fn total_control_tags(&self) -> u32 {
        self.total_control_tags
    }

    /// Provides read-only access to `total_treatment_tags`.
    pub fn total_treatment_tags(&self) -> u32 {
        self.total_treatment_tags
    }

    /// Computes the Poisson p-value for a given region.
    pub fn poisson_pval(
        &self,
        control_tags: u32,
        treatment_tags: u32,
    ) -> Result<f64, PoissonPvalError> {
        if self.total_control_tags == 0 || self.total_treatment_tags == 0 {
            return Err(PoissonPvalError::DivisionByZero);
        }

        // cast to f64 to avoid overflow. Additionally, the poisson function requires
        // f64 as input for the poisson parameter mu. treatment_tags is set to a u64
        // and the result is f64 in precision.
        let hop_ratio = self.total_treatment_tags as f64 / self.total_control_tags as f64;
        let mu = (control_tags as f64 + self.pseudocount) * hop_ratio;

        if mu <= 0.0 {
            return Err(PoissonPvalError::InvalidPoissonParameter);
        }

        let poisson = Poisson::new(mu).map_err(|_| PoissonPvalError::InvalidPoissonParameter)?;

        let x = treatment_tags as u64;
        let p_val = (1.0 - poisson.cdf(x)) + poisson.pmf(x);

        // Check if it's out of bounds before clamping
        if p_val > 1.0 {
            eprintln!(
                "⚠️  Warning: computed Poisson p-value exceeds 1.0 ({:.17}) for region (control_tags={}, treatment_tags={}, mu={:.6})",
                p_val, control_tags, treatment_tags, mu
            );
        } else if p_val < 0.0 {
            eprintln!(
                "⚠️  Warning: computed Poisson p-value is negative ({:.17}) for region (control_tags={}, treatment_tags={}, mu={:.6})",
                p_val, control_tags, treatment_tags, mu
            );
        }

        // Note that this clamps the p-value to the range [0, 1]. This is not strictly necessary
        Ok(p_val.clamp(0.0, 1.0))
    }

    /// Computes the enrichment score for a given region.
    pub fn enrichment(
        &self,
        control_tags: u32,
        treatment_tags: u32,
    ) -> Result<f64, EnrichmentError> {
        // see poisson_pval() for a discussion on precision
        if self.total_control_tags == 0 || self.total_treatment_tags == 0 {
            return Err(EnrichmentError::DivisionByZero);
        }

        let numerator = treatment_tags as f64 / self.total_treatment_tags as f64;
        let denominator = (control_tags as f64 + self.pseudocount) / self.total_control_tags as f64;

        let enrichment = numerator / denominator;

        if enrichment.is_nan() || enrichment.is_infinite() || enrichment < 0.0 {
            return Err(EnrichmentError::InvalidComputation);
        }

        Ok(enrichment)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_enrichment_valid() {
        let tag_counts = RegionStats {
            total_control_tags: 1000,
            total_treatment_tags: 2000,
            pseudocount: 0.1,
        };

        let result = tag_counts.enrichment(50, 100);
        assert!(result.is_ok());
        println!("Enrichment: {:?}", result);
        assert!((result.unwrap() - 0.998004).abs() < 1e-6);
    }

    #[test]
    fn test_enrichment_zero_control_tags() {
        let tag_counts = RegionStats {
            total_control_tags: 0,
            total_treatment_tags: 2000,
            pseudocount: 0.1,
        };

        let result = tag_counts.enrichment(50, 100);
        assert!(matches!(result, Err(EnrichmentError::DivisionByZero)));
    }

    #[test]
    fn test_enrichment_nan_or_infinite() {
        let tag_counts = RegionStats {
            total_control_tags: 1000,
            total_treatment_tags: 2000,
            pseudocount: 0.0,
        };

        let result = tag_counts.enrichment(0, 0); // This may cause division by zero
        assert!(matches!(result, Err(EnrichmentError::InvalidComputation)));
    }

    #[test]
    fn test_poisson_pval_valid_1() {
        let tag_counts = RegionStats {
            total_control_tags: 100,
            total_treatment_tags: 10,
            pseudocount: 0.1,
        };

        let result = tag_counts.poisson_pval(5, 2);
        assert!(result.is_ok());
        let p_val = result.unwrap();
        println!("Test 1 - Poisson p-value: {:.6}", p_val);
        assert!((p_val - 0.093252).abs() < 1e-6);
    }

    #[test]
    fn test_poisson_pval_valid_2() {
        let tag_counts = RegionStats {
            total_control_tags: 200,
            total_treatment_tags: 20,
            pseudocount: 0.1,
        };

        let result = tag_counts.poisson_pval(10, 4);
        assert!(result.is_ok());
        let p_val = result.unwrap();
        println!("Test 2 - Poisson p-value: {:.6}", p_val);
        assert!((p_val - 0.019607).abs() < 1e-6);
    }

    #[test]
    fn test_poisson_pval_zero_control_tags() {
        let tag_counts = RegionStats {
            total_control_tags: 0,
            total_treatment_tags: 10,
            pseudocount: 0.1,
        };

        let result = tag_counts.poisson_pval(5, 2);
        assert!(matches!(result, Err(PoissonPvalError::DivisionByZero)));
    }

    #[test]
    fn test_poisson_pval_zero_treatment_tags() {
        let tag_counts = RegionStats {
            total_control_tags: 100,
            total_treatment_tags: 0,
            pseudocount: 0.1,
        };

        let result = tag_counts.poisson_pval(5, 2);
        assert!(matches!(result, Err(PoissonPvalError::DivisionByZero)));
    }

    #[test]
    fn test_poisson_pval_invalid_mu() {
        let tag_counts = RegionStats {
            total_control_tags: 100,
            total_treatment_tags: 10,
            pseudocount: 0.0,
        };

        let result = tag_counts.poisson_pval(0, 2);
        assert!(matches!(
            result,
            Err(PoissonPvalError::InvalidPoissonParameter)
        ));
    }
}
