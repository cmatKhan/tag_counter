use statrs::distribution::{DiscreteCDF, Poisson};

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
    total_control_tags: u32,
    total_treatment_tags: u32,
    pub pseudocount: f32,
}

impl RegionStats {
    pub fn new(total_control_tags: u32, total_treatment_tags: u32, pseudocount: f32) -> Self {
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

        let hop_ratio = self.total_treatment_tags as f32 / self.total_control_tags as f32;
        let mu = (control_tags as f32 + self.pseudocount) * hop_ratio;

        if mu <= 0.0 {
            return Err(PoissonPvalError::InvalidPoissonParameter);
        }

        let poisson =
            Poisson::new(mu.into()).map_err(|_| PoissonPvalError::InvalidPoissonParameter)?;
        let p_val = 1.0 - poisson.cdf(treatment_tags as u64);

        Ok(p_val)
    }

    /// Computes the enrichment score for a given region.
    pub fn enrichment(
        &self,
        control_tags: u32,
        treatment_tags: u32,
    ) -> Result<f32, EnrichmentError> {
        if self.total_control_tags == 0 || self.total_treatment_tags == 0 {
            return Err(EnrichmentError::DivisionByZero);
        }

        let numerator = treatment_tags as f32 / self.total_treatment_tags as f32;
        let denominator = (control_tags as f32 + self.pseudocount) / self.total_control_tags as f32;

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
        let result = enrichment(1000, 2000, 50, 100, 0.1);
        assert!(result.is_ok());
        println!("Enrichment: {:?}", result);
        assert!((result.unwrap() - 0.998004) < 1e-6);
    }

    #[test]
    fn test_enrichment_zero_control_tags() {
        let result = enrichment(0, 2000, 50, 100, 0.1);
        assert!(matches!(result, Err(EnrichmentError::DivisionByZero)));
    }

    #[test]
    fn test_enrichment_nan_or_infinite() {
        let result = enrichment(1000, 2000, 0, 0, 0.0); // This may cause division by zero
        assert!(matches!(result, Err(EnrichmentError::InvalidComputation)));
    }
}
