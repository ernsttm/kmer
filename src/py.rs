use pyo3::prelude::*;

use crate::kmer::ktuple_distance;

#[pymodule]
fn libkmer(_py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "ktuple_distance")]
    fn ktuple_distance(_py: Python, first_sequence: &str, second_sequence: &str, kmer_size: usize) -> f64 {
        crate::kmer::ktuple_distance(first_sequence, second_sequence, kmer_size)
    }

    #[pyfn(m, "neeleman_wunch_distance")]
    fn needleman_wunsh_distance(_py: Python, first_sequence: &str, second_sequence: &str) -> f64 {
        crate::needlman_wunsch::needleman_wunsh_distance(first_sequence, second_sequence)
    }

    Ok(())
}