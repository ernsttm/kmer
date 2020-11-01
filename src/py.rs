use pyo3::prelude::*;

use crate::kmer::find_alignment;

#[pymodule]
fn libkmer(_py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "find_alignment")]
    fn find_alignment_py(_py: Python, first_sequence: &str, second_sequence: &str, kmer_size: usize) ->  f64 {
        find_alignment(first_sequence, second_sequence, kmer_size)
    }

    Ok(())
}