use std::cmp::{min, max};

const MATCH_SCORE: i32 = 2;
const MISMATCH_SCORE: i32 = -3;
const GAP_OPEN: i32 = -11;
const GAP_EXTEND: i32 = -1;

pub fn needleman_wunsh_distance(first_sequence: &str, second_sequence: &str) -> f64 {
    let m = first_sequence.len() + 1;
    let n = second_sequence.len() + 1;
    let mut scoring_matrix: Vec<Vec<i32>> = Vec::new();
    let mut traceback_matrix: Vec<Vec<char>> = Vec::new();

    let mut index = 0;
    scoring_matrix.resize(m, Vec::new());
    for row in &mut scoring_matrix {
        row.resize(n, 0);
        row[0] = -7 * index;

        if index == 0 {
            for row_index in 0..row.len() {
                row[row_index] = -7 * row_index as i32;
            }
        }
    }

    let mut initial = true;
    traceback_matrix.resize(m, Vec::new());
    for row in &mut traceback_matrix{
        row.resize(n, ' ');

        row[0] = 'V';
        if initial {
            for row_index in 0..row.len() {
                row[row_index] = 'H';
            }
        }
    }

    for i in 1..m {
        for j in 1..n {
            let match_score = if first_sequence == second_sequence { MATCH_SCORE } else { MISMATCH_SCORE };

            let diag = scoring_matrix[i - 1][j - 1] + match_score;
            let horizontal = scoring_matrix[i][j - 1] - 7;
            let vertical = scoring_matrix[i - 1][j] - 7;
            let minimum_score = min(diag, min(horizontal, vertical));

            scoring_matrix[i][j] = minimum_score;
            if minimum_score == diag {
                traceback_matrix[i][j] = 'D';
            } else if minimum_score == vertical {
                traceback_matrix[i][j] = 'V';
            } else {
                traceback_matrix[i][j] = 'H';
            }
        }
    }

    let mut dist = 0;
    let mut steps = 0;
    let mut i = m - 1;
    let mut j = n - 1;
    while i != 0 && j != 0 {
        steps = steps + 1;

        if traceback_matrix[i][j] == 'D' {
            if first_sequence.as_bytes()[i - 1] != first_sequence.as_bytes()[j - 1] {
                dist = dist + 1;
            }

            i = i - 1;
            j = j - 1;
        } else if traceback_matrix[i][j] == 'H' {
            dist = dist + 1;
            i = i - 1;
        } else {
            dist = dist + 1;
            j = j - 1;
        }
    }

    dist as f64 / max(m, n) as f64
}