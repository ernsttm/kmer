use std::collections::HashMap;
use std::cmp::max;

/// A 2d reference into the query and subject sequences where a kmer match occurred.
struct KmerMatchIndex {
    /// The index within the query sequence the match occurred at.
    query_index: usize,

    /// The index within the subject sequence the match occurred at.
    subject_index: usize,
}

/// A kmer match found between two sequences.
struct MatchState {
    /// The indices within the query and subject sequences the match was found at.
    index: (usize, usize),

    /// The
    match_index: i32,
    score: i64,
}

fn parse_kmers(genome: &str, kmer_size: usize) -> HashMap<&str, Vec<usize>> {
    let mut kmers: HashMap<&str, Vec<usize>> = HashMap::new();
    for index in 0..genome.len() - kmer_size + 1 {
        let kmer = &genome[index..index + kmer_size];

        match kmers.get_mut(&kmer) {
            Some(indices) => indices.push(index),
            None => {
                let mut indices = Vec::new();
                indices.push(index);
                kmers.insert(kmer, indices);
            }
        }
    }

    kmers
}

fn cross_product(a: &Vec<usize>, b: &Vec<usize>, output: &mut Vec<(usize, usize)>) {
    for a_index in a {
        for b_index in b {
            output.push((a_index.clone(), b_index.clone()));
        }
    }
}

fn calculate_and_sort_matches(
    first_kmers: &HashMap<&str, Vec<usize>>,
    second_kmers: &HashMap<&str, Vec<usize>>,
) -> Vec<(usize, usize)> {
    let mut matches = Vec::new();

    for (kmer, first_indices) in first_kmers {
        match second_kmers.get(kmer) {
            Some(second_indices) => {
                cross_product(first_indices, second_indices, &mut matches);
            }
            None => { /* No op */ }
        }
    }

    matches.sort();
    matches
}

fn distance_function(i: usize, j: usize) -> i64 {
    // This is a version of the Kronecker-Delta function, with a c value of .75
    if i == j {
        0
    } else {
        1
    }
}

fn find_p_alignment(kmer_matches: &Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    let mut backtrace_map: HashMap<usize, i32> = HashMap::new();
    let mut p_alignment: Vec<MatchState> = Vec::new();

    p_alignment.push(MatchState {
        index: (0, 0),
        match_index: -1,
        score: 0,
    });

    for (match_index, (i, j)) in kmer_matches.iter().enumerate() {
        let mut kmer_match_state = MatchState {
            index: (i.clone(), j.clone()),
            match_index: match_index as i32,
            score: 1,
        };

        loop {
            let mut updated_score = false;
            for p in &mut p_alignment {
                if i <= &p.index.0 || j <= &p.index.1 {
                    continue;
                }

                let i_val = i - p.index.0;
                let j_val = j - p.index.1;
                if p.score - distance_function(i_val, j_val) + 1 > kmer_match_state.score {
                    kmer_match_state.score =
                        p.score - distance_function(i_val, j_val) + kmer_match_state.score;
                    backtrace_map.insert(match_index, p.match_index);
                    updated_score = true;
                    break;
                }
            }

            if !updated_score {
                match p_alignment
                    .binary_search_by(|p| kmer_match_state.score.partial_cmp(&p.score).unwrap())
                {
                    Ok(pos) => p_alignment.insert(pos, kmer_match_state),
                    Err(pos) => p_alignment.insert(pos, kmer_match_state),
                }
                break;
            }
        }
    }

    let mut pairs = Vec::new();
    let mut trace_p = p_alignment.first().unwrap();
    while backtrace_map.contains_key(&(trace_p.match_index as usize)) {
        pairs.push(trace_p.index);
        let index_of_interest =
            kmer_matches[backtrace_map[&(trace_p.match_index as usize)] as usize];
        trace_p = p_alignment
            .iter()
            .find(|&x| x.index == index_of_interest)
            .unwrap();
    }

    pairs.push(trace_p.index);
    pairs
}

fn calculate_distances(first_sequence_len: usize, second_sequence_len: usize, match_pairs: &Vec<(usize, usize)>) -> f64 {
    let mut total_index: usize = 0;
    let mut first_index: usize = 0;
    let mut second_index: usize = 0;
    let mut first_gap: usize = 0;
    let mut second_gap: usize = 0;
    for (query_pair_index, subject_pair_index) in match_pairs {
        while total_index < max(query_pair_index.clone(), subject_pair_index.clone()) {
            if first_index < query_pair_index.clone() {
                first_index = first_index + 1;
            } else {
                first_gap = first_gap + 1;
            }

            if second_index < subject_pair_index.clone() {
                second_index = second_index + 1;
            } else {
                second_gap = second_gap + 1;
            }

            total_index = total_index + 1;
        }
    }

    let alignment_length = max(first_sequence_len + first_gap, second_sequence_len + second_gap) as f64;
    println!("Pair length : {} | alignment length : {}", match_pairs.len(), alignment_length);
    1.0 - (match_pairs.len() as f64 / alignment_length)
}

pub fn ktuple_distance(first_sequence: &str, second_sequence: &str, kmer_size: usize) -> f64{
    let first_kmers = parse_kmers(&first_sequence, kmer_size);
    let second_kmers = parse_kmers(&second_sequence, kmer_size);
    let matches = calculate_and_sort_matches(&first_kmers, &second_kmers);

    let pairs = find_p_alignment(&matches);
    calculate_distances(first_sequence.len(), second_sequence.len(), &pairs)
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
