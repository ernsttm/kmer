use std::collections::HashMap;
use std::time::Instant;

struct MatchState {
    index: (usize, usize),
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
            output.push((a_index.clone() + 1, b_index.clone() + 1));
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

fn find_p_alignment(kmer_matches: &Vec<(usize, usize)>) {
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

    let mut trace_p = p_alignment.first().unwrap();
    while backtrace_map.contains_key(&(trace_p.match_index as usize)) {
        println!("{:?}", trace_p.index);
        let index_of_interest =
            kmer_matches[backtrace_map[&(trace_p.match_index as usize)] as usize];
        trace_p = p_alignment
            .iter()
            .find(|&x| x.index == index_of_interest)
            .unwrap();
    }

    println!("{:?}", trace_p.index);
}

pub fn find_alignment(first_sequence: &str, second_sequence: &str, kmer_size: usize) {
    let now = Instant::now();
    let first_kmers = parse_kmers(&first_sequence, kmer_size);
    let second_kmers = parse_kmers(&second_sequence, kmer_size);
    let matches = calculate_and_sort_matches(&first_kmers, &second_kmers);
    println!("4 step time : {}", now.elapsed().as_micros());

    find_p_alignment(&matches);
    println!("5 step time : {}", now.elapsed().as_micros());
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
