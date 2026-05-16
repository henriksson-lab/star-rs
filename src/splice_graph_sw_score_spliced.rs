#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SpliceGraph::swScoreSpliced` at STAR/source/SpliceGraph_swScoreSpliced.cpp:8. Args: readSeq: char, readLen: uint32, superTr: SuperTranscript, cigar: vector<array<uint32,2>>"]
pub fn splicegraph_swscorespliced_l8_splicegraph_swscorespliced(
    splice_graph: &mut crate::splice_graph::SpliceGraph,
    read_seq: &[u8],
    read_len: u32,
    super_tr: &crate::super_transcriptome::SuperTranscript,
    cigar: &mut Vec<[u32; 2]>,
) -> i32 {
    let super_tr_len = super_tr.length as usize;
    let read_len_usize = read_len as usize;
    let mut score_max_global = 0_i32;
    let mut i_donor = 0_i32;
    let mut i_acceptor = 0_i32;
    let sj_yes = !super_tr.sj_c.is_empty();

    splice_graph
        .direction_matrix
        .resize(read_len_usize * super_tr_len, 0);
    splice_graph.direction_matrix.fill(0);
    if splice_graph.score_two_columns[0].len() < read_len_usize + 1 {
        splice_graph.score_two_columns[0].resize(read_len_usize + 1, 0);
    }
    if splice_graph.score_two_columns[1].len() < read_len_usize + 1 {
        splice_graph.score_two_columns[1].resize(read_len_usize + 1, 0);
    }

    let mut score_column = vec![0_i32; read_len_usize + 1];
    let mut i_two_column = 0_usize;
    splice_graph.score_two_columns[i_two_column][..=read_len_usize].fill(0);

    for col in 0..super_tr_len {
        let score_column_prev = score_column;
        i_two_column = if i_two_column == 0 { 1 } else { 0 };
        score_column = splice_graph.score_two_columns[i_two_column][..=read_len_usize].to_vec();

        let mut sj_n = 0_u32;
        let mut donor_column: Option<usize> = None;
        if sj_yes {
            while i_acceptor < super_tr.sj_c.len() as i32
                && col as u32 == super_tr.sj_c[i_acceptor as usize][1]
            {
                if splice_graph.sj_dindex.len() <= sj_n as usize {
                    splice_graph.sj_dindex.resize(sj_n as usize + 1, 0);
                }
                splice_graph.sj_dindex[sj_n as usize] = super_tr.sj_c[i_acceptor as usize][2];
                sj_n += 1;
                i_acceptor += 1;
            }

            if i_donor < super_tr.sj_donor.len() as i32
                && col as u32 == super_tr.sj_donor[i_donor as usize]
            {
                donor_column = Some(i_donor as usize);
                if splice_graph.scoring_matrix.len() <= i_donor as usize {
                    splice_graph
                        .scoring_matrix
                        .resize(i_donor as usize + 1, vec![0; read_len_usize + 1]);
                }
                if splice_graph.scoring_matrix[i_donor as usize].len() < read_len_usize + 1 {
                    splice_graph.scoring_matrix[i_donor as usize].resize(read_len_usize + 1, 0);
                }
                score_column =
                    splice_graph.scoring_matrix[i_donor as usize][..=read_len_usize].to_vec();
                i_donor += 1;
            }
        }

        score_column[0] = 0;
        for row in 1..=read_len_usize {
            let mut dir_ind_max = 0_u8;
            let mut score_max = 0_i32;

            let score1 = score_column[row - 1] + splice_graph.gap_penalty as i32;
            if score1 > score_max {
                dir_ind_max = 1;
                score_max = score1;
            }

            let score1 = score_column_prev[row] + splice_graph.gap_penalty as i32;
            if score1 > score_max {
                dir_ind_max = 2;
                score_max = score1;
            }

            let score1 = score_column_prev[row - 1]
                + if read_seq[row - 1] == super_tr.seq[col] {
                    splice_graph.match_score as i32
                } else {
                    splice_graph.mismatch_penalty as i32
                };
            if score1 > score_max {
                dir_ind_max = 3;
                score_max = score1;
            }

            for ii in 0..sj_n as usize {
                let sj_dindex = splice_graph.sj_dindex[ii] as usize;
                let score1 =
                    splice_graph.scoring_matrix[sj_dindex][row] + splice_graph.gap_penalty as i32;
                if score1 > score_max {
                    dir_ind_max = (4 + ii * 2) as u8;
                    score_max = score1;
                }

                let score1 = splice_graph.scoring_matrix[sj_dindex][row - 1]
                    + if read_seq[row - 1] == super_tr.seq[col] {
                        splice_graph.match_score as i32
                    } else {
                        splice_graph.mismatch_penalty as i32
                    };
                if score1 > score_max {
                    dir_ind_max = (5 + ii * 2) as u8;
                    score_max = score1;
                }
            }

            splice_graph.direction_matrix[(row - 1) + col * read_len_usize] = dir_ind_max;
            score_column[row] = score_max;
            if score_max_global < score_max {
                score_max_global = score_max;
                splice_graph.align_info.a_end[0] = row as u32;
                splice_graph.align_info.a_end[1] = col as u32;
            }
        }

        splice_graph.score_two_columns[i_two_column][..=read_len_usize]
            .copy_from_slice(&score_column[..=read_len_usize]);
        if let Some(donor) = donor_column {
            splice_graph.scoring_matrix[donor][..=read_len_usize]
                .copy_from_slice(&score_column[..=read_len_usize]);
        }
    }
    splice_graph.align_info.a_end[0] = splice_graph.align_info.a_end[0].wrapping_sub(1);

    cigar.clear();
    cigar.reserve(read_len_usize);
    let mut row = splice_graph.align_info.a_end[0] as i32;
    let mut col = splice_graph.align_info.a_end[1] as i32;

    splice_graph.align_info.n_map = 0;
    splice_graph.align_info.n_mm = 0;
    splice_graph.align_info.n_i = 0;
    splice_graph.align_info.n_d = 0;
    splice_graph.align_info.n_sj = 0;
    let mut i_acceptor = super_tr.sj_c.len() as i32 - 1;

    if row != read_len as i32 - 1 {
        cigar.push([BAM_CIGAR_S, read_len - 1 - row as u32]);
    }

    let mut cigar_op;
    let mut cigar_len = 0_u32;
    let mut cigar_op_prev = u32::MAX;
    let mut sj_gap = 0_u32;
    while col >= 0 && row >= 0 {
        let dir1 =
            splice_graph.direction_matrix[row as usize + col as usize * read_len_usize] as u32;

        if dir1 == 0 {
            break;
        }

        match dir1 {
            1 => {
                row -= 1;
                splice_graph.align_info.n_i += 1;
                cigar_op = BAM_CIGAR_I;
            }
            2 => {
                col -= 1;
                splice_graph.align_info.n_d += 1;
                cigar_op = BAM_CIGAR_D;
            }
            3 => {
                splice_graph.align_info.n_map += 1;
                splice_graph.align_info.n_mm +=
                    (read_seq[row as usize] != super_tr.seq[col as usize]) as u32;
                cigar_op = BAM_CIGAR_M;
                row -= 1;
                col -= 1;
            }
            _ => {
                splice_graph.align_info.n_sj += 1;
                while i_acceptor + 1 != 0 && col as u32 <= super_tr.sj_c[i_acceptor as usize][1] {
                    i_acceptor -= 1;
                }

                if (dir1 - 4) % 2 == 1 {
                    splice_graph.align_info.n_map += 1;
                    splice_graph.align_info.n_mm +=
                        (read_seq[row as usize] != super_tr.seq[col as usize]) as u32;
                    row -= 1;
                    cigar_op = BAM_CIGAR_M;
                } else {
                    splice_graph.align_info.n_d += 1;
                    cigar_op = BAM_CIGAR_D;
                }

                sj_gap = col as u32;
                let sj_index = (i_acceptor + 1 + ((dir1 - 4) / 2) as i32) as usize;
                col = super_tr.sj_donor[super_tr.sj_c[sj_index][2] as usize] as i32;
                sj_gap = sj_gap - col as u32 - 1;
            }
        }

        if cigar_op != cigar_op_prev {
            if cigar_len > 0 {
                cigar.push([cigar_op_prev, cigar_len]);
            }
            cigar_len = 0;
            cigar_op_prev = cigar_op;
        }
        cigar_len += 1;
        if sj_gap > 0 {
            cigar.push([cigar_op, cigar_len]);
            cigar.push([BAM_CIGAR_N, sj_gap]);
            cigar_len = 0;
            cigar_op_prev = u32::MAX;
            sj_gap = 0;
        }
    }

    if cigar_len > 0 {
        cigar.push([cigar_op_prev, cigar_len]);
    }

    row += 1;
    col += 1;
    splice_graph.align_info.a_start[0] = row as u32;
    splice_graph.align_info.a_start[1] = col as u32;
    if row > 0 {
        cigar.push([BAM_CIGAR_S, row as u32]);
    }

    cigar.reverse();
    score_max_global
}
