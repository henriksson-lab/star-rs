#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original `SoloFeature::umiArrayCorrect_Graph` at STAR/source/SoloFeature_collapseUMI_Graph.cpp:16. Args: nU0: uint32, umiArr: uintUMI, readInfoRec: bool, nUMIyes: bool, umiCorr: unordered_map <uintUMI,uintUMI>"]
pub fn solofeature_collapseumi_graph_l16_solofeature_umiarraycorrect_graph(
    p_solo: &crate::parameters_solo::ParametersSolo,
    n_u0: u32,
    umi_arr: &mut [u32],
    umi_array_stride: u32,
    read_info_rec: bool,
    n_umi_yes: bool,
    umi_corr: &mut std::collections::BTreeMap<u32, u32>,
) -> u32 {
    let stride = umi_array_stride as usize;
    let n_total = n_u0 as usize * stride;
    let n_records = n_u0 as usize;
    let mut n_u1 = n_u0;
    let mut n_u2 = n_u0;
    let mut graph_n = 0_u32;
    let mut graph_conn = Vec::<[u32; 2]>::new();
    let mut graph_comp = Vec::<u32>::new();

    for iu in (0..n_total).step_by(stride) {
        umi_arr[iu + 2] = u32::MAX;
    }

    for i in 0..n_records {
        for j in (i + 1)..n_records {
            let a = i * stride;
            let b = j * stride;
            if umi_arr[a] > umi_arr[b] {
                for k in 0..stride {
                    umi_arr.swap(a + k, b + k);
                }
            }
        }
    }

    solofeature_collapseumi_graph_l80_collapseumiwith1mmlowhalf(
        umi_arr,
        umi_array_stride,
        p_solo.umi_mask_low,
        n_u0,
        &mut n_u1,
        &mut n_u2,
        &mut graph_n,
        &mut graph_conn,
    );

    for iu in (0..n_total).step_by(stride) {
        parameterssolo_l496_parameterssolo_umiswaphalves(p_solo, &mut umi_arr[iu]);
    }

    for i in 0..n_records {
        for j in (i + 1)..n_records {
            let a = i * stride;
            let b = j * stride;
            if umi_arr[a] > umi_arr[b] {
                for k in 0..stride {
                    umi_arr.swap(a + k, b + k);
                }
            }
        }
    }

    solofeature_collapseumi_graph_l80_collapseumiwith1mmlowhalf(
        umi_arr,
        umi_array_stride,
        p_solo.umi_mask_low,
        n_u0,
        &mut n_u1,
        &mut n_u2,
        &mut graph_n,
        &mut graph_conn,
    );

    let n_conn_comp = solofeature_collapseumi_graph_l142_graphnumberofconnectedcomponents(
        graph_n,
        &graph_conn,
        &mut graph_comp,
    );
    n_u1 = n_u1.wrapping_add(n_conn_comp);

    if read_info_rec {
        for ii in 0..graph_comp.len() {
            if graph_comp[ii] == u32::MAX {
                graph_comp[ii] = ii as u32;
            }
        }

        let bit_top_mask = !(1_u32 << 31);
        let mut umi_best = vec![[0_u32, 0_u32]; graph_n as usize];
        let mut umi_corr_color = std::collections::BTreeMap::<u32, u32>::new();
        for iu in (0..n_total).step_by(stride) {
            parameterssolo_l496_parameterssolo_umiswaphalves(p_solo, &mut umi_arr[iu]);
            if umi_arr[iu + 2] == u32::MAX {
                continue;
            }
            let color1 = graph_comp[umi_arr[iu + 2] as usize];
            let count1 = umi_arr[iu + 1] & bit_top_mask;
            if umi_best[color1 as usize][0] < count1 {
                umi_best[color1 as usize][0] = count1;
                umi_best[color1 as usize][1] = umi_arr[iu];
            }
            umi_corr_color.insert(umi_arr[iu], color1);
        }

        for iu in (0..n_total).step_by(stride) {
            let umi = umi_arr[iu];
            if let Some(&color) = umi_corr_color.get(&umi) {
                umi_corr.insert(umi, umi_best[color as usize][1]);
            }
        }
    }

    if n_umi_yes { n_u1 } else { 0 }
}

#[doc = "Original `collapseUMIwith1MMlowHalf` at STAR/source/SoloFeature_collapseUMI_Graph.cpp:80. Args: umiArr: uint32, umiArrayStride: uint32, umiMaskLow: uint32, nU0: uint32, nU1: uint32, nU2: uint32, nC: uint32, vC: vector<array<uint32,2>>"]
pub fn solofeature_collapseumi_graph_l80_collapseumiwith1mmlowhalf(
    umi_arr: &mut [u32],
    umi_array_stride: u32,
    umi_mask_low: u32,
    n_u0: u32,
    n_u1: &mut u32,
    n_u2: &mut u32,
    n_c: &mut u32,
    v_c: &mut Vec<[u32; 2]>,
) {
    let bit_top = 1_u32 << 31;
    let bit_top_mask = !bit_top;
    let stride = umi_array_stride as usize;
    let n_total = stride * n_u0 as usize;

    for iu in (0..n_total).step_by(stride) {
        for iuu in ((iu + stride)..n_total).step_by(stride) {
            let uu_xor = umi_arr[iu] ^ umi_arr[iuu];
            if uu_xor > umi_mask_low {
                break;
            }

            if uu_xor != 0 {
                let shift = (uu_xor.trailing_zeros() / 2) * 2;
                if (uu_xor >> shift) > 3 {
                    continue;
                }
            }

            if umi_arr[iu + 2] == u32::MAX && umi_arr[iuu + 2] == u32::MAX {
                umi_arr[iu + 2] = *n_c;
                umi_arr[iuu + 2] = *n_c;
                *n_c = n_c.wrapping_add(1);
                *n_u1 = n_u1.wrapping_sub(2);
            } else if umi_arr[iu + 2] == u32::MAX {
                umi_arr[iu + 2] = umi_arr[iuu + 2];
                *n_u1 = n_u1.wrapping_sub(1);
            } else if umi_arr[iuu + 2] == u32::MAX {
                umi_arr[iuu + 2] = umi_arr[iu + 2];
                *n_u1 = n_u1.wrapping_sub(1);
            } else if umi_arr[iuu + 2] != umi_arr[iu + 2] {
                v_c.push([umi_arr[iu + 2], umi_arr[iuu + 2]]);
            }

            if (umi_arr[iuu + 1] & bit_top) == 0
                && (umi_arr[iu + 1] & bit_top_mask)
                    > 2_u32
                        .wrapping_mul(umi_arr[iuu + 1] & bit_top_mask)
                        .wrapping_sub(1)
            {
                umi_arr[iuu + 1] |= bit_top;
                *n_u2 = n_u2.wrapping_sub(1);
            } else if (umi_arr[iu + 1] & bit_top) == 0
                && (umi_arr[iuu + 1] & bit_top_mask)
                    > 2_u32
                        .wrapping_mul(umi_arr[iu + 1] & bit_top_mask)
                        .wrapping_sub(1)
            {
                umi_arr[iu + 1] |= bit_top;
                *n_u2 = n_u2.wrapping_sub(1);
            }
        }
    }
}

#[doc = "Original `graphDepthFirstSearch` at STAR/source/SoloFeature_collapseUMI_Graph.cpp:132. Args: n: uint32, nodeEdges: vector<vector<uint32>>, nodeColor: vector <uint32>"]
pub fn solofeature_collapseumi_graph_l132_graphdepthfirstsearch(
    n: u32,
    node_edges: &[Vec<u32>],
    node_color: &mut [u32],
) {
    for &nn in &node_edges[n as usize] {
        if node_color[nn as usize] == u32::MAX {
            node_color[nn as usize] = node_color[n as usize];
            solofeature_collapseumi_graph_l132_graphdepthfirstsearch(nn, node_edges, node_color);
        }
    }
}

#[doc = "Original `graphNumberOfConnectedComponents` at STAR/source/SoloFeature_collapseUMI_Graph.cpp:142. Args: N: uint32, V: vector<array<uint32,2>>, nodeColor: vector<uint32>"]
pub fn solofeature_collapseumi_graph_l142_graphnumberofconnectedcomponents(
    n: u32,
    v: &[[u32; 2]],
    node_color: &mut Vec<u32>,
) -> u32 {
    node_color.clear();
    node_color.resize(n as usize, u32::MAX);

    if v.is_empty() {
        return n;
    }

    let mut node_edges = vec![Vec::<u32>::new(); n as usize];
    for edge in v {
        node_edges[edge[0] as usize].push(edge[1]);
        node_edges[edge[1] as usize].push(edge[0]);
    }

    let mut n_conn_comp = 0_u32;
    for ii in 0..n {
        if node_edges[ii as usize].is_empty() {
            n_conn_comp += 1;
        } else if node_color[ii as usize] == u32::MAX {
            n_conn_comp += 1;
            node_color[ii as usize] = ii;
            solofeature_collapseumi_graph_l132_graphdepthfirstsearch(ii, &node_edges, node_color);
        }
    }
    n_conn_comp
}
