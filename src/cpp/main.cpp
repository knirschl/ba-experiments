//
// Created by knirschl on 18.04.23.
//

#include <sstream>
#include <iomanip>
#include "io/parse_cli.h"
#include "io/parse_file.h"
#include "io/write_file.h"
#include "nj/tree.h"
#include "nj/NJSimple.h"
#include "misc/meta.h"

static const std::string CORRECTION_IDENT{"S~G"};
static const std::string MATRX_PHY_FILE{".matrix.phy"};
static const std::string GTREE_NWK_FILE{".geneTree.newick"};

/**
 * Resets tree to a set of leafs.
 *
 * @param alignment_ids locus identifier of every leaf
 * @return the forest of leafs
 */
auto reset(const std::vector<std::string> &alignment_ids) {
    // create tree
    std::shared_ptr<Tree> tree = std::make_shared<Tree>();
    tree->make_leafs(alignment_ids);
    return tree;
}

dist_matrix_t correct_matrix(const double scale, const dist_matrix_t &species_tree_mat,
                             const dist_matrix_t &alignment_mat,
                             const vector_t<std::pair<int, int>> &speciation_pairs) {
    dist_matrix_t corrected_matrix{alignment_mat};
    for (auto &pair : speciation_pairs) {
        // lowest common ancestor was not dup
        std::string locus1 = idx2leafname[pair.first];
        std::string locus2 = idx2leafname[pair.second];
        int ai = leafname2matidx[locus1];
        int aj = leafname2matidx[locus2];
        int si = groupname2matidx[leafname2groupname[locus1]];
        int sj = groupname2matidx[leafname2groupname[locus2]];

        corrected_matrix[ai][aj] += scale * species_tree_mat[si][sj];
        corrected_matrix[ai][aj] /= scale + 1;
        corrected_matrix[aj][ai] = corrected_matrix[ai][aj];
    }
    return corrected_matrix;
}

/*
-s
/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_DL_s100_f100_sites200_GTR_bl1.0_d1.0_l1.0_t0.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed1650734/species_trees/speciesTree.matrix.phy
-a
/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_DL_s100_f100_sites200_GTR_bl1.0_d1.0_l1.0_t0.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed1650734/families/family_042/alignment.msa.matrix.phy
-p
/home/fili/Documents/KIT/2023/BA/code/src/test/output/ba.exp.
-m
/home/fili/Documents/KIT/2023/BA/code/output/families/ssim_DL_s100_f100_sites200_GTR_bl1.0_d1.0_l1.0_t0.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed1650734/families/family_042/mappings/mapping.link
 */
/**
 * Main entry point of this project. Reads the input in and computes either a corrected
 * gene alignment distance matrix or a neighbor-joined gene tree.
 *
 * @param argc  argument count
 * @param argv  arguments
 * @return 0
 */
int main(int argc, char *argv[]) {
    // --- parse ---
    // get cli inputs
    auto cli_parser = build_parser("thesis", "0.2");
    parse(cli_parser, argc, argv);

    // read species tree
    auto species_tree_pair = parse_phylip_mat_from_file<dist_t>(get_species_matrix(cli_parser));
    auto species_tree_mat = species_tree_pair.first;
    auto species_tree_ids = species_tree_pair.second;
    int s_cnt{};
    for_each(species_tree_ids.begin(), species_tree_ids.end(),
             [&s_cnt](auto &s) { groupname2matidx.emplace(s, s_cnt++); });

    // read alignment
    auto alignment_pair = parse_phylip_mat_from_file<dist_t>(get_alignment_matrix(cli_parser));
    auto alignment_mat = alignment_pair.first;
    auto alignment_ids = alignment_pair.second;
    int a_cnt{};
    for_each(alignment_ids.begin(), alignment_ids.end(),
             [&a_cnt](auto &s) { leafname2matidx.emplace(s, a_cnt++); });

    // read mapping
    auto map_config{get_mapping_config(cli_parser)};

    // --- pre-calculate ---
    // get starting tree
    std::shared_ptr<Tree> tree;
    std::vector<int> active;
    if (has_user_specified_tree(cli_parser)) {
        tree = parse_newick_from_file(get_starting_tree(cli_parser));
        active = leaf_indices;
    } else {
        tree = reset(alignment_ids);
        active = leaf_indices;
        // NJ gene tree with only alignment matrix (0S+G)
        neighborJoining<>(alignment_mat, tree, active);
    }
    //std::cout << "Start tree : " << tree->to_newick() << "\n" << tree->node_info() << "\n";
    tree->reroot_APro();
    auto speciation_pairs{tree->get_speciation_pairs()};
    //std::cout << "A-Pro tree: " << tree->to_newick() << "\n" << tree->node_info() << "\n";

    // --- calculate ---
    // NJ gene tree with corrected values
    double scales[] = { /*
            0, 0.5, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,1.9, 2.0, 2.05,
            2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65,
            2.7, 2.75, 2.85, 2.95, 3.1, 3.2, 3.3, 3.4, 3.5, 4, 4.5, 5, 7.5, 10
            */
            0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55,
            1.65, 1.75, 1.85, .95, 2.8, 2.9, 3.0, 3.05, 3.15, 3.25, 3.35, 3.45, 3.55, 3.6, 3.65,
            3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.55,
            4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95};
    // parallelize?: #pragma omp parallel for
    for (double scale: scales) {
        // correct with scaling
        dist_matrix_t corrected_matrix{
                correct_matrix(scale, species_tree_mat, alignment_mat, speciation_pairs)};

        // convert double to string without trailing zeros and with a specified precision
        std::ostringstream oss;
        oss << std::setprecision(8) << std::noshowpoint << scale << CORRECTION_IDENT;
        auto scale_id{oss.str()};

        // output/compute depending on "-c"
        bool success{true};;
        if (int c{get_c(cli_parser)}) {
            // output matrix
            success = write_phylip(corrected_matrix, alignment_ids,
                                   get_output_prefix(cli_parser).append(scale_id).append(
                                           MATRX_PHY_FILE));
            if (c == 2) {
                // no tree
                continue;
            }
        }
        // compute NJ tree and output
        tree = reset(alignment_ids);
        active = leaf_indices;
        neighborJoining(corrected_matrix, tree, active);
        std::cout << "Neighbor-joined tree (" << scale_id << "): " << tree->to_newick() << '\n';
        success &= write_newick(*tree, get_output_prefix(cli_parser).append(scale_id).append(
                GTREE_NWK_FILE));
        if (!success) {
            std::cout << "Failed to complete " << scale_id << " task\n";
        }
    }

    return 0;
}


