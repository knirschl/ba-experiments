//
// Created by knirschl on 18.04.23.
//

#include <sstream>
#include <iomanip>
#include "io/parse_cli.h"
#include "io/parse_file.h"
#include "io/write_file.h"
#include "io/exec_cmd.h"
#include "nj/tree.h"
#include "nj/NJSimple.h"
#include "nj/metadata.h"

static const std::string CORRECTION_IDENT{"S~G"};
static const std::string MATRX_PHY_FILE{".matrix.phy"};
static const std::string GTREE_NWK_FILE{".geneTree.newick"};
static const std::string START_TREE_NWK_FILE{".startGeneTree.newick"};
static const std::string MADROOT_BIN{"/hits/basement/cme/knirsch/github/"
                                     "MADroot/bin/madRoot"};

/**
 * Resets tree to a set of leafs.
 *
 * @param alignment_ids locus identifier of every leaf
 * @return the forest of leafs
 */
auto reset(const std::vector<std::string> &alignment_ids, const metadata &mdata) {
    // create tree
    std::shared_ptr<Tree> tree = std::make_shared<Tree>(mdata);
    tree->make_leafs(alignment_ids);
    return tree;
}

/**
 * Correct the distances between every two entries in {@param alignment_mat} with the corresponding
 * scaled entry in {@param species_tree_mat}.
 *
 * @param scale
 * @param species_tree_mat
 * @param alignment_mat
 * @param mdata
 * @return the corrected matrix
 */
dist_matrix_t correct_matrix(const double scale, const dist_matrix_t &species_tree_mat,
                             const dist_matrix_t &alignment_mat, metadata &mdata) {
    const size_t n{alignment_mat.size()};
    dist_matrix_t corrected_matrix(n, dist_vector_t(n));
    for (int i{}; i < n; i++) {
        for (int j{}; j < n; j++) {
            std::string locus1 = mdata.idx2nodename[i];
            std::string locus2 = mdata.idx2nodename[j];
            int ai = mdata.leafname2matidx[locus1]; // probably always same as i, same for aj
            int aj = mdata.leafname2matidx[locus2];
            int si = glob_mdata.groupname2matidx[mdata.leafname2groupname[locus1]];
            int sj = glob_mdata.groupname2matidx[mdata.leafname2groupname[locus2]];

            corrected_matrix[ai][aj] += scale * species_tree_mat[si][sj];
            corrected_matrix[ai][aj] /= scale + 1;
            corrected_matrix[aj][ai] = corrected_matrix[ai][aj];
        }
    }
    return corrected_matrix;
}

/**
 * Correct the distances between two entries in {@param alignment_mat} with the corresponding scaled
 * entry in {@param species_tree_mat} for every pair in {@param speciation_pairs}.
 *
 * @param scale
 * @param species_tree_mat
 * @param alignment_mat
 * @param speciation_pairs
 * @param mdata
 * @return the corrected matrix
 */
dist_matrix_t correct_matrix(const double scale, const dist_matrix_t &species_tree_mat,
                             const dist_matrix_t &alignment_mat,
                             const vector_t<std::pair<int, int>> &speciation_pairs,
                             metadata &mdata) {
    dist_matrix_t corrected_matrix{alignment_mat};
    for (auto &pair: speciation_pairs) {
        // lowest common ancestor was not dup
        std::string locus1 = mdata.idx2nodename[pair.first];
        std::string locus2 = mdata.idx2nodename[pair.second];
        int ai = mdata.leafname2matidx[locus1];
        int aj = mdata.leafname2matidx[locus2];
        int si = glob_mdata.groupname2matidx[mdata.leafname2groupname[locus1]];
        int sj = glob_mdata.groupname2matidx[mdata.leafname2groupname[locus2]];

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
             [&s_cnt](auto &s) { glob_mdata.groupname2matidx.emplace(s, s_cnt++); });

    // read alignment
    auto alignment_pair = parse_phylip_mat_from_file<dist_t>(get_alignment_matrix(cli_parser));
    auto alignment_mat = alignment_pair.first;
    auto alignment_ids = alignment_pair.second;
    if (alignment_ids.size() < 4) {
        std::cout << "Too few loci: " << std::to_string(alignment_ids.size()) << "\n";
        return 0;
    }
    int a_cnt{};
    metadata starting_tree_mdata;
    for_each(alignment_ids.begin(), alignment_ids.end(),
             [&a_cnt, &starting_tree_mdata](auto &s) {
                 starting_tree_mdata.leafname2matidx.emplace(s, a_cnt++);
             });

    // read mapping
    auto map_config{get_mapping_config(cli_parser)};
    std::string out_prefix{get_output_prefix(cli_parser)};
    // fill map
    if (get<0>(map_config).empty()) {
        // no mapping provided
        for (int i{}; i < alignment_ids.size(); i++) {
            starting_tree_mdata.leafname2groupname.emplace(alignment_ids[i], species_tree_ids[i]);
        }
    } else {
        starting_tree_mdata.leafname2groupname = parse_mapping_from_cfg(map_config);
    }

    // --- pre-calculate ---
    // get starting tree
    std::shared_ptr<Tree> tree_tagged;
    if (has_user_specified_tree(cli_parser)) {
        tree_tagged = parse_newick_from_file(get_starting_tree(cli_parser), starting_tree_mdata);
        out_prefix.append("u");
    } else {
        tree_tagged = reset(alignment_ids, starting_tree_mdata);
        double spec_mat_scale{0.0}; // TODO change to 2.4
        dist_matrix_t start_mat{correct_matrix(spec_mat_scale, species_tree_mat, alignment_mat,
                                               tree_tagged->mdata)};
        // NJ gene tree with only alignment matrix (0S+G)
        neighborJoining<>(start_mat, tree_tagged, tree_tagged->mdata.leaf_indices);
        std::ostringstream oss;
        oss << std::setprecision(4) << std::noshowpoint << spec_mat_scale;
        out_prefix.append(oss.str());
    }
    //std::cout << "Start tree : " << tree_tagged->to_newick() << "\n" << tree_tagged->node_info() << "\n";
    switch (get_algo(cli_parser)) {
        case 0:
            tree_tagged->reroot_APro();
            out_prefix.append("a.");
            break;
        case 1:
            //tree_tagged->tag_APro(tree_tagged->reroot_MAD());
            write_newick(*tree_tagged, out_prefix + START_TREE_NWK_FILE);
            std::make_shared<Tree>(
                    exec((MADROOT_BIN + " " + out_prefix + START_TREE_NWK_FILE).c_str()),
                    tree_tagged->mdata);
            out_prefix.append("m.");
            break;
        default:
            // as dup is initialized with false this is the same as S+G (just a bit slower)
            out_prefix.append("+.");
            break;
    }
    auto speciation_pairs{tree_tagged->get_speciation_pairs()};
    //std::cout << "Tagged tree: " << tree_tagged->to_newick() << "\n" << tree_tagged->node_info() << "\n";

    // --- calculate ---
    // NJ gene tree with corrected values
    double scales[] = {
            0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
            0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55,
            1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35,
            2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
            3.8, 4, 4.25, 4.5, 4.75, 5, 5.5, 6, 7, 8, 9, 10, 25, 50, 100};
    int c{get_c(cli_parser)};
//#pragma omp parallel for default(shared)
    for (double scale: scales) {
        // correct with scaling
        dist_matrix_t corrected_matrix{
                correct_matrix(scale, species_tree_mat, alignment_mat, speciation_pairs,
                               tree_tagged->mdata)};

        // convert double to string without trailing zeros and with a specified precision
        std::ostringstream oss;
        oss << std::setprecision(8) << std::noshowpoint << scale << CORRECTION_IDENT;
        auto scale_id{oss.str()};

        // output/compute depending on "-c"
        bool success{true};
        if (c) {
            // output matrix
            success = write_phylip(corrected_matrix, alignment_ids,
                                   out_prefix + scale_id.append(
                                           MATRX_PHY_FILE));
            if (c == 2) {
                // no tree
                continue;
            }
        }
        // TODO fix. prob because of the global maps -> prob fixed
        // compute NJ tree_tagged and output
        std::shared_ptr<Tree> tree_corrected{reset(alignment_ids, tree_tagged->mdata)};
        neighborJoining(corrected_matrix, tree_corrected, tree_corrected->mdata.leaf_indices);
        //std::cout << "Neighbor-joined tree_tagged (" << scale_id << "): " << tree_tagged->to_newick() << '\n';
        success &= write_newick(*tree_corrected, out_prefix + scale_id.append(
                GTREE_NWK_FILE));
        if (!success) {
            std::cout << "Failed to complete " << scale_id << " task\n";
        }
    }

    return 0;
}
