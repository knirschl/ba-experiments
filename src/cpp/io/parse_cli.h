//
// Created by knirschl on 05.06.23.
//

#ifndef BA_PARSE_CLI_H
#define BA_PARSE_CLI_H

#include <argparse/argparse.hpp>

// ./thesis --speciestree "species_tree.matrix.phy" --alignmat "alignment.msa.matrix.phy" -p "output-dir/"

auto build_parser(std::string const& name, std::string const& version) {
    argparse::ArgumentParser program{name, version};

    program.add_description("");
    program.add_epilog("");

    program.add_argument("-s", "--species-mat")
            .required()
            .help("specify the species tree distance matrix in PHYLIP format");

    program.add_argument("-a", "--align-mat")
            .required()
            .help("specify the gene alignment distance matrix in PHYLIP format");

    program.add_argument("-t", "--starting-tree")
            .help("specify the gene tree to compute dup/loss on (NEWICK format");

    // r := tag and reroot algorithm
    program.add_argument("-r", "--tag-algo")
            .default_value(0)
            .scan<'i', int>()
            .help("0: APro, 1: MAD, 2: none (everything gets corrected (S+G))");

    // c := compute
    program.add_argument("-c", "--compute-tree")
            .default_value(2)
            .scan<'i', int>()
            .help("0: only tree, 1: both, 2: only matrix");

    program.add_argument("-p", "--prefix")
            .default_value(std::string{"/"})
            .help("specify the output prefix");

    program.add_argument("-m", "--map")
            .default_value(std::string{""})
            .help("specify the mapping between species names and locus names");

    // direction s->l or l->s
    program.add_argument("-l")
            .default_value(false)
            .implicit_value(true)
            .help("if set, mapping file maps locus names to species names");

    // delimiter between species and locus names
    program.add_argument("-d", "--delimiter")
            .default_value(std::string{":"})
            .help("mapping separator");

    return program;
}

auto parse(argparse::ArgumentParser& program, int argc, char *argv[]) {
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    return 0;
}

std::string get_species_matrix(argparse::ArgumentParser& program) {
    return program.get("-s");
}

std::string get_alignment_matrix(argparse::ArgumentParser &program) {
    return program.get("-a");
}

bool has_user_specified_tree(argparse::ArgumentParser &program) {
    return program.is_used("-t");
}

std::string get_starting_tree(argparse::ArgumentParser &program) {
    return program.get("-t");
}

int get_algo(argparse::ArgumentParser &program) {
    return program.get<int>("-r");
}

int get_c(argparse::ArgumentParser &program) {
    return program.get<int>("-c");
}

std::string get_output_prefix(argparse::ArgumentParser &program) {
    return program.get("-p");
}

std::tuple<std::string, bool, std::string> get_mapping_config(argparse::ArgumentParser &program) {
    return {program.get("-m"), program.get<bool>("-l"), program.get("-d")};
}

#endif //BA_PARSE_CLI_H
