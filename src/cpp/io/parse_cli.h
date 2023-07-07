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

    program.add_argument("-s", "--speciesmat")
            .required()
            .help("specify the species tree distance matrix in PHYLIP format");

    program.add_argument("-a", "--alignmat")
            .required()
            .help("specify the gene alignment distance matrix in PHYLIP format");

    //program.add_argument("-f", "--families")
    //        .help("specify the families file");

    program.add_argument("-m", "--mapping")
            .default_value(std::string{""})
            .help("specify the mapping between species names and locus names");

    program.add_argument("-l")
            .default_value(false)
            .implicit_value(true)
            .help("if set, mapping file maps locus names to species names");

    program.add_argument("-d", "--delimiter")
            .default_value(std::string{":"})
            .help("mapping separator");

    program.add_argument("-p", "--prefix")
            .default_value(std::string{"/"})
            .help("specify the output prefix");

    return program;
}

auto parse_and_get(argparse::ArgumentParser& program, int argc, char *argv[]) {
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        exit(1);
    }
    std::vector<std::string> res{};
    res.push_back(program.get("-s"));
    res.push_back(program.get("-a"));
    res.push_back(program.get("-p"));

    return res;
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

template<typename T>
T get(argparse::ArgumentParser& program, std::string const& opt) {
    return program.get<T>(opt);
}

std::string getS(argparse::ArgumentParser& program) {
    return program.get("-s");
}

std::string getA(argparse::ArgumentParser &program) {
    return program.get("-a");
}

std::string getP(argparse::ArgumentParser &program) {
    return program.get("-p");
}

std::tuple<std::string, bool, std::string> getMappingConfig(argparse::ArgumentParser &program) {
    return {program.get("-m"), program.get<bool>("-l"), program.get("-d")};
}

#endif //BA_PARSE_CLI_H
