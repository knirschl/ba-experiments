//
// Created by fili on 23.08.23.
//

#ifndef BA_EXEC_CMD_H
#define BA_EXEC_CMD_H

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

static const int BUFFER_SIZE{128};

int exec(const char *cmd, std::string &output) {
    std::array<char, BUFFER_SIZE> buffer{};

    auto pipe(popen(cmd, "r"));
    if (!pipe) {
        throw std::runtime_error("call to popen() failed!");
    }
    size_t count;
    while (!feof(pipe)) {
        while ((count = fread(buffer.data(), 1, BUFFER_SIZE, pipe))) {
            output.insert(output.end(), std::begin(buffer), std::next(std::begin(buffer), count));
        }
    }

    return pclose(pipe);
}

int exec(const std::string &cmd, std::string &output) {
    return exec(cmd.c_str(), output);
}

#endif //BA_EXEC_CMD_H
