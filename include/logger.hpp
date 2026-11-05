#pragma once

#include <iostream>
#include <string>
#include <format>

using namespace std;

class Logger
{
private:
    static constexpr const char *RESET = "\033[0m";
    static constexpr const char *RED = "\033[31m";
    static constexpr const char *GREEN = "\033[32m";
    static constexpr const char *YELLOW = "\033[33m";
    static constexpr const char *BLUE = "\033[34m";

public:
    static void info(const string &msg)
    {
        cout << format("{}[INFO] {}{}", BLUE, msg, RESET) << endl;
    }

    static void success(const string &msg)
    {
        cout << format("{}[SUCCESS] {}{}", GREEN, msg, RESET) << endl;
    }

    static void warn(const string &msg)
    {
        cout << format("{}[WARN] {}{}", YELLOW, msg, RESET) << endl;
    }

    static void error(const string &msg)
    {
        cerr << format("{}[ERROR] {}{}", RED, msg, RESET) << endl;
    }
};